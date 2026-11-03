"""
compare_juice_facsimile.py
--------------------------
Compare JUICE vs FACSIMILE output for a given scenario.

Usage:
    python compare_juice_facsimile.py scene_03 [--tol 1e-4] [--all]

Options:
    --tol      Max relative deviation threshold (default: 1e-4)
    --floor    Concentration floor in molec/cm³ (default: 100 = the solver's atol_floor): deviations
               are measured against max(max|FAC|, floor), so species the solver does not resolve
               (a few molecules per cm³) are not reported as 100 % off
    --all      Write all columns to log (default: top 100, with count warning)
    --no-stop  Continue through all groups after the first one that exceeds --tol

Paths come from FR_GUI.toml [paths] (+ FR_GUI.local.toml); without them the repository folders
(easy_root = MCM_tested_scenarios/IDL Input, output_root = …/JUICE archive/JUICE) are used.
Exit code 0 = every compared group within --tol, 1 = not (or a file is missing).
compare() returns the same verdict as a summary dict (used by run_benchmark.py).
Files compared (in order):
  1. s-values  : settings (KMT, J-values, COSX, T, M, …)
  2. k-values  : reaction rate constants  (needs *.kv)
  3. c-values  : concentrations (ODE species)
  4. v-values  : reaction turnover rates  (needs *.kv)
  5. pd-values : production/destruction   (needs *.pd)
Stops and reports the first column whose max relative deviation exceeds --tol.
"""

import sys
import io
import os
import re
import argparse
import tomllib
import datetime as _dt
import numpy as np
import pandas as pd
from column_file_reader import load_column_file

# ── Reaction-string canonicaliser (mirrors easy_parser._canon_rxn_key) ───────
def _canon_rxn_str(rxn_str: str) -> str:
    """Normalise a reaction string to canonical form.

    Handles: reactant/product order, merged stoichiometry (0.5+0.5 → 1),
    coefficient notation (NA+NA → 2*NA).  Mirrors EasyParser._canon_rxn_key
    so that JUICE and FACSIMILE reaction labels can be matched regardless of
    the order in which reactants or products were written.

    Examples:
        'O3+NO-->NO2'         →  'NO+O3-->NO2'
        'NO+O3-->NO2'         →  'NO+O3-->NO2'
        'N2O5-->NA+NA'        →  'N2O5-->2*NA'
        'N2O5-->2*NA'         →  'N2O5-->2*NA'
    """
    def _side(s: str) -> str:
        coeff_map: dict = {}
        for term in s.upper().replace(' ', '').split('+'):
            if not term:
                continue
            m = re.match(r'^([\d.]+)\*(.+)$', term)
            if m:
                try:
                    c = float(m.group(1))
                except ValueError:
                    c = 1.0
                sp = m.group(2)
            else:
                c, sp = 1.0, term
            coeff_map[sp] = coeff_map.get(sp, 0.0) + c
        parts = []
        for sp in sorted(coeff_map):
            c = coeff_map[sp]
            c_int = int(round(c))
            if abs(c - c_int) < 1e-9:
                parts.append(sp if c_int == 1 else f'{c_int}*{sp}')
            else:
                parts.append(f'{c:.6g}*{sp}')
        return '+'.join(parts)

    if '-->' in rxn_str:
        lhs, rhs = rxn_str.split('-->', 1)
        rhs_norm = rhs.strip().upper().replace(' ', '')
        # Treat explicit "(NONE)", "NONE", and empty string as identical (no products)
        if rhs_norm in ('(NONE)', 'NONE', ''):
            return _side(lhs) + '-->'
        return _side(lhs) + '-->' + _side(rhs)
    return rxn_str.upper().replace(' ', '')


def _canon_kv_label(label: str) -> str:
    """Normalise a kv-column label such as 'k_R00006:O3+NO-->NO2'.

    Extracts the reaction part (after optional 'k_'/'v_' prefix and
    optional 'RNNNnnn:' index), canonicalises it, and returns a
    prefix-free canonical reaction string for cross-model matching.
    """
    # Strip k_/v_ prefix and optional Rnnnnn: index
    s = re.sub(r'^[kv]_', '', label)
    s = re.sub(r'^R\d+:', '', s)
    return _canon_rxn_str(s)


# ── Settings heuristic (mirrors results_viewer._is_setting_col) ────────────
_SETTING_PREFIXES = ('KMT', 'KRO', 'KAP', 'KBP', 'KFP', 'KCH', 'KNO', 'KDE',
                     'KROP', 'KFPAN', 'KBPAN', 'KDEC', 'KAPHO', 'KAPNO')
_SETTING_NAMES    = {'T', 'TEMP', 'P', 'PRESS', 'M', 'HV', 'COSX',
                     'LNO2', 'MNO2', 'NNO2', 'SZ', 'SZA', 'RH', 'FLOW',
                     'DILUTE', 'JSCALE', 'DIL', 'N2', 'O2'}

def _is_setting(col):
    cu = col.upper()
    if cu in _SETTING_NAMES:
        return True
    if any(cu.startswith(p) for p in _SETTING_PREFIXES):
        return True
    if cu.startswith('J') and len(cu) > 1 and (cu[1].isdigit() or '_CS' in cu):
        return True
    return False


def _trim_boundary_rows(df):
    """Remove the duplicate last-timestamp row written by FACSIMILE at end-of-run.

    FACSIMILE always writes a duplicate final row (same timestamp, wrong s-values)
    because the WHENEVER fertig=1 block calls CALL print after the last ENZ row
    has already been read for the next (non-existent) step.

    The initial conditions row (t=t_0, c=INITIAL[*]) is NOT removed here —
    it is handled automatically by the time-range alignment in compare():
    JUICE skips its own initial row (iloc[1:]), so t_min = t_1, and the
    FACSIMILE t_0 row falls outside the overlapping range.

    New FACSIMILE files (after the LABEL 63/64 fix in easy_model2facsimile.pro)
    no longer write the initial conditions row; the duplicate last row remains
    and is still removed here.
    """
    if len(df) < 2:
        return df
    return df[~df.index.duplicated(keep='first')]  # drop duplicate last timestamp


# ── FACSIMILE .cs parser ───────────────────────────────────────────────────
def load_facsimile_cs(path):
    cols = {}
    data_lines = []
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('COLUMN'):
                m = re.match(r'COLUMN\s+(\d+)=(.+)', s)
                if m:
                    idx  = int(m.group(1)) - 1
                    name = m.group(2).split('<')[0].strip()
                    cols[idx] = 'TIME_s' if idx == 0 else name
            elif s and s[0].isdigit():
                fixed = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', s)
                data_lines.append(fixed)
    if not cols or not data_lines:
        raise ValueError(f"Could not parse FACSIMILE file: {path}")
    headers = [cols.get(i, f'COL{i}') for i in range(max(cols)+1)]
    rows = [[float(v) for v in r.split()] for r in data_lines]
    df = pd.DataFrame(rows, columns=headers[:len(rows[0])])
    df = df.set_index('TIME_s')
    return _trim_boundary_rows(df)


# ── FACSIMILE .kv parser ──────────────────────────────────────────────────
def load_facsimile_kv(path):
    """Load a FACSIMILE .kv file into two DataFrames (k-values, v-values).

    Column headers have the form:
        COLUMN N=k<i,0> !U{reactants}-->{products}!N   (rate constants)
        COLUMN N=v<i,0> !U{reactants}-->{products}!N   (turnover rates)
    Reaction names are extracted from !U...!N and canonicalised so they can
    be matched against JUICE kv columns regardless of reactant order.
    """
    k_cols, v_cols = {}, {}
    data_lines = []
    in_data = False
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('COLUMN'):
                m = re.match(r'COLUMN\s+(\d+)=([kv])<\d+[^>]*>\s*!U(.*?)!N', s, re.I)
                if m:
                    idx   = int(m.group(1)) - 1
                    kv    = m.group(2).lower()
                    rxn   = _canon_rxn_str(m.group(3).strip())
                    if kv == 'k':
                        k_cols[idx] = rxn
                    else:
                        v_cols[idx] = rxn
            elif '&' in s:
                in_data = True
            elif in_data and s and s[0].isdigit():
                fixed = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', s)
                data_lines.append(fixed)
    if not data_lines:
        raise ValueError(f"No data in FACSIMILE kv file: {path}")
    all_cols = {0: 'TIME_s'}
    all_cols.update({i: f'k_{n}' for i, n in k_cols.items()})
    all_cols.update({i: f'v_{n}' for i, n in v_cols.items()})
    headers = [all_cols.get(i, f'COL{i}') for i in range(max(all_cols)+1)]
    rows = [[float(v) for v in r.split()] for r in data_lines]
    df = pd.DataFrame(rows, columns=headers[:len(rows[0])]).set_index('TIME_s')
    df = _trim_boundary_rows(df)
    k_df = df[[c for c in df.columns if c.startswith('k_')]]
    v_df = df[[c for c in df.columns if c.startswith('v_')]]
    k_df.columns = [c[2:] for c in k_df.columns]
    v_df.columns = [c[2:] for c in v_df.columns]
    return k_df, v_df


# ── FACSIMILE .pd parser ──────────────────────────────────────────────────
def load_facsimile_pd(path):
    """Load a FACSIMILE .pd file → DataFrame with P_{sp} and D_{sp} columns."""
    p_cols, d_cols = {}, {}
    data_lines = []
    in_data = False
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('COLUMN'):
                m = re.match(r'COLUMN\s+(\d+)=([pd])<\d+>\s+(\S+)', s, re.I)
                if m:
                    idx  = int(m.group(1)) - 1
                    pd_  = m.group(2).lower()
                    sp   = m.group(3).upper()
                    if pd_ == 'p':
                        p_cols[idx] = sp
                    else:
                        d_cols[idx] = sp
            elif '&' in s:
                in_data = True
            elif in_data and s and s[0].isdigit():
                fixed = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', s)
                data_lines.append(fixed)
    if not data_lines:
        raise ValueError(f"No data in FACSIMILE pd file: {path}")
    all_cols = {0: 'TIME_s'}
    all_cols.update({i: f'P_{n}' for i, n in p_cols.items()})
    all_cols.update({i: f'D_{n}' for i, n in d_cols.items()})
    headers = [all_cols.get(i, f'COL{i}') for i in range(max(all_cols)+1)]
    rows = [[float(v) for v in r.split()] for r in data_lines]
    df = pd.DataFrame(rows, columns=headers[:len(rows[0])]).set_index('TIME_s')
    return _trim_boundary_rows(df)


def _show_juice_only_kv(kv_jui, k_jui, v_jui):
    for label, jui_map in (('3. k-values', k_jui), ('4. v-values', v_jui)):
        print(f"\n{'─'*60}")
        print(f"  {label} ({len(jui_map)} reactions)")
        print(f"{'-'*60}")
        print("  (no FACSIMILE kv reference — JUICE values shown)")
        for rxn, col in list(jui_map.items())[:5]:
            vals = kv_jui[col].values
            print(f"  {rxn[:50]:50}  max={np.nanmax(vals):.3e}")
        if len(jui_map) > 5:
            print(f"  ... ({len(jui_map)-5} more)")


def _show_juice_only_pd(pd_jui, p_jui):
    print(f"\n{'─'*60}")
    print(f"  5. pd-values ({len(p_jui)} species)")
    print(f"{'-'*60}")
    print("  (no FACSIMILE pd reference — JUICE values shown)")
    for sp in [c[2:] for c in p_jui[:5]]:
        p = pd_jui[f'P_{sp}'].values
        d = pd_jui[f'D_{sp}'].values
        print(f"  {sp:15}  P_max={np.nanmax(p):.3e}  D_max={np.nanmax(d):.3e}")


# ── Main comparison ────────────────────────────────────────────────────────
# molec/cm³ — the solver's absolute tolerance (FR_GUI.toml [solver] atol_floor = 100, ≈ 4e-6 ppt):
# differences below it are not resolved by the solver in either model, so they are not compared.
DEFAULT_FLOOR = 100.0


def _toml_path() -> str:
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "FR_GUI.toml")


def compare(scenario: str, tol: float, show_all: bool = False, no_stop: bool = False,
            juice_dir_override: str = None, fac_dir_override: str = None,
            floor: float = DEFAULT_FLOOR) -> dict:
    """Compare one scenario and print the report.  Returns a summary:

        {"ok": bool, "groups": {label: {"n": int, "n_exceed": int,
                                        "worst": [(max_rel, column), ...]}}}

    `floor` (molec/cm³) is the smallest scale a concentration's deviation is measured against:
    max|JUICE − FAC| / max(max|FAC|, floor).  Species that never rise above the solver's absolute
    tolerance (SO3, Cl at 1e-10 … 1e-170, late MCM products at a few molec/cm³) then count as
    agreeing instead of as 100 % off.
    Settings (J values, rate constants) are compared without a floor."""
    # ── Paths: [paths] of the TOML (+ .local.toml), defaulting to the repository folders ──
    from scenario_prep import load_paths
    paths = load_paths(_toml_path())
    easy_root   = paths.get("easy_root",   "")
    output_root = paths.get("output_root", "")
    summary = {"ok": True, "groups": {}}

    juice_dir = juice_dir_override or os.path.join(output_root, scenario)
    easy_dir  = fac_dir_override   or os.path.join(easy_root,  scenario)

    cs_fac_path = os.path.join(easy_dir,  f"{scenario}.cs")
    # Support both naming conventions.
    # Prefer the newer file when both exist (elog-based runs use _juice suffix).
    def _find_juice_file(suffix):
        plain = os.path.join(juice_dir, f"{scenario}{suffix}")
        juice = os.path.join(juice_dir, f"{scenario}_juice{suffix}")
        if os.path.isfile(plain) and os.path.isfile(juice):
            return juice if os.path.getmtime(juice) >= os.path.getmtime(plain) else plain
        if os.path.isfile(plain):
            return plain
        if os.path.isfile(juice):
            return juice
        return plain  # return plain so the "not found" error names the expected file

    cs_jui_path = _find_juice_file(".cs")
    kv_jui_path = _find_juice_file(".kv")
    pd_jui_path = _find_juice_file(".pd")
    kv_fac_path = os.path.join(easy_dir, f"{scenario}.kv")
    pd_fac_path = os.path.join(easy_dir, f"{scenario}.pd")

    # ── Load files ───────────────────────────────────────────────────────
    def _mtime(p):
        try:
            return _dt.datetime.fromtimestamp(os.path.getmtime(p)).strftime('%Y-%m-%d %H:%M:%S')
        except Exception:
            return '(unknown)'

    print(f"\n{'='*60}")
    print(f"Scenario : {scenario}")
    print(f"Tolerance: {tol:g} (max relative deviation)")
    print(f"{'='*60}")

    if not os.path.isfile(cs_fac_path):
        raise FileNotFoundError(f"FACSIMILE .cs not found: {cs_fac_path}")
    if not os.path.isfile(cs_jui_path):
        raise FileNotFoundError(f"JUICE .cs not found: {cs_jui_path}")

    fac = load_facsimile_cs(cs_fac_path)   # duplicate last row removed by _trim_boundary_rows
    jui = load_column_file(cs_jui_path)
    jui = jui.iloc[1:]   # JUICE row 0 = initial conditions (c=0); skip so t_ref starts at t_1

    # Align: interpolate FACSIMILE to JUICE time points (JUICE is the reference)
    t_jui = jui.index.values.astype(float)
    t_fac = fac.index.values.astype(float)
    t_min = max(t_jui.min(), t_fac.min())
    t_max = min(t_jui.max(), t_fac.max())
    mask  = (t_jui >= t_min) & (t_jui <= t_max)
    t_ref = t_jui[mask]
    jui   = jui.iloc[mask]

    # Interpolate each FACSIMILE column onto t_ref
    fac_aligned = {}
    for col in fac.columns:
        fac_aligned[col] = np.interp(t_ref, t_fac, fac[col].values.astype(float))
    fac = pd.DataFrame(fac_aligned, index=t_ref)

    if len(t_ref) == 0:
        raise ValueError("No overlapping time range between JUICE and FACSIMILE output.")

    # Case-insensitive column matching: FACSIMILE uses lowercase ENZ names (jno2),
    # JUICE uses uppercase (JNO2). Rename JUICE columns to FAC case where they
    # differ only in capitalisation so they are treated as the same column.
    _fac_lower = {c.lower(): c for c in fac.columns}
    _jui_lower = {c.lower(): c for c in jui.columns}
    _jui_rename = {_jui_lower[k]: _fac_lower[k]
                   for k in set(_fac_lower) & set(_jui_lower)
                   if _fac_lower[k] != _jui_lower[k]}
    if _jui_rename:
        jui = jui.rename(columns=_jui_rename)

    common_cols = sorted(set(fac.columns) & set(jui.columns))

    # Classify s vs c: a setting is anything defined by a function in the EASY
    # script (s_entries) that is not a k-value. Load from the JUICE elog if available.
    _s_names = set()
    try:
        _jui_stem = os.path.splitext(cs_jui_path)[0].replace('_cs', '')
        _jui_json  = _jui_stem + '.json'
        if not os.path.isfile(_jui_json):  # try without stem replacement
            _jui_json = os.path.splitext(cs_jui_path)[0] + '.json'
        if os.path.isfile(_jui_json):
            import json as _json
            with open(_jui_json) as _fh:
                _jcfg = _json.load(_fh)
            _elog = _jcfg.get('elog_path', '')
            if _elog and os.path.isfile(_elog):
                from easy_parser import EasyParser as _EP
                _mech = _EP().parse(_elog)
                _s_names = {e.name.lower() for e in _mech.s_entries}
    except Exception:
        pass

    def _is_s_col(c):
        if _s_names:
            return c.lower() in _s_names    # from actual EASY script s_entries
        return _is_setting(c)               # fallback: name-based heuristic

    s_cols = [c for c in common_cols if _is_s_col(c)]
    c_cols = [c for c in common_cols if not _is_s_col(c)]

    fac_only = sorted(set(fac.columns) - set(jui.columns))
    jui_only = sorted(set(jui.columns) - set(fac.columns))

    print(f"\nReference     : FACSIMILE  ({cs_fac_path})")
    print(f"               Date        {_mtime(cs_fac_path)}")
    print(f"Simulation    : JUICE      ({cs_jui_path})")
    print(f"               Date        {_mtime(cs_jui_path)}")
    print(f"Deviation sign: positive = JUICE > FACSIMILE")
    print(f"\nCommon columns: {len(common_cols)}  "
          f"(s={len(s_cols)}, c={len(c_cols)})")
    if fac_only:
        print(f"  *** In FACSIMILE but MISSING from JUICE: {', '.join(fac_only)}")
    if jui_only:
        print(f"  (Only in JUICE, not in FACSIMILE: {', '.join(jui_only[:10])}"
              f"{'...' if len(jui_only) > 10 else ''})")
    print(f"Time points   : {len(t_ref)}\n")

    def _compare_group(label, cols, jui_df, fac_df, floor_=0.0):
        print(f"{'-'*60}")
        print(f"  {label} ({len(cols)} columns" + (f"; floor {floor_:g} molec/cm³" if floor_ > 0 else "") + ")")
        print(f"{'-'*60}")
        worst = []
        for col in cols:
            j = jui_df[col].values.astype(float)
            f = fac_df[col].values.astype(float)
            diff = j - f
            scale = max(np.nanmax(np.abs(f)), floor_)
            rel = np.nanmax(np.abs(j)) if scale == 0 else np.nanmax(np.abs(diff)) / scale
            # Direction at the point of maximum absolute deviation
            idx_max = np.argmax(np.abs(diff))
            direction = "JUICE>REF" if diff[idx_max] > 0 else "JUICE<REF"
            mean_j = np.nanmean(j);  mean_f = np.nanmean(f)
            worst.append((rel, col, direction, mean_j, mean_f))
        # Critical variables always appear first regardless of deviation
        _CRITICAL = {'T', 'TEMP', 'P', 'PRESS', 'HV', 'RH', 'H2O', 'DIL',
                     'M', 'N2', 'O2', 'DILUTE', 'FLOW', 'SZ', 'SZA'}
        critical = sorted([(r, c, d, mj, mf) for r, c, d, mj, mf in worst
                           if c.upper() in _CRITICAL], reverse=True)
        others   = sorted([(r, c, d, mj, mf) for r, c, d, mj, mf in worst
                           if c.upper() not in _CRITICAL], reverse=True)
        worst    = critical + others
        n_total  = len(worst)
        n_exceed = sum(1 for r, *_ in worst if r > tol)
        limit    = n_total if show_all else 100
        if not show_all and n_total > limit:
            print(f"  (showing top {limit} of {n_total} columns; use --all for full list)")
        for rel, col, direction, mean_j, mean_f in worst[:limit]:
            flag = "  *** EXCEEDS TOLERANCE ***" if rel > tol else ""
            print(f"  {col:20}  max_rel={rel:.2e}  {direction}  "
                  f"mean(JUICE)={mean_j:.3e}  mean(REF)={mean_f:.3e}{flag}")
        if not show_all and n_total > limit:
            remaining_exceed = sum(1 for r, *_ in worst[limit:] if r > tol)
            if remaining_exceed:
                print(f"  ... {remaining_exceed} more columns exceed tolerance (not shown; use --all)")
        _ranked = [(float(r), c) for r, c, *_ in sorted(worst, reverse=True)]
        summary["groups"][label] = {"n": n_total, "n_exceed": n_exceed, "worst": _ranked[:10], "all": _ranked}
        exceeded = [(rel, col) for rel, col, *_ in worst if rel > tol]
        if exceeded:
            summary["ok"] = False
            rel, col = exceeded[0]
            return False, col, rel
        return True, None, None

    groups = [
        ("1. s-values (settings)", s_cols, jui, fac, 0.0),
        ("2. c-values (concentrations)", c_cols, jui, fac, floor),
    ]

    for label, cols, jui_df, fac_df, floor_ in groups:
        if not cols:
            print(f"\n  {label}: no common columns — skipped")
            continue
        ok, bad_col, bad_rel = _compare_group(label, cols, jui_df, fac_df, floor_)
        if not ok:
            print(f"\n*** {'CONTINUES' if no_stop else 'STOP'}: {label} — '{bad_col}' max_rel={bad_rel:.2e} exceeds {tol:g} ***")
            if not no_stop:
                return summary

    # k-values and v-values: compare against FACSIMILE .kv if available
    if os.path.isfile(kv_jui_path):
        kv_jui = load_column_file(kv_jui_path)
        # Build canonical-key → JUICE-column mapping
        k_jui = {_canon_kv_label(c): c for c in kv_jui.columns if c.startswith('k_')}
        v_jui = {_canon_kv_label(c): c for c in kv_jui.columns if c.startswith('v_')}

        if os.path.isfile(kv_fac_path):
            try:
                k_fac_df, v_fac_df = load_facsimile_kv(kv_fac_path)
                # Align time grids: interpolate FACSIMILE onto JUICE time points
                t_jui = kv_jui.index.values.astype(float)
                t_fac_k = k_fac_df.index.values.astype(float)
                t_fac_v = v_fac_df.index.values.astype(float)

                for label, jui_map, fac_df, t_fac in (
                        ('3. k-values (rate constants)', k_jui, k_fac_df, t_fac_k),
                        ('4. v-values (turnover rates)',  v_jui, v_fac_df, t_fac_v)):
                    common = sorted(set(jui_map) & set(fac_df.columns))
                    print(f"\n{'─'*60}")
                    print(f"  {label} — {len(common)} common reactions "
                          f"(JUICE: {len(jui_map)}, FAC: {len(fac_df.columns)})")
                    print(f"{'-'*60}")
                    worst = []
                    for rxn in common:
                        j = np.interp(t_jui, t_fac, fac_df[rxn].values.astype(float))
                        jv = kv_jui[jui_map[rxn]].values.astype(float)
                        scale = np.nanmax(np.abs(j))
                        rel = (np.nanmax(np.abs(jv)) if scale == 0
                               else np.nanmax(np.abs(jv - j)) / scale)
                        direction = 'JUICE>REF' if (jv - j)[np.argmax(np.abs(jv - j))] > 0 else 'JUICE<REF'
                        worst.append((rel, rxn, direction, np.nanmean(jv), np.nanmean(j)))
                    worst.sort(key=lambda x: x[0], reverse=True)
                    n_exceed = sum(1 for r, *_ in worst if r > tol)
                    summary["groups"][label] = {"n": len(worst), "n_exceed": n_exceed,
                                                "worst": [(float(r), x) for r, x, *_ in worst[:10]]}
                    if n_exceed:
                        summary["ok"] = False
                    for rel, rxn, direction, mj, mf in worst[:10]:
                        flag = '  *** EXCEEDS TOLERANCE ***' if rel > tol else ''
                        print(f"  {rxn:45}  max_rel={rel:.2e}  {direction}  "
                              f"mean(J)={mj:.3e}  mean(F)={mf:.3e}{flag}")
                    if len(worst) > 10:
                        print(f"  ... ({len(worst)-10} more; {n_exceed} total exceed tolerance)")
            except Exception as e:
                print(f"\n  [WARNING] Could not load FACSIMILE kv: {e}")
                _show_juice_only_kv(kv_jui, k_jui, v_jui)
        else:
            _show_juice_only_kv(kv_jui, k_jui, v_jui)
    else:
        print(f"\n  3+4. kv-values: {kv_jui_path} not found — run with kv/pd checkbox enabled")

    # pd values: compare against FACSIMILE .pd if available
    if os.path.isfile(pd_jui_path):
        pd_jui = load_column_file(pd_jui_path)
        p_jui = [c for c in pd_jui.columns if c.startswith('P_')]
        d_jui = [c for c in pd_jui.columns if c.startswith('D_')]

        if os.path.isfile(pd_fac_path):
            try:
                pd_fac = load_facsimile_pd(pd_fac_path)
                t_jui = pd_jui.index.values.astype(float)
                t_fac = pd_fac.index.values.astype(float)
                common_sp = sorted(set(c[2:] for c in p_jui) &
                                   set(c[2:] for c in pd_fac.columns if c.startswith('P_')))
                print(f"\n{'─'*60}")
                print(f"  5. pd-values — {len(common_sp)} species")
                print(f"{'-'*60}")
                for sp in common_sp[:5]:
                    pj = pd_jui[f'P_{sp}'].values.astype(float)
                    dj = pd_jui[f'D_{sp}'].values.astype(float)
                    pf = np.interp(t_jui, t_fac, pd_fac[f'P_{sp}'].values.astype(float))
                    df_ = np.interp(t_jui, t_fac, pd_fac[f'D_{sp}'].values.astype(float))
                    # Scale = max of both models; floor = 1 molec/cm³/s to avoid
                    # division-by-zero when one model has zero production/destruction.
                    p_scale = max(np.nanmax(np.abs(pf)), np.nanmax(np.abs(pj)), 1.0)
                    d_scale = max(np.nanmax(np.abs(df_)), np.nanmax(np.abs(dj)), 1.0)
                    p_rel = np.nanmax(np.abs(pj - pf)) / p_scale
                    d_rel = np.nanmax(np.abs(dj - df_)) / d_scale
                    print(f"  {sp:15}  P_rel={p_rel:.2e}  D_rel={d_rel:.2e}")
                if len(common_sp) > 5:
                    print(f"  ... ({len(common_sp)-5} more)")
            except Exception as e:
                print(f"\n  [WARNING] Could not load FACSIMILE pd: {e}")
                _show_juice_only_pd(pd_jui, p_jui)
        else:
            _show_juice_only_pd(pd_jui, p_jui)
    else:
        print(f"\n  5. pd-values: {pd_jui_path} not found — run with kv/pd checkbox enabled")

    # ── 6. kOH comparison (D_OH / [OH]) ─────────────────────────────────────
    print(f"\n{'─'*60}")
    print(f"  6. kOH  (JUICE: Σk·[X]  |  FACSIMILE: D_OH/([OH]+1))")
    print(f"{'─'*60}")
    _koh_done = False
    if os.path.isfile(pd_jui_path):
        _pd_jui = load_column_file(pd_jui_path)
        if 'kOH' in _pd_jui.columns:
            _koh_j = np.interp(t_ref, _pd_jui.index.values.astype(float),
                               _pd_jui['kOH'].values.astype(float))
        elif 'D_OH' in _pd_jui.columns and 'OH' in jui.columns:
            _d_oh_j = np.interp(t_ref, _pd_jui.index.values.astype(float),
                                _pd_jui['D_OH'].values.astype(float))
            _oh_j = jui['OH'].values.astype(float)
            with np.errstate(divide='ignore', invalid='ignore'):
                _koh_j = np.where(_oh_j > 1e2, _d_oh_j / _oh_j, np.nan)
        else:
            _koh_j = None
        if _koh_j is not None:
            _valid_j = np.sum(np.isfinite(_koh_j))
            _mean_j  = np.nanmean(_koh_j)
            print(f"  JUICE   : mean={_mean_j:.3e} s⁻¹  ({_valid_j} valid points)")

            if os.path.isfile(pd_fac_path):
                try:
                    _pd_fac2 = load_facsimile_pd(pd_fac_path)
                    if 'D_OH' in _pd_fac2.columns and 'OH' in fac.columns:
                        _d_oh_f = np.interp(t_ref, _pd_fac2.index.values.astype(float),
                                            _pd_fac2['D_OH'].values.astype(float))
                        _oh_f = fac['OH'].values.astype(float)
                        # D_OH/([OH]+1): avoids 0/0 at night, correct when [OH]>>1
                        _koh_f = _d_oh_f / (_oh_f + 1.0)
                        _valid_f = np.sum(np.isfinite(_koh_f))
                        _mean_f  = np.nanmean(_koh_f)
                        print(f"  FACSIMILE: mean={_mean_f:.3e} s⁻¹  ({_valid_f} valid points)")
                        _mask = np.isfinite(_koh_j) & np.isfinite(_koh_f)
                        if _mask.any():
                            _sc = np.maximum(np.abs(_koh_j[_mask]), np.abs(_koh_f[_mask]))
                            _sc = np.maximum(_sc, 1e-30)
                            _rel = np.max(np.abs(_koh_j[_mask] - _koh_f[_mask]) / _sc)
                            _tag = "*** EXCEEDS TOLERANCE ***" if _rel > tol else ""
                            print(f"  max_rel  = {_rel:.2e}  {_tag}")
                            summary["groups"]["6. kOH"] = {"n": 1, "n_exceed": int(_rel > tol),
                                                           "worst": [(float(_rel), "kOH")]}
                            if _rel > tol:
                                summary["ok"] = False
                        _koh_done = True
                except Exception as _e:
                    print(f"  [WARNING] FACSIMILE pd: {_e}")
            if not _koh_done:
                print("  (no FACSIMILE .pd — JUICE only)")
        else:
            print("  D_OH column missing in JUICE pd file.")
    else:
        print("  kOH not computed (pd file or OH missing).")

    print(f"\n{'='*60}")
    if summary["ok"]:
        print(f"  All compared groups within tolerance ({tol:g}).")
    else:
        for label, g in summary["groups"].items():
            if g["n_exceed"]:
                print(f"  {label}: {g['n_exceed']} of {g['n']} exceed {tol:g}"
                      f" (worst: {g['worst'][0][1]} {g['worst'][0][0]:.2e})")
    print(f"{'='*60}\n")
    return summary


if __name__ == "__main__":
    # µ, ³, ─ on a Windows console; only when run as a script, never when imported
    if sys.platform == "win32" and hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    parser = argparse.ArgumentParser(
        description="Compare JUICE vs FACSIMILE output for a scenario.")
    parser.add_argument("scenario", help="Scenario name (e.g. scene_03)")
    parser.add_argument("--tol", type=float, default=1e-4,
                        help="Max relative deviation threshold (default: 1e-4)")
    parser.add_argument("--all", action="store_true",
                        help="Write all columns to log (default: top 100)")
    parser.add_argument("--no-stop", action="store_true",
                        help="Continue through all groups even if tolerance is exceeded")
    parser.add_argument("--juice-dir", default=None,
                        help="Override JUICE output directory (default: output_root/scenario)")
    parser.add_argument("--fac-dir", default=None,
                        help="Override FACSIMILE reference directory (default: easy_root/scenario)")
    parser.add_argument("--floor", type=float, default=DEFAULT_FLOOR,
                        help=f"Concentration floor in molec/cm³ for the relative deviation "
                             f"(default: {DEFAULT_FLOOR:g}; 0 = none)")
    args = parser.parse_args()

    # Log next to the JUICE output (output_root from the TOML, repository default without [paths])
    from scenario_prep import load_paths
    _juice_dir = args.juice_dir or os.path.join(load_paths(_toml_path())["output_root"], args.scenario)
    os.makedirs(_juice_dir, exist_ok=True)
    log_path = os.path.join(_juice_dir,
                            f"{args.scenario}_difference_JUICE_Facsimile.log")

    # Tee stdout to both console and log file
    import datetime as _dt
    class _Tee:
        def __init__(self, console, logfile):
            self.console = console
            self.logfile = logfile
        def write(self, s):
            self.console.write(s)
            self.logfile.write(s)
        def flush(self):
            self.console.flush()
            self.logfile.flush()

    _saved_stdout = sys.stdout
    log_fh = open(log_path, "w", encoding="utf-8")
    log_fh.write(f"# JUICE vs FACSIMILE comparison  --  "
                 f"{_dt.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
    sys.stdout = _Tee(_saved_stdout, log_fh)

    _summary = {"ok": False}
    try:
        _summary = compare(args.scenario, args.tol, show_all=args.all, no_stop=args.no_stop,
                           juice_dir_override=args.juice_dir, fac_dir_override=args.fac_dir,
                           floor=args.floor)
    except (FileNotFoundError, ValueError) as _e:
        print(f"\n*** {_e}")
    finally:
        sys.stdout = _saved_stdout
        log_fh.close()
    print(f"\nLog written to: {log_path}")
    sys.exit(0 if _summary.get("ok") else 1)
