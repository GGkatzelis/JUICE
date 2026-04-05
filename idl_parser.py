"""
idl_parser.py
=============
Parse an IDL-generated SAPHIR EASY master file (.easy) and its ENZ file,
producing a scenario_config dict that faithfully mirrors the IDL/FACSIMILE
model.

Key design:
- Injection reactions from the EASY file are preserved as explicit kinetic
  reactions in ``injection_reactions``.  Each "qXXX = input(COLUMN)" switch
  becomes a zeroth-order source (empty reactants) whose rate expression
  includes the ENZ switch column name (0/1 signal already in ctx via
  load_saphir_environment).  TwoReaction chain intermediates are kept as
  normal first-order reactions.
- Chemistry files are resolved recursively from the chemistry/ subfolder.
- Burst-injection summaries (timing + estimated ppb) are generated for the
  Burst Injections tab (informational, marked "IDL", unchecked in GUI).

Usage (API):
    from idl_parser import parse_idl_easy
    cfg = parse_idl_easy("path/to/scenario.easy", "path/to/kammer.enz")

Usage (CLI):
    python idl_parser.py  scenario.easy  [kammer.enz]  [output.json]
"""

import re
import os
import json
import datetime
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# ENZ epoch
# ---------------------------------------------------------------------------
ENZ_EPOCH = datetime.datetime(2000, 1, 1, 0, 0, 0)


# ---------------------------------------------------------------------------
# Standard ENZ column -> (env_map key, canonical column name) mapping
# ---------------------------------------------------------------------------
ENZ_COL_TO_ENV = {
    "t":        ("TEMP",    "t"),
    "jno2":     ("JNO2",    "jNO2"),
    "jo1d":     ("JO1D",    "jO1D"),
    "dilute":   ("FLOW",    "dilute"),
    "hv":       ("HV",      "hv"),
    "h2o":      ("H2O",     "H2O"),
    "sz":       ("SZ",      "sz"),
    "jhono":    ("JHONO",   "jHONO"),
    "jh2o2":    ("JH2O2",   "jH2O2"),
    "jhcho_m":  ("JHCHO_m", "jHCHO_m"),
    "jhcho_r":  ("JHCHO_r", "jHCHO_r"),
}

# QH2O reactions are handled separately in the H2O impurities block (as injection_reactions
# gated by QH2O). QCO2 is always-on and goes to continuous_sources. Neither goes through
# the generic injection_source_rxns loop.
_NON_INJECTION_QNAMES = {"QCO2", "QH2O"}


# ---------------------------------------------------------------------------
# ENZ helpers
# ---------------------------------------------------------------------------

def _load_enz(enz_path: str) -> pd.DataFrame:
    """Return DataFrame with original-case column names, TIME in SAPHIR-epoch s."""
    col_map = {}
    data_start = 0
    with open(enz_path, "r", encoding="latin-1") as fh:
        for i, line in enumerate(fh):
            if line.startswith("COLUMN"):
                parts = line.split("=", 1)
                name = parts[1].split("/")[0].split("(")[0].strip()
                idx = int(parts[0].split()[1]) - 1
                col_map[idx] = name
            if "&&" in line:
                data_start = i + 1
                break
    headers = [col_map[i] for i in sorted(col_map)]
    df = pd.read_csv(enz_path, sep=r",\s*", engine="python",
                     skiprows=data_start, header=None, names=headers)
    return df


def _s_to_utc(s: float) -> str:
    dt = ENZ_EPOCH + datetime.timedelta(seconds=float(s))
    return dt.strftime("%d.%m.%Y %H:%M:%S")


# ---------------------------------------------------------------------------
# File-system helpers
# ---------------------------------------------------------------------------

def _find_file_recursive(name: str, base_dir: str):
    """Search base_dir and all subdirectories for *name*.
    Returns the path relative to base_dir (forward slashes), or None."""
    for root, _dirs, files in os.walk(base_dir):
        if name in files:
            rel = os.path.relpath(os.path.join(root, name), base_dir)
            return rel.replace("\\", "/")
    return None


# ---------------------------------------------------------------------------
# EASY parsing helpers
# ---------------------------------------------------------------------------

def _strip_inline_comment(line: str) -> str:
    return line.split(";")[0].rstrip()


# Substitutions: FACSIMILE lowercase/mixed variable names → Python ctx names.
# ctx keys are the uppercase ENZ column names (set by runner via ctx[k.upper()]).
# Order matters: longer patterns first to avoid partial matches.
_EASY_VAR_SUBS = [
    # J-value columns (EASY uses lowercase, ctx stores uppercase)
    (r'(?<![A-Za-z0-9_])jhcho_m(?![A-Za-z0-9_])', 'JHCHO_m'),
    (r'(?<![A-Za-z0-9_])jhcho_r(?![A-Za-z0-9_])', 'JHCHO_r'),
    (r'(?<![A-Za-z0-9_])jhono(?![A-Za-z0-9_])',   'JHONO'),
    (r'(?<![A-Za-z0-9_])jh2o2(?![A-Za-z0-9_])',   'JH2O2'),
    (r'(?<![A-Za-z0-9_])jno2(?![A-Za-z0-9_])',    'JNO2'),
    (r'(?<![A-Za-z0-9_])jo1d(?![A-Za-z0-9_])',    'JO1D'),
    # HV_in is the FACSIMILE alias for the HV (roof open/closed) ENZ column
    (r'(?<![A-Za-z0-9_])HV_in(?![A-Za-z0-9_])',   'HV'),
    (r'(?<![A-Za-z0-9_])hv_in(?![A-Za-z0-9_])',   'HV'),
]


def _easy_to_py(expr: str) -> str:
    """Convert a FACSIMILE rate expression to a Python-evaluable string.

    - @ → **  (FACSIMILE power operator)
    - lowercase/mixed J-value names → uppercase ctx names (JNO2, JO1D, …)
    - HV_in → HV  (ENZ column name used in ctx)
    """
    out = expr.replace("@", "**")
    for pattern, replacement in _EASY_VAR_SUBS:
        out = re.sub(pattern, replacement, out, flags=re.IGNORECASE)
    return out


def _eval_const(expr: str, M: float = 2.5e19):
    """Evaluate a simple numeric FACSIMILE constant. Returns float or None."""
    safe = expr.replace("@", "**").upper()
    try:
        return float(eval(safe, {"__builtins__": {}}, {"M": M, "E": 2.71828}))
    except Exception:
        return None


def _ppb_from_rate(rate_expr: str):
    """Extract ppb/min from pattern M*1e-9/60*X. Returns float or None."""
    m = re.search(r"M\s*\*\s*1[eE]-9\s*/\s*60\s*\*\s*([0-9eE.+\-*/]+)",
                  rate_expr, re.I)
    if m:
        try:
            return float(eval(m.group(1)))
        except Exception:
            pass
    return None


def _read_easy_lines(easy_path: str):
    """Yield (clean_line, is_commented) for each line, with ; comments stripped."""
    with open(easy_path, "r", encoding="latin-1") as fh:
        for line in fh:
            stripped = line.strip()
            is_commented = stripped.startswith((";", "*"))
            clean = _strip_inline_comment(stripped)
            yield clean, is_commented


def _rxn_key(lhs_up: str, rhs_str: str) -> str:
    return lhs_up.replace(" ", "") + "-->" + rhs_str.replace(" ", "").upper()


def _products_of(rhs_str: str) -> list:
    """Return list of UPPERCASE species names from a reaction RHS string."""
    out = []
    for term in re.split(r"\+", rhs_str):
        sp = re.sub(r"^[\d.]+\s*\*\s*", "", term.strip()).upper()
        if sp and not sp[0].isdigit():
            out.append(sp)
    return out


# ---------------------------------------------------------------------------
# Main parser
# ---------------------------------------------------------------------------

def parse_idl_easy(easy_path: str, enz_path: str = None,
                   chemistry_folder: str = None,
                   chamber_vol: float = 276.0) -> dict:
    """
    Parse an IDL-generated EASY master file and optional ENZ file.

    Parameters
    ----------
    easy_path         : path to the .easy file
    enz_path          : path to the .enz file (auto-discovered if None)
    chemistry_folder  : path to chemistry/ folder (auto-discovered if None)
    chamber_vol       : chamber volume in m3 (276 for SAPHIR)

    Returns
    -------
    dict -- scenario_config ready for use with the Python MCM tool.
    The new key ``injection_reactions`` contains the IDL-faithful kinetic
    injection reactions.  ``burst_injections`` contains informational
    summaries (marked with ``"source": "IDL"``).
    """

    easy_dir = os.path.dirname(os.path.abspath(easy_path))

    if chemistry_folder is None:
        chemistry_folder = os.path.join(easy_dir, "chemistry")

    # Auto-discover ENZ file
    if enz_path is None:
        for f in os.listdir(easy_dir):
            if f.lower().endswith(".enz"):
                enz_path = os.path.join(easy_dir, f)
                break

    enz_df = None
    if enz_path and os.path.isfile(enz_path):
        try:
            enz_df = _load_enz(enz_path)
            print(f"[idl_parser] ENZ loaded: {len(enz_df.columns)} columns, "
                  f"{len(enz_df)} rows")
        except Exception as e:
            print(f"[idl_parser] WARNING: could not load ENZ: {e}")

    # ---- parse EASY file line by line ----------------------------------------
    mechanism_files_raw = []
    env_map             = {}
    pressure_hpa        = None
    initials            = {}
    j_scalings          = {}

    all_k_assignments   = {}   # normalised_key -> rate_expr
    all_reactions       = []   # (lhs_upper, rhs_str) active only
    injection_switches  = {}   # var_upper -> ENZ_col (original case) or None
    source_consts       = {}   # var_upper -> py_expr  (SHONO/SHCHO/WALLDEP etc.)

    y_active    = False;  y_eqn     = ""
    scale       = "1.0"

    co2_bg_rate = None
    h2o_imps    = []
    skipped_notes = []   # human-readable notes about skipped/inactive lines

    for clean, is_commented in _read_easy_lines(easy_path):
        if not clean:
            continue

        # #include
        inc = re.match(r"#include\s+(\S+)", clean, re.I)
        if inc and not is_commented:
            mechanism_files_raw.append(inc.group(1).strip())
            continue

        # k[lhs --> rhs] = CONST(rate)  or  k[...] = rate
        km = re.match(r"k\s*\[(.+?)\]\s*=\s*(.*)", clean, re.I)
        if km and not is_commented:
            rxn_part = km.group(1).strip()
            val_part = km.group(2).strip()
            # Unwrap outer CONST(...)
            inner_m = re.match(r"CONST\s*\(\s*(.*?)\s*\)\s*$", val_part, re.I)
            inner = inner_m.group(1) if inner_m else val_part
            if "-->" in rxn_part:
                lhs_r, rhs_r = rxn_part.split("-->", 1)
                key = _rxn_key(lhs_r.strip().upper(), rhs_r.strip())
            else:
                key = rxn_part.strip().replace(" ", "").upper() + "-->"
            all_k_assignments[key] = inner
            continue

        # reaction: A --> B  (no = sign, no k[)
        if "-->" in clean and "=" not in clean and "k[" not in clean \
                and not is_commented:
            lhs, rhs = clean.split("-->", 1)
            lhs = lhs.strip().upper()
            rhs = rhs.strip()
            if lhs:
                all_reactions.append((lhs, rhs))
            continue

        # FACSIMILE-only constructs — note and skip
        if re.match(r"\w+\s*=\s*FUNCTION\s*\[", clean, re.I) and not is_commented:
            fname = clean.split("=")[0].strip()
            skipped_notes.append(f"FUNCTION definition '{fname}': FACSIMILE-only macro, "
                                 f"not needed in Python ODE")
            continue
        if re.match(r"FILES\s*\[", clean, re.I) and not is_commented:
            skipped_notes.append(f"FILES[] directive: FACSIMILE file config, not needed")
            continue
        if re.match(r"FACS\s*\[", clean, re.I) and not is_commented:
            skipped_notes.append(f"FACS[] directive: FACSIMILE solver config, not needed")
            continue

        # X = INPUT(col)
        inp = re.match(r"(\w+)\s*=\s*INPUT\s*\(\s*(\w+)\s*\)", clean, re.I)
        if inp and not is_commented:
            var_up  = inp.group(1).upper()
            enz_col = inp.group(2).strip()
            mapping = (ENZ_COL_TO_ENV.get(enz_col.lower()) or
                       ENZ_COL_TO_ENV.get(var_up.lower()))
            if mapping:
                env_key, canonical = mapping
                env_map[env_key] = {"type": "SAPHIR File", "value": canonical}
            else:
                injection_switches[var_up] = enz_col
            continue

        # X = CONST(expr)
        cst = re.match(r"(\w+)\s*=\s*CONST\s*\(\s*(.*?)\s*\)\s*$", clean, re.I)
        if cst and not is_commented:
            var_up = cst.group(1).upper()
            expr   = cst.group(2).strip()

            if var_up in ("P", "PRESS"):
                try:
                    pressure_hpa = float(expr)
                except Exception:
                    pressure_hpa = 1013.25

            elif var_up == "Y":
                y_active = True
                ym = (re.search(r"([0-9eE.+\-]+)[eE]-9\s*\*\s*M", expr, re.I) or
                      re.search(r"M\s*\*\s*([0-9eE.+\-]+)[eE]-9", expr, re.I))
                y_eqn = ym.group(1) if ym else expr

            elif var_up in ("QY",):
                y_active = True
                ym = (re.search(r"([0-9eE.+\-]+)[eE]-9\s*\*\s*M", expr, re.I) or
                      re.search(r"M\s*\*\s*([0-9eE.+\-]+)[eE]-9", expr, re.I))
                y_eqn = ym.group(1) if ym else expr

            elif var_up == "QCO2":
                # Always-on CO2 source
                injection_switches[var_up] = None

            elif (var_up.lower().startswith("j") and
                  var_up not in ("JHONO", "JNO2", "JO1D", "JH2O2",
                                 "JHCHO_M", "JHCHO_R", "JSCALE")):
                j_scalings[var_up] = _easy_to_py(expr)

            else:
                # Collect as potential chamber-source constant (SHONO, SHCHO, WALLDEP, …)
                source_consts[var_up] = _easy_to_py(expr)
            continue

        # INITIAL[X] = CONST(expr)
        init = re.match(r"INITIAL\s*\[\s*(\w+)\s*\]\s*=\s*CONST\s*\((.*?)\)",
                        clean, re.I)
        if init and not is_commented:
            spec = init.group(1).upper()
            expr = init.group(2).strip()
            if spec not in ("*", "ALL"):
                val = _eval_const(expr)
                if val is not None and val > 0:
                    ppb = val / 2.5e19 * 1e9
                    if ppb > 0:
                        initials[spec] = round(ppb, 6)
            continue

    # ---- build reaction lookup -----------------------------------------------
    # Rebuild lookup from collected reactions for accurate key matching
    k_lookup = {}
    for lhs_up, rhs_str in all_reactions:
        key = _rxn_key(lhs_up, rhs_str)
        if key in all_k_assignments:
            k_lookup[key] = all_k_assignments[key]

    def _k(lhs_up, rhs_str):
        return k_lookup.get(_rxn_key(lhs_up, rhs_str), "")

    # ---- classify reactions --------------------------------------------------
    injection_source_rxns = []         # (lhs_up, rhs_str, rate, switch_col|None)
    injection_intermediate_species = set()
    regular_reactions = []

    for lhs_up, rhs_str in all_reactions:
        if lhs_up in injection_switches:
            rate = _k(lhs_up, rhs_str)
            injection_source_rxns.append(
                (lhs_up, rhs_str, rate, injection_switches[lhs_up]))
            for sp in _products_of(rhs_str):
                injection_intermediate_species.add(sp)
        else:
            regular_reactions.append((lhs_up, rhs_str))

    # Extend intermediates iteratively (chain reactions)
    injection_chain_rxns = []
    changed = True
    while changed:
        changed = False
        still_regular = []
        for lhs_up, rhs_str in regular_reactions:
            if lhs_up in injection_intermediate_species:
                rate = _k(lhs_up, rhs_str)
                injection_chain_rxns.append((lhs_up, rhs_str, rate))
                for sp in _products_of(rhs_str):
                    if sp not in injection_intermediate_species:
                        injection_intermediate_species.add(sp)
                        changed = True
            else:
                still_regular.append((lhs_up, rhs_str))
        regular_reactions = still_regular

    # ---- build injection_reactions -------------------------------------------
    def _is_zero_rate(rate_expr: str) -> bool:
        norm = rate_expr.replace(" ", "")
        if norm in ("0", "0.", "0.0", ""):
            return True
        # Match *0 that is genuinely zero: not followed by any nonzero digit sequence
        # i.e. *0, *0., *0.0, *0.00... but NOT *0.01, *0.5, *0.12, *0e3
        return bool(re.search(r"\*0(?!\.?[0-9]*[1-9])(?!\.?[0-9]*[eE])", norm))

    injection_reactions = []

    for lhs_up, rhs_str, rate_expr, switch_col in injection_source_rxns:
        if lhs_up in _NON_INJECTION_QNAMES:
            continue
        if _is_zero_rate(rate_expr):
            skipped_notes.append(f"zero-rate injection: {lhs_up} --> {rhs_str}  (k={rate_expr})")
            continue

        rate_py = _easy_to_py(rate_expr)
        # Gate by ENZ switch column (uppercase matches ctx keys from ENZ loader)
        if switch_col is not None:
            full_rate = f"({rate_py}) * {switch_col.upper()}"
        else:
            full_rate = rate_py   # always-on

        # Use the switch column name as the reactant label (EASY file consistency).
        # The ODE does NOT consume it — QNO-type names are ENZ constrained inputs
        # that get discarded from all_species by the runner.  The actual gating is
        # already embedded in full_rate via "* QNO".
        reactant_label = lhs_up if switch_col is None else switch_col.upper()
        injection_reactions.append({
            "reactants": reactant_label,
            "products":  rhs_str,
            "rate":      full_rate,
            "label":     f"{lhs_up} injection",
            "group":     "IDL Injection",
            "comment":   f"k[{lhs_up} --> {rhs_str}] = {rate_expr}",
        })

    for lhs_up, rhs_str, rate_expr in injection_chain_rxns:
        if _is_zero_rate(rate_expr):
            continue
        injection_reactions.append({
            "reactants": lhs_up,
            "products":  rhs_str,
            "rate":      _easy_to_py(rate_expr),
            "label":     f"{lhs_up} chain",
            "group":     "IDL Injection",
            "comment":   f"k[{lhs_up} --> {rhs_str}] = {rate_expr}",
        })

    # ---- categorise remaining (chamber physics) reactions --------------------
    # These are added as explicit reactions (consistent with the EASY file approach)
    # rather than as artifact-checkbox entries.  The CONST value of each source
    # species (SHONO, SHCHO) is multiplied into the rate so the full production
    # rate is self-contained: rate_full = k_rate * source_const.
    chamber_source_reactions = []

    for lhs_up, rhs_str in regular_reactions:
        rhs_up = rhs_str.upper()
        rate   = _k(lhs_up, rhs_str)

        if lhs_up in ("SHONO", "SHONO+HV") and "HONO" in rhs_up:
            if not rate:
                skipped_notes.append(f"no k-rate found for {lhs_up} --> {rhs_str}")
                continue
            k_py = _easy_to_py(rate)
            const_py = source_consts.get("SHONO", "")
            # SHONO = CONST(...) in FACSIMILE — it is a fixed parameter, not an ODE species.
            # Its value is already folded into full_rate, so reactants are empty (zeroth-order).
            # SPECIAL CASE: if lhs is "SHONO+HV", HV is an explicit reactant in the EASY file
            # and must be multiplied in.  For ethanol-style "SHONO=CONST(HV*...)", HV is
            # already inside const_py so this branch is not taken.
            if lhs_up == "SHONO+HV" and const_py and "HV" not in const_py.upper():
                full_rate = f"({k_py}) * ({const_py}) * HV"
            else:
                full_rate = f"({k_py}) * ({const_py})" if const_py else k_py
            chamber_source_reactions.append({
                "reactants": "",
                "products":  rhs_str,
                "rate":      full_rate,
                "label":     f"SHONO={const_py or '1'} -> HONO (chamber photolytic source)",
                "group":     "Chamber Source",
                "comment":   f"k={rate}  *  SHONO={const_py or '1'}",
            })

        elif lhs_up in ("SHCHO", "SHCHO+HV", "QSHCHO+HV") and "HCHO" in rhs_up:
            if not rate:
                skipped_notes.append(f"no k-rate found for {lhs_up} --> {rhs_str}")
                continue
            k_py = _easy_to_py(rate)
            const_py = source_consts.get("SHCHO", source_consts.get(lhs_up, ""))
            # SHCHO = CONST(...) — same as SHONO: fold const into rate, empty reactants.
            full_rate = f"({k_py}) * ({const_py})" if const_py else k_py
            chamber_source_reactions.append({
                "reactants": "",
                "products":  rhs_str,
                "rate":      full_rate,
                "label":     f"{lhs_up}={const_py or '1'} -> HCHO (chamber photolytic source)",
                "group":     "Chamber Source",
                "comment":   f"k={rate}  *  {lhs_up}={const_py or '1'}",
            })

        elif lhs_up == "SNO2" and "NO2" in rhs_up:
            # SNO2 is an ODE species (not a CONST).  Active only if INITIAL > 0
            # or SNO2=CONST(1) is present.  Check source_consts and initials.
            sno2_const = source_consts.get("SNO2", "")
            sno2_init  = initials.get("SNO2", 0.0)
            if not sno2_const and sno2_init == 0.0:
                skipped_notes.append(
                    f"SNO2 source inactive: no CONST def and INITIAL[SNO2]=0  "
                    f"(k={rate or 'none'})")
                continue
            if not rate:
                skipped_notes.append(f"no k-rate found for SNO2 --> {rhs_str}")
                continue
            k_py = _easy_to_py(rate)
            if sno2_const:
                full_rate = f"({k_py}) * ({sno2_const})"
            else:
                full_rate = k_py   # SNO2 will be in ODE state at its initial value
            chamber_source_reactions.append({
                "reactants": "SNO2",
                "products":  rhs_str,
                "rate":      full_rate,
                "label":     "SNO2 → NO2 (chamber NO2 source)",
                "group":     "Chamber Source",
                "comment":   f"k={rate}  SNO2_const={sno2_const or 'ODE species'}",
            })

        elif lhs_up == "H2O" and "HONO" in rhs_up:
            # Wall dark-HONO: H2O is an env species → multiplied via reactants_env
            if not rate:
                skipped_notes.append(f"no k-rate found for H2O --> {rhs_str}")
                continue
            chamber_source_reactions.append({
                "reactants": "H2O",
                "products":  "HONO",
                "rate":      _easy_to_py(rate),
                "label":     "H2O → HONO (wall dark source)",
                "group":     "Chamber Source",
                "comment":   f"k={rate}  [H2O from ENZ via env_map]",
            })

        elif "Y" in lhs_up.upper() and "OH" in lhs_up.upper():
            # Y + OH background reactivity — handled via y_active artifact, not as reaction
            skipped_notes.append(f"Y + OH --> HO2: handled as Y-background artifact "
                                 f"(y_eqn={y_eqn} ppb, k=2.4e-13)")
        else:
            skipped_notes.append(f"unclassified regular reaction: {lhs_up} --> {rhs_str}")

    # ---- CO2 background (always-on QCO2 switch) -----------------------------
    if "QCO2" in injection_switches:
        for lhs_up, rhs_str, rate_expr, sc in injection_source_rxns:
            if lhs_up == "QCO2" and "CO2" in rhs_str.upper():
                if not _is_zero_rate(rate_expr):
                    co2_bg_rate = _easy_to_py(rate_expr)
                break

    # ---- H2O injection impurities -------------------------------------------
    # qH2O --> X reactions (non-H2O, non-HONO) are gated by the QH2O ENZ switch.
    # They are treated as injection_reactions (same mechanism as qNO --> NO),
    # NOT as h2o_injection_impurities, because in the IDL model there is no separate
    # H2O burst injection to act as the timing trigger — the QH2O column IS the gate.
    if "QH2O" in injection_switches:
        qh2o_col = injection_switches["QH2O"]  # original ENZ column name (e.g. "QH2O")
        _seen_qh2o_products = set()  # deduplicate: same rhs can appear in multiple #include files
        for lhs_up, rhs_str, rate_expr, sc in injection_source_rxns:
            if lhs_up != "QH2O":
                continue
            rhs_up = rhs_str.upper()
            if "H2O" in rhs_up:   # skip the main H2O refill itself
                skipped_notes.append(f"qH2O --> H2O (main H2O injection): not modelled "
                                     f"(H2O is a constrained ENV variable, k={rate_expr})")
                continue
            if _is_zero_rate(rate_expr):
                skipped_notes.append(f"zero-rate qH2O impurity: QH2O --> {rhs_str}  "
                                     f"(k={rate_expr})")
                continue
            dedup_key = (rhs_up, rate_expr.replace(" ", ""))
            if dedup_key in _seen_qh2o_products:
                skipped_notes.append(f"duplicate qH2O impurity skipped: QH2O --> {rhs_str}")
                continue
            _seen_qh2o_products.add(dedup_key)
            rate_py  = _easy_to_py(rate_expr)
            full_rate = f"({rate_py}) * {qh2o_col.upper()}"
            tgt = _products_of(rhs_str)
            tgt_name = tgt[0] if tgt else rhs_str.strip()
            injection_reactions.append({
                "reactants": qh2o_col.upper(),
                "products":  rhs_str,
                "rate":      full_rate,
                "label":     f"QH2O impurity: {tgt_name}",
                "group":     "IDL Injection",
                "comment":   f"k[QH2O --> {rhs_str}] = {rate_expr}",
            })

    # ---- burst injection summaries (informational) ---------------------------
    burst_injections = []
    if enz_df is not None:
        time_col = enz_df.columns[0]
        t_arr = enz_df[time_col].values.astype(float)
        enz_cols_lc = {c.lower(): c for c in enz_df.columns}

        processed_switches = set()
        for lhs_up, rhs_str, rate_expr, switch_col in injection_source_rxns:
            if lhs_up in _NON_INJECTION_QNAMES or switch_col is None:
                continue
            sc_lc = switch_col.lower()
            if sc_lc in processed_switches:
                continue
            processed_switches.add(sc_lc)

            enz_col_orig = enz_cols_lc.get(sc_lc)
            if enz_col_orig is None:
                continue
            sig = enz_df[enz_col_orig].values.astype(float)
            # ENZ switch columns are not always binary 0/1 — they can be fractional
            # (e.g. QNO=0.2 meaning 20% of timestep active). Use sig > 0 to detect any injection.
            active_idx = np.where(sig > 0)[0]
            if len(active_idx) == 0:
                continue

            start_s   = float(t_arr[active_idx[0]])
            end_s     = float(t_arr[active_idx[-1]])
            # Add the duration of the last active row so that N active rows
            # → dur = N × dt (not (N-1) × dt which misses the final row's contribution).
            last_idx = int(active_idx[-1])
            dt_last = float(t_arr[last_idx + 1] - t_arr[last_idx]) \
                      if last_idx + 1 < len(t_arr) else 60.0
            dur_s     = max(60.0, end_s - start_s + dt_last)
            start_utc = _s_to_utc(start_s)

            src_products = _products_of(rhs_str)
            # True TwoReaction: source products are CONSUMED by chain reactions
            # (i.e. they appear as LHS of injection_chain_rxns), not just listed
            # in injection_intermediate_species (which includes ALL source products).
            chain_lhs = {c_lhs for c_lhs, _, _ in injection_chain_rxns}
            is_two_rxn = any(sp in chain_lhs for sp in src_products)

            # Determine ultimate target
            ultimate_target = "UNKNOWN"
            if is_two_rxn:
                for _il, ir, _cr in injection_chain_rxns:
                    for sp in _products_of(ir):
                        if sp not in injection_intermediate_species:
                            ultimate_target = sp
                            break
                    if ultimate_target != "UNKNOWN":
                        break
                if ultimate_target == "UNKNOWN" and src_products:
                    ultimate_target = src_products[0]
            else:
                all_src_rxns = [(r, p, k, sc)
                                for r, p, k, sc in injection_source_rxns
                                if r == lhs_up and sc and sc.lower() == sc_lc]
                best = max(all_src_rxns,
                           key=lambda x: (_ppb_from_rate(x[2]) or 0),
                           default=(lhs_up, rhs_str, rate_expr, switch_col))
                prods = _products_of(best[1])
                ultimate_target = prods[0] if prods else "UNKNOWN"

            ppb_per_min = _ppb_from_rate(rate_expr)
            # Scale by the mean ENZ switch value during active window (may be < 1.0)
            mean_switch = float(np.mean(sig[active_idx])) if len(active_idx) else 1.0
            ppb_est = str(round(ppb_per_min * mean_switch * dur_s / 60.0, 2)) \
                      if ppb_per_min else "?"

            entry = {
                "name":       f"{ultimate_target} (IDL)",
                "target":     ultimate_target,
                "ppb":        ppb_est,
                "start_utc":  start_utc,
                "duration_s": str(int(round(dur_s))),
                "profile":    "TwoReaction" if is_two_rxn else "Gaussian",
                "source":     "IDL",
                "comment":    f"k-rate: {rate_expr}",
                "t90_frac":   "0.5",
            }

            if is_two_rxn:
                k_fast = k_slow = k_delivery = None
                for _il, _ir, cr in injection_chain_rxns:
                    v = _eval_const(cr)
                    v_str = str(round(v, 6)) if v else cr
                    prods_il = _products_of(_ir)
                    if any(p in injection_intermediate_species for p in prods_il):
                        if k_fast is None:
                            k_fast = v_str
                        elif k_slow is None:
                            k_slow = v_str
                    else:
                        k_delivery = v_str
                entry["k_fast"]     = k_fast     or "0.016667"
                entry["k_slow"]     = k_slow     or "0.000556"
                entry["k_delivery"] = k_delivery or "0.016667"
                entry["fast_frac"]  = "0.9"

            burst_injections.append(entry)

    # ---- resolve mechanism files recursively --------------------------------
    resolved_mech_files = []
    if os.path.isdir(chemistry_folder):
        for fname in mechanism_files_raw:
            rel = _find_file_recursive(fname, chemistry_folder)
            resolved_mech_files.append(rel if rel else fname)
    else:
        resolved_mech_files = list(mechanism_files_raw)

    # ---- pressure -> env_map ------------------------------------------------
    if pressure_hpa is not None:
        env_map["PRESS"] = {"type": "Fixed Constant", "value": str(pressure_hpa)}

    # ---- assemble config ----------------------------------------------------
    config = {
        "folder":                    easy_dir,
        "chemistry_folder":          chemistry_folder,
        "chamber_vol":               chamber_vol,
        "mechanism_files":           resolved_mech_files,
        "env_map":                   env_map,
        "j_scalings":                j_scalings,
        "artifacts": {
            "active":      True,
            "scale":       scale,
            "hono_active": False,   # handled as chamber_source_reaction
            "hono_eqn":    "",
            "hcho_active": False,   # handled as chamber_source_reaction
            "hcho_eqn":    "",
            "no2_active":  False,   # handled as chamber_source_reaction (if active)
            "no2_eqn":     "",
            "y_active":    y_active,
            "y_eqn":       y_eqn,
            "wall_active": False,   # handled as chamber_source_reaction
            "dark_hono":   "",
        },
        "initials":                  initials,
        "continuous_sources":        ([{"target": "CO2", "rate": co2_bg_rate}]
                                      if co2_bg_rate else []),
        "burst_injections":          burst_injections,
        "injection_reactions":       injection_reactions,
        "chamber_source_reactions":  chamber_source_reactions,
        "custom_reactions":          [],
        "excluded_reactions":        [],
        "blocked_parameters":        [],
        "h2o_injection_impurities":  h2o_imps,
        "_skipped":                  skipped_notes,
    }

    _print_summary(config)
    return config


def _print_summary(cfg: dict):
    W = 60
    print(f"\n[idl_parser] {'=' * W}")
    print(f"[idl_parser]  EASY FILE INTEGRATION REPORT")
    print(f"[idl_parser] {'=' * W}")

    print(f"[idl_parser] MECHANISM FILES ({len(cfg['mechanism_files'])}):")
    for f in cfg["mechanism_files"]:
        print(f"[idl_parser]   [OK]{f}")

    print(f"[idl_parser] CONSTRAINED VARIABLES ({len(cfg['env_map'])}):")
    for k, v in cfg["env_map"].items():
        tag = f"SAPHIR File ->{v['value']}" if v['type'] == "SAPHIR File" else f"Fixed = {v['value']}"
        print(f"[idl_parser]   [OK]{k:12s}  [{tag}]")

    print(f"[idl_parser] J-SCALINGS ({len(cfg['j_scalings'])}):")
    for k, v in cfg["j_scalings"].items():
        print(f"[idl_parser]   [OK]{k} = {v}")

    print(f"[idl_parser] INJECTION REACTIONS ({len(cfg['injection_reactions'])}):")
    for r in cfg["injection_reactions"]:
        print(f"[idl_parser]   [OK][{r['reactants']}] --> {r['products']}  "
              f"rate={r['rate'][:60]}")

    csr = cfg.get("chamber_source_reactions", [])
    print(f"[idl_parser] CHAMBER SOURCE REACTIONS ({len(csr)}):")
    for r in csr:
        print(f"[idl_parser]   [OK][{r['reactants']}] --> {r['products']}  "
              f"rate={r['rate'][:60]}")

    arts = cfg["artifacts"]
    y_on = arts.get("y_active")
    print(f"[idl_parser] Y BACKGROUND REACTIVITY: "
          f"{'[OK] active  y_ppb=' + str(arts.get('y_eqn')) if y_on else '[--] inactive'}")

    print(f"[idl_parser] INITIAL CONDITIONS ({len(cfg['initials'])}):")
    for sp, val in cfg["initials"].items():
        print(f"[idl_parser]   [OK] {sp} = {val} ppb")

    cs = cfg["continuous_sources"]
    print(f"[idl_parser] CO2 BACKGROUND SOURCE: "
          f"{'[OK] ' + cs[0]['rate'] if cs else '[--] none'}")

    imps = cfg["h2o_injection_impurities"]
    print(f"[idl_parser] H2O INJECTION IMPURITIES ({len(imps)}):")
    for imp in imps:
        print(f"[idl_parser]   [OK]{imp['target']} @ {imp['rate_ppb_per_min']} ppb/min")

    bi = cfg["burst_injections"]
    print(f"[idl_parser] BURST INJECTION SUMMARIES ({len(bi)}):")
    for b in bi:
        print(f"[idl_parser]   [OK]{b['name']}  target={b['target']}  "
              f"~{b['ppb']} ppb  start={b['start_utc']}")

    sk = cfg.get("_skipped", [])
    print(f"[idl_parser] SKIPPED / INACTIVE ({len(sk)}):")
    for note in sk:
        print(f"[idl_parser]   [--]{note}")

    print(f"[idl_parser] {'=' * W}\n")


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: idl_parser.py  <scenario.easy>  [kammer.enz]  [output.json]")
        sys.exit(1)

    easy_p = sys.argv[1]
    enz_p  = sys.argv[2] if len(sys.argv) > 2 else None
    out_p  = (sys.argv[3] if len(sys.argv) > 3
              else os.path.splitext(easy_p)[0] + "_idl.json")

    cfg = parse_idl_easy(easy_p, enz_p)

    with open(out_p, "w", encoding="utf-8") as fh:
        json.dump(cfg, fh, indent=4)
    print(f"[idl_parser] Written -> {out_p}")
