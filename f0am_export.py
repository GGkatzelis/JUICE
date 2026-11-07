"""
f0am_export.py — write a JUICE scenario as an F0AM run with identical inputs.

Purpose: a model intercomparison in which only the ODE integration differs.  F0AM gets

  * the same reactions and stoichiometry JUICE integrates (after include/override handling);
  * for every reaction and every output step, the rate constant JUICE used in that step
    (JUICE, like FACSIMILE, freezes all settings at the start of each output step), so photolysis,
    T/p/H2O dependence, chamber sources, wall loss, dilution and injection switches are identical;
  * the same initial concentrations and the same output grid (one F0AM step = one output step).

How the pieces map onto F0AM (format: F0AM mechanism scripts, read from F0AM's documentation and
example files; no F0AM code is used or copied):

  JUICE                                        F0AM
  reaction, rate constant per step             k(:,i) = constant  or  column of a k table (.mat)
  reactants that are settings (DIL, M, H2O,    folded into k (they are frozen per step in JUICE)
    injection switches Q..., ENZ inputs)
  injection  Q --> X                           zero-order reaction (no reactant, F0AM's 'ONE')
  RO2 sum (setting, frozen per step)           RO2ToAdd + 'RO2' as second reactant (F0AM sums it
                                                 at every RHS call: a deliberate numerical difference)
  O1D in pseudo-steady state                   O1D as an ordinary species
  rescaling y *= M_new/M_old at step bounds    first-order term per species, spread over the step
  F0AM's own M (6.022e23 P/(8.314e4 T))        only converts ppb <-> molec/cm3; rate constants
                                                 already contain JUICE's M, initial values and output
                                                 are converted with F0AM's own NumberDensity

Files per scenario (in <out>/<stem>/):
  <stem>_JUICE_mech.m     F0AM mechanism script (species, RO2 list, reactions)
  <stem>_JUICE_k.mat      time-varying rate constants, KT(step, column)
  <stem>_JUICE_inputs.mat P, T, H2O per step, initial concentrations, time grid
  run_<stem>.m            run script (F0AM v3/v4 on the MATLAB path); writes <stem>_F0AM.csv
  <stem>_JUICE_J.m        placeholder J function (one zero): F0AM's InitializeChemistry needs one
  <stem>_JUICE.csv        the JUICE result on the same grid (molec/cm3), for the comparison
  <stem>_export.txt       what was exported, checks, notes

After the MATLAB run:  --compare  reports F0AM vs JUICE and vs FACSIMILE with the benchmark metric.

Check without MATLAB:  --check  integrates the written files with an F0AM-like step scheme in Python
(RO2 summed dynamically, one stiff solve per step) and compares that with JUICE.

    python f0am_export.py scene_04 scene_05 uv_lights_on_off --out paper/F0AM_scenarios --check
"""
from __future__ import annotations

import argparse
import contextlib
import datetime as _dt
import io
import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
TOML = os.path.join(HERE, "FR_GUI.toml")
_NAME_OK = re.compile(r"^[A-Za-z][A-Za-z0-9_]{0,61}$")
_RESERVED = {"ONE", "RO2"}
COMPARE_FLOOR = 2.5e4          # molec/cm3, the floor used by the FACSIMILE benchmark


# ─────────────────────────────────────────────────────────────────────────────
# Recording JUICE's run
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Recording:
    stem: str
    species: List[str]                     # F0AM species (JUICE ODE species + PSS species)
    ro2: List[str]
    rxn_names: List[str]
    rxn_reactants: List[List[str]]         # <= 2 species names each ([] = zero order)
    rxn_stoich: List[Dict[str, float]]     # species -> net coefficient
    k: np.ndarray                          # (n_steps, n_rxn) rate constants incl. folded settings
    t0: np.ndarray                         # start of each step (s since 2000-01-01, JUICE clock)
    dt: float
    P: np.ndarray
    T: np.ndarray
    H2O: np.ndarray
    M: np.ndarray                          # JUICE's M per step (reported only)
    y0: Dict[str, float]                   # molec/cm3
    k_loss: Optional[np.ndarray]           # per-species first-order term (expansion/dilution) or None
    juice_t: np.ndarray
    juice_y: Dict[str, np.ndarray]
    notes: List[str] = field(default_factory=list)


def _parse_side(tokens) -> List[tuple]:
    out = []
    for tok in tokens:
        c, n = 1.0, tok.strip().upper()
        if not n:
            continue
        if "*" in n:
            a, b = n.split("*", 1)
            c, n = float(a), b.strip()
        out.append((c, n))
    return out


def record(json_path: str, stem: str) -> Recording:
    """Run JUICE on a prepared scenario and record what the solver used in every output step."""
    import scenario_runner as sr

    r = sr.ScenarioRunner(json_path)
    snaps: Dict[float, tuple] = {}

    def hook(t_abs, ctx_s, y_state, active_s):
        if not hasattr(hook, "buf"):
            hook.buf = np.zeros(len(r.mechanism.reactions))
            hook.buf2 = np.zeros(len(r.mechanism.reactions))
            hook.ro2_rx = np.array([i for i, x in enumerate(r.mechanism.reactions)
                                    if "RO2" in {n.upper() for n in x.safe_expr.co_names}], dtype=int)
            hook.count = 0
            hook.ro2_bad = 0.0
        buf = hook.buf
        r._compute_rates_fn(ctx_s, buf)
        np.nan_to_num(buf, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
        row = buf.copy()
        if hook.ro2_rx.size:
            # rate constants of the RO2 reactions per unit RO2 (F0AM multiplies by its own RO2 sum)
            ro2_now = float(ctx_s["RO2"]) if "RO2" in ctx_s else 0.0
            ctx2 = sr.SafeDict(ctx_s); ctx2["RO2"] = 1.0
            r._compute_rates_fn(ctx2, hook.buf2)
            np.nan_to_num(hook.buf2, nan=0.0, posinf=0.0, neginf=0.0, copy=False)
            if hook.count % 50 == 0:
                # proportional to RO2, and nothing else depends on it
                want = hook.buf2.copy(); want[hook.ro2_rx] *= ro2_now
                den = np.maximum(np.abs(buf), 1e-300)
                dev = np.abs(buf - want) / den
                dev[np.abs(buf) < 1e-300] = 0.0
                hook.ro2_bad = max(hook.ro2_bad, float(dev.max()))
            row[hook.ro2_rx] = hook.buf2[hook.ro2_rx]
        hook.count += 1
        env = {n: float(ctx_s[n]) if n in ctx_s else 0.0 for n in hook.env_names}
        T = float(ctx_s["T"]); P = float(ctx_s["P"]) if "P" in ctx_s else float(ctx_s.get("PRESS", 1013.25))
        snaps[float(t_abs)] = (row, env, float(ctx_s["RO2"]) if "RO2" in ctx_s else 0.0,
                               np.array(y_state, dtype=float), T, P, float(ctx_s.get("H2O", 0.0)),
                               float(r.get_M(t_abs)),
                               [float(eval(c, sr.__dict__, ctx_s)) for c in hook.cs_codes])

    with contextlib.redirect_stdout(io.StringIO()):
        r._setup_environment()
        r._build_mechanism()

    cfg = r.cfg
    problems = []
    # continuous source terms (chamber set-up, e.g. a CO2 background source): zero-order reactions,
    # evaluated per step with the frozen settings, exactly as the solver evaluates them
    cs = [x for x in cfg.get("continuous_sources", []) if x.get("target") and x.get("rate")]
    hook.cs_codes = [compile(x["rate"], "<string>", "eval") for x in cs]
    if any(b.get("ppb") and b.get("source") != "IDL" for b in cfg.get("burst_injections", [])):
        problems.append("burst injections outside the EASY/ENZ chain are not exported")
    if cfg.get("h2o_injection_impurities"):
        problems.append("H2O-injection impurities are not exported")
    if problems:
        raise ValueError("; ".join(problems))

    names_env = set()
    for rxn in r.mechanism.reactions:
        for c, n in _parse_side(rxn.reactants):
            names_env.add(n)
    hook.env_names = sorted(names_env)          # superset; species are filtered later
    r.snap_hook = hook

    with contextlib.redirect_stdout(io.StringIO()):
        r._solve()

    sp_names = [s.upper() for s in r.solver.species]
    sp_set = set(sp_names)
    pss = {s.upper() for s in getattr(r, "_pss_species", set())}
    co_names = set()
    for rxn in r.mechanism.reactions:
        co_names.update(n.upper() for n in rxn.safe_expr.co_names)
    dep = co_names & (sp_set | pss)
    if dep:
        raise ValueError(f"rate expressions depend on species {sorted(dep)[:8]} — not row-constant, "
                         "cannot be exported as a k table")
    if getattr(hook, "ro2_bad", 0.0) > 1e-9:
        raise ValueError(f"rate constants are not proportional to RO2 (deviation {hook.ro2_bad:.2e})")
    if not snaps:
        raise ValueError("no settings snapshots recorded — run the scenario in a per-ENZ-row solver mode")

    t_span = np.asarray(r.t_span, dtype=float)
    steps = t_span[:-1]
    dts = np.diff(t_span)
    if not np.allclose(dts, dts[0], rtol=0, atol=1e-6):
        raise ValueError("output steps are not uniform; F0AM needs one IntTime for all steps")
    dt = float(dts[0])
    keys = np.array(sorted(snaps))
    rows = []
    for t in steps:
        j = np.searchsorted(keys, t - 1e-6)
        if j >= len(keys) or abs(keys[j] - t) > 1e-3:
            raise ValueError(f"no snapshot recorded for the step starting at {t}")
        rows.append(snaps[keys[j]])
    n_steps = len(rows)

    # RO2 members (last definition wins, as in EASY)
    ro2_def = ""
    for e in r.mechanism.s_entries:
        if e.name.upper() == "RO2":
            ro2_def = e.definition
    ro2 = [tok for tok in dict.fromkeys(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", ro2_def.upper()))
           if tok in sp_set]

    f_species = list(sp_names) + sorted(pss - sp_set)
    bad = [s for s in f_species if not _NAME_OK.match(s) or s in _RESERVED]
    if bad:
        raise ValueError(f"species names not usable in F0AM: {bad[:10]}")
    f_set = set(f_species)

    rxn_names, rxn_reac, rxn_st, env_entries, ro2_rxn = [], [], [], [], []
    for i, rxn in enumerate(r.mechanism.reactions):
        reac, st = [], {}
        for c, n in _parse_side(rxn.reactants):
            if n in f_set:
                if abs(c - round(c)) > 1e-12 or c < 1:
                    raise ValueError(f"reaction {i}: non-integer reactant coefficient {c}*{n}")
                reac += [n] * int(round(c))
                st[n] = st.get(n, 0.0) - c
            else:
                env_entries.append((i, n, c))
        for c, n in _parse_side(rxn.products):
            if n in f_set:
                st[n] = st.get(n, 0.0) + c
        if len(reac) > 2:
            raise ValueError(f"reaction {i} has {len(reac)} species reactants; F0AM allows two")
        if "RO2" in {x.upper() for x in rxn.safe_expr.co_names}:
            if len(reac) > 1:
                raise ValueError(f"reaction {i}: RO2 reaction with two species reactants")
            ro2_rxn.append(i)
        rxn_names.append("+".join(x for x in (rxn.reactants or [])) + "=" + "+".join(rxn.products or []))
        rxn_reac.append(reac)
        rxn_st.append({k: v for k, v in st.items() if v != 0.0})

    K = np.vstack([row[0] for row in rows])                          # (n_steps, n_rxn)
    fac = np.ones(K.shape)
    for i, n, c in env_entries:
        vals = np.array([row[1].get(n, 0.0) for row in rows])
        fac[:, i] *= vals ** c
    K = K * fac
    cs_notes = []
    for j, x in enumerate(cs):
        tgt = x["target"].upper()
        used = {n.upper() for n in hook.cs_codes[j].co_names}
        if used & (sp_set | pss):
            raise ValueError(f"continuous source of {tgt} depends on species {sorted(used & sp_set)}")
        if tgt not in f_set:
            continue                                     # JUICE ignores it as well
        K = np.column_stack([K, np.array([row[8][j] for row in rows])])
        rxn_names.append(f"={tgt}"); rxn_reac.append([]); rxn_st.append({tgt: 1.0})
        cs_notes.append(f"continuous source of {tgt} ({x['rate']}) as a zero-order reaction"
                        + (" - zero in every step, as in JUICE" if not np.any(K[:, -1]) else ""))
    n_rxn = len(rxn_names)
    ro2_snap = np.array([row[2] for row in rows])
    notes = list(cs_notes)
    if ro2_rxn:
        for i in ro2_rxn:
            rxn_reac[i] = rxn_reac[i] + ["RO2"]
        notes.append(f"{len(ro2_rxn)} reactions use RO2 ({len(ro2)} members); JUICE freezes RO2 per step, "
                     "F0AM sums it at every RHS call")

    # rescaling at step boundaries (M changes) + analytical dilution (only without explicit DIL reactions)
    Ms = np.array([row[7] for row in rows])
    k_loss = np.zeros(n_steps)
    k_loss[1:] = -np.log(Ms[1:] / Ms[:-1]) / dt
    if np.any(k_loss != 0):
        notes.append(f"air density changes between steps (max {np.abs(np.log(Ms[1:]/Ms[:-1])).max():.2e} "
                     "per step): JUICE rescales all concentrations at the boundary, F0AM gets the same "
                     "factor as a first-order term over the following step")
    if cfg.get("analytical_dilution", True) and not r._has_explicit_dil:
        flow = r.ctx_data["FLOW"] if "FLOW" in r.ctx_data else (lambda t: 0.0)
        dil_sw = np.array([row[1].get("DIL", 1.0) for row in rows])
        kd = np.array([(float(flow(t)) + 0.6) / (cfg["chamber_vol"] * 3600.0) for t in steps]) * dil_sw
        k_loss = k_loss + kd
        notes.append("analytical dilution (no explicit DIL reactions) exported as first-order loss")
    if not np.any(k_loss != 0):
        k_loss = None

    y_first = rows[0][3]
    y0 = {s: float(max(y_first[j], 0.0)) for j, s in enumerate(sp_names)}
    for s in pss - sp_set:
        y0[s] = 0.0
    sol_t = np.asarray(r.sol.t, dtype=float)
    juice_y = {s: np.maximum(r.sol.y[j], 0.0) for j, s in enumerate(sp_names)}
    if "O1D" in pss:
        notes.append("O1D: pseudo-steady state in JUICE, an ordinary species in F0AM")
    return Recording(stem=stem, species=f_species, ro2=ro2, rxn_names=rxn_names, rxn_reactants=rxn_reac,
                     rxn_stoich=rxn_st, k=K, t0=steps, dt=dt,
                     P=np.array([row[5] for row in rows]), T=np.array([row[4] for row in rows]),
                     H2O=np.array([row[6] for row in rows]), M=Ms, y0=y0, k_loss=k_loss,
                     juice_t=sol_t, juice_y=juice_y, notes=notes)


# ─────────────────────────────────────────────────────────────────────────────
# Writing the F0AM files
# ─────────────────────────────────────────────────────────────────────────────

def _g(x: float) -> str:
    return repr(float(x)).replace("e+", "e")


def write(rec: Recording, out_root: str) -> str:
    from scipy.io import savemat
    stem = rec.stem
    d = os.path.join(out_root, stem)
    os.makedirs(d, exist_ok=True)
    n_steps, n_rxn = rec.k.shape

    const = np.all(rec.k == rec.k[0:1, :], axis=0)
    var_idx = np.where(~const)[0]
    col_of = {int(i): c + 1 for c, i in enumerate(var_idx)}      # MATLAB 1-based columns
    KT = rec.k[:, var_idx]
    loss_col = None
    if rec.k_loss is not None:
        KT = np.column_stack([KT, rec.k_loss]) if KT.size else rec.k_loss[:, None]
        loss_col = KT.shape[1]
    savemat(os.path.join(d, f"{stem}_JUICE_k.mat"), {"KT": KT if KT.size else np.zeros((n_steps, 0))},
            do_compression=True)

    names = rec.species
    init = np.array([rec.y0.get(s, 0.0) for s in names])
    savemat(os.path.join(d, f"{stem}_JUICE_inputs.mat"), {
        "P": rec.P[:, None], "T": rec.T[:, None], "H2O": rec.H2O[:, None], "M_JUICE": rec.M[:, None],
        "t0": float(rec.t0[0]), "dt": rec.dt, "nsteps": float(n_steps),
        "init_names": np.array(names, dtype=object)[:, None], "init_molec": init[:, None]},
        do_compression=True)

    today = _dt.date.today().isoformat()
    L = [f"% {stem}_JUICE_mech.m - F0AM mechanism script generated by JUICE f0am_export.py ({today})",
         f"% {len(names)} species, {n_rxn} reactions"
         + (f" + {len(names)} first-order terms (density change / dilution)" if loss_col else ""),
         "% Rate constants are the values JUICE used in each output step (settings frozen per step):",
         f"% constants are written inline, time-varying ones come from {stem}_JUICE_k.mat (KT).",
         "", "SpeciesToAdd = {..."]
    L += _cell_lines(names)
    L += ["};", "RO2ToAdd = {..."] + _cell_lines(rec.ro2) + ["};", "AddSpecies", ""]
    if KT.size:
        L += [f"JUICEk_ = load('{stem}_JUICE_k.mat','KT'); JUICEk_ = JUICEk_.KT;", ""]
    for i in range(n_rxn):
        kexpr = f"JUICEk_(:,{col_of[i]})" if i in col_of else _g(rec.k[0, i])
        L.append("i=i+1;")
        L.append(f"Rnames{{i}} = '{rec.rxn_names[i]}';")
        L.append(f"k(:,i) = {kexpr};")
        for j, g in enumerate(rec.rxn_reactants[i]):
            L.append(f"Gstr{{i,{j + 1}}} = '{g}';")
        fl = [f"f{s}(i)=f{s}(i){'+' if c > 0 else '-'}{_g(abs(c))};" for s, c in rec.rxn_stoich[i].items()]
        if fl:
            L.append(" ".join(fl))
        L.append("")
    if loss_col:
        L.append("% density change between steps (JUICE: y*M_new/M_old at the step boundary) and dilution")
        for s in names:
            L += ["i=i+1;", f"Rnames{{i}} = '{s}=(density change/dilution)';", f"k(:,i) = JUICEk_(:,{loss_col});",
                  f"Gstr{{i,1}} = '{s}';", f"f{s}(i)=f{s}(i)-1;", ""]
    with open(os.path.join(d, f"{stem}_JUICE_mech.m"), "w", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join(L) + "\n")

    write_scripts(d, stem, today)

    # JUICE result on the F0AM output grid (end of every step), molec/cm3
    t_end = rec.t0 + rec.dt
    jt = rec.juice_t
    idx = [int(np.argmin(np.abs(jt - t))) for t in t_end]
    cols = [s for s in names if s in rec.juice_y]
    with open(os.path.join(d, f"{stem}_JUICE.csv"), "w", encoding="ascii") as fh:
        fh.write("TIME," + ",".join(cols) + "\n")
        Y = np.column_stack([rec.juice_y[s][idx] for s in cols])
        for t, row in zip(t_end, Y):
            fh.write(f"{t:.1f}," + ",".join(f"{v:.10e}" for v in row) + "\n")

    with open(os.path.join(d, f"{stem}_export.txt"), "w", encoding="utf-8") as fh:
        fh.write(f"{stem}: F0AM export of the JUICE scenario ({today})\n")
        fh.write(f"species {len(names)}, reactions {n_rxn}, RO2 members {len(rec.ro2)}, steps {n_steps} x {rec.dt:g} s\n")
        fh.write(f"time-varying rate constants: {len(var_idx)} of {n_rxn} (in {stem}_JUICE_k.mat)\n")
        fh.write(f"T {rec.T.min():.2f}-{rec.T.max():.2f} K, p {rec.P.min():.1f}-{rec.P.max():.1f} hPa\n")
        for n in rec.notes:
            fh.write(f"note: {n}\n")
        if np.any(rec.H2O < 0):
            fh.write(f"note: H2O input is negative in {int(np.sum(rec.H2O < 0))} steps (min {rec.H2O.min():.3g}); "
                     "JUICE used these values in the rate constants, F0AM's Met.H2O is clipped at 0 (inert there)\n")
    return d


def write_scripts(d: str, stem: str, today: Optional[str] = None) -> None:
    """run_<stem>.m and the placeholder J function. Separate from write() so that a change to the run
    template can be applied to an existing export (--scripts-only) without running JUICE again."""
    today = today or _dt.date.today().isoformat()
    with open(os.path.join(d, f"run_{stem}.m"), "w", encoding="ascii", errors="replace") as fh:
        fh.write(RUN_TEMPLATE.replace("{STEM}", stem).replace("{DATE}", today))
    # F0AM's InitializeChemistry (v3 core, also in F0AM4GG) evaluates ChemFiles{2} to build its jcorr_all
    # output and leaves it unassigned when the entry is empty ("Output argument jcorr_all not assigned").
    # Photolysis is already inside the rate constants, so the J function returns a single unused zero.
    with open(os.path.join(d, f"{stem}_JUICE_J.m"), "w", encoding="ascii", errors="replace") as fh:
        fh.write("\n".join([
            f"function J = {stem}_JUICE_J(Met)",
            f"% {stem}_JUICE_J.m - placeholder J function generated by JUICE f0am_export.py ({today})",
            "% All photolysis frequencies are inside the rate constants from JUICE. F0AM's InitializeChemistry",
            "% still needs a J function (it sets jcorr_all from it), so this returns one unused J value of 0.",
            "J.J0 = zeros(size(Met.T(:)));",
            ""]))


def _cell_lines(names, per=8):
    out = []
    for i in range(0, len(names), per):
        out.append(" ".join(f"'{n}';" for n in names[i:i + per]) + " ...")
    return out


RUN_TEMPLATE = r"""%% run_{STEM}.m - F0AM run of the JUICE scenario {STEM} (generated {DATE} by JUICE f0am_export.py)
% Identical inputs: every reaction and every rate constant per output step is what JUICE used;
% one F0AM step = one JUICE/FACSIMILE output step. Only the ODE integration differs.
% Needs F0AM (v3 or v4) on the MATLAB path. Writes {STEM}_F0AM.csv (molec/cm3, end of each step,
% time in s since 2000-01-01 as in the JUICE and FACSIMILE output).
%
% TIGHT = true repeats the run with ode15s RelTol 1e-6 (F0AM's default is MATLAB's 1e-3). F0AM has
% no option for this, so a patched copy of F0AM_ModelCore is written into this folder for that run.

clear
TIGHT = false;

here = fileparts(mfilename('fullpath'));
cd(here)
D = load('{STEM}_JUICE_inputs.mat');
n = round(D.nsteps);

Met = {...
    'P'     D.P;...                 % hPa (= mbar), as JUICE
    'T'     D.T;...                 % K
    'H2O'   max(D.H2O, 0);...       % molec/cm3, carried along only: every rate already contains JUICE's
    ...                             %   H2O (negative measured values included); F0AM rejects negatives
    'kdil'  zeros(n,1);...          % dilution is in the mechanism, as in JUICE
    };

M1 = NumberDensity(D.P(1), D.T(1));         % F0AM's own conversion, so the initial molec/cm3 are exact
names = cellstr(D.init_names(:));
InitConc = [names num2cell(D.init_molec(:)./M1.*1e9) num2cell(zeros(numel(names),1))];

ChemFiles = {...
    '';...                          % no generic K: rate constants come from JUICE
    '{STEM}_JUICE_J(Met)';...        % placeholder (one J = 0): photolysis is inside the rate constants,
    '{STEM}_JUICE_mech';...         %   but F0AM's InitializeChemistry needs a J function for jcorr_all
    };

BkgdConc = {'DEFAULT' 0};

ModelOptions.Verbose       = 1;
ModelOptions.EndPointsOnly = 1;
ModelOptions.LinkSteps     = 1;
ModelOptions.Repeat        = 1;
ModelOptions.IntTime       = D.dt;

core = @F0AM_ModelCore;
tag = '';
if TIGHT
    src = fileread(which('F0AM_ModelCore'));
    src = strrep(src, 'function S = F0AM_ModelCore(', 'function S = F0AM_ModelCore_tight(');
    nrep = numel(strfind(src, 'odeset(''Jacobian'''));
    if nrep == 0, error('could not find the odeset call in F0AM_ModelCore'); end
    src = strrep(src, 'odeset(''Jacobian''', 'odeset(''RelTol'',1e-6,''Jacobian''');
    fid = fopen(fullfile(here, 'F0AM_ModelCore_tight.m'), 'w'); fwrite(fid, src); fclose(fid);
    rehash
    core = @F0AM_ModelCore_tight;
    tag = '_tight';
end

timer_ = tic;                       % own handle: F0AM calls tic itself and would reset a bare toc
S = core(Met, InitConc, ChemFiles, BkgdConc, ModelOptions);
runtime = toc(timer_);

t = D.t0 + D.dt*(1:n)';
sp = setdiff(fieldnames(S.Conc), {'ONE'}, 'stable');
X = zeros(n, numel(sp));
for j = 1:numel(sp)
    X(:,j) = S.Conc.(sp{j}) .* S.Met.M ./ 1e9;          % ppb -> molec/cm3 with F0AM's M
end
out = ['{STEM}_F0AM' tag '.csv'];
fid = fopen(out, 'w');
fprintf(fid, 'TIME'); fprintf(fid, ',%s', sp{:}); fprintf(fid, '\n');
fprintf(fid, [repmat('%.10e,', 1, numel(sp)) '%.10e\n'], [t X]');
fclose(fid);
fid = fopen(['{STEM}_F0AM' tag '_runtime.txt'], 'w'); fprintf(fid, '%.3f s\n', runtime); fclose(fid);
fprintf('F0AM run time %.1f s, written %s\n', runtime, out)
"""


# ─────────────────────────────────────────────────────────────────────────────
# Check without MATLAB: integrate the written files the way F0AM does
# ─────────────────────────────────────────────────────────────────────────────

def read_mech(path: str) -> dict:
    """Parse a mechanism script written by write() (this format only)."""
    txt = open(path, encoding="ascii").read()
    def cell(name):
        m = re.search(name + r" = \{\.\.\.\n(.*?)\n\};", txt, re.S)
        return re.findall(r"'([^']+)'", m.group(1)) if m else []
    species, ro2 = cell("SpeciesToAdd"), cell("RO2ToAdd")
    rx = []
    for blk in txt.split("i=i+1;\n")[1:]:
        kk = re.search(r"k\(:,i\) = (.+?);", blk).group(1)
        g = re.findall(r"Gstr\{i,\d\} = '([^']+)';", blk)
        st = {s: float(c) for s, c in re.findall(r"f(\w+)\(i\)=f\w+\(i\)([+-][0-9.eE+-]+);", blk)}
        rx.append((kk, g, st))
    return {"species": species, "ro2": ro2, "rx": rx}


def emulate(d: str, stem: str, rtol: float = 1e-7, ro2_frozen: bool = False) -> tuple:
    """Integrate the exported files step by step (F0AM scheme: k per step, RO2 summed at every call).

    ro2_frozen=True holds RO2 at its step-start value instead (the JUICE/FACSIMILE scheme): then the
    emulation must reproduce JUICE to solver precision, which tests the export itself."""
    from scipy.io import loadmat
    from scipy.integrate import solve_ivp
    from scipy import sparse
    mech = read_mech(os.path.join(d, f"{stem}_JUICE_mech.m"))
    inp = loadmat(os.path.join(d, f"{stem}_JUICE_inputs.mat"))
    KT = loadmat(os.path.join(d, f"{stem}_JUICE_k.mat"))["KT"]
    n = int(round(inp["nsteps"].item()))
    sp = mech["species"]; ix = {s: j for j, s in enumerate(sp)}; ns = len(sp); nr = len(mech["rx"])
    ro2_idx = np.array([ix[s] for s in mech["ro2"]], dtype=int)
    kconst = np.zeros(nr); kcol = np.full(nr, -1)
    g1 = np.full(nr, -1); g2 = np.full(nr, -1)        # -1 = ONE, -2 = RO2
    rows, cols, vals = [], [], []
    for r, (kk, g, st) in enumerate(mech["rx"]):
        m = re.match(r"JUICEk_\(:,(\d+)\)", kk)
        if m: kcol[r] = int(m.group(1)) - 1
        else: kconst[r] = float(kk)
        gi = [(-2 if x == "RO2" else ix[x]) for x in g]
        if len(gi) > 0: g1[r] = gi[0]
        if len(gi) > 1: g2[r] = gi[1]
        for s, c in st.items():
            rows.append(r); cols.append(ix[s]); vals.append(c)
    F = sparse.csr_matrix((vals, (rows, cols)), shape=(nr, ns))
    FT = F.T.tocsr()
    is_ro2_2 = g2 == -2
    Rm = sparse.csr_matrix((np.ones(int(is_ro2_2.sum()) * len(ro2_idx)),
                            (np.repeat(np.where(is_ro2_2)[0], len(ro2_idx)),
                             np.tile(ro2_idx, int(is_ro2_2.sum())))), shape=(nr, ns))

    def conc_of(y, g, ro2):
        out = np.ones(nr)
        sel = g >= 0; out[sel] = y[g[sel]]
        out[g == -2] = ro2
        return out

    y = np.array([float(v) for v in inp["init_molec"].ravel()])
    names0 = [str(x[0]) if hasattr(x, "__len__") else str(x) for x in inp["init_names"].ravel()]
    y = np.array([y[names0.index(s)] for s in sp])
    out = np.zeros((n, ns))
    dt = float(inp["dt"].item())
    r1 = np.arange(nr)
    for i in range(n):
        k = kconst.copy()
        sel = kcol >= 0
        k[sel] = KT[i, kcol[sel]]

        ro2_0 = y[ro2_idx].sum() if len(ro2_idx) else 0.0

        def f(t, yy):
            yy = np.maximum(yy, 0.0)
            ro2 = ro2_0 if ro2_frozen else (yy[ro2_idx].sum() if len(ro2_idx) else 0.0)
            v = k * conc_of(yy, g1, ro2) * conc_of(yy, g2, ro2)
            return FT @ v

        def jac(t, yy):
            yy = np.maximum(yy, 0.0)
            ro2 = ro2_0 if ro2_frozen else (yy[ro2_idx].sum() if len(ro2_idx) else 0.0)
            c1, c2 = conc_of(yy, g1, ro2), conc_of(yy, g2, ro2)
            rr, cc, vv = [], [], []
            s1 = g1 >= 0; rr.append(r1[s1]); cc.append(g1[s1]); vv.append((k * c2)[s1])
            s2 = g2 >= 0; rr.append(r1[s2]); cc.append(g2[s2]); vv.append((k * c1)[s2])
            D = sparse.csr_matrix((np.concatenate(vv), (np.concatenate(rr), np.concatenate(cc))), shape=(nr, ns))
            if is_ro2_2.any() and not ro2_frozen:
                D = D + sparse.diags(np.where(is_ro2_2, k * c1, 0.0)) @ Rm
            return (FT @ D).tocsc()

        s = solve_ivp(f, (0.0, dt), y, method="BDF", rtol=rtol, atol=1e-3, jac=jac, t_eval=[dt])
        if not s.success:
            raise RuntimeError(f"step {i}: {s.message}")
        y = np.maximum(s.y[:, -1], 0.0)
        out[i] = y
    return sp, out


def read_csv(path: str) -> tuple:
    with open(path) as fh:
        head = fh.readline().strip().split(",")
    data = np.loadtxt(path, delimiter=",", skiprows=1, ndmin=2)
    return head[1:], data[:, 0], data[:, 1:]


def deviation_table(names_a, A, names_b, B, floor=COMPARE_FLOOR) -> list:
    """Max |a-b| / max(max|b|, floor) per species (the benchmark metric); B is the reference."""
    ib = {s: j for j, s in enumerate(names_b)}
    out = []
    for j, s in enumerate(names_a):
        if s not in ib:
            continue
        a, b = A[:, j], B[:, ib[s]]
        den = max(np.abs(b).max(), floor)
        out.append((float(np.abs(a - b).max() / den), s, float(np.abs(b).max())))
    return sorted(out, reverse=True)


def check(d: str, stem: str, ro2_frozen: bool = False) -> str:
    sp, Y = emulate(d, stem, ro2_frozen=ro2_frozen)
    names_j, tj, J = read_csv(os.path.join(d, f"{stem}_JUICE.csv"))
    dev = deviation_table(sp, Y, names_j, J)
    big = [x for x in dev if x[2] > COMPARE_FLOOR]
    n2 = sum(1 for x in big if x[0] <= 0.02); n01 = sum(1 for x in big if x[0] <= 1e-3)
    main = {s: v for v, s, _ in dev if s in ("O3", "NO", "NO2", "OH", "HO2")}
    lines = [f"{stem}: {'JUICE-scheme (RO2 frozen per step)' if ro2_frozen else 'F0AM-scheme'} emulation of the written files vs JUICE "
             f"(species above {COMPARE_FLOOR:g} molec/cm3: {len(big)})",
             f"  within 2 %: {n2}/{len(big)}   within 0.1 %: {n01}/{len(big)}",
             "  O3/NO/NO2/OH/HO2: " + ", ".join(f"{s} {main[s]*100:.3g} %" for s in main),
             "  worst: " + ", ".join(f"{s} {v*100:.3g} %" for v, s, _ in big[:5])]
    with open(os.path.join(d, f"{stem}_export.txt"), "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    return "\n".join(lines)


def compare_outputs(d: str, stem: str) -> str:
    """After the MATLAB run: F0AM vs JUICE and, where a FACSIMILE result exists, vs FACSIMILE."""
    out = []
    names_j, tj, J = read_csv(os.path.join(d, f"{stem}_JUICE.csv"))
    refs = [("JUICE", names_j, tj, J)]
    fac = os.path.join(HERE, "MCM_tested_scenarios", "IDL Input", stem, f"{stem}.cs")
    if os.path.isfile(fac):
        from run_benchmark import parse_cs_file
        df = parse_cs_file(fac)
        refs.append(("FACSIMILE", [c.upper() for c in df.columns[1:]], df["TIME"].values, df.iloc[:, 1:].values))
    for tag in ("", "_tight"):
        p = os.path.join(d, f"{stem}_F0AM{tag}.csv")
        if not os.path.isfile(p):
            continue
        names_f, tf, Fv = read_csv(p)
        for label, nr, tr, R in refs:
            idx = [int(np.argmin(np.abs(tr - t))) for t in tf]
            ok = np.abs(tr[idx] - tf) < 1.0
            dev = deviation_table(names_f, Fv[ok], nr, R[idx][ok])
            big = [x for x in dev if x[2] > COMPARE_FLOOR]
            main = {s_: v for v, s_, _ in dev if s_ in ("O3", "NO", "NO2", "OH", "HO2")}
            out += [f"{stem}: F0AM{tag} vs {label} ({int(ok.sum())} common output times, "
                    f"{len(big)} species above {COMPARE_FLOOR:g} molec/cm3)",
                    f"  within 2 %: {sum(1 for x in big if x[0] <= 0.02)}/{len(big)}   "
                    f"within 0.1 %: {sum(1 for x in big if x[0] <= 1e-3)}/{len(big)}",
                    "  O3/NO/NO2/OH/HO2: " + ", ".join(f"{k_} {v*100:.3g} %" for k_, v in main.items()),
                    "  worst: " + ", ".join(f"{s_} {v*100:.3g} %" for v, s_, _ in big[:5])]
    return "\n".join(out) if out else f"{stem}: no {stem}_F0AM.csv yet (run run_{stem}.m in MATLAB)"


# ─────────────────────────────────────────────────────────────────────────────
# Scenario preparation (benchmark scenes and the UV-lamp case)
# ─────────────────────────────────────────────────────────────────────────────

def prepare(scene: str, work: str, profile: str = "precise") -> str:
    import scenario_prep as sp
    if scene == "uv_lights_on_off":
        import check_lamp_chamber as cl
        import experiment_spec as es
        spec = cl.lamp_spec()
        with contextlib.redirect_stdout(io.StringIO()):
            res = es.write_experiment(spec, os.path.join(work, scene))
            return sp.prepare_scenario(res, sp.load_paths(TOML)["include_dir"], toml_path=TOML,
                                       solver=sp.solver_from_toml(TOML, profile))
    import run_benchmark as rb
    rb.PROFILE = profile
    with contextlib.redirect_stdout(io.StringIO()):
        return rb.prepare(scene, work, kv_pd=False)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description="Write JUICE scenarios as F0AM runs with identical inputs.")
    ap.add_argument("scenes", nargs="+")
    ap.add_argument("--out", default=os.path.join(HERE, "paper", "F0AM_scenarios"))
    ap.add_argument("--profile", default="precise", choices=("fast", "precise"))
    ap.add_argument("--check", action="store_true", help="integrate the written files in Python and compare")
    ap.add_argument("--compare", action="store_true",
                    help="only compare existing F0AM results (<stem>_F0AM.csv) with JUICE and FACSIMILE")
    ap.add_argument("--scripts-only", action="store_true",
                    help="only rewrite run_<stem>.m and <stem>_JUICE_J.m of an existing export (no JUICE run)")
    a = ap.parse_args(argv)
    if a.compare:
        for scene in a.scenes:
            print(compare_outputs(os.path.join(a.out, scene), scene))
        return 0
    if a.scripts_only:
        for scene in a.scenes:
            d = os.path.join(a.out, scene)
            if not os.path.isfile(os.path.join(d, f"{scene}_JUICE_inputs.mat")):
                print(f"{scene}: no export in {d}")
                return 1
            write_scripts(d, scene)
            print(f"{scene}: run_{scene}.m and {scene}_JUICE_J.m rewritten in {d}")
        return 0
    work = os.path.join(a.out, "_work")
    rc = 0
    for scene in a.scenes:
        try:
            js = prepare(scene, work, a.profile)
            rec = record(js, scene)
            d = write(rec, a.out)
            print(f"{scene}: written to {d} ({len(rec.species)} species, {len(rec.rxn_names)} reactions, "
                  f"{rec.k.shape[0]} steps)")
            if a.check:
                print(check(d, scene, ro2_frozen=True))
                print(check(d, scene))
        except Exception as e:
            rc = 1
            print(f"{scene}: FAILED - {type(e).__name__}: {e}")
    return rc


if __name__ == "__main__":
    sys.exit(main())
