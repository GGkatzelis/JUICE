"""
check_lamp_chamber.py — the lamp-chamber counterpart of the SAPHIR golden tests

SAPHIR runs are checked against FACSIMILE (run_benchmark.py) and the IDL-written EASY files
(experiment_spec.py --selftest).  A lamp chamber has no FACSIMILE reference, so this script pins
one fixed experiment in a generic 10 m³ UV-lamp chamber and checks two things every time:

  1. physics that must hold whatever the mechanism does —
       dark start:      every J = 0, no O3, NO2 falls at exactly the dilution rate
       lights on, NOx:  extended photostationary state  jNO2·[NO2] = [NO]·(k[O3] + k[HO2] + k[RO2])  (ratio 1.00)
                        H2O2 is lost at J3 + k[OH]
                        J4 = jNO2 = jNO2(full power) × lights, J<n> = jNO2 × ratio of the lamp spectrum
       + ethanol:       OH appears (1e5–3e7 cm⁻³), ethanol is lost faster than by dilution, Ox is produced
       lights off:      every J = 0 again, OH collapses, Ox is no longer produced, N2O5 builds up
  2. regression against the pinned result of the first accepted run (reference .cs, EASY, ENZ)
     in validation/UV lamps/uv_lights_on_off/ — every species within ±2 % (median ratio).

    python -m juice.tools.check_lamp_chamber            # regenerate, run (~1 min), check physics, compare with the reference
    python -m juice.tools.check_lamp_chamber --update   # accept the current result as the new reference (after a deliberate change)

Timeline (UTC, 2024-07-09): 08:00 start with 20 ppb NO2 and 1 ppm H2O2 (the usual OH source of UV-chamber
experiments) in the dark · 09:00 lights on · 10:30 200 ppb ethanol injected (10 min) · 13:00 lights off · 15:00 end.  298 K, 1013 hPa, 1 % H2O, 0.5 m³/h flow
(dilution 5 %/h).  Chamber profile generic_lamps: jNO2 = 2e-3 1/s at full power, synthetic 350 nm
blacklight spectrum — change the profile and the reference changes with it (then run --update).
"""
from __future__ import annotations
from juice._root import ROOT as _JUICE_ROOT

import contextlib
import io
import math
import os
import shutil
import sys
import tempfile

import numpy as np

HERE = _JUICE_ROOT
NAME = "uv_lights_on_off"
REF_DIR = os.path.join(HERE, "validation", "UV lamps", NAME)
CHAMBER = "generic_lamps"


def lamp_spec():
    from juice.core import experiment_spec as es
    s = es.ExperimentSpec(name=NAME, title="UV-lamp chamber check: NO2 dark -> lights on -> ethanol -> lights off")
    es.apply_chamber(s, CHAMBER)
    c = s.conditions
    c.source, c.date, c.start, c.end, c.step_s = "manual", "2024-07-09", "08:00", "15:00", 60
    c.temperature_K, c.temperature_mode, c.pressure_hPa, c.h2o_percent, c.flow_m3h = 298.0, "const", 1013.0, 1.0, 0.5
    c.roof_open = [("09:00", "13:00")]
    s.precursors = ["C2H5OH"]
    s.initials = [("NO2", "20e-9*M"), ("H2O2", "1000e-9*M")]
    s.injections = [es.InjectionEvent(name="ethanol", products="C2H5OH", total_ppb="200", start="10:30", duration_min=10,
                                      comment="200 ppb ethanol under the lamps")]
    return s


def run(spec, out_dir: str) -> dict:
    from juice.core import experiment_spec as es
    from juice.core.scenario_prep import prepare_scenario, load_paths
    res = es.write_experiment(spec, out_dir)
    js = prepare_scenario(res, load_paths("")["include_dir"])
    from juice.core import scenario_runner as sr
    with contextlib.redirect_stdout(io.StringIO()):
        r = sr.ScenarioRunner(js); r._setup_environment(); r._build_mechanism(); r._solve(); r._export_results()
    res["cs"] = os.path.join(out_dir, f"{spec.name}.cs")
    return res


def physics_checks(cs_path: str, spec) -> list:
    """Returns the list of failures (empty = all good); prints one line per check."""
    from juice.core import chambers
    from juice.chemistry import lamp_photolysis as lp
    from juice.tools.run_benchmark import parse_cs_file
    fails = []

    def ok(cond, text):
        print(("OK   " if cond else "FAIL ") + text)
        if not cond:
            fails.append(text)

    d = parse_cs_file(cs_path)
    p = chambers.chamber(CHAMBER)
    jfull = float(p.light["jno2_full_s1"])
    lj = lp.lamp_j(p)
    t = (d["TIME"].values - d["TIME"].values[0]) / 3600.0 + 8.0          # hour of day
    M = d["M"].values
    hv = d["HV"].values
    ppb = lambda sp: d[sp].values / M * 1e9
    lit = hv > 0.5
    T = spec.conditions.temperature_K
    k_dil = (spec.conditions.flow_m3h + p.leak_flow_m3h) / (p.volume_m3 * 3600.0)
    k_no_o3 = 1.4e-12 * math.exp(-1310.0 / T)
    win = lambda a, b: (t >= a) & (t < b)

    # 1) the light: J columns
    ok(np.allclose(d["J4"].values, jfull * hv, rtol=1e-6, atol=1e-12) and np.allclose(d["JNO2"].values, jfull * hv, rtol=1e-6, atol=1e-12),
       f"J4 = jNO2 = {jfull:g} 1/s x lights (0 when dark)")
    for j in (1, 3, 7, 11, 12, 13):
        col = f"J{j}"
        if col in d.columns:
            ok(np.allclose(d[col].values, lj.ratio(j) * jfull * hv, rtol=1e-3, atol=1e-15),
               f"J{j} = jNO2 x {lj.ratio(j):.3g} ({lj.method}); 0 when dark")

    # 2) dark start 08:00–09:00
    w = win(8.1, 8.95)
    ok(ppb("O3")[w].max() < 0.05, f"dark start: O3 stays < 0.05 ppb (max {ppb('O3')[w].max():.3g})")
    ok(d["OH"].values[w].max() < 1e4, f"dark start: OH < 1e4 cm-3 (max {d['OH'].values[w].max():.2g})")
    i0, i1 = np.argmax(t >= 8.1), np.argmax(t >= 8.9)
    k_obs = -math.log(d["NO2"].values[i1] / d["NO2"].values[i0]) / ((t[i1] - t[i0]) * 3600.0)
    ok(abs(k_obs / k_dil - 1) < 0.03, f"dark start: NO2 falls at the dilution rate ({k_obs:.3e} vs {k_dil:.3e} 1/s)")

    # 3) lights on, NOx only 09:20–10:25: Leighton ratio
    k_ho2_no = 3.45e-12 * math.exp(270.0 / T); k_ro2_no = 2.3e-12 * math.exp(360.0 / T)      # MCM: HO2+NO, CH3O2+NO (RO2 sum)
    def leighton(w):
        conv = k_no_o3 * d["O3"].values[w] + k_ho2_no * d["HO2"].values[w] + (k_ro2_no * d["RO2"].values[w] if "RO2" in d.columns else 0)
        return (d["J4"].values[w] * d["NO2"].values[w]) / (d["NO"].values[w] * conv)
    w = win(9.33, 10.42); phi = leighton(w)
    ok(0.97 < np.median(phi) < 1.03 and phi.min() > 0.95 and phi.max() < 1.05,
       f"extended photostationary state (NOx + H2O2): jNO2[NO2]/([NO](k[O3]+k[HO2]+k[RO2])) median {np.median(phi):.3f} ({phi.min():.3f}-{phi.max():.3f})")
    w2 = win(10.75, 12.92); phi2 = leighton(w2)
    ok(0.95 < np.median(phi2) < 1.05, f"extended photostationary state with ethanol: median {np.median(phi2):.3f} ({phi2.min():.3f}-{phi2.max():.3f})")
    i0, i1 = np.argmax(t >= 9.33), np.argmax(t >= 12.9)
    k_h2o2 = -math.log(d["H2O2"].values[i1] / d["H2O2"].values[i0]) / ((t[i1] - t[i0]) * 3600.0) - k_dil
    j3 = lj.ratio(3) * jfull
    oh_w = np.median(d["OH"].values[win(9.33, 12.92)])
    k_exp = j3 + 2.9e-12 * math.exp(-160.0 / T) * oh_w                     # photolysis + OH; HO2+HO2 re-forms a little
    ok(0.7 * k_exp < k_h2o2 < 1.5 * k_exp, f"H2O2 lost under the lamps at {k_h2o2:.2e} 1/s (J3 {j3:.2e} + k[OH] {k_exp - j3:.2e} = {k_exp:.2e})")
    ok(ppb("O3")[w].min() > 1.0, f"O3 formed under the lamps ({ppb('O3')[w].min():.1f}-{ppb('O3')[w].max():.1f} ppb)")

    # 4) ethanol under the lamps 10:45–12:55
    w = win(10.75, 12.92)
    oh = d["OH"].values[w]
    ok(5e5 < np.median(oh) < 1e8, f"ethanol + lamps: OH median {np.median(oh):.2e} cm-3 (5e5-1e8 expected with 1 ppm H2O2)")
    i0, i1 = np.argmax(t >= 10.75), np.argmax(t >= 12.9)
    k_etoh = -math.log(d["C2H5OH"].values[i1] / d["C2H5OH"].values[i0]) / ((t[i1] - t[i0]) * 3600.0)
    ok(k_etoh > 1.2 * k_dil, f"ethanol lost faster than by dilution ({k_etoh:.2e} vs dilution {k_dil:.2e} 1/s -> OH ~ {(k_etoh - k_dil) / 3.2e-12:.1e} cm-3)")
    ox = ppb("O3") + ppb("NO2")
    ox_corr = ox * np.exp(k_dil * (t - 8.0) * 3600.0)                     # dilution-corrected
    ok(ox_corr[i1] > ox_corr[i0] * 1.02, f"Ox (O3+NO2, dilution-corrected) produced: {ox_corr[i0]:.1f} -> {ox_corr[i1]:.1f} ppb")

    # 5) lights off 13:10–14:55
    w = win(13.17, 14.92)
    ok(all(d[f"J{j}"].values[w].max() == 0 for j in (1, 3, 4, 7, 11, 12) if f"J{j}" in d.columns), "lights off: every J = 0")
    ok(d["OH"].values[w].max() < 1e5, f"lights off: OH collapses (max {d['OH'].values[w].max():.2g} cm-3)")
    i0, i1 = np.argmax(t >= 13.17), np.argmax(t >= 14.9)
    ok(ox_corr[i1] <= ox_corr[i0] * 1.01, f"lights off: no Ox production ({ox_corr[i0]:.1f} -> {ox_corr[i1]:.1f} ppb, dilution-corrected)")
    i_off = np.argmax(t >= 13.0); i_30 = np.argmax(t >= 13.5)
    if ppb("NO")[i_off] > ppb("O3")[i_off]:                                  # NO in excess: it titrates the ozone
        ok(ppb("O3")[i_30] < 0.05 * ppb("O3")[i_off], f"lights off: excess NO titrates O3 ({ppb('O3')[i_off]:.1f} -> {ppb('O3')[i_30]:.2f} ppb in 30 min)")
    else:                                                                     # O3 in excess: NO3/N2O5 night chemistry
        ok(ppb("N2O5")[i1] > 0.01, f"lights off: N2O5 builds up from O3 + NO2 ({ppb('N2O5')[i1] * 1000:.0f} ppt at 14:54)")
    return fails


def compare_reference(cs_path: str, ref_cs: str, tol: float = 0.02) -> list:
    """Median ratio run/reference per species (> 1 ppt in both) must be within ±tol."""
    from juice.tools.run_benchmark import parse_cs_file
    a, b = parse_cs_file(cs_path), parse_cs_file(ref_cs)
    bad, n = [], 0
    for sp in b.columns:
        if sp in ("TIME",) or sp not in a.columns or not sp[0].isalpha() or sp.startswith(("J", "KMT")):
            continue
        x, y = a[sp].values, b[sp].values
        if len(x) != len(y):
            return [f"different length ({len(x)} vs {len(y)} rows)"]
        thr = 1e-12 * b["M"].values if sp not in ("M", "T", "HV", "DILUTE", "P", "H2O") else 0
        mask = (np.abs(x) > thr) & (np.abs(y) > thr)
        if mask.sum() < 10:
            continue
        n += 1
        r = np.median(x[mask] / y[mask])
        if abs(r - 1) > tol:
            bad.append(f"{sp}: median ratio {r:.4f}")
    print(f"{'OK  ' if not bad else 'FAIL'} reference comparison: {n} species, {len(bad)} outside ±{tol:.0%}" + (": " + ", ".join(bad[:8]) if bad else ""))
    return bad


def main(argv) -> int:
    update = "--update" in argv
    from juice.core import experiment_spec as es
    from juice.core.easy_scene import statements_equivalent
    spec = lamp_spec()
    tmp = tempfile.mkdtemp(prefix="lamp_check_")
    print(f"UV-lamp chamber check: {CHAMBER}, work dir {tmp}")
    res = run(spec, os.path.join(tmp, NAME))
    fails = physics_checks(res["cs"], spec)
    ref_easy, ref_enz, ref_cs = (os.path.join(REF_DIR, NAME + ext) for ext in (".easy", ".enz", ".reference.cs"))
    if update:
        if fails:
            print(f"NOT updated: {len(fails)} physics check(s) failed"); return 1
        os.makedirs(REF_DIR, exist_ok=True)
        es.save_spec(spec, os.path.join(REF_DIR, NAME + ".experiment.json"))
        shutil.copy2(res["easy"], ref_easy); shutil.copy2(res["enz"], ref_enz); shutil.copy2(res["cs"], ref_cs)
        print(f"reference written to {REF_DIR}"); return 0
    if not os.path.isfile(ref_cs):
        print("no reference yet — run with --update"); return 1 if fails else 0
    same, diffs = statements_equivalent(ref_easy, res["easy"])
    print(f"{'OK  ' if same else 'FAIL'} EASY statement-equivalent to the reference" + ("" if same else f": {diffs[:3]}"))
    enz_same = open(res["enz"], "rb").read() == open(ref_enz, "rb").read()
    print(f"{'OK  ' if enz_same else 'FAIL'} ENZ identical to the reference")
    bad = compare_reference(res["cs"], ref_cs)
    n_fail = len(fails) + len(bad) + (0 if same else 1) + (0 if enz_same else 1)
    print(f"{'LAMP CHECK PASSED' if n_fail == 0 else f'LAMP CHECK FAILED ({n_fail})'}")
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
