"""Run all test scenarios in both interpolation modes and compare against FACSIMILE.
Saves results as simulation_results_<mode>.csv and prints comparison tables.

Usage: python run_benchmark.py [--skip-mcm-full]
"""
import os, sys, json, time, shutil, re
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

BASE = os.path.dirname(os.path.abspath(__file__))
SCENARIOS = [
    ("scene_01", os.path.join(BASE, "MCM_tested_scenarios", "IDL Input", "scene_01")),
    ("scene_02", os.path.join(BASE, "MCM_tested_scenarios", "IDL Input", "scene_02")),
    ("ethanol",  os.path.join(BASE, "MCM_tested_scenarios", "IDL Input", "ethanol")),
    ("mcm_full", os.path.join(BASE, "MCM_tested_scenarios", "IDL Input", "mcm_full")),
]
MODES = ["linear", "piecewise"]


def parse_cs_file(cs_path):
    """Parse FACSIMILE .cs output file."""
    with open(cs_path, encoding='latin-1') as f:
        lines = f.readlines()
    col_names = ['TIME']
    data_start = 0
    for i, line in enumerate(lines):
        m = re.match(r'\s*COLUMN\s+\d+=(\S+)', line)
        if m:
            name = m.group(1).split('<')[0]
            if name != 'TIME':
                col_names.append(name)
        if i > 5 and line.strip():
            try:
                vals = line.strip().split()
                if len(vals) > 5:
                    float(vals[0]); float(vals[1])
                    data_start = i
                    break
            except:
                pass
    rows = []
    for line in lines[data_start:]:
        line = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', line.strip())
        if not line:
            continue
        try:
            vals = [float(x) for x in line.split()]
            if len(vals) == len(col_names):
                rows.append(vals)
        except:
            pass
    return pd.DataFrame(rows, columns=col_names)


def compare_vs_facsimile(scenario_dir, sim_csv_path, ppb_threshold=0.001):
    """Compare a simulation CSV against the FACSIMILE .cs in the same dir."""
    cs_files = [f for f in os.listdir(scenario_dir)
                if f.endswith('.cs') and not f.endswith('.cs.ps')]
    if not cs_files:
        return None

    fac = parse_cs_file(os.path.join(scenario_dir, cs_files[0]))
    sim = pd.read_csv(sim_csv_path)
    sim_t = sim['TIME_s'].values
    fac_t = fac['TIME'].values
    common = [sp for sp in fac.columns if sp != 'TIME' and sp in sim.columns]

    results = {}
    for sp in common:
        fi = interp1d(fac_t, fac[sp].values, fill_value='extrapolate')(sim_t)
        mid = len(sim_t) // 2
        s, f, M = sim[sp].values[mid:], fi[mid:], sim['M'].values[mid:]
        s_ppb, f_ppb = s / M * 1e9, f / M * 1e9
        mask = (np.abs(f_ppb) > ppb_threshold) & (np.abs(s_ppb) > ppb_threshold)
        if mask.sum() < 3:
            continue
        pct = (np.median(s_ppb[mask] / f_ppb[mask]) - 1) * 100
        results[sp] = pct
    return results


def run_scenario(name, scenario_dir, mode):
    """Run a scenario with the given interpolation mode, save results."""
    config_path = os.path.join(scenario_dir, "scenario_config.json")
    with open(config_path, 'r') as f:
        cfg = json.load(f)

    # Set interpolation mode
    cfg['interpolation'] = mode

    # Write temporary config
    tmp_config = os.path.join(scenario_dir, f"_tmp_config_{mode}.json")
    with open(tmp_config, 'w') as f:
        json.dump(cfg, f, indent=4)

    print(f"\n{'='*65}")
    print(f" Running: {name} [{mode}]")
    print(f"{'='*65}")

    t0 = time.time()
    try:
        # Import fresh each time to avoid stale state
        import scenario_runner as sr
        runner = sr.ScenarioRunner(tmp_config)
        runner._setup_environment()
        runner._build_mechanism()
        runner._solve()
        runner._export_results()
        elapsed = time.time() - t0
        print(f">>> Completed in {elapsed//60:.0f}m {elapsed%60:.1f}s")
    except Exception as e:
        elapsed = time.time() - t0
        print(f">>> FAILED after {elapsed:.1f}s: {e}")
        import traceback
        traceback.print_exc()
    finally:
        # Clean up temp config
        if os.path.exists(tmp_config):
            os.remove(tmp_config)

    # Rename output to mode-specific filename
    src = os.path.join(scenario_dir, "simulation_results.csv")
    dst = os.path.join(scenario_dir, f"simulation_results_{mode}.csv")
    if os.path.exists(src):
        shutil.copy2(src, dst)
        print(f"  Saved: {dst}")
        return dst, elapsed
    return None, elapsed


def print_comparison_table(all_results):
    """Print formatted comparison of all scenarios and modes."""
    KEY_SPECIES = ['OH','O3','HO2','NO','NO2','HCHO','H2O2','HNO3',
                   'CH3CHO','PAN','CO','C2H5OH','CH3COCH3']

    print(f"\n\n{'='*80}")
    print(f" COMPREHENSIVE COMPARISON: JUICE vs FACSIMILE")
    print(f"{'='*80}\n")

    for name in ["scene_01", "scene_02", "ethanol", "mcm_full"]:
        if name not in all_results:
            continue
        data = all_results[name]
        print(f"\n--- {name} ---")

        # Print timing
        for mode in MODES:
            if mode in data and data[mode].get('time'):
                t = data[mode]['time']
                print(f"  {mode:12s}: {t//60:.0f}m {t%60:.1f}s")

        # Print accuracy summary
        for mode in MODES:
            if mode in data and data[mode].get('results'):
                r = data[mode]['results']
                n2 = sum(1 for v in r.values() if abs(v) <= 2.0)
                n5 = sum(1 for v in r.values() if abs(v) <= 5.0)
                print(f"  {mode:12s}: {n2}/{len(r)} within +/-2%, "
                      f"{n5}/{len(r)} within +/-5%")

        # Print key species side-by-side
        has_linear = 'linear' in data and data['linear'].get('results')
        has_pw = 'piecewise' in data and data['piecewise'].get('results')
        if has_linear or has_pw:
            print(f"\n  {'Species':14s} {'Linear':>10s}  {'Piecewise':>10s}")
            print(f"  {'':14s} {'':>10s}  {'':>10s}")
            for sp in KEY_SPECIES:
                l_val = data.get('linear', {}).get('results', {}).get(sp)
                p_val = data.get('piecewise', {}).get('results', {}).get(sp)
                l_str = f"{l_val:+.2f}%" if l_val is not None else "  ---"
                p_str = f"{p_val:+.2f}%" if p_val is not None else "  ---"
                if l_val is not None or p_val is not None:
                    print(f"  {sp:14s} {l_str:>10s}  {p_str:>10s}")


def main():
    skip_mcm = '--skip-mcm-full' in sys.argv

    all_results = {}

    for name, scenario_dir in SCENARIOS:
        if skip_mcm and name == 'mcm_full':
            print(f"\n>>> Skipping mcm_full (--skip-mcm-full)")
            continue

        all_results[name] = {}

        for mode in MODES:
            csv_path, elapsed = run_scenario(name, scenario_dir, mode)
            all_results[name][mode] = {'time': elapsed, 'csv': csv_path}

            if csv_path:
                results = compare_vs_facsimile(scenario_dir, csv_path)
                if results:
                    all_results[name][mode]['results'] = results

    print_comparison_table(all_results)


if __name__ == '__main__':
    main()
