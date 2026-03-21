import os, json, re, numpy as np, pandas as pd
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from mechanism_solver import MechanismSolver
from easy_parser import EasyParser
from mechanism import Mechanism, Reaction

def run_simulation(config_path):
    with open(config_path, 'r') as f:
        cfg = json.load(f)
    
    # 1. Setup Data Interpolators
    root = cfg['folder']
    enz_file = next(f for f in os.listdir(root) if f.endswith('.enz'))
    data = pd.read_csv(os.path.join(root, enz_file), skiprows=35, header=None)
    t_span = data.iloc[:, 0].values
    
    # Header mapping
    with open(os.path.join(root, enz_file), 'r', encoding='latin-1') as f:
        h_lines = [next(f) for _ in range(35)]
    headers = [re.findall(r'COLUMN\s*\d+=(.*)', l)[0].split('/')[0].strip().lower() for l in h_lines if "COLUMN" in l]

    ctx_data = {}
    for var, m in cfg['env_map'].items():
        if m['type'] == "SAPHIR File":
            idx = headers.index(m['value'].lower())
            ctx_data[var] = interp1d(t_span, data.iloc[:, idx], fill_value="extrapolate")
        else:
            ctx_data[var] = lambda t, v=float(m['value']): v

    def get_M(t):
        return (ctx_data['PRESS'](t) * 100) / (1.3806e-23 * ctx_data['TEMP'](t)) * 1e-6

    # 2. Build Mechanism + Artifacts
    parser = EasyParser(base_path=os.path.join(root, "chemistry"))
    mech = Mechanism()
    for f in cfg['mechanism_files']:
        m_part = parser.parse(os.path.join(root, "chemistry", f))
        mech.species.extend(m_part.species)
        mech.reactions.extend(m_part.reactions)

    # Add SAPHIR Artifacts if Active
    if cfg['artifacts']['active']:
        scale = eval(cfg['artifacts']['scale']) # Handles "2.6 * 1.1 * 0.7"
        # HONO Source
        mech.reactions.append(Reaction([], ["HONO"], f"({cfg['artifacts']['hono_eqn']}) * {scale}"))
        # HCHO Source
        mech.reactions.append(Reaction([], ["HCHO"], f"{cfg['artifacts']['hcho_eqn']}"))
        # Dark Wall
        mech.reactions.append(Reaction(["H2O"], ["HONO"], cfg['artifacts']['dark_hono']))

    # 3. Solve
    solver = MechanismSolver(mech, cfg['env_map'])
    y0 = np.zeros(solver.num_species)
    for spec, ppb in cfg['initials'].items():
        if spec in solver.species_to_idx:
            y0[solver.species_to_idx[spec]] = (ppb * 1e-9) * get_M(t_span[0])

    def system_ode(t, y):
        dydt = solver.calculate_rhs(t, y, ctx_data, get_M)
        # Apply Global Dilution
        k_dil = (ctx_data['FLOW'](t) + 0.6) / (cfg['chamber_vol'] * 3600)
        return dydt - y * k_dil

    sol = solve_ivp(system_ode, [t_span[0], t_span[-1]], y0, method='BDF', t_eval=t_span)
    
    # 4. Export
    out = pd.DataFrame(sol.y.T, index=sol.t, columns=solver.species)
    out.to_csv(os.path.join(root, "model_output.csv"))
    print("✅ Simulation complete. Results saved.")

if __name__ == "__main__":
    run_simulation("path/to/your/scenario_config.json")