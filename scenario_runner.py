import os
import json
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from mechanism_solver import MechanismSolver
from easy_parser import EasyParser
from mechanism import Mechanism, Reaction
from load_saphir_environment import load_saphir_environment
import sys

class ScenarioRunner:
    def __init__(self, config_path):
        with open(config_path, 'r') as f:
            self.cfg = json.load(f)
        self.root = self.cfg['folder']
        self.mechanism = None
        self.solver = None
        self.t_span = None
        self.ctx_data = {}

    def run(self):
        self._setup_environment()
        self._build_mechanism()
        self._solve()
        self._export_results()

    def _setup_environment(self):
        enz_file = next(f for f in os.listdir(self.root) if 'kammer' in f.lower() and f.endswith('.enz'))
        providers, t_min, t_max, headers = load_saphir_environment(os.path.join(self.root, enz_file))
        self.t_span = np.linspace(t_min, t_max, num=1000)

        self.ctx_data.update(providers)

        for var, m in self.cfg['env_map'].items():
            if m['type'] != "SAPHIR File":
                self.ctx_data[var] = lambda t, v=float(m['value']): np.full_like(t, v) if isinstance(t, np.ndarray) else v


    def get_M(self, t):
        temp = self.ctx_data['TEMP'](t)
        press = self.ctx_data['PRESS'](t)
        return (press * 100) / (1.3806e-23 * temp) * 1e-6


    def _build_mechanism(self):
        parser = EasyParser(base_path=os.path.join(self.root, "chemistry"))
        self.mechanism = Mechanism()
        for f in self.cfg['mechanism_files']:
            m_part = parser.parse(f)
            self.mechanism.species.extend(m_part.species)
            self.mechanism.reactions.extend(m_part.reactions)
        self._add_artifacts()

    def _add_artifacts(self):
        if self.cfg['artifacts']['active']:
            scale = eval(self.cfg['artifacts']['scale'])
            if self.cfg['artifacts']['hono_active']:
                self.mechanism.reactions.append(Reaction([], ["HONO"], f"({self.cfg['artifacts']['hono_eqn']}) * {scale}"))
            if self.cfg['artifacts']['hcho_active']:
                self.mechanism.reactions.append(Reaction([], ["HCHO"], f"{self.cfg['artifacts']['hcho_eqn']}"))
            if self.cfg['artifacts']['no2_active']:
                self.mechanism.reactions.append(Reaction([], ["NO2"], f"{self.cfg['artifacts']['no2_eqn']}"))
            if self.cfg['artifacts']['y_active']:
                self.mechanism.species.append('Y')
                self.mechanism.reactions.append(Reaction(['Y', 'OH'], ["HO2"], f"{self.cfg['artifacts']['y_eqn']}"))
            if self.cfg['artifacts']['wall_active']:
                self.mechanism.reactions.append(Reaction(["H2O"], ["HONO"], self.cfg['artifacts']['dark_hono']))

    def _solve(self):
        # Custom Reactions must be added before solver initialization
        for r in self.cfg.get('custom_reactions', []):
            self.mechanism.reactions.append(Reaction(r['reactants'].split(' + '), r['products'].split(' + '), r['rate']))
            for s in r['reactants'].split(' + '):
                if s.strip() not in self.mechanism.species: self.mechanism.species.append(s.strip())
            for s in r['products'].split(' + '):
                if s.strip() not in self.mechanism.species: self.mechanism.species.append(s.strip())

        self.solver = MechanismSolver(self.mechanism, self.cfg['env_map'])
        y0 = np.zeros(self.solver.num_species)
        for spec, ppb in self.cfg['initials'].items():
            if spec in self.solver.species_to_idx:
                y0[self.solver.species_to_idx[spec]] = (ppb * 1e-9) * self.get_M(self.t_span[0])

        # Pre-compile injection logic
        cont_sources_eval = []
        for s in self.cfg.get('continuous_sources', []):
            if s['target'] in self.solver.species_to_idx:
                idx = self.solver.species_to_idx[s['target']]
                cont_sources_eval.append({'idx': idx, 'rate_str': s['rate']})

        burst_sources_eval = []
        for b in self.cfg.get('burst_injections', []):
            if b['target'] in self.solver.species_to_idx:
                idx = self.solver.species_to_idx[b['target']]
                start_time = pd.to_datetime(b['start_utc'], format='%d.%m.%Y %H:%M:%S').timestamp()
                end_time = start_time + float(b['duration_s'])
                rate_per_s = (float(b['ppb']) * 1e-9) * self.get_M(start_time) / float(b['duration_s'])
                burst_sources_eval.append({'idx': idx, 'start': start_time, 'end': end_time, 'rate': rate_per_s})

        def system_ode(t, y):
            dydt = self.solver.calculate_rhs(t, y, self.ctx_data, self.get_M)
            
            # Dilution
            k_dil = (self.ctx_data['FLOW'](t) + 0.6) / (self.cfg['chamber_vol'] * 3600)
            dydt -= y * k_dil

            # Continuous Sources
            M_val = self.get_M(t)
            for source in cont_sources_eval:
                rate = eval(source['rate_str'], {'M': M_val})
                dydt[source['idx']] += rate

            # Burst Injections
            for burst in burst_sources_eval:
                if burst['start'] <= t <= burst['end']:
                    dydt[burst['idx']] += burst['rate']
            
            return dydt

        self.sol = solve_ivp(system_ode, [self.t_span[0], self.t_span[-1]], y0, method='BDF', t_eval=self.t_span)

    def _export_results(self):
        out = pd.DataFrame(self.sol.y.T, index=self.sol.t, columns=self.solver.species)
        out.to_csv(os.path.join(self.root, "model_output.csv"))
        print("✅ Simulation complete. Results saved.")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        config_path = sys.argv[1]
    else:
        config_path = "C:/Users/g.gkatzelis/Desktop/My Folders/My Coding/Python/GitHub/MCM_SAPHIR_Tool/MCM_tested_scenarios/ethanol_scenario/scenario_config.json"
    runner = ScenarioRunner(config_path)
    runner.run()