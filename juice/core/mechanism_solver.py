import numpy as np
from juice.core import mcm_kinetics # Ensure this is imported for KMT/Troe functions

class MechanismSolver:
    def __init__(self, mechanism, env_map):
        self.mechanism = mechanism
        self.species = sorted(list(set(mechanism.species)))
        self.species_to_idx = {s: i for i, s in enumerate(self.species)}
        self.num_species = len(self.species)
        self.env_map = env_map 

    def calculate_rhs(self, t, y, ctx_data, get_M):
        dydt = np.zeros(self.num_species)
        M_val = get_M(t)
        temp = ctx_data['TEMP'](t)
        
        # Build evaluation context with MCM kinetic functions
        eval_context = {
            "T": temp, "M": M_val, "exp": np.exp, "LOG10": np.log10,
            "KMT06": mcm_kinetics.KMT06, "KMT07": mcm_kinetics.KMT07,
            "KRO2": mcm_kinetics.KRO2, "JNO2": ctx_data['JNO2'](t),
            "JO1D": ctx_data['JO1D'](t), "HV": ctx_data['HV'](t),
            "RH": self.calculate_rh(t, ctx_data, M_val, temp)
        }

        for rxn in self.mechanism.reactions:
            try:
                # Evaluate rate string within the scientific context
                k = eval(rxn.rate_expression, {"__builtins__": None}, eval_context)
                
                # Reaction Velocity: v = k * [Reactant1] * [Reactant2]...
                v = k
                for r in rxn.reactants:
                    # Handle stoichiometric coefficients (e.g., 2*OH)
                    coeff = 1.0
                    name = r
                    if '*' in r:
                        coeff_str, name = r.split('*')
                        coeff = float(coeff_str)
                    
                    if name in self.species_to_idx:
                        v *= (y[self.species_to_idx[name]] ** coeff)

                # Distribute production and loss
                for r in rxn.reactants:
                    name = r.split('*')[-1]
                    if name in self.species_to_idx: dydt[self.species_to_idx[name]] -= v
                for p in rxn.products:
                    name = p.split('*')[-1]
                    if name in self.species_to_idx: dydt[self.species_to_idx[name]] += v
            except:
                continue 

        return dydt

    def calculate_rh(self, t, ctx_data, M, temp):
        """Standard SAPHIR WVC2RH conversion"""
        h2o_mixing = ctx_data['H2O'](t) 
        p_sat = np.exp(21.36469 - 5339.66 / temp)
        rh = (h2o_mixing * 100) / (p_sat / (1.38e-16 * temp) * M)
        return np.clip(rh, 0, 100)