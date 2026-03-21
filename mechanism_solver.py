import numpy as np
import sympy as sp

class MechanismSolver:
    def __init__(self, mechanism, env_map):
        self.mechanism = mechanism
        self.species = sorted(list(set(mechanism.species)))
        self.species_to_idx = {s: i for i, s in enumerate(self.species)}
        self.num_species = len(self.species)
        self.env_map = env_map # From JSON

    def calculate_rhs(self, t, y, ctx_data, get_M):
        dydt = np.zeros(self.num_species)
        M = get_M(t)
        temp = ctx_data['TEMP'](t)
        
        # Build evaluation context for rate constants
        # In MCM, 'T' is used for Temperature and 'M' for air density
        eval_context = {"T": temp, "M": M, "exp": np.exp, "LOG10": np.log10}
        
        # Add J-values from context
        eval_context["JNO2"] = ctx_data['JNO2'](t)
        eval_context["JO1D"] = ctx_data['JO1D'](t)
        eval_context["HV"] = ctx_data['HV'](t)
        eval_context["RH"] = self.calculate_rh(t, ctx_data, M, temp)

        for rxn in self.mechanism.reactions:
            # Calculate Rate Constant k
            try:
                k = eval(rxn.rate_str, {"__builtins__": None}, eval_context)
            except:
                k = 0.0 # Safety for malformed rates

            # Reaction Velocity: v = k * [A] * [B]
            v = k
            for r in rxn.reactants:
                if r in self.species_to_idx:
                    v *= y[self.species_to_idx[r]]
            
            # Update Derivatives
            for r in rxn.reactants:
                if r in self.species_to_idx: dydt[self.species_to_idx[r]] -= v
            for p in rxn.products:
                if p in self.species_to_idx: dydt[self.species_to_idx[p]] += v

        return dydt

    def calculate_rh(self, t, ctx_data, M, temp):
        # WVC2RH conversion logic from .easy
        h2o_mixing = ctx_data['H2O'](t) # molec/cm3
        # Simple RH approximation for the console/artifact logic
        p_sat = np.exp(21.36469 - 5339.66 / temp)
        rh = (h2o_mixing * 100) / (p_sat / (1.38e-16 * temp) * M)
        return np.clip(rh, 0, 100)