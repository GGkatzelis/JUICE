"""mechanism_solver.py — the species list and index of a mechanism, as the runner uses them.

scenario_runner builds one MechanismSolver per run and reads `species`, `species_to_idx` and `num_species` from
it; the right-hand side, the Jacobian and the solver itself live in scenario_runner.py.  The class once also
evaluated the right-hand side in Python (calculate_rhs, with the KMT functions of mcm_kinetics); nothing called
that any more, and it was removed on 25 Nov 2026.
"""


class MechanismSolver:
    def __init__(self, mechanism, env_map):
        self.mechanism = mechanism
        self.species = sorted(list(set(mechanism.species)))
        self.species_to_idx = {s: i for i, s in enumerate(self.species)}
        self.num_species = len(self.species)
        self.env_map = env_map
