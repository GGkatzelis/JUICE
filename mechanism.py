from dataclasses import dataclass, field
from typing import List, Dict, Tuple

@dataclass(order=True)
class Species:
    """Represents a chemical species."""
    name: str

@dataclass
class Reaction:
    """Represents a chemical reaction."""
    reactants: Tuple[str, ...]
    products: Tuple[str, ...]
    rate_constant: str

@dataclass
class Mechanism:
    """Represents the entire chemical mechanism."""
    enz_files: Dict[str, str] = field(default_factory=dict)
    facsimile_params: Dict[str, str] = field(default_factory=dict)
    initial_conditions: Dict[str, str] = field(default_factory=dict)
    functions: Dict[str, str] = field(default_factory=dict)
    variables: Dict[str, str] = field(default_factory=dict)
    aliases: Dict[str, str] = field(default_factory=dict)
    reactions: List[Reaction] = field(default_factory=list)
    species: List[Species] = field(default_factory=list)

    def __post_init__(self):
        """Extract species from reactions after initialization."""
        self._extract_species()

    def _extract_species(self):
        """Extract unique species from reactions."""
        species_names = set()
        for reaction in self.reactions:
            for reactant in reaction.reactants:
                if reactant not in self.variables:
                    species_names.add(reactant)
            for product in reaction.products:
                if product not in self.variables:
                    species_names.add(product)
        
        # Also add species from initial conditions
        for species_name in self.initial_conditions.keys():
            if species_name != '*' and species_name not in self.variables:
                species_names.add(species_name)

        self.species = sorted([Species(name=name) for name in species_names])
