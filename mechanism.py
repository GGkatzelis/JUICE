from dataclasses import dataclass, field
from typing import List, Dict

@dataclass
class Reaction:
    reactants: List[str]
    products: List[str]
    rate_expression: str

@dataclass
class Mechanism:
    species: List[str] = field(default_factory=list)
    parameters: Dict[str, str] = field(default_factory=dict)
    reactions: List[Reaction] = field(default_factory=list)
    initial_conditions: Dict[str, str] = field(default_factory=dict)
    input_vars: List[str] = field(default_factory=list)
    functions: Dict[str, str] = field(default_factory=dict) 
    aliases: Dict[str, str] = field(default_factory=dict)