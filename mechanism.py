from dataclasses import dataclass, field
from typing import List, Dict


@dataclass
class Reaction:
    reactants: List[str]
    products: List[str]
    rate_expression: str


# ── New structured entries (Phase 1) ─────────────────────────────────────────
# These are populated by easy_parser alongside the existing dict-based fields
# and will replace them in Phase 2.

@dataclass
class FunctionEntry:
    """A user-defined FUNCTION macro from the EASY script.
    Example: DENSITY = FUNCTION[(%1)/((%2)*1.379E-19)]
    """
    name: str        # "DENSITY"
    definition: str  # "(%1)/((%2)*1.379E-19)"


@dataclass
class KEntry:
    """Rate constant for one chemical reaction.
    Example: k[C2H5OH+OH-->C2H5O] = CONST(3.0D-12*EXP(20/(T))*0.05)
    """
    name: str        # "C2H5OH+OH-->C2H5O"
    definition: str  # "CONST(3.0D-12*EXP(20/(T))*0.05)"  — full call incl. wrapper
    reactants: List[str] = field(default_factory=list)
    products:  List[str] = field(default_factory=list)
    value: float = 0.0


@dataclass
class SEntry:
    """A setting/parameter defined by a function call (not integrated by ODE).
    Examples: T=INPUT(t), KMT06=CONST(f(T,H2O)), M=DENSITY(P,T), J4=PHTFKT(...)
    Also includes prescribed species: NO=CONST(M*1e-9*0.5*HV)
    """
    name: str        # "T", "KMT06", "M", "J4", "NO" (if prescribed)
    definition: str  # "INPUT(t)", "CONST(f(T,H2O))", "DENSITY(P,T)"
    value: float = 0.0
    is_synthetic: bool = False  # True = synthesized from reaction rate, not in EASY script


@dataclass
class CEntry:
    """A chemical species integrated by the ODE solver.
    Appears in reaction equations but has no function-call definition.
    """
    name: str        # "OH", "NO2", "CH3O2"
    value: float = 0.0


@dataclass
class EnzEntry:
    """One column from the ENZ file at the current time step.
    Name preserves original ENZ case (e.g. 't', 'hv', 'jNO2').
    """
    name: str        # original ENZ column name (case-sensitive)
    value: float = 0.0


@dataclass
class Mechanism:
    # ── Existing fields (Phase 1: kept for backward compatibility) ────────────
    species: List[str] = field(default_factory=list)
    parameters: Dict[str, str] = field(default_factory=dict)
    commented_parameters: Dict[str, str] = field(default_factory=dict)
    reactions: List[Reaction] = field(default_factory=list)
    initial_conditions: Dict[str, str] = field(default_factory=dict)
    input_vars: List[str] = field(default_factory=list)
    functions: Dict[str, str] = field(default_factory=dict)
    aliases: Dict[str, str] = field(default_factory=dict)

    # ── New structured entries (Phase 1) ─────────────────────────────────────
    function_entries: List[FunctionEntry] = field(default_factory=list)
    k_entries:        List[KEntry]        = field(default_factory=list)
    s_entries:        List[SEntry]        = field(default_factory=list)
