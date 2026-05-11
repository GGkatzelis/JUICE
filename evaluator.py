"""
evaluator.py — unified evaluation engine for k, s, c, enz, function vectors.

Public API:
  evaluate(call_string, function_entries, s_vec, c_vec, enz_vec)  → scalar
  build_ctx(s_vec, c_vec, enz_vec)                                → dict
  topological_sort(s_entries)                                     → ordered list
  update_s_vec(s_entries, function_entries, c_vec, enz_vec)       → in-place
"""

import re
import numpy as np
from easy_parser import expand_macros, _split_top_level_args


# ── Math functions available in all eval contexts ────────────────────────────
_MATH_CTX = {
    'EXP': np.exp,   'LOG': np.log,   'LOG10': np.log10,
    'ABS': np.abs,   'SQRT': np.sqrt,
    'AMAX': np.maximum, 'AMIN': np.minimum,
    'COS': np.cos,   'SIN': np.sin,   'TAN': np.tan,
    'np': np,
}


def build_ctx(s_vec, c_vec, enz_vec):
    """Build a flat evaluation context dict from all three vectors.

    Priority (higher overrides lower): s > c > enz.
    All names are stored uppercase.
    """
    ctx = dict(_MATH_CTX)
    for e in enz_vec:
        ctx[e.name.upper()] = e.value
    for e in c_vec:
        ctx[e.name.upper()] = e.value
    for e in s_vec:
        ctx[e.name.upper()] = e.value
    return ctx


def evaluate(call_string, function_entries, s_vec, c_vec, enz_vec):
    """Evaluate a definition string using the k/s/c/enz/function vectors.

    call_string must be a function call: "CONST(expr)", "INPUT(col)",
    "DENSITY(P,T)", "PHTFKT(l,m,n)", etc.

    Returns a scalar float (or numpy scalar for scalar inputs).
    """
    call_string = call_string.strip()

    # ── Parse: function_name(args_str) ───────────────────────────────────────
    m = re.match(r'^(\w+)\s*\((.*)\)\s*$', call_string, re.DOTALL | re.IGNORECASE)
    if not m:
        # Bare expression with no function wrapper — evaluate directly.
        ctx = build_ctx(s_vec, c_vec, enz_vec)
        expr = call_string.replace('@', '**')
        expr = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expr)
        return float(eval(expr, {"__builtins__": {}}, ctx))

    func_name = m.group(1).upper()
    args_str  = m.group(2)

    # ── Special case: INPUT ───────────────────────────────────────────────────
    # INPUT(col) reads the raw ENZ value. ENZ column names are case-sensitive,
    # but we match case-insensitively for robustness.
    if func_name == 'INPUT':
        col = args_str.strip()
        for e in enz_vec:
            if e.name.upper() == col.upper():
                return e.value
        raise KeyError(
            f"evaluate(): INPUT('{col}') — column not found in enz_vec.\n"
            f"  Available: {[e.name for e in enz_vec]}"
        )

    # ── Find function definition ──────────────────────────────────────────────
    func_def = None
    func_dict = {f.name: f.definition for f in function_entries}
    if func_name in func_dict:
        func_def = func_dict[func_name]
    if func_def is None:
        raise KeyError(
            f"evaluate(): function '{func_name}' not found in function_entries.\n"
            f"  Available: {list(func_dict.keys())}"
        )

    # ── Expand template: %1 → first arg, %2 → second arg, … ─────────────────
    args = _split_top_level_args(args_str)
    expanded = func_def
    for i, arg in enumerate(args, 1):
        expanded = re.sub(rf'%{i}\b', f'({arg.strip()})', expanded)

    # ── Recursively expand nested user-defined function calls ─────────────────
    expanded = expand_macros(expanded, func_dict)

    # ── Normalise and convert to Python syntax ────────────────────────────────
    expanded = expanded.upper()                                  # case-normalise (FUNCTION bodies may be mixed-case)
    expanded = expanded.replace('@', '**')                       # power operator
    expanded = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expanded) # Fortran D→E

    # ── Evaluate ─────────────────────────────────────────────────────────────
    ctx = build_ctx(s_vec, c_vec, enz_vec)
    try:
        return float(eval(expanded, {"__builtins__": {}}, ctx))
    except Exception as ex:
        raise RuntimeError(
            f"evaluate(): error evaluating '{func_name}({args_str})'.\n"
            f"  Expanded: {expanded}\n"
            f"  Error: {ex}"
        ) from ex


# ── Topological sort for s_entries ───────────────────────────────────────────

def _deps_of(definition, func_dict):
    """Return the set of names (uppercase) referenced in a definition string.

    After expanding all user-defined function templates, the remaining
    names in the expression are the dependencies on other s/c/enz entries.
    INPUT(col) contributes the col name as an ENZ dependency (not an s dep).
    """
    m = re.match(r'^(\w+)\s*\((.*)\)\s*$', definition.strip(),
                 re.DOTALL | re.IGNORECASE)
    if not m:
        # Bare expression
        try:
            code = compile(definition.replace('@', '**'), '<string>', 'eval')
            return {n.upper() for n in code.co_names}
        except Exception:
            return set()

    fname = m.group(1).upper()
    if fname == 'INPUT':
        return set()  # ENZ read — no s/c dependency

    args_str = m.group(2)
    # Expand template
    template = func_dict.get(fname, '%1')
    args = _split_top_level_args(args_str)
    expanded = template
    for i, arg in enumerate(args, 1):
        expanded = re.sub(rf'%{i}\b', f'({arg.strip()})', expanded)
    expanded = expand_macros(expanded, func_dict)
    expanded = expanded.upper().replace('@', '**')
    expanded = re.sub(r'(?<=\d)[dD](?=[+-]?\d)', 'E', expanded)
    try:
        code = compile(expanded, '<string>', 'eval')
        return {n for n in code.co_names if not n.startswith('__')}
    except Exception:
        return set()


def topological_sort(s_entries, function_entries):
    """Return s_entries in evaluation order (dependencies before dependents).

    Cycles (e.g. mutual references) are broken by processing remaining nodes
    in original order — typically not present in valid EASY scripts.
    """
    func_dict = {f.name: f.definition for f in function_entries}
    s_names   = {e.name.upper() for e in s_entries}
    deps      = {e.name.upper(): _deps_of(e.definition, func_dict) & s_names
                 for e in s_entries}

    ordered, visited = [], set()

    def visit(name):
        if name in visited:
            return
        visited.add(name)
        for dep in deps.get(name, set()):
            visit(dep)
        # Find the entry with this name (last-wins: take the last one)
        entry = next((e for e in reversed(s_entries)
                      if e.name.upper() == name), None)
        if entry is not None:
            ordered.append(entry)

    for e in s_entries:
        visit(e.name.upper())

    return ordered


def update_s_vec(s_entries, function_entries, c_vec, enz_vec):
    """Evaluate all s entries in topological order and update their .value.

    This replaces the 3-pass var_assignment loop.
    c_vec and enz_vec must already hold current values.
    """
    ordered = topological_sort(s_entries, function_entries)
    for entry in ordered:
        try:
            entry.value = evaluate(
                entry.definition, function_entries, s_entries, c_vec, enz_vec)
        except Exception as ex:
            # Keep previous value (0.0) — report for debugging
            import warnings
            warnings.warn(f"update_s_vec: error evaluating '{entry.name}' "
                          f"({entry.definition[:60]}): {ex}")


def update_enz_vec(enz_vec, providers, t):
    """Update enz_vec values from ENZ interpolators at time t."""
    for e in enz_vec:
        fn = providers.get(e.name)
        if fn is not None:
            try:
                e.value = float(fn(t))
            except Exception:
                e.value = 0.0
