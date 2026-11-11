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
from juice.core.easy_parser import expand_macros, _split_top_level_args


# ── Math functions available in all eval contexts ────────────────────────────
_MATH_CTX = {
    'EXP': np.exp,   'LOG': np.log,   'LOG10': np.log10,
    'ABS': np.abs,   'SQRT': np.sqrt,
    'AMAX': np.maximum, 'AMIN': np.minimum,
    'COS': np.cos,   'SIN': np.sin,   'TAN': np.tan,
    'np': np,
}


# Fortran/FACSIMILE double exponents (1.0D-31 -> 1.0E-31), only inside NUMERIC LITERALS: the number must
# start at a token boundary, so a '<digit>D<digit>' inside an identifier (GECKO-A species G2D5004) is left
# alone.  From F. Rohrer's branch (2026-11-08); the old lookbehind regex renamed such species.
_D_EXP = re.compile(r'(?<![A-Za-z_0-9.])(\d+\.?\d*|\.\d+)[dD]([+-]?\d+)')


def fortran_d_to_e(expr: str) -> str:
    return _D_EXP.sub(r'\1E\2', expr)


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
        expr = fortran_d_to_e(expr)
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
    expanded = fortran_d_to_e(expanded) # Fortran D→E

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
    expanded = fortran_d_to_e(expanded)
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


def precompile_s_vec(s_entries, function_entries):
    """Pre-compile each s entry's definition to a code object.

    Call once after topological_sort. The hot loop then uses
    eval(entry._code, ctx) directly, avoiding string parsing overhead.
    INPUT entries get _code=None; _input_col holds the ENZ column name.
    """
    func_dict = {f.name: f.definition for f in function_entries}
    for entry in s_entries:
        definition = entry.definition.strip()
        m = re.match(r'^(\w+)\s*\((.*)\)\s*$', definition, re.DOTALL | re.IGNORECASE)
        if m and m.group(1).upper() == 'INPUT':
            entry._code      = None
            entry._input_col = m.group(2).strip().upper()
            continue
        # Expand template + nested macros, then compile
        expr = definition
        if m:
            fname    = m.group(1).upper()
            args_str = m.group(2)
            template = func_dict.get(fname, '%1')
            args     = _split_top_level_args(args_str)
            expr     = template
            for i, arg in enumerate(args, 1):
                expr = re.sub(rf'%{i}\b', f'({arg.strip()})', expr)
            expr = expand_macros(expr, func_dict)
        expr = expr.upper().replace('@', '**')
        expr = fortran_d_to_e(expr)
        try:
            entry._code      = compile(expr.strip(), '<string>', 'eval')
            entry._expr_str  = expr.strip()   # kept for batch-function generation
            entry._input_col = None
        except Exception:
            entry._code      = None
            entry._expr_str  = None
            entry._input_col = None


def update_s_vec(s_entries, function_entries, c_vec, enz_vec, _ordered=None):
    """Evaluate all s entries in topological order and update their .value.

    Requires precompile_s_vec() to have been called first (entry._code set).
    The hot loop uses eval(code, ctx) — same speed as the old var_assignments.
    """
    ordered = _ordered if _ordered is not None else topological_sort(s_entries, function_entries)
    ctx = build_ctx(s_entries, c_vec, enz_vec)
    for entry in ordered:
        try:
            if getattr(entry, '_input_col', None) is not None:
                for e in enz_vec:
                    if e.name.upper() == entry._input_col:
                        entry.value = e.value
                        ctx[entry.name] = e.value
                        break
            elif getattr(entry, '_code', None) is not None:
                entry.value = float(eval(entry._code, {"__builtins__": {}}, ctx))
                ctx[entry.name] = entry.value
            else:
                entry.value = evaluate(
                    entry.definition, function_entries, s_entries, c_vec, enz_vec)
                ctx[entry.name] = entry.value
        except Exception as ex:
            import warnings
            warnings.warn(f"update_s_vec: '{entry.name}' ({entry.definition[:50]}): {ex}")


def update_enz_vec(enz_vec, providers, t):
    """Update enz_vec values from ENZ interpolators at time t."""
    for e in enz_vec:
        fn = providers.get(e.name)
        if fn is not None:
            try:
                e.value = float(fn(t))
            except Exception:
                e.value = 0.0
