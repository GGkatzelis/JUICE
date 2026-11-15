#!/usr/bin/env python
"""live_forms.py -- the two pure expression analysers of the live-settings block of scenario_runner._solve
(moved out unchanged 2026-11-15; no numpy, no JUICE import).

    monomial_form(src, live, other)   rate expression = F * product(LIVE**e), F without live settings
                                      -> {LIVE: e}; else None
    affine_form(src)                  setting expression = c0 + sum(coef * NAME) with numerical coefficients
                                      -> ({NAME: coef}, c0); else None

scenario_runner uses monomial_form to evaluate the rate constants that depend on live settings (KCOND*C/COA,
K*RO2 ...) as array statements, and affine_form to evaluate settings that are linear in the species (RO2, COA,
sums) as a dot product with the state vector.

    python live_forms.py --selftest
"""
import sys


def monomial_form(src, live, other):
    """expression = F * product(LIVE**e), F without live settings -> {LIVE: e}; else None.
    `other` = names that also change inside a row (species): any of them -> None."""
    import ast
    try:
        tree = ast.parse(src.strip(), mode='eval').body
    except Exception:
        return None

    def names(n):
        return {x.id for x in ast.walk(n) if isinstance(x, ast.Name)}

    def walk(n):
        nm = names(n)
        if nm & other:
            return None
        if any(x != x.upper() for x in nm if x.upper() in live):
            return None
        if not (nm & live):
            return {}
        if isinstance(n, ast.Name):
            return {n.id: 1.0}
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            return walk(n.operand)
        if isinstance(n, ast.BinOp):
            if isinstance(n.op, (ast.Mult, ast.Div)):
                a, b = walk(n.left), walk(n.right)
                if a is None or b is None:
                    return None
                f = 1.0 if isinstance(n.op, ast.Mult) else -1.0
                d = dict(a)
                for k, v in b.items():
                    d[k] = d.get(k, 0.0) + f * v
                return d
            if isinstance(n.op, ast.Pow):
                e = n.right
                if isinstance(e, ast.UnaryOp) and isinstance(e.op, ast.USub) and isinstance(e.operand, ast.Constant):
                    ev = -e.operand.value
                elif isinstance(e, ast.Constant):
                    ev = e.value
                else:
                    return None
                if not isinstance(ev, (int, float)) or isinstance(ev, bool):
                    return None
                a = walk(n.left)
                return None if a is None else {k: v * float(ev) for k, v in a.items()}
        return None
    try:
        d = walk(tree)
    except Exception:
        return None
    if not d:
        return None
    d = {k: v for k, v in d.items() if v != 0.0}
    return d or None


def affine_form(src):
    """expression = c0 + sum(coef * NAME) with numerical coefficients -> ({NAME: coef}, c0), else None"""
    import ast
    if not src:
        return None

    def scale(a, f):
        return {k: v * f for k, v in a[0].items()}, a[1] * f

    def walk(n):
        if isinstance(n, ast.Constant) and isinstance(n.value, (int, float)) and not isinstance(n.value, bool):
            return {}, float(n.value)
        if isinstance(n, ast.Name):
            return {n.id.upper(): 1.0}, 0.0
        if isinstance(n, ast.UnaryOp) and isinstance(n.op, (ast.USub, ast.UAdd)):
            a = walk(n.operand)
            return None if a is None else (scale(a, -1.0) if isinstance(n.op, ast.USub) else a)
        if isinstance(n, ast.BinOp):
            a = walk(n.left)
            b = walk(n.right) if a is not None else None
            if a is None or b is None:
                return None
            if isinstance(n.op, (ast.Add, ast.Sub)):
                f = 1.0 if isinstance(n.op, ast.Add) else -1.0
                d = a[0]
                for k, v in b[0].items():
                    d[k] = d.get(k, 0.0) + f * v
                return d, a[1] + f * b[1]
            if isinstance(n.op, ast.Mult):
                if not a[0]: return scale(b, a[1])
                if not b[0]: return scale(a, b[1])
            if isinstance(n.op, ast.Div) and not b[0] and b[1] != 0.0:
                return scale(a, 1.0 / b[1])
            if isinstance(n.op, ast.Pow) and not a[0] and not b[0]:
                return {}, a[1] ** b[1]
        return None
    _rl = sys.getrecursionlimit()
    try:
        sys.setrecursionlimit(max(_rl, 200000))
        return walk(ast.parse(src.strip(), mode='eval').body)
    except Exception:
        return None
    finally:
        sys.setrecursionlimit(_rl)


def _selftest():
    mf, af = monomial_form, affine_form
    none = set()
    # --- monomial_form: rate expression = (factor fixed in the row) x product of live settings ^ exponent
    assert mf("KCOND*4.2E+01/COA", {'COA'}, none) == {'COA': -1.0}
    assert mf(" KCOND*4.2E+01/COA ", {'COA'}, none) == {'COA': -1.0}                # blanks around the text
    assert mf("KCOND*4.2E+01/COA", {'COA', 'KCOND'}, none) == {'COA': -1.0, 'KCOND': 1.0}
    assert mf("2.7E-12*EXP(360/T)*RO2", {'RO2'}, none) == {'RO2': 1.0}
    assert mf("A*B", {'A', 'B'}, none) == {'A': 1.0, 'B': 1.0}                      # two live names
    assert mf("A*B/A", {'A', 'B'}, none) == {'B': 1.0}                              # exponent 0 is dropped
    assert mf("K*COA/COA", {'COA'}, none) is None                                   # nothing live is left
    assert mf("K*RO2**2", {'RO2'}, none) == {'RO2': 2.0}
    assert mf("K*COA**-1", {'COA'}, none) == {'COA': -1.0}
    assert mf("K/(COA*RO2)", {'COA', 'RO2'}, none) == {'COA': -1.0, 'RO2': -1.0}
    assert mf("-K*RO2", {'RO2'}, none) == {'RO2': 1.0}                              # the sign belongs to the factor
    assert mf("EXP(RO2)", {'RO2'}, none) is None                                    # not a monomial
    assert mf("K*(1+RO2)", {'RO2'}, none) is None                                   # a sum with a live name
    assert mf("K*RO2**N", {'RO2'}, none) is None                                    # exponent is not a number
    assert mf("K*RO2*NO", {'RO2'}, {'NO'}) is None                                  # a species in the expression
    assert mf("K*ro2", {'RO2'}, none) is None                                       # live name in another spelling
    assert mf("0.5+1.6606E-12*(280.2*SOA_A+136.2*SOA_B)", {'COA'}, none) is None    # no live name at all
    assert mf("2.7E-12*EXP(360/T)", {'RO2'}, none) is None
    assert mf("K*(RO2", {'RO2'}, none) is None                                      # not an expression
    # --- affine_form: setting = c0 + sum(coef * NAME)
    assert af("0.5+1.6606E-12*(280.2*SOA_A+136.2*SOA_B)") == (
        {'SOA_A': 280.2 * 1.6606E-12, 'SOA_B': 136.2 * 1.6606E-12}, 0.5)
    assert af("A+B") == ({'A': 1.0, 'B': 1.0}, 0.0)
    assert af("a+B") == ({'A': 1.0, 'B': 1.0}, 0.0)                                 # names are returned upper-case
    assert af("A+A") == ({'A': 2.0}, 0.0)
    assert af("(A+B)/2") == ({'A': 0.5, 'B': 0.5}, 0.0)
    assert af("1-A") == ({'A': -1.0}, 1.0)
    assert af("-A+1") == ({'A': -1.0}, 1.0)
    assert af("2**3+A") == ({'A': 1.0}, 8.0)
    assert af("3") == ({}, 3.0)
    assert af("A*B") is None                                                        # product of two names
    assert af("EXP(RO2)") is None                                                   # a function call
    assert af("KCOND*4.2E+01/COA") is None                                          # division by a name
    assert af("2.7E-12*EXP(360/T)*RO2") is None
    assert af("A/0") is None
    assert af("A+") is None                                                         # not an expression
    assert af("") is None and af(None) is None
    # a long sum (RO2, the SOA mass): the recursion limit is raised inside and restored.  (The size is limited
    # by the interpreter: Python 3.13 on Windows builds the syntax tree of a flat sum up to 2994 terms; a longer
    # one gives None = not analysed.)
    rl = sys.getrecursionlimit()
    n = 1000
    got = af("0.5+" + "+".join(f"{i + 1}*S{i}" for i in range(n)))
    assert got is not None and len(got[0]) == n and got[0]['S0'] == 1.0 and got[0][f'S{n - 1}'] == float(n) and got[1] == 0.5
    assert sys.getrecursionlimit() == rl
    print("live_forms selftest OK")


if __name__ == "__main__":
    if "--selftest" in sys.argv[1:]:
        _selftest()
    else:
        print(__doc__)
