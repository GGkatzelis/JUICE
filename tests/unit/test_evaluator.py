"""evaluator.py, function by function: the evaluation context, one definition string, the order in which
settings are evaluated, and the update of settings and ENZ values.  Every input is small enough to check by hand."""
import warnings
from types import SimpleNamespace as NS

import pytest

from juice.core.easy_parser import expand_macros
from juice.core.evaluator import (build_ctx, build_ctx_fast, evaluate, precompile_s_vec, topological_sort,
                                  update_enz_vec, update_s_vec)


def entry(name, value=0.0, definition=""):
    return NS(name=name, value=value, definition=definition)


def function(name, definition):
    return NS(name=name, definition=definition)


FUNCTIONS = [
    function("CONST", "%1"),
    function("DOUBLE", "2*%1"),
    function("RATIO", "%1/%2"),
    function("QUAD", "DOUBLE(DOUBLE(%1))"),          # a function that calls another one
    function("ARRH", "%1*exp(-%2/TEMP)"),            # lower-case body, as hand-written EASY files have
]


# ── the context ──────────────────────────────────────────────────────────────

def test_context_priority_is_s_over_c_over_enz():
    ctx = build_ctx([entry("x", 3.0)], [entry("X", 2.0), entry("y", 20.0)], [entry("X", 1.0), entry("Y", 10.0),
                                                                             entry("z", 100.0)])
    assert (ctx["X"], ctx["Y"], ctx["Z"]) == (3.0, 20.0, 100.0)


def test_context_names_are_upper_case_and_math_is_there():
    ctx = build_ctx([entry("temp", 298.0)], [], [])
    assert "TEMP" in ctx and "temp" not in ctx
    assert ctx["EXP"](0.0) == 1.0 and ctx["AMAX"](2.0, 5.0) == 5.0 and ctx["AMIN"](2.0, 5.0) == 2.0


def test_fast_context_equals_the_plain_one():
    s, c, enz = [entry("A", 1.0)], [entry("B", 2.0), entry("A", 9.0)], [entry("C", 3.0), entry("B", 9.0)]
    fast = build_ctx_fast(s, ["B", "A"], [2.0, 9.0], enz)
    assert fast == build_ctx(s, c, enz)


def test_fast_context_extra_overrides_c_but_not_s():
    fast = build_ctx_fast([entry("A", 1.0)], ["A", "O1D"], [9.0, 9.0], [], extra={"O1D": 5.0, "A": 7.0})
    assert fast["O1D"] == 5.0 and fast["A"] == 1.0


# ── one definition string ────────────────────────────────────────────────────

def test_bare_expression_with_power_and_fortran_exponent():
    assert evaluate("2@3", [], [], [], []) == 8.0
    assert evaluate("1.0D-3*T", [], [entry("T", 300.0)], [], []) == pytest.approx(0.3)


def test_input_reads_the_enz_column_whatever_the_case():
    enz = [entry("jNO2", 0.005)]
    assert evaluate("INPUT(JNO2)", [], [], [], enz) == 0.005
    assert evaluate("input( jno2 )", [], [], [], enz) == 0.005


def test_input_of_a_missing_column_names_the_column():
    with pytest.raises(KeyError, match="NOPE"):
        evaluate("INPUT(NOPE)", [], [], [], [entry("T", 1.0)])


def test_function_arguments_are_substituted_as_a_whole():
    # RATIO(1+1, 2+2) must be (1+1)/(2+2) = 0.5, not 1+1/2+2 = 3.5
    assert evaluate("RATIO(1+1, 2+2)", FUNCTIONS, [], [], []) == 0.5


def test_argument_that_holds_a_comma_in_brackets_stays_one_argument():
    assert evaluate("RATIO(AMAX(1, 8), 2)", FUNCTIONS, [], [], []) == 4.0


def test_nested_functions_expand():
    assert evaluate("QUAD(3)", FUNCTIONS, [], [], []) == 12.0


def test_a_function_inside_its_own_argument():
    """Bug of 24 Nov 2026: DOUBLE(DOUBLE(3)) expanded to '(2*((2*(3))))(2*(3)))', the inner call written twice."""
    assert expand_macros("DOUBLE(DOUBLE(3))", {"DOUBLE": "2*%1"}) == "(2*((2*(3))))"
    assert evaluate("DOUBLE(DOUBLE(3))", FUNCTIONS, [], [], []) == 12.0
    assert evaluate("RATIO(RATIO(8, 2), RATIO(4, 2))", FUNCTIONS, [], [], []) == 2.0


def test_a_function_defined_by_its_own_name_ends():
    # the SAPHIR includes hold AMAX = FUNCTION[AMAX(%1,%2)]; the expansion must stop and still evaluate
    self_named = [function("AMAX", "AMAX(%1,%2)")]
    assert evaluate("AMAX(AMAX(1, 5), 3)", self_named, [], [], []) == 5.0


def test_function_body_is_case_normalised():
    value = evaluate("ARRH(2.0, 0)", FUNCTIONS, [entry("TEMP", 298.0)], [], [])
    assert value == 2.0


def test_function_name_may_be_lower_case_in_the_call():
    assert evaluate("double(4)", FUNCTIONS, [], [], []) == 8.0


def test_unknown_function_is_an_error_that_lists_the_known_ones():
    with pytest.raises(KeyError, match="NOSUCH"):
        evaluate("NOSUCH(1)", FUNCTIONS, [], [], [])


def test_unknown_name_inside_a_function_reports_the_expanded_expression():
    with pytest.raises(RuntimeError, match=r"Expanded: 2\*\(MISSING\)"):
        evaluate("DOUBLE(MISSING)", FUNCTIONS, [], [], [])


# ── the order of evaluation ──────────────────────────────────────────────────

def names(entries):
    return [e.name for e in entries]


def test_dependencies_come_first():
    s = [entry("C", definition="CONST(B*2)"), entry("B", definition="CONST(A+1)"), entry("A", definition="CONST(5)")]
    assert names(topological_sort(s, FUNCTIONS)) == ["A", "B", "C"]


def test_independent_entries_keep_their_order():
    s = [entry("Z", definition="CONST(1)"), entry("A", definition="CONST(2)"), entry("M", definition="INPUT(M)")]
    assert names(topological_sort(s, FUNCTIONS)) == ["Z", "A", "M"]


def test_dependency_through_a_function_body_is_seen():
    s = [entry("K", definition="ARRH(1, 100)"), entry("TEMP", definition="INPUT(T)")]
    assert names(topological_sort(s, FUNCTIONS)) == ["TEMP", "K"]


def test_a_name_defined_twice_appears_once_with_its_last_definition():
    s = [entry("A", definition="CONST(1)"), entry("B", definition="CONST(A)"), entry("A", definition="CONST(2)")]
    ordered = topological_sort(s, FUNCTIONS)
    assert names(ordered) == ["A", "B"] and ordered[0].definition == "CONST(2)"


def test_a_cycle_ends_and_keeps_every_entry():
    s = [entry("A", definition="CONST(B)"), entry("B", definition="CONST(A)")]
    assert sorted(names(topological_sort(s, FUNCTIONS))) == ["A", "B"]


# ── settings and ENZ values ──────────────────────────────────────────────────

def prepared(s):
    ordered = topological_sort(s, FUNCTIONS)
    precompile_s_vec(ordered, FUNCTIONS)
    return ordered


def test_settings_are_evaluated_in_order_from_enz_and_concentrations():
    s = [entry("K", definition="CONST(ARRH(2, 0)*M)"), entry("M", definition="CONST(P/T)"),
         entry("T", definition="INPUT(TEMP_IN)"), entry("TEMP", definition="CONST(T)"),
         entry("R", definition="DOUBLE(OH)")]
    ordered = prepared(s)
    update_s_vec(s, FUNCTIONS, [entry("OH", 1.5e6)], [entry("TEMP_IN", 300.0), entry("P", 600.0)], _ordered=ordered)
    value = {e.name: e.value for e in s}
    assert value == {"T": 300.0, "TEMP": 300.0, "M": 2.0, "K": 4.0, "R": 3.0e6}


def test_fast_replacement_is_used_instead_of_the_definition():
    s = [entry("RO2", definition="CONST(A+B)"), entry("TWICE", definition="DOUBLE(RO2)")]
    ordered = prepared(s)
    update_s_vec(s, FUNCTIONS, [entry("A", 1.0), entry("B", 2.0)], [], _ordered=ordered,
                 _fast={"RO2": lambda ctx: 10.0})
    assert (s[0].value, s[1].value) == (10.0, 20.0)


def test_a_failing_setting_warns_and_the_others_are_still_evaluated():
    s = [entry("BAD", value=-1.0, definition="CONST(1/ZERO)"), entry("GOOD", definition="CONST(7)")]
    ordered = prepared(s)
    with pytest.warns(UserWarning, match="BAD"):
        update_s_vec(s, FUNCTIONS, [entry("ZERO", 0.0)], [], _ordered=ordered)
    assert s[1].value == 7.0
    assert s[0].value == -1.0            # the previous value stays


def test_enz_values_follow_their_providers():
    enz = [entry("T", 0.0), entry("P", 5.0), entry("BROKEN", 3.0)]

    def broken(t):
        raise ValueError("outside the table")

    update_enz_vec(enz, {"T": lambda t: 280.0 + t, "BROKEN": broken}, 10.0)
    assert enz[0].value == 290.0         # from its provider
    assert enz[1].value == 5.0           # no provider: untouched
    assert enz[2].value == 0.0           # a provider that fails gives zero
