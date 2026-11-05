"""A FUNCTION macro defined more than once: EASY 6.00 applies the last definition everywhere,
EASY 5.96b applies a definition from its own line on.  JUICE follows the EASY that made the run."""
import os

import pytest

import elog_generator as eg


def test_easy_version_is_read_from_the_fac_header(tmp_path):
    (tmp_path / "x.fac").write_text("* Facsimile-Script / EASY Ver. 5.96b 2010-5-19;\n* Sun Jul 27 15:51:15 2025;\n", encoding="latin-1")
    assert eg.easy_version(str(tmp_path), stem="x") == "5.96b"
    assert eg.easy_version(str(tmp_path / "nothing")) == ""


@pytest.mark.parametrize("version,scope", [("5.96b", "use"), ("6.00", "last"), ("", "last"), ("7.1", "last")])
def test_function_scope_for(version, scope):
    assert eg.function_scope_for(version) == scope


RAW = ("F = FUNCTION[(%1)*2]\n"
       "A = F(1)\n"
       "F = FUNCTION[(%1)*3]\n"
       "B = F(1)\n"
       "X = CONST(1)\n"
       "X = CONST(2)\n")


def _canonical(tmp_path, scope):
    raw = tmp_path / "t_juice.elog"; raw.write_text(RAW, encoding="utf-8")
    out = eg.make_canonical_elog(str(raw), function_scope=scope)
    return [l.strip() for l in open(out, encoding="utf-8") if l.strip() and not l.startswith(";")]


def test_canonical_last_keeps_only_the_last_function(tmp_path):
    lines = _canonical(tmp_path, eg.FUNCTION_SCOPE_LAST)
    assert lines == ["A = F(1)", "F = FUNCTION[(%1)*3]", "B = F(1)", "X = CONST(2)"]


def test_canonical_use_keeps_every_function_definition(tmp_path):
    lines = _canonical(tmp_path, eg.FUNCTION_SCOPE_USE)
    assert lines == ["F = FUNCTION[(%1)*2]", "A = F(1)", "F = FUNCTION[(%1)*3]", "B = F(1)", "X = CONST(2)"]


def test_parser_versions_a_redefined_function(tmp_path):
    """With both definitions present, A uses the first body and B the second — in the inline
    expansion and in the evaluator's entries."""
    from easy_parser import EasyParser
    from evaluator import precompile_s_vec, topological_sort, update_s_vec
    p = tmp_path / "t.canonical.elog"; p.write_text(RAW, encoding="utf-8")
    parser = EasyParser(); mech = parser.parse(str(p))
    assert "*2" in mech.parameters["A"].replace(" ", "") and "*3" in mech.parameters["B"].replace(" ", "")
    defs = {e.name: e.definition for e in mech.s_entries}
    assert defs["A"] == "F(1)" and defs["B"] == "F__2(1)"
    assert {e.name for e in mech.function_entries} >= {"F", "F__2"}
    ordered = topological_sort(mech.s_entries, mech.function_entries)
    precompile_s_vec(ordered, mech.function_entries)
    update_s_vec(mech.s_entries, mech.function_entries, {}, {}, _ordered=ordered)
    vals = {e.name: e.value for e in mech.s_entries}
    assert vals["A"] == pytest.approx(2.0) and vals["B"] == pytest.approx(3.0)
