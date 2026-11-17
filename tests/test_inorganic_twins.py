"""CO + OH once: a complete MCM export holds the inorganic reactions that JUICE's inorganic block adds again.

Reactions with the same name are one reaction, so the overlap is harmless except where the names differ: the
block writes 'OH + CO --> HO2 + CO2', the MCM 'CO + OH = HO2'.  Until 17 Nov 2026 both ran in every experiment
built from the complete MCM file (OH loss to CO doubled); the FACSIMILE scenes use IDL-written files without
inorganic reactions and never showed it."""
import re

from juice.chemistry import mcm_converter as mc
from juice.chemistry import mechanism_registry as mreg


def test_only_the_co2_twin_is_dropped():
    kept = mc._drop_inorganic_twins([(" CO + OH  -->  HO2 ", "k1"), (" HCOCO  -->  CO + OH ", "k2"),
                                     (" O3 + OH  -->  HO2 ", "k3"), (" C2H5OH + OH  -->  CH3CHO + HO2 ", "k4")])
    assert [r for _m, r in kept] == ["k2", "k3", "k4"]


def test_subset_of_the_complete_mcm_has_co_oh_once(tmp_path):
    base = mreg.load_bases()["MCM_v3.3.1"]
    import os
    res = mc.convert_mcm_file(os.path.join(base.folder, base.files["fac"]), output_dir=tmp_path, precursors=["C2H5OH"])
    mech = next(p for p in tmp_path.iterdir() if p.name.endswith(".init.mech")).read_text(encoding="latin-1")
    co_oh = [ln for ln in mech.splitlines() if re.match(r"^\s*(CO \+ OH|OH \+ CO)\s+-->", ln)]
    assert len(co_oh) == 1 and "CO2" in co_oh[0], co_oh
    assert res


def test_no_mcm_reaction_runs_beside_the_inorganic_block_under_another_name():
    """General guard, for every installed MCM download and whatever the block holds: an MCM reaction with the reactants
    of a block reaction must be one the block has itself (same products, so the two merge by name).  Anything else would
    run twice, as CO + OH did."""
    import os
    block = [mc._reaction_sides(m) for m, _r in mc.INORGANIC_DARK + mc.INORGANIC_PHOT]
    block_reactants = {s[0] for s in block}
    checked = 0
    for bid, base in mreg.load_bases().items():
        path = os.path.join(base.folder, base.files.get("fac", "")) if base.kind == "base" else ""
        if not os.path.isfile(path):
            continue
        lines = mc._combine_continuations(mc._read_lines(path))
        voc = []
        for rate, mech in mc._extract_reactions(lines):
            voc.extend(mc._convert_reaction(rate, mech))
        stray = sorted({m.strip() for m, _r in mc._drop_inorganic_twins(voc)
                        if mc._reaction_sides(m)[0] in block_reactants and mc._reaction_sides(m) not in block})
        assert not stray, (bid, stray)
        checked += 1
    assert checked
