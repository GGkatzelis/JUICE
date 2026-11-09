"""Fortran D exponents become E only inside numbers; GECKO-A species names such as G2D5004 stay as they are."""
from evaluator import fortran_d_to_e


def test_numbers_are_converted():
    assert fortran_d_to_e("1.0D-31*EXP(-1000/T)") == "1.0E-31*EXP(-1000/T)"
    assert fortran_d_to_e("4.00D-32 + 2d+3 - .5D0") == "4.00E-32 + 2E+3 - .5E0"
    assert fortran_d_to_e("K = 1.2e-15*H2O") == "K = 1.2e-15*H2O"           # already E: untouched


def test_identifiers_are_not():
    for name in ("G2D5004", "C2D5", "GD2000", "X1D2", "G2D5004*1.0D-9"):
        out = fortran_d_to_e(name)
        assert out.startswith(name.split("*")[0]), out
    assert fortran_d_to_e("G2D5004*1.0D-9") == "G2D5004*1.0E-9"
