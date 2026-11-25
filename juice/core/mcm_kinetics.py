"""
mcm_kinetics.py — two MCM rate expressions as Python functions, for the results viewer's budget engine.

The solver does not use this module: every KMT, COMPRC and PHTFKT of a mechanism is a parameter or a FUNCTION of
its EASY files and is expanded where it is used (easy_parser.expand_macros).  The viewer puts these two functions
into its evaluation context as a fallback for an expression that names them without defining them.

Until 25 Nov 2026 the module also held KMT01–KMT16 as Troe fall-off functions with a fixed broadening factor and
no N term, KMT08–16 all as one generic form.  They were removed: a missing function leaves a visible gap in a
budget, an approximate one draws a wrong number that looks right.  The mechanisms define their own KMT values.
"""
import numpy as np


def KRO2(TEMP):
    """Generic RO2 + NO rate constant [cm3 molecule-1 s-1]."""
    return 2.9E-12 * np.exp(350 / TEMP)


def COMPRC(k0_full, ki, Fc):
    """
    Lindemann-Hinshelwood Troe rate constant — Python equivalent of the EASY/FACSIMILE
    COMPRC FUNCTION[].

    k0_full : low-pressure limit ALREADY multiplied by M and temperature factor [s^-1 or cm3/molec/s]
    ki      : high-pressure limit [cm3/molec/s]
    Fc      : broadening factor (dimensionless, typically 0.3–0.85)

    Formula (IUPAC): k = k0*ki/(k0+ki) * F
      where F = Fc ^ (1 / (1 + (log10(k0/ki) / N)^2))
      and   N = 0.75 - 1.27*log10(Fc)
    """
    ki_safe = np.clip(ki, 1e-300, None)
    ratio = np.clip(k0_full / ki_safe, 1e-30, None)
    N = 0.75 - 1.27 * np.log10(np.clip(Fc, 1e-30, None))
    F = np.float_power(np.clip(Fc, 1e-30, None), 1.0 / (1.0 + (np.log10(ratio) / N) ** 2))
    return k0_full * ki / (k0_full + ki_safe) * F
