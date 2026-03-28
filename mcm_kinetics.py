import numpy as np

def _troe(k0, kinf, M, Fc=0.6):
    """Universal Troe fall-off equation engine"""
    k0_m = k0 * M
    rel = k0_m / kinf
    # Add a small epsilon to prevent log10(0) domain errors if concentrations drop
    return (k0_m / (1.0 + rel)) * Fc ** (1.0 / (1.0 + (np.log10(np.clip(rel, 1e-30, None)))**2))

def KMT01(A, B, C, D, TEMP, M):
    return _troe(A * np.exp(B / TEMP), C * np.exp(D / TEMP), M)

def KMT02(A, B, C, D, E, TEMP, M):
    return _troe(A * (TEMP/300)**B * np.exp(C / TEMP), D * np.exp(E / TEMP), M)

def KMT03(A, B, C, D, E, TEMP, M):
    return _troe(A * np.exp(B / TEMP), C * (TEMP/300)**D * np.exp(E / TEMP), M)

def KMT04(A, B, C, D, E, F, TEMP, M):
    return _troe(A * (TEMP/300)**B * np.exp(C / TEMP), D * (TEMP/300)**E * np.exp(F / TEMP), M)

def KMT05(A, B, C, D, E, TEMP, M):
    # KMT05 explicitly passes Fc as the E variable
    return _troe(A * np.exp(B / TEMP), C * np.exp(D / TEMP), M, Fc=E)

def KMT06(A, B, C, D, TEMP, M):
    return KMT01(A, B, C, D, TEMP, M)

def KMT07(A, B, C, D, E, F, TEMP, M):
    return KMT04(A, B, C, D, E, F, TEMP, M)

def fallback_kmt(*args):
    """A bulletproof fallback for unspecified KMT functions to prevent solver crashing."""
    A, B, C, D = args[0], args[1], args[2], args[3]
    TEMP, M = args[-2], args[-1]
    return _troe(A * np.exp(B / TEMP), C * np.exp(D / TEMP), M)

# Map the remaining KMTs to the generic fallback so they never evaluate to 0.0
KMT08 = KMT09 = KMT10 = KMT11 = KMT12 = KMT13 = KMT14 = KMT15 = KMT16 = fallback_kmt

def KRO2(TEMP):
    """Generic RO2 + NO reaction rate"""
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