import numpy as np

def KMT06(A, B, C, D, TEMP, M):
    """MCM standard Three-body rate law"""
    k0 = A * np.exp(B / TEMP)
    kinf = C * np.exp(D / TEMP)
    k0_m = k0 * M
    rel = k0_m / kinf
    return (k0_m / (1 + rel)) * 0.6**(1 + (np.log10(rel))**2)**-1

def KRO2(TEMP):
    """Generic RO2 + NO reaction rate"""
    return 2.9E-12 * np.exp(350 / TEMP)

def KMT07(A, B, C, D, E, F, TEMP, M):
    """MCM standard Troe fall-off"""
    k0 = A * (TEMP/300)**B * np.exp(C/TEMP)
    kinf = D * (TEMP/300)**E * np.exp(F/TEMP)
    k0_m = k0 * M
    rel = k0_m / kinf
    return (k0_m / (1 + rel)) * 0.6**(1 + (np.log10(rel))**2)**-1

def KMT08(A, B, C, D, E, TEMP, M):
    """Pressure dependent with specific exponent"""
    k0 = A * np.exp(B / TEMP)
    kinf = C * np.exp(D / TEMP)
    k_inf_m = kinf / M
    # ... implementation specific to the MCM definition ...
    return A * np.exp(B/TEMP) # Placeholder for specific KMT08 logic