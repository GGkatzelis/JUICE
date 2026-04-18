import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import re

def load_saphir_environment(file_path, interpolation='linear'):
    """Load ENZ environment file and build interpolation providers.

    interpolation: 'linear'    — linear interpolation for smooth variables (T, J-rates,
                                 H2O, FLOW, ...) and piecewise for Q-switches/HV. Physically
                                 more realistic (chamber variables evolve continuously between
                                 ENZ rows). Fast: minimal BDF restarts. DEFAULT.
                   'piecewise' — zero-order hold on ALL columns, matching FACSIMILE/EASY
                                 convention. Slower (BDF restart at every ENZ row) and
                                 physically less realistic, but reproduces legacy FACSIMILE
                                 output exactly.
    """
    col_map = {}
    data_start_line = 0
    headers = []

    with open(file_path, 'r', encoding='latin-1') as f:
        for i, line in enumerate(f):
            if line.startswith('COLUMN'):
                parts = line.split('=')
                name = parts[1].split('/')[0].split('(')[0].strip().lower()
                idx = int(parts[0].split()[1]) - 1
                col_map[name] = idx
                headers.append(name)
            if '&&' in line:
                data_start_line = i + 1
                break

    df = pd.read_csv(file_path, sep=r',\s*', engine='python', skiprows=data_start_line, header=None)
    rev_map = {v: k for k, v in col_map.items()}
    df = df.rename(columns=rev_map)

    # Q-switch columns (binary on/off) are ALWAYS piecewise-constant
    _q_cols = {c for c in df.columns if c.upper().startswith('Q') or c.upper() == 'HV'}

    providers = {}
    for col in df.columns:
        if interpolation == 'piecewise' or col in _q_cols:
            # Zero-order hold: value stays constant until next ENZ row
            providers[col] = interp1d(df['time'], df[col], kind='previous',
                                      fill_value="extrapolate", bounds_error=False)
        else:
            # Linear interpolation between ENZ rows (smooth transition)
            providers[col] = interp1d(df['time'], df[col], kind='linear',
                                      fill_value="extrapolate", bounds_error=False)

    return providers, df['time'].min(), df['time'].max(), headers