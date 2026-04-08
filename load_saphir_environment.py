import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import re

def load_saphir_environment(file_path):
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

    providers = {}
    for col in df.columns:
        # Q-prefixed switch columns (QNO, QH2O, qVCP2, …) are discrete duty-cycle
        # signals that FACSIMILE treats as piecewise-constant (each row's value
        # applies for its full duration). Linear interpolation produces an
        # unphysical triangular ramp that delivers the same integrated mass over
        # a wider time window, causing peak under-shoot for fast-titrating species
        # like NO. Use 'previous' (zero-order hold) for these to match FACSIMILE.
        if col.lower().startswith('q') and col.lower() != 'q_seed':
            providers[col] = interp1d(df['time'], df[col], kind='previous',
                                      fill_value="extrapolate", bounds_error=False)
        else:
            providers[col] = interp1d(df['time'], df[col], fill_value="extrapolate")

    return providers, df['time'].min(), df['time'].max(), headers