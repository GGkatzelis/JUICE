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
        # FACSIMILE treats ALL ENZ columns as piecewise-constant: each row's
        # value applies until the next row (zero-order hold). Use kind='previous'
        # for all columns to match this convention exactly.
        providers[col] = interp1d(df['time'], df[col], kind='previous',
                                  fill_value="extrapolate", bounds_error=False)

    return providers, df['time'].min(), df['time'].max(), headers