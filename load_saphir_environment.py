import pandas as pd
import numpy as np
from scipy.interpolate import interp1d

def load_saphir_environment(file_path):
    col_map = {}
    data_start_line = 0
    
    with open(file_path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('COLUMN'):
                parts = line.split('=')
                # Clean name: 'TIME / SEC...' -> 'time'
                name = parts[1].split('/')[0].split()[0].strip().lower()
                idx = int(parts[0].split()[1]) - 1
                col_map[name] = idx
            if '&&' in line:
                data_start_line = i + 1
                break

    df = pd.read_csv(file_path, sep=r',\s*', engine='python', skiprows=data_start_line, header=None)
    rev_map = {v: k for k, v in col_map.items()}
    df = df.rename(columns=rev_map)

    env_providers = {
        'TEMP': interp1d(df['time'], df['t'], fill_value="extrapolate"),
        'M':    interp1d(df['time'], (1006.5 / (df['t'] * 1.379E-19)), fill_value="extrapolate"),
        'DIL':  interp1d(df['time'], df['dilute'], fill_value="extrapolate")
    }

    j_providers = {col: interp1d(df['time'], df[col], fill_value="extrapolate") 
                   for col in df.columns if col.startswith('j')}
    
    q_providers = {col: interp1d(df['time'], df[col], fill_value="extrapolate") 
                   for col in df.columns if col.startswith('q')}

    return env_providers, j_providers, q_providers, df['time'].min(), df['time'].max()