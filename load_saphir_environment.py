import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import fnmatch
import re

def load_saphir_environment(file_path, interpolation='linear', step_col_patterns=None):
    """Load ENZ environment file and build interpolation providers.

    interpolation: 'linear'    — linear interpolation for smooth variables (T, J-rates,
                                 H2O, FLOW, ...) and piecewise for Q-switches/HV. Physically
                                 more realistic (chamber variables evolve continuously between
                                 ENZ rows). Fast: minimal BDF restarts. DEFAULT.
                   'piecewise' — zero-order hold on ALL columns, matching FACSIMILE/EASY
                                 convention. Slower (BDF restart at every ENZ row) and
                                 physically less realistic, but reproduces legacy FACSIMILE
                                 output exactly.

    step_col_patterns: list of name prefixes/exact names (case-insensitive) that identify
                       step-function (zero-order hold) control columns, e.g. ['Q', 'HV'].
                       A column matches if its name starts with or equals any pattern.
                       Defaults to ['Q', 'HV'] if None.
    """
    if step_col_patterns is None:
        step_col_patterns = ['Q*', 'HV']

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

    # Identify step-function (zero-order hold) columns using fnmatch-style patterns.
    # Patterns are matched case-insensitively against column names.
    # Examples: "Q*" matches Q, QNO, Q_NO2; "HV" matches only HV; "Q*_NO" matches Q_NO, Q25_NO.
    def _is_step_col(name):
        nu = name.upper()
        for p in step_col_patterns:
            if fnmatch.fnmatch(nu, p.upper()):
                return True
        return False

    _q_cols = {c for c in df.columns if _is_step_col(c)}

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

    # Compute step_change_times: ENZ row times where at least one step-function column
    # actually changes its value.  These are the times requiring BDF solver restarts in
    # linear mode.  In piecewise mode all ENZ rows are restart points.
    step_change_times = []
    if interpolation == 'piecewise':
        # Every internal ENZ row is a potential step change
        step_change_times = list(df['time'].values[1:-1].astype(float))
    else:
        for i in range(len(df) - 1):
            for col in _q_cols:
                if col in df.columns and df[col].iloc[i] != df[col].iloc[i + 1]:
                    step_change_times.append(float(df['time'].iloc[i + 1]))
                    break

    return (providers, float(df['time'].min()), float(df['time'].max()),
            headers, df['time'].values.astype(float), step_change_times)
