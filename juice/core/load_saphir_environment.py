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

    # Step-function (piecewise) columns: an ENZ line holds from its own time to the next
    # line; a query at exactly the time of a line returns THAT line (scipy kind='previous').
    # Until 2026-11-17 all breakpoints except the first were shifted by +1 ns.  With times
    # in seconds since 2000 (about 7.7E8, float64 spacing 119 ns) that shift rounds away,
    # so every such scenario already ran with the unshifted times.  With small times
    # (an ENZ file starting at 0) the shift was real: the row starting at a line still read
    # the line before, e.g. a T step at 3600 s was never applied (F. Rohrer's scene_T1;
    # his fix 9e4fc04 in MCM_SAPHIR_Tool).  No shift now: the same for every time base.
    _t_piecewise = df['time'].values.copy().astype(float)

    providers = {}
    for col in df.columns:
        if interpolation == 'piecewise' or col in _q_cols:
            # Zero-order hold with left-hand boundary convention (see above)
            providers[col] = interp1d(_t_piecewise, df[col], kind='previous',
                                      fill_value="extrapolate", bounds_error=False)
        else:
            # Linear interpolation between ENZ rows (smooth transition)
            providers[col] = interp1d(df['time'], df[col], kind='linear',
                                      fill_value="extrapolate", bounds_error=False)

    # Columns that require a BDF restart when they change:
    #   - Step-function columns (HV, Q*): photolysis J-values change
    #   - T / P (piecewise mode only): k-values + M change discontinuously
    # In linear mode T/P change smoothly — no restart needed for T/P there.
    _temp_cols = {c for c in df.columns if c.upper() in ('T', 'TEMP', 'TEMPERATURE')}
    _pres_cols = {c for c in df.columns if c.upper() in ('P', 'PRESS', 'PRESSURE')}
    if interpolation == 'piecewise':
        _restart_cols = _q_cols | _temp_cols | _pres_cols   # T/P jump → k + M change
    else:
        _restart_cols = _q_cols                              # T/P smooth → no restart

    # Compute step_change_times: internal ENZ row times where at least one
    # restart-triggering column actually changes its value.
    # NOTE: This uses raw ENZ values.  scenario_runner.py filters further
    # using effective ctx_data values (after CONST overrides are applied)
    # so that e.g. hv=0→1 in ENZ but HV=CONST(1) does not trigger a restart.
    step_change_times = []
    for i in range(len(df) - 1):
        t_next = float(df['time'].iloc[i + 1])
        if t_next <= df['time'].iloc[0] or t_next >= df['time'].iloc[-1]:
            continue  # skip first/last
        for col in _restart_cols:
            if col in df.columns and df[col].iloc[i] != df[col].iloc[i + 1]:
                step_change_times.append(t_next)
                break

    return (providers, float(df['time'].min()), float(df['time'].max()),
            headers, df['time'].values.astype(float), step_change_times)
