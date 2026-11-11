"""
column_file_reader.py
---------------------
Universal reader for column-format files used in JUICE and FACSIMILE:

  Format A — New JUICE (default output since JUICE 2.2):
      # comment
      COLUMN 1=TIME_s
      COLUMN 2=OH
      ...
      NUMBER OF COLUMNS=N
      TIME_s,OH,...           ← CSV header
      value,value,...         ← comma-separated data

  Format B — FACSIMILE cs / pd (space-separated):
      COLUMN 1=TIME_s
      COLUMN 2=OH
      ...
       773812800 1.5e5 ...    ← space-separated data (no explicit CSV header)

  Format C — ENZ (space-separated, with &&&&& separator):
      COLUMN 1=TIME / SEC SINCE 1.1.2000 00:00:00
      COLUMN 2=sz
      ...
      NUMBER OF COLUMNS=N
      &&&&&&&&&&&&&&&&&&&&&
       773812800 35.0 ...

  Format D — Old JUICE CSV (backward compat):
      # comment
      TIME_s,OH,...
      value,value,...

Public API
----------
  load_column_file(path, index_col=0, comment='#') -> pd.DataFrame
      Returns a DataFrame with TIME_s as the index.

  detect_format(path) -> str   ('juice_new' | 'facsimile' | 'enz' | 'juice_old')
"""

import re
import io
import pandas as pd
import numpy as np


# ── Format detection ──────────────────────────────────────────────────────────

def detect_format(path: str) -> str:
    """Return one of: 'juice_new', 'facsimile', 'enz', 'juice_old'."""
    has_column_defs = False
    has_number_of_cols = False
    has_ampersand = False
    has_csv_header = False   # comma-separated line after NUMBER OF COLUMNS

    with open(path, encoding='latin-1') as fh:
        after_number = False
        for line in fh:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith(';'):
                continue
            if s.upper().startswith('COLUMN ') or s.upper().startswith('COLUMN\t'):
                has_column_defs = True
                continue
            if s.upper().startswith('NUMBER OF COLUMNS'):
                has_number_of_cols = True
                after_number = True
                continue
            if after_number:
                if s.startswith('&&&'):
                    has_ampersand = True
                elif ',' in s and not s[0].isdigit():
                    has_csv_header = True
                break

    if has_column_defs and has_number_of_cols and has_ampersand:
        return 'enz'
    if has_column_defs and has_number_of_cols and has_csv_header:
        return 'juice_new'
    if has_column_defs:
        return 'facsimile'
    return 'juice_old'


# ── Column-definition parser ──────────────────────────────────────────────────

def _normalize_col_name(raw: str, col_idx: int) -> str:
    """Normalize a FACSIMILE-style COLUMN header string to a plain identifier.

    Rules (applied in order):
      col_idx 0 (TIME)                                   → 'TIME_s'
      'p<0> HONO'  / 'd<0> NO'                          → 'P_HONO' / 'D_NO'
      'k<5,0> !UCH4+OH-->CH3!N' / 'v<…>'               → 'k_R00005:CH4+OH-->CH3'
      'OH<0>' / 'JSCALE<0>' (any NAME<...> suffix)      → 'OH' / 'JSCALE'
    """
    if col_idx == 0:
        return 'TIME_s'
    # p<n> SP or d<n> SP  (pd-file columns)
    m = re.match(r'^([pd])<\d+>\s+(\S+)', raw, re.I)
    if m:
        return m.group(1).upper() + '_' + m.group(2)
    # k<n,0> !U{rxn}!N or v<n,0> !U{rxn}!N  (kv-file columns)
    m = re.match(r'^([kv])<(\d+),\d+>\s+!U(.+?)!N', raw, re.I)
    if m:
        kv = m.group(1).lower()
        n  = int(m.group(2))
        rxn = m.group(3).strip()
        return f'{kv}_R{n:05d}:{rxn}'
    # NAME<…> — strip angle-bracket suffix
    return raw.split('<')[0].strip()


def _parse_column_defs(path: str):
    """Return ordered list of column names from COLUMN N=name lines."""
    cols = {}
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            m = re.match(r'^COLUMN\s+(\d+)\s*=\s*(.+)$', s, re.I)
            if m:
                idx  = int(m.group(1)) - 1
                name = _normalize_col_name(m.group(2).strip(), idx)
                cols[idx] = name
    if not cols:
        return []
    return [cols[i] for i in range(max(cols) + 1) if i in cols]


# ── Per-format readers ────────────────────────────────────────────────────────

def _read_juice_new(path: str, index_col: int, comment: str) -> pd.DataFrame:
    """New JUICE format: COLUMN headers + CSV data."""
    # The CSV header and data follow NUMBER OF COLUMNS=N.
    # pd.read_csv with comment='#' handles the leading comment line.
    # The COLUMN N= lines are not valid CSV → skip them.
    # Strategy: collect all lines from after NUMBER OF COLUMNS.
    buf = []
    past_header = False
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if s.upper().startswith('NUMBER OF COLUMNS'):
                past_header = True
                continue
            if past_header:
                buf.append(line)
    text = ''.join(buf)
    df = pd.read_csv(io.StringIO(text), index_col=index_col, comment=comment)
    return df


def _read_facsimile(path: str, index_col: int) -> pd.DataFrame:
    """FACSIMILE cs/pd format: COLUMN headers + space-separated data."""
    col_names = _parse_column_defs(path)
    data_lines = []
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if not s or re.match(r'^COLUMN', s, re.I):
                continue
            if s.upper().startswith('NUMBER OF COLUMNS'):
                continue
            if s and s[0].isdigit():
                # Fix missing-E scientific notation: 3.84-116 → 3.84E-116
                fixed = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', s)
                data_lines.append(fixed)
    if not data_lines:
        raise ValueError(f"No data found in FACSIMILE file: {path}")
    rows = [[float(v) for v in r.split()] for r in data_lines]
    n = max(len(row) for row in rows)
    headers = col_names[:n] if col_names else [f'COL{i}' for i in range(n)]
    # Pad columns if shorter than header
    if len(headers) < n:
        headers += [f'COL{i}' for i in range(len(headers), n)]
    df = pd.DataFrame(rows, columns=headers[:n])
    # Remove duplicate last row (FACSIMILE WHENEVER artefact)
    df = df[~df.iloc[:, 0].duplicated(keep='first')]
    if index_col is not None and 0 <= index_col < len(df.columns):
        df = df.set_index(df.columns[index_col])
    return df


def _read_enz(path: str, index_col: int) -> pd.DataFrame:
    """ENZ format: COLUMN headers + &&&& + space-separated data."""
    col_names = _parse_column_defs(path)
    data_lines = []
    in_data = False
    with open(path, encoding='latin-1') as fh:
        for line in fh:
            s = line.strip()
            if s.startswith('&&&'):
                in_data = True
                continue
            if in_data and s:
                # ENZ data may use Fortran comma-decimal or E notation quirks
                clean = s.replace(',', '.') if ',' in s and '.' not in s else s
                # Fix missing-E: 1.23-04 → 1.23E-04
                clean = re.sub(r'(\d)([-+])(\d)', r'\1E\2\3', clean)
                try:
                    data_lines.append([float(v) for v in clean.split()])
                except ValueError:
                    pass
    if not data_lines:
        raise ValueError(f"No data found in ENZ file: {path}")
    n = max(len(row) for row in data_lines)
    headers = col_names[:n] if col_names else [f'COL{i}' for i in range(n)]
    if len(headers) < n:
        headers += [f'COL{i}' for i in range(len(headers), n)]
    df = pd.DataFrame(data_lines, columns=headers[:n])
    # Remove duplicate time rows (segment-boundary artefact in JUICE ENZ output)
    df = df[~df.iloc[:, 0].duplicated(keep='last')]
    if index_col is not None and 0 <= index_col < len(df.columns):
        df = df.set_index(df.columns[index_col])
    return df


def _read_juice_old(path: str, index_col: int, comment: str) -> pd.DataFrame:
    """Old JUICE CSV format: # comment + CSV header + CSV data."""
    return pd.read_csv(path, index_col=index_col, comment=comment)


# ── Public API ────────────────────────────────────────────────────────────────

def load_column_file(path: str, index_col: int = 0,
                     comment: str = '#') -> pd.DataFrame:
    """Load any JUICE/FACSIMILE/ENZ column file into a DataFrame.

    Auto-detects format; returns DataFrame with TIME_s (or first column) as index.
    Handles:
      - New JUICE format (COLUMN headers + CSV data)
      - FACSIMILE cs/pd (COLUMN headers + space-separated data)
      - ENZ (COLUMN headers + &&&& + space-separated data)
      - Old JUICE CSV (# comment + CSV data, no COLUMN headers)
    """
    fmt = detect_format(path)
    if fmt == 'juice_new':
        return _read_juice_new(path, index_col, comment)
    elif fmt == 'facsimile':
        return _read_facsimile(path, index_col)
    elif fmt == 'enz':
        return _read_enz(path, index_col)
    else:  # juice_old
        return _read_juice_old(path, index_col, comment)
