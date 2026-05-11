"""
elog_generator.py
-----------------
JUICE preprocessor: expands all #include directives recursively
(like FACSIMILE preprocessor) and produces two output files:

  <scenario>_juice.elog           — raw expanded file, FACSIMILE-compatible
                                    annotation format, all comments preserved
  <scenario>_juice.canonical.elog — deduplicated canonical form:
                                    for every parameter, rate-constant, initial
                                    condition and ALIAS definition, only the
                                    LAST occurrence is kept ("last wins").
                                    This is the authoritative input for JUICE.

Raw elog annotation format (mirrors FACSIMILE):
    ;#include filename ! /path/to/containing/dir/   <- opening
    ... lines from that file ...
    ;#include filename                              <- closing

"""

import os
import re
import datetime


# ---------------------------------------------------------------------------
# Step 1: expand includes → raw elog
# ---------------------------------------------------------------------------

def generate_juice_elog(easy_path: str,
                        search_dirs: list = None,
                        output_path: str = None) -> str:
    """Expand all #include directives recursively and save as _juice.elog.

    Annotation format mirrors FACSIMILE:
      ;#include name ! /containing/dir/
      ... content ...
      ;#include name

    All comment lines from source files are preserved.
    $ line-continuation is resolved before output.
    Returns the path to the generated .elog file.
    """
    search_dirs = [d for d in (search_dirs or []) if d and os.path.isdir(d)]

    def find_file(name: str, base_dir: str) -> str:
        candidate = os.path.join(base_dir, name)
        if os.path.isfile(candidate):
            return candidate
        basename = os.path.basename(name)
        for d in search_dirs:
            for root, _, files in os.walk(d):
                if basename in files:
                    return os.path.join(root, basename)
        return None

    def expand(file_path: str, visited: set = None) -> list:
        if visited is None:
            visited = set()
        real = os.path.realpath(file_path)
        if real in visited:
            return [f";#include {os.path.basename(file_path)} ! CIRCULAR INCLUDE SKIPPED\n"]
        visited.add(real)

        base_dir = os.path.dirname(os.path.abspath(file_path))
        lines_out = []

        with open(file_path, encoding='latin-1') as fh:
            raw = fh.read()

        # Resolve $ line-continuation (FACSIMILE convention)
        raw = re.sub(r'\$\s*\n', ' ', raw)

        for line in raw.splitlines(keepends=True):
            stripped = line.strip()
            is_commented = stripped.startswith((';', '*'))

            if not is_commented:
                m = re.match(r'^[^;]*#include\s+["\']?(\S+?)["\']?\s*(?:;.*)?$',
                             line, re.I)
                if m:
                    inc_name = m.group(1).strip()
                    inc_path = find_file(inc_name, base_dir)
                    if inc_path:
                        inc_dir = os.path.dirname(os.path.abspath(inc_path))
                        # FACSIMILE-compatible opening annotation
                        lines_out.append(
                            f";#include {inc_name} ! {inc_dir}{os.sep}\n")
                        lines_out.extend(expand(inc_path, visited))
                        # FACSIMILE-compatible closing annotation
                        lines_out.append(f";#include {inc_name}\n")
                    else:
                        lines_out.append(
                            f";#include {inc_name} ! NOT FOUND\n")
                        lines_out.append(line)
                    continue

            # All other lines (content and comments) are kept verbatim
            lines_out.append(line)

        return lines_out

    all_lines = expand(easy_path)

    if output_path is None:
        stem = os.path.splitext(easy_path)[0]
        output_path = stem + '_juice.elog'

    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(output_path, 'w', encoding='utf-8') as fh:
        fh.write(f"; JUICE elog generated from: {easy_path}\n")
        fh.write(f"; Date: {ts}\n")
        fh.write(f";\n")
        fh.writelines(all_lines)

    return output_path


# ---------------------------------------------------------------------------
# Step 2: deduplicate → canonical elog  ("last definition wins")
# ---------------------------------------------------------------------------

def _definition_key(line: str):
    """Return a hashable key if the line is a definition, else None.

    Definitions whose key can appear more than once and must be deduplicated:
      PARAM:X      for  X = anything
      INITIAL:X    for  INITIAL[X] = anything
      ALIAS:X      for  ALIAS[X] = anything
      RATE:raw     for  k[reaction] = anything
      FACS:param   for  FACS[param] = anything

    Comments, blank lines and bare reaction lines (A+B-->C without k[])
    are not definitions — they return None and are always kept.
    """
    stripped = line.strip()

    # Skip blank lines and comment lines
    if not stripped or stripped.startswith((';', '*', '#')):
        return None

    # Skip lines that are purely a reaction equation (contain --> but NOT k[...]=)
    if '-->' in stripped and not re.match(r'^k\s*\[', stripped, re.I):
        return None

    # k[reaction] = rate
    m = re.match(r'^k\s*\[(.+?)\]\s*=', stripped, re.I)
    if m:
        return ('RATE', m.group(1).strip().upper())

    # INITIAL[X] = value
    m = re.match(r'^INITIAL\s*\[([^\]]+)\]\s*=', stripped, re.I)
    if m:
        return ('INITIAL', m.group(1).strip().upper())

    # ALIAS[X] = ...
    m = re.match(r'^ALIAS\s*\[([^\]]+)\]\s*=', stripped, re.I)
    if m:
        return ('ALIAS', m.group(1).strip().upper())

    # FACS[param] = value
    m = re.match(r'^FACS\s*\[([^\]]+)\]\s*=', stripped, re.I)
    if m:
        return ('FACS', m.group(1).strip().upper())

    # X = anything  (parameter / function / constant assignment)
    m = re.match(r'^([A-Za-z_]\w*)\s*=', stripped)
    if m:
        return ('PARAM', m.group(1).strip().upper())

    return None


def make_canonical_elog(raw_elog_path: str,
                         output_path: str = None) -> str:
    """Deduplicate the raw elog so every definition key appears only once.

    Algorithm:
      Forward pass  — collect {key → last line-index} for all definitions.
      Backward pass — write only lines where this is the last occurrence of
                      their key (all non-definition lines are always written).

    The result is the canonical elog:
      - Parameter X defined N times → only the LAST definition is kept.
      - Comment and blank lines are preserved unchanged.
      - Bare reaction equations (A+B-->C) are always kept.

    Returns the path to the written .canonical.elog file.
    """
    with open(raw_elog_path, encoding='utf-8', errors='replace') as fh:
        lines = fh.readlines()

    # Forward pass: key → index of last occurrence
    key_last: dict = {}
    for i, line in enumerate(lines):
        k = _definition_key(line)
        if k is not None:
            key_last[k] = i   # overwrite → last index wins

    # Filter: keep line if it has no definition key, or if it IS the last occurrence
    kept = []
    removed = 0
    for i, line in enumerate(lines):
        k = _definition_key(line)
        if k is None or key_last[k] == i:
            kept.append(line)
        else:
            removed += 1

    if output_path is None:
        stem = os.path.splitext(raw_elog_path)[0]
        # Remove _juice suffix if present so we get <scenario>_juice.canonical.elog
        output_path = stem + '.canonical.elog'

    ts = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    with open(output_path, 'w', encoding='utf-8') as fh:
        fh.write(f"; JUICE canonical elog — last definition wins\n")
        fh.write(f"; Erzeugt aus: {raw_elog_path}\n")
        fh.write(f"; Datum: {ts}\n")
        fh.write(f"; Entfernte Mehrfach-Definitionen: {removed}\n")
        fh.write(f";\n")
        fh.writelines(kept)

    return output_path


# ---------------------------------------------------------------------------
# Combined entry point: generate both files
# ---------------------------------------------------------------------------

def generate_both(easy_path: str,
                  search_dirs: list = None,
                  output_dir: str = None) -> tuple:
    """Generate raw elog and canonical elog from an EASY file.

    output_dir: directory for the output files (default: same dir as easy_path).
                Set to the JUICE scenario folder so the elogs land there.

    Returns (raw_elog_path, canonical_elog_path).
    """
    stem = os.path.splitext(os.path.basename(easy_path))[0]
    if output_dir and os.path.isdir(output_dir):
        raw_path   = os.path.join(output_dir, stem + '_juice.elog')
        canon_path = os.path.join(output_dir, stem + '_juice.canonical.elog')
    else:
        raw_path   = None   # default: next to easy_path
        canon_path = None

    raw_path   = generate_juice_elog(easy_path,
                                     search_dirs=search_dirs,
                                     output_path=raw_path)
    canon_path = make_canonical_elog(raw_path, output_path=canon_path)

    # Remove stale copies next to the .easy file (from runs before output_dir was set).
    if output_dir and os.path.isdir(output_dir):
        _easy_stem = os.path.splitext(easy_path)[0]
        for _old in (_easy_stem + '_juice.elog',
                     _easy_stem + '_juice.canonical.elog'):
            if os.path.isfile(_old) and os.path.abspath(_old) not in (
                    os.path.abspath(raw_path), os.path.abspath(canon_path)):
                try:
                    os.remove(_old)
                except OSError:
                    pass  # ignore if locked or already gone

    return raw_path, canon_path


if __name__ == "__main__":
    import sys, argparse
    ap = argparse.ArgumentParser(
        description="Generate JUICE elog from an EASY file.")
    ap.add_argument("easy", help="Path to the .easy file")
    ap.add_argument("search_dirs", nargs="*",
                    help="Additional #include search directories")
    ap.add_argument("--output-dir", default=None,
                    help="Directory for output elogs (default: same dir as .easy file). "
                         "Set to the JUICE scenario folder to keep EASY directory clean.")
    args = ap.parse_args()
    raw, canon = generate_both(args.easy,
                               search_dirs=args.search_dirs,
                               output_dir=args.output_dir)
    print(f"Raw elog:       {raw}")
    print(f"Canonical elog: {canon}")
