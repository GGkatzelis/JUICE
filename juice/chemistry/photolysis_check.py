"""Photolysis reactions that the light switch does not close.

In a chamber under the sun the measured J values come from outside the chamber (SAPHIR: the building roof), so they
are non-zero while the roof is closed.  A photolysis reaction is closed in one of two ways:

    X + hv --> ...                  hv (HV) as a reactant: the MCM, its add-ons, GECKO-A
    J4 = CONST(jno2*HV)             the switch inside the definition of the J the rate uses: lumped mechanisms;
    J41 = PHTFKT(...)               the clear-sky PHTFKT of j_clear_sky.incl carries HV as well

A reaction with neither runs in the dark (SAPHFIRE 2026-07-02: furfural of the NOAA_BB add-on, rated by a PHTFKT
that saphir.jmeasured_02.incl redefines without HV).  `ungated(spec)` finds such reactions in what a person may
have written: the add-ons, the includes of a lumped base, extra includes and the experiment's own lines.  The
MCM subset and a GECKO-A base are written by JUICE's converters with hv and are not read here.
"""
from __future__ import annotations

import os
import re
from functools import lru_cache
from typing import Dict, List, Tuple

from juice.chemistry import mechanism_registry as mreg

_DEF = re.compile(r"^\s*([A-Za-z_]\w*)\s*=\s*(.+)$", re.S)
_RATE = re.compile(r"^\s*k\[(.+?)\]\s*=\s*(.+)$", re.I | re.S)
_ID = re.compile(r"(?<![\w.])[A-Za-z_]\w*")
_J_NAME = re.compile(r"^(J\w*|KHV\w*)$", re.I)                 # J41, JN52, JNO2_RACM2, jno2, jhcho_r, KHV_G1700
_J_FUNCTIONS = ("PHTFKT", "PHOTMMSC", "GPHOT", "PHTFKT_CS")
_SWITCH = ("HV", "HV_IN")


def _logical_lines(text: str):
    buf = ""
    for raw in text.splitlines():
        code = raw.split(";", 1)[0].rstrip()
        if code.endswith("$"):
            buf += code[:-1] + " "
            continue
        line, buf = buf + code, ""
        if line.strip():
            yield line


def _parse(text: str) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    """(definitions NAME -> right-hand side, [(reaction, rate expression)]) of EASY text; last definition wins."""
    defs: Dict[str, str] = {}
    rates: List[Tuple[str, str]] = []
    for line in _logical_lines(text):
        m = _RATE.match(line)
        if m:
            rates.append((re.sub(r"\s+", " ", m.group(1).strip()), m.group(2).strip()))
            continue
        m = _DEF.match(line)
        if m and "-->" not in line:
            defs[m.group(1).upper()] = m.group(2).strip()
    return defs, rates


@lru_cache(maxsize=64)
def _parse_file(path: str, mtime: float) -> Tuple[Dict[str, str], List[Tuple[str, str]]]:
    with open(path, encoding="latin-1") as fh:
        return _parse(fh.read())


def _files(spec) -> List[str]:
    """The EASY files of this experiment that are not written by a converter with hv."""
    out: List[str] = []
    try:
        b = mreg.base(spec.base_mechanism)
        if not spec.mechanism_stem and mreg.is_lumped(spec.base_mechanism):
            out += [os.path.join(b.folder, f) for f in sorted(os.listdir(b.folder)) if f.lower().endswith(".incl")]
        if spec.addons:
            for a in mreg.addons_for(b.id, list(spec.addons)):
                out += a.paths()
    except Exception:
        pass
    for name in getattr(spec, "extra_includes", []) or []:
        for d in [""] + list(getattr(spec, "include_dirs", []) or []):
            p = os.path.join(d, name) if d else name
            if os.path.isfile(p):
                out.append(p)
                break
    return list(dict.fromkeys(out))


def find(easy_text: str, files: List[str], inputs_follow_switch: bool = False) -> List[Tuple[str, str, str]]:
    """[(reaction, rate expression, where)] of the photolysis reactions without a light switch.

    `easy_text` is the experiment's EASY file (its J definitions, overrides and #include lines), `files` the
    mechanism files to read with it.  `inputs_follow_switch`: the measured J are zero without light anyway
    (lamps: jNO2 is the lamps' own time series)."""
    main_defs, main_rates = _parse(easy_text)
    defs: Dict[str, str] = {}
    rates: List[Tuple[str, str, str]] = []
    for p in files:
        try:
            d, r = _parse_file(p, os.path.getmtime(p))
        except OSError:
            continue
        defs.update(d)
        rates += [(rx, ex, os.path.basename(p)) for rx, ex in r]
    defs.update(main_defs)                                   # the experiment's lines come last
    rates += [(rx, ex, "the experiment's own lines") for rx, ex in main_rates]
    # saphir.jmeasured_*.incl redefines PHTFKT without HV (its own reactions carry '+ hv')
    phtfkt_gated = not re.search(r"^\s*#include\s+\S*jmeasured", easy_text, re.I | re.M)

    state: Dict[str, Tuple[bool, bool]] = {}                 # NAME -> (is a photolysis frequency, switch applied)

    def look(name: str, depth: int = 0) -> Tuple[bool, bool]:
        if name in state:
            return state[name]
        state[name] = (bool(_J_NAME.match(name)), False)     # guard against circular definitions
        rhs = defs.get(name)
        if rhs is None:                                      # not defined here: a J<n> of the mechanism's own PHTFKT list
            state[name] = (bool(_J_NAME.match(name)), phtfkt_gated)
            return state[name]
        state[name] = expr(rhs, depth + 1, named_j=bool(_J_NAME.match(name)))
        return state[name]

    def expr(text: str, depth: int = 0, named_j: bool = False) -> Tuple[bool, bool]:
        """(the expression is a photolysis frequency, the light switch is applied to it)"""
        if re.match(r"\s*INPUT\s*\(", text, re.I):           # a measured column
            return named_j, bool(inputs_follow_switch)
        ids = [i.upper() for i in _ID.findall(text)]
        switched = any(i in _SWITCH for i in ids)
        found, gated = False, True                           # found: a J function or another J inside
        for i in ids:
            if i in _J_FUNCTIONS:
                found, gated = True, gated and i == "PHTFKT" and phtfkt_gated
            elif i not in _SWITCH and depth < 12 and (i in defs or _J_NAME.match(i)):
                j, g = look(i, depth)
                if j:
                    found, gated = True, gated and g
        if not (found or named_j):
            return False, False
        if switched:
            return True, True
        if not found:                                        # a J given as a number: only zero needs no switch
            body = re.sub(r"^CONST\((.*)\)$", r"\1", text.strip(), flags=re.I | re.S)
            try:
                return True, float(body.replace("D", "E").replace("d", "e")) == 0.0
            except ValueError:
                return True, False
        return True, gated

    out = []
    for rx, ex, where in rates:
        lhs = rx.split("-->", 1)[0]
        if re.search(r"(?<![\w])HV(?![\w])", lhs, re.I):
            continue
        # a reactant that is itself a constant carrying the switch (SAPHIR's HONO source: SHONO = CONST(HV_in*...))
        if any(t.upper() in defs and any(i.upper() in _SWITCH for i in _ID.findall(defs[t.upper()]))
               for t in _ID.findall(lhs)):
            continue
        is_j, gated = expr(ex)
        if is_j and not gated:
            out.append((rx, re.sub(r"^CONST\((.*)\)$", r"\1", ex, flags=re.I | re.S), where))
    return out


_RESULTS: Dict[tuple, List[Tuple[str, str, str]]] = {}


def ungated(spec) -> List[Tuple[str, str, str]]:
    """`find` for an experiment (ExperimentSpec); [] when the experiment cannot be rendered yet."""
    from juice.core import experiment_spec as es
    try:
        text = es.render_easy(spec)
        lamps = es._profile(spec.conditions).light_type != "sun"
    except Exception:
        return []
    files = _files(spec)
    key = (text, tuple((p, os.path.getmtime(p)) for p in files if os.path.isfile(p)), lamps)
    if key not in _RESULTS:                                  # the wizard checks after every change: same text, same answer
        if len(_RESULTS) > 32:
            _RESULTS.clear()
        _RESULTS[key] = find(text, files, inputs_follow_switch=lamps)
    return _RESULTS[key]


LIGHTS = (("measured", "SAPHIR"), ("clear_sky", "SAPHIR"), ("clear_sky", "generic_lamps"))


def _test_experiment(base_id: str, addons: List[str], j: str, chamber: str):
    import juice
    exp = juice.new_experiment("photolysis_check", chamber=chamber)
    c = exp.conditions
    c.source, c.date, c.start, c.end = "manual", "2024-07-09", "08:00", "10:00"
    c.temperature_mode, c.j_source, c.roof_open = "const", j, [("09:00", "10:00")]
    if mreg.is_complete(base_id):
        exp.base_mechanism, exp.precursors = base_id, []
    else:
        exp.base_mechanism, exp.precursors = base_id, ["TOLUENE"]
    exp.addons = list(addons)
    exp.injections = [juice.Injection(name="NO", products="NO", total_ppb="5", start="08:05")]
    return exp


def check_installed(lights=LIGHTS) -> List[str]:
    """Every lumped base and every add-on installed under chemistry/mechanisms, under sunlight with measured and
    clear-sky J and under lamps: the photolysis reactions the light switch does not close, as one line each.
    [] when all are clean.  An add-on is checked on each installed base it is compatible with.  Unlike `ungated`,
    an experiment that cannot be written is a problem here, not a pass.  tests/test_photolysis_check.py runs this,
    so a new add-on or converted mechanism is checked by the test suite as soon as its folder is in place."""
    from juice.core import experiment_spec as es
    bases = mreg.load_bases()
    cases = [(b, []) for b, m in bases.items() if m.is_lumped]
    for a in mreg.load_addons().values():
        on = [b for b in bases if a.compatible(b) and not bases[b].is_gecko] or ([mreg.DEFAULT_BASE] if not a.compatible_with else [])
        if not on:
            continue
        cases.append((mreg.DEFAULT_BASE if mreg.DEFAULT_BASE in on else on[0], [a.id]))
    out = []
    for base_id, addons in cases:
        label = " + ".join([base_id] + addons)
        for j, chamber in lights:
            spec = _test_experiment(base_id, addons, j, chamber)
            try:
                text = es.render_easy(spec)
                lamps = es._profile(spec.conditions).light_type != "sun"
            except Exception as e:
                out.append(f"{label} ({j}, {chamber}): the experiment cannot be written: {e}")
                continue
            for rx, ex, where in find(text, _files(spec), inputs_follow_switch=lamps):
                out.append(f"{label} ({j}, {chamber}): {rx}  rate {ex}  in {where}")
    return out


def messages(spec, limit: int = 4) -> List[str]:
    hits = ungated(spec)
    if not hits:
        return []
    shown = "; ".join(f"{rx} (rate {ex}, {where})" for rx, ex, where in hits[:limit])
    more = f"; and {len(hits) - limit} more" if len(hits) > limit else ""
    return [f"Photolysis without the light switch: {shown}{more}. "
            f"{'These reactions run' if len(hits) > 1 else 'This reaction runs'} with the light off as well: "
            "the photolysis frequency in the rate does not carry the switch. Write the reaction with '+ hv' "
            "(as the MCM does) or multiply its rate by HV."]
