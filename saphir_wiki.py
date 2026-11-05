"""
saphir_wiki.py — read a SAPHIR experiment from the saphirwiki action log.

The wiki (https://saphir.fz-juelich.de/saphirwiki/YYYY-MM-DD) is a MoinMoin instance
inside the FZJ network; reading a day page needs a wiki login.  This module accepts
the page in any of these forms and turns the action log into wizard content:

    * live fetch with a login            fetch_wiki_text(date, user, password)
    * a page saved from the browser      load_wiki_file("2024-07-09.html")
    * the raw wiki text (?action=raw)    load_wiki_file("2024-07-09.txt")
    * pasted text                        parse_action_log(text)

Action-log lines look like   HH:MM <code> <person> <TYPE>: description
Codes (same convention as enz_builder):
    R / r   roof open / close                 H / h   humidification start / stop
    S / s   seed start / stop                 X       VOC injection
    O       ozone        M   NO     D   NO2   c   CO      C   CO2 tracer
    $ and anything else = status / note (kept, not modelled)

`events_to_import()` builds roof times, InjectionEvents with their action type,
MCM precursors and a list of checks for the user ("wiki comments can be wrong, but
they are the best starting point").  Nothing here needs Qt.
"""
from __future__ import annotations

import html
import os
import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

from experiment_spec import InjectionEvent, ACTIONS, WIKI_CODE_ACTION, H2O_IMPURITY_DEFAULTS
import injection_calc as ic

WIKI_BASE = "https://saphir.fz-juelich.de/saphirwiki"
ENV_USER, ENV_PASS = "SAPHIR_WIKI_USER", "SAPHIR_WIKI_PASS"     # login defaults (never stored by JUICE)


def wiki_credentials(user: str = "", password: str = "") -> Tuple[str, str]:
    """Given values win; empty ones fall back to the SAPHIR_WIKI_USER / SAPHIR_WIKI_PASS environment variables."""
    return (user or os.environ.get(ENV_USER, "")).strip(), password or os.environ.get(ENV_PASS, "")
FAST_DURATION_MIN = 1.0          # single-event injections (NO, O3, VOC …) in the wiki are "instant"
SAPHIR_VOLUME_M3 = 276.0         # chamber volume used by the IDL dilution formula

# emission mixtures that have no single MCM species (recorded only unless the user maps them)
MIXTURE_WORDS = {"gasoline": "GASOLINE", "diesel": "DIESEL", "cooking": "COOKING", "biogenic": "BIOGENIC",
                 "vcp": "VCP", "mixture": "MIXTURE", "mix": "MIXTURE"}

_UNITS = [                                   # a stated mixing ratio beats a liquid volume
    (r"(\d+(?:[.,]\d+)?)\s*ppbv?\b", "ppb"),
    (r"(\d+(?:[.,]\d+)?)\s*ppmv?\b", "ppm"),
    (r"(\d+(?:[.,]\d+)?)\s*(?:µl|μl|ul|microliter|mikroliter|microl|mikrol)\b", "ul"),
    (r"(\d+(?:[.,]\d+)?)\s*(?:ml|milliliter)\b", "ml"),
    (r"(\d+(?:[.,]\d+)?)\s*(?:nl|nanoliter)\b", "nl"),
]
_DURATION = re.compile(r"(\d+(?:[.,]\d+)?)\s*(sec|seconds|s|min|minutes)\b(?!\s*-1)")
_KOH = re.compile(r"(\d+(?:[.,]\d+)?)\s*s\s*-1\b")
STOP_WORDS = ("stop", "stopp", "end", "ende", "finish", "switch to experiment")


def parse_duration_min(desc: str) -> Optional[float]:
    """'50 sec' -> 0.83, '5 min' -> 5; None when the line gives no duration."""
    m = _DURATION.search(desc.lower())
    if not m:
        return None
    v = _num(m.group(1))
    return v / 60.0 if m.group(2).startswith("s") else v


def is_stop_line(desc: str) -> bool:
    return any(w in desc.lower() for w in STOP_WORDS)


@dataclass
class WikiEvent:
    time: str                 # HH:MM
    code: str                 # single character
    description: str          # rest of the line (person / TYPE: text)
    line: str = ""            # original line
    use: bool = True          # user may untick in the import dialog

    @property
    def action(self) -> str:
        if self.code == "R":
            return "roof open"
        if self.code == "r":
            return "roof close"
        if self.code in WIKI_CODE_ACTION:
            a = WIKI_CODE_ACTION[self.code]
            if self.code in ("h", "s"):
                return ACTIONS[a]["label"] + " — stop"
            return ACTIONS[a]["label"]
        return "note"


@dataclass
class WikiImport:
    date: str = ""
    roof_open: List[Tuple[str, str]] = field(default_factory=list)
    injections: List[InjectionEvent] = field(default_factory=list)
    precursors: List[str] = field(default_factory=list)
    notes: List[str] = field(default_factory=list)             # status / free-text lines
    checks: List[Tuple[str, str]] = field(default_factory=list)   # (level, message) for the user


# ─────────────────────────────────────────────────────────────────────────────
# Fetching / loading
# ─────────────────────────────────────────────────────────────────────────────

def strip_html(text: str) -> str:
    if "<" not in text:
        return text
    text = re.sub(r"(?is)<(script|style).*?</\1>", " ", text)
    text = re.sub(r"(?i)<br\s*/?>|</p>|</div>|</li>|</tr>|</h\d>", "\n", text)
    text = re.sub(r"<[^>]+>", " ", text)
    text = html.unescape(text)
    return re.sub(r"[ \t\xa0]+", " ", text)


def wiki_url(date: str) -> str:
    return f"{WIKI_BASE}/{date}"


def fetch_wiki_text(date: str, user: str = "", password: str = "", timeout: int = 20) -> str:
    """Return the raw wiki text of the day page.  Logs in with user/password if given.
    Raises PermissionError when the page needs a login we do not have, and
    requests exceptions when the server is unreachable (VPN)."""
    import requests
    import urllib3
    urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)   # FZJ internal CA
    sess = requests.Session()
    url = wiki_url(date)
    user, password = wiki_credentials(user, password)

    def _raw() -> Tuple[int, str]:
        r = sess.get(url, params={"action": "raw"}, timeout=timeout, verify=False)
        try:
            return r.status_code, r.content.decode("utf-8")
        except UnicodeDecodeError:
            return r.status_code, r.text

    status, text = _raw()
    needs_login = status == 403 or "not allowed to view" in text.lower()
    if needs_login and user:
        sess.post(url, params={"action": "login"},
                  data={"action": "login", "name": user, "password": password, "login": "Login"},
                  timeout=timeout, verify=False)
        status, text = _raw()
        needs_login = status == 403 or "not allowed to view" in text.lower()
    if needs_login:
        raise PermissionError("The wiki page needs a login (MoinMoin: 'not allowed to view this page'). "
                              + ("The login with user '%s' was not accepted; check the credentials, " % user if user else
                                 "Enter your wiki user name and password (or set %s / %s), " % (ENV_USER, ENV_PASS))
                              + "or save the page from the browser and load the file.")
    if status == 404:
        raise FileNotFoundError(f"No wiki page for {date} ({url})")
    if status != 200:
        raise RuntimeError(f"HTTP {status} from {url}")
    return text


def load_wiki_file(path: str) -> str:
    """Text of a wiki page saved from the browser (.html/.htm) or as raw text (.txt/.md)."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        return fh.read()


# ─────────────────────────────────────────────────────────────────────────────
# Parsing
# ─────────────────────────────────────────────────────────────────────────────

_LINE = re.compile(r"^\s*(\d{1,2}):(\d{2})(?::(\d{2}))?\s+([A-Za-z$])(?=\s|$)(.*)$")
# 'HH:MM H.F. CHANGE: Roof open' — no one-letter code; the action is read from the words
_FREE = re.compile(r"^\s*(\d{1,2}):(\d{2})(?::(\d{2}))?\s+(.+)$")
_INITIALS_TAG = re.compile(r"^\s*(?:(?:[A-ZÄÖÜ][a-zäöü]?\.\s*){1,3})?\s*(?:(CHANGE|ZUSTAND|STATUS|STATE|INFO|NOTE|NOTIZ|COMMENT|KOMMENTAR)\s*:)?\s*(.*)$", re.IGNORECASE)
_STATE_TAGS = ("zustand", "status", "state", "info", "note", "notiz", "comment", "kommentar")
_ADD_WORDS = r"(inject|injection|injected|addition|added|add |zugabe|dosier|eingeleitet|filled|fill )"


def infer_code(text: str) -> str:
    """The one-letter action code for a log line written without one ('' = a note).
    Reads the words: roof open/closed, humidification on/off, seed on/off, ozone / NO / NO2 / CO / CO2
    additions, and any other injection or addition → X (VOC or mixture).  State reports
    ('ZUSTAND:', 'STATUS:') and flushing / valve / instrument lines stay notes."""
    m = _INITIALS_TAG.match(text)
    tag, body = (m.group(1) or "").lower(), (m.group(2) if m else text)
    if tag in _STATE_TAGS:
        return ""
    t = " " + body.lower() + " "
    off = re.search(r"\b(off|stop|stopp|end|ende|aus|closed|close|zu|geschlossen|finished)\b", t) is not None
    if re.search(r"\b(roof|dach)\b", t):
        return "r" if re.search(r"\b(clos|zu\b|geschlossen|schlie)", t) else ("R" if re.search(r"\b(open|opened|auf|offen|geöffnet|öffn)", t) else "")
    if re.search(r"(humidif|befeucht)", t):
        return "h" if off else "H"
    if re.search(r"\b(seed|aerosol|ammonium ?sulfat)", t):
        return "s" if off else "S"
    adding = re.search(_ADD_WORDS, t) is not None or re.search(r"\d\s*(ppb|ppm|µl|μl|ul|microl|ml|sccm)\b", t) is not None
    if not adding:
        return ""
    if re.search(r"\b(o3|ozone|ozon)\b", t):
        return "O"
    if re.search(r"\bno2\b", t):
        return "D"
    if re.search(r"\b(no|nitric oxide|stickstoffmonoxid)\b", t):
        return "M"
    if re.search(r"\bco2\b", t):
        return "C"
    if re.search(r"\b(co|carbon monoxide|kohlenmonoxid)\b", t):
        return "c"
    return "X"


def parse_action_log(text: str) -> List[WikiEvent]:
    """All action-log lines of a page (HTML or raw wiki text), in page order.
    'HH:MM <code> …' lines keep their code; 'HH:MM <initials> CHANGE: …' lines without a code get
    one from infer_code (the original line is kept, so the dialog shows what was read)."""
    text = strip_html(text)
    events: List[WikiEvent] = []
    for line in text.splitlines():
        m = _LINE.match(line)
        if m:
            hh, mm, _ss, code, rest = m.groups()
            desc = re.sub(r"\s{2,}", " ", rest).strip(" :|")
            events.append(WikiEvent(f"{int(hh):02d}:{mm}", code, desc, line.strip()))
            continue
        m = _FREE.match(line)
        if not m or not re.match(r"^\s*[A-ZÄÖÜa-zäöü]", m.group(4)):      # a time followed by words, not by numbers / another time
            continue
        hh, mm, _ss, rest = m.groups()
        if int(hh) > 23 or int(mm) > 59:
            continue
        code = infer_code(rest) or "-"
        desc = re.sub(r"\s{2,}", " ", rest).strip(" :|")
        events.append(WikiEvent(f"{int(hh):02d}:{mm}", code, desc, line.strip()))
    return events


def _hhmm_to_min(t: str) -> int:
    h, m = t.split(":")
    return int(h) * 60 + int(m)


def _num(s: str) -> float:
    return float(s.replace(",", "."))


_MULT_VOL = re.compile(r"(?:(\d+)\s*[x×]\s*)?(\d+(?:[.,]\d+)?)\s*(µl|μl|ul|microliter|mikroliter|microl|mikrol|ml|milliliter|nl)\b", re.IGNORECASE)
_FLOW = re.compile(r"(\d+(?:[.,]\d+)?)\s*(sccm|slpm|slm|l/min|ml/min)\b", re.IGNORECASE)
_FC = re.compile(r"\bFC\s*(\d+(?:[.,]\d+)?)\s*l\b\s*(\d+(?:[.,]\d+)?)\s*%", re.IGNORECASE)      # FC 0.5l 62.5% → 0.5 L/min full scale at 62.5 %
_BOTTLE = re.compile(r"\((\d+(?:[.,]\d+)?)\s*(ppm|%)\)", re.IGNORECASE)                            # Propene (1%), NO (3835ppm)
_PPBC = re.compile(r"(\d+(?:[.,]\d+)?)\s*ppbc\b", re.IGNORECASE)
_PPM_BOTTLE_WORDS = ("bottle", "cylinder", "standard", "flasche", "gas mix")


def parse_amounts(desc: str) -> Dict[str, Optional[float]]:
    """Every quantity the operator wrote: stated ppb, liquid µl (all terms summed, '2x3.2 ml' honoured),
    gas ml, MFC flow (sccm), duration (s), bottle ppm.  Missing → None."""
    low = desc.lower()
    out: Dict[str, Optional[float]] = {"ppb": None, "ul": None, "ml": None, "sccm": None, "duration_s": None, "bottle_ppm": None, "ppbc": None}
    m = re.search(r"(\d+(?:[.,]\d+)?)\s*ppbv?\b", low)
    if m:
        out["ppb"] = _num(m.group(1))
    mb = _BOTTLE.search(desc)
    if mb:                                      # "(1%)" or "(3835ppm)" right after the compound = bottle concentration
        v = _num(mb.group(1)); out["bottle_ppm"] = v * 1e4 if mb.group(2) == "%" else v
    m = re.search(r"(\d+(?:[.,]\d+)?)\s*ppmv?\b", low)
    if m and not mb:
        if any(w in low for w in _PPM_BOTTLE_WORDS):
            out["bottle_ppm"] = _num(m.group(1))
        elif out["ppb"] is None:
            out["ppb"] = _num(m.group(1)) * 1000.0
    mp = _PPBC.search(desc)
    out["ppbc"] = _num(mp.group(1)) if mp else None
    ul = ml = 0.0
    for n, v, unit in _MULT_VOL.findall(desc):
        q = (int(n) if n else 1) * _num(v)
        u = unit.lower()
        if u in ("ml", "milliliter"):
            ml += q
        elif u == "nl":
            ul += q * 1e-3
        else:
            ul += q
    out["ul"] = ul or None
    out["ml"] = ml or None
    m = _FLOW.search(desc)
    if m:
        v, u = _num(m.group(1)), m.group(2).lower()
        out["sccm"] = v * (1000.0 if u in ("slpm", "slm", "l/min") else 1.0)
    mf = _FC.search(desc)
    if mf:                                      # full scale [L/min] × setting [%]
        out["sccm"] = _num(mf.group(1)) * 1000.0 * _num(mf.group(2)) / 100.0
    d = parse_duration_min(desc)
    out["duration_s"] = d * 60.0 if d is not None else None
    return out


def amount_for_compound(desc: str, rc: "ic.ResolvedCompound", T_K: float, P_hPa: float, ch: "ic.ChamberProfile"
                        ) -> Tuple[Optional[float], Dict[str, str], str, List[str]]:
    """(ppb or None, delivery dict, explanation, alerts) from the wiki text for one compound.
    A stated ppb wins; otherwise liquid µl, gas ml (pure bottle assumed), or MFC × time."""
    a = parse_amounts(desc)
    alerts: List[str] = []
    calc: Optional[Tuple[float, str]] = None
    delivery: Dict[str, str] = {}
    if a["ul"] and rc.mw and rc.density:
        calc = ic.ppb_liquid(a["ul"], rc.mw, rc.density, T_K, P_hPa, ch)
        delivery = {"mode": "liquid", "vol_ul": f"{a['ul']:g}", "mw": f"{rc.mw:g}", "density": f"{rc.density:g}"}
    elif a["ml"] and rc.mw:
        if rc.density:
            calc = ic.ppb_liquid(a["ml"] * 1000.0, rc.mw, rc.density, T_K, P_hPa, ch)
            delivery = {"mode": "liquid", "vol_ul": f"{a['ml'] * 1000.0:g}", "mw": f"{rc.mw:g}", "density": f"{rc.density:g}"}
        else:
            frac = (a["bottle_ppm"] or 0) * 1e-6 or 1.0
            calc = ic.ppb_gas_syringe(a["ml"], frac, T_K, P_hPa, ch)
            delivery = {"mode": "gas_syringe", "vol_ml": f"{a['ml']:g}", "bottle_value": f"{a['bottle_ppm'] or 1:g}",
                        "bottle_unit": "ppm" if a["bottle_ppm"] else "pure"}
            if not a["bottle_ppm"]:
                alerts.append("gas volume given without bottle concentration: pure gas assumed")
    elif a["sccm"] and a["duration_s"]:
        if a["bottle_ppm"]:
            calc = ic.ppb_mfc(a["sccm"], a["duration_s"], a["bottle_ppm"] * 1e-6, T_K, P_hPa, ch)
        delivery = {"mode": "mfc", "flow_sccm": f"{a['sccm']:g}", "duration_s": f"{a['duration_s']:g}",
                    "bottle_value": f"{a['bottle_ppm'] or ''}", "bottle_unit": "ppm"}
        if not a["bottle_ppm"]:
            alerts.append("MFC flow and time given but no bottle concentration: enter it in the amount calculator")
    elif (a["ul"] or a["ml"]) and not (rc.mw and rc.density):
        alerts.append("volume given but molar mass/density unknown: enter them in the amount calculator")
    if a["ppb"] is not None:                 # a stated value wins; the calculation is only a cross-check
        expl = f"{a['ppb']:g} ppb as stated in the wiki"
        alerts = []
        if calc:
            expl += f" (calculated {calc[0]:.3g} ppb: {calc[1]})"
            if abs(calc[0] - a["ppb"]) > 0.15 * max(a["ppb"], 1e-9):
                alerts.append(f"stated {a['ppb']:g} ppb vs calculated {calc[0]:.3g} ppb differ by more than 15 %")
        return a["ppb"], {"mode": "ppb", "value": f"{a['ppb']:g}", **({"calc": calc[1]} if calc else {})}, expl, alerts
    if calc:
        return calc[0], delivery, calc[1], alerts
    return None, delivery, "no usable amount in the wiki line: enter ppb", alerts + ["no amount found"]


def mixture_label(desc: str) -> str:
    """Name of an unknown mixture from the line: the words between 'Inject' and the first
    amount/'for', e.g. 'Phenolic + Furanoid Mix 43ul for 234ppbc' → 'Phenolic_Furanoid_Mix'."""
    txt = desc.split("CHANGE:")[-1]
    txt = re.sub(r"^\s*[A-ZÄÖÜ][a-z]?\.\s*[A-ZÄÖÜ][a-z]?\.?\s*", "", txt)      # operator initials "R.W." / "Ö.E."
    txt = re.sub(r"(?i)\binject(?:ion|ed)?\b(?:\s+of)?", " ", txt)
    txt = re.split(r"\s+for\s+|\d+(?:[.,]\d+)?\s*(?:µl|μl|ul|ml)\b", txt, maxsplit=1)[0]
    words = re.findall(r"[A-Za-z][A-Za-z0-9]*", txt)
    return "_".join(words)[:28] if words else ""


def identify_mixture_word(desc: str) -> str:
    low = desc.lower()
    for word, label in MIXTURE_WORDS.items():
        if re.search(rf"\b{word}\b", low):
            return label
    return ""


def events_to_import(events: List[WikiEvent], date: str = "", T_K: float = 298.0, P_hPa: float = 1013.0,
                     chamber_name: str = "SAPHIR", mixtures: Optional[Dict[str, "ic.Mixture"]] = None) -> WikiImport:
    """Roof schedule, injections (with action types, amounts and delivery info) and MCM
    precursors from the action log.  Mixture lines ("VCP solution (180 uL)") are expanded
    with the recipe from mixtures/*.mixture.json when the name or an alias matches."""
    imp = WikiImport(date=date)
    ch = ic.chamber(chamber_name)
    if mixtures is None:
        mixtures = ic.load_mixtures()
    opens: List[str] = []
    pending: Dict[str, Tuple[str, WikiEvent]] = {}     # slow action (humidification/seed) → (start, event)
    counts: Dict[str, int] = {}

    def unique(base: str) -> Tuple[str, int]:
        counts[base] = counts.get(base, 0) + 1
        k = counts[base]
        return (base if k == 1 else f"{base}_{k}"), k

    def add_event(action: str, name: str, products: str, ppb: str, start: str, dur_min: float, comment: str,
                  extra: Optional[List[Tuple[str, str]]] = None, delivery: Optional[Dict[str, str]] = None):
        d = ACTIONS[action]
        uname, k = unique(name)
        col = d["column"] or f"Q_{name}"
        sw = d["switch"] or ""
        if k > 1:                                  # repeated injection of the same thing → own column/switch
            col = f"{col}_{k}"
            sw = f"{sw}_{k}" if sw else ""
        comment = comment.encode("ascii", "replace").decode()        # EASY comments must stay ASCII
        ev = InjectionEvent(uname, products, ppb, start, dur_min, col, sw,
                            extra_channels=list(extra or []), comment=comment, action=action,
                            delivery=dict(delivery or {}))
        imp.injections.append(ev)
        return ev

    for ev in events:
        if not ev.use:
            continue
        code, t, desc = ev.code, ev.time, ev.description
        if code == "R":
            opens.append(t)
        elif code == "r":
            if opens:
                imp.roof_open.append((opens.pop(0), t))
            else:
                imp.checks.append(("warning", f"{t}: roof closed without a preceding open"))
        elif code in ("H", "S") and is_stop_line(desc) and WIKI_CODE_ACTION[code] in pending:
            a = WIKI_CODE_ACTION[code]                 # operators often log the stop with the same code
            st, sev = pending.pop(a)
            _close_slow(imp, add_event, a, st, t, f"wiki {st}-{t}: {sev.description[:60]}")
        elif code in ("H", "S"):
            a = WIKI_CODE_ACTION[code]
            if a in pending:
                st, _ = pending.pop(a)
                imp.checks.append(("warning", f"{t}: {ACTIONS[a]['label']} started twice ({st} and {t}) — first one has no stop"))
                _close_slow(imp, add_event, a, st, "", f"wiki: {desc[:60]} (no stop found)")
            pending[a] = (t, ev)
        elif code in ("h", "s"):
            a = WIKI_CODE_ACTION[code]
            if a in pending:
                st, sev = pending.pop(a)
                _close_slow(imp, add_event, a, st, t, f"wiki {st}-{t}: {sev.description[:60]}")
            else:
                imp.checks.append(("warning", f"{t}: {ACTIONS[a]['label']} stopped without a start"))
        elif code == "X":
            dur = parse_duration_min(desc) or FAST_DURATION_MIN
            mix = ic.find_mixture(desc, mixtures)
            if mix is not None:
                a = parse_amounts(desc)
                mev, precs, alerts = ic.mixture_event(mix, a["ul"], t, dur, T_K, P_hPa, ch)
                mev.comment = f"wiki: {desc[:60]} | " + mev.comment
                add_event("mixture", mev.name, mev.products, mev.total_ppb, t, dur, mev.comment, delivery=mev.delivery)
                for sp in precs:
                    if sp not in imp.precursors:
                        imp.precursors.append(sp)
                for lvl, msg in alerts:
                    imp.checks.append((lvl, f"{t}: {msg}"))
                if not a["ul"]:
                    imp.checks.append(("warning", f"{t}: mixture {mix.name} without injected volume: recipe amounts used as is"))
                continue
            a0 = parse_amounts(desc)
            if a0["ppbc"] and ("mix" in desc.lower() or "solution" in desc.lower()):
                label = mixture_label(desc) or identify_mixture_word(desc) or "MIX"
                mk = _KOH.search(desc)
                extra = f"; OH reactivity {mk.group(1)} s-1" if mk else ""
                add_event("mixture", label, "", "", t, dur,
                          f"wiki: {desc[:60]} | {a0['ppbc']:g} ppbC stated{extra}; no recipe '{label}' in mixtures/ — define it to model the components")
                imp.checks.append(("warning", f"{t}: '{label}' ({a0['ppbc']:g} ppbC) has no recipe file — recorded only until the mixture is defined"))
                continue
            rc = ic.resolve_compound(desc)
            if rc.status == "mcm":
                ppb, delivery, expl, alerts = amount_for_compound(desc, rc, T_K, P_hPa, ch)
                if rc.mcm not in imp.precursors:
                    imp.precursors.append(rc.mcm)
                add_event("voc", rc.mcm, rc.mcm, f"{ppb:.4g}" if ppb is not None else "0", t, dur,
                          f"wiki: {desc[:60]} | {expl}", delivery=delivery)
                for msg in alerts:
                    imp.checks.append(("warning", f"{t}: {rc.mcm} — {msg}"))
            elif rc.status == "no_mcm":
                ppb, delivery, expl, alerts = amount_for_compound(desc, rc, T_K, P_hPa, ch)
                label = ic.safe_name(rc.note.split("'")[1] if "'" in rc.note else "compound")
                add_event("voc", label, "", f"{ppb:.4g}" if ppb is not None else "0", t, dur,
                          f"wiki: {desc[:60]} | {expl} | NOT in MCM: give a surrogate species as product", delivery=delivery)
                imp.checks.append(("warning", f"{t}: '{label}' has no MCM species — recorded only until a surrogate is set as product"))
            else:
                label = identify_mixture_word(desc) or "VOC"
                add_event("mixture", label, "", "", t, dur,
                          f"wiki: {desc[:60]} (recorded only; define the mixture or map to MCM surrogates)")
                imp.checks.append(("warning", f"{t}: '{desc[:40]}' not recognised (no MCM species, no mixture file) — recorded only"))
        elif code in WIKI_CODE_ACTION:                  # O M D c C
            a = WIKI_CODE_ACTION[code]
            d = ACTIONS[a]
            rc = ic.resolve_compound(d["products"].split("+")[0])
            ppb, delivery, expl, alerts = amount_for_compound(desc, rc, T_K, P_hPa, ch)
            dur = parse_duration_min(desc) or FAST_DURATION_MIN
            add_event(a, d["name"], d["products"], f"{ppb:.4g}" if ppb is not None else d["total_ppb"], t, dur,
                      f"wiki: {desc[:60]} | {expl}", delivery=delivery)
            for msg in alerts:
                imp.checks.append(("warning", f"{t}: {d['label']} — {msg}" + ("" if ppb is not None else f" (default {d['total_ppb']} ppb used)")))
        else:
            imp.notes.append(f"{t} {code} {desc}")
    for a, (st, sev) in pending.items():
        imp.checks.append(("warning", f"{st}: {ACTIONS[a]['label']} has no stop entry — duration set to default"))
        _close_slow(imp, add_event, a, st, "", f"wiki: {sev.description[:60]} (no stop found)")
    for o in opens:
        imp.checks.append(("warning", f"{o}: roof opened without a close entry — closes at the end of the window"))
        imp.roof_open.append((o, "23:59"))
    # no roof entries = a dark experiment: not a problem (the experiment check notes it as such)
    if not imp.injections:
        imp.checks.append(("warning", "no injection entries found"))
    imp.injections.sort(key=lambda e: _hhmm_to_min(e.start) if e.start else 0)
    return imp


def _close_slow(imp: WikiImport, add_event, action: str, start: str, stop: str, comment: str):
    d = ACTIONS[action]
    dur = (_hhmm_to_min(stop) - _hhmm_to_min(start)) if (start and stop) else d["duration_min"]
    if dur <= 0:
        imp.checks.append(("warning", f"{start}: {d['label']} stop is not after start — default duration used"))
        dur = d["duration_min"]
    extra = list(H2O_IMPURITY_DEFAULTS[1:]) if action == "humidification" else []
    add_event(action, d["name"], d["products"], d["total_ppb"], start, float(dur), comment, extra)


def archive_day_folder(archive_root: str, date: str) -> str:
    """Z:\\IEK8-SAPHIR\\data\\YYYY\\MM\\DD for the date ('' if it does not exist)."""
    try:
        y, m, d = date.split("-")
        p = os.path.join(archive_root, y, m, d)
        return p if os.path.isdir(p) else ""
    except Exception:
        return ""


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("usage: python saphir_wiki.py <YYYY-MM-DD | saved page .html/.txt> [user password]")
        print("       login defaults: environment variables SAPHIR_WIKI_USER / SAPHIR_WIKI_PASS"); sys.exit(1)
    src = sys.argv[1]
    text = load_wiki_file(src) if os.path.isfile(src) else fetch_wiki_text(src, *(sys.argv[2:4]))
    evs = parse_action_log(text)
    for e in evs:
        print(f"{e.time} {e.code} [{e.action}] {e.description}")
    imp = events_to_import(evs, src if not os.path.isfile(src) else "")
    print("roof:", imp.roof_open)
    print("precursors:", imp.precursors)
    for inj in imp.injections:
        print(f"  {inj.action:14s} {inj.name:10s} {inj.products[:28]:28s} {inj.total_ppb:>8s} ppb  {inj.start} +{inj.duration_min:g} min  {inj.column}/{inj.switch}  [{inj.delivery.get('mode','')}] {inj.comment[:70]}")
    for lvl, msg in imp.checks:
        print(f"  [{lvl}] {msg}")
