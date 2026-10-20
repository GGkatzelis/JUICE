"""
easy_scene.py — read / model / write SAPHIR EASY master files (scene.easy)
==========================================================================

Goal 2 of JUICE: generate the scene.easy exactly as the IDL workflow does, so
that FACSIMILE and JUICE receive identical input.

Two levels:

1. **Document level** (`EasyDocument`, `parse_document`, `render_document`):
   every physical line is kept verbatim and classified.  Rendering a parsed
   document is byte-identical to the source — nothing is ever lost.

2. **Scene level** (`SceneSpec`, `extract_scene`, `write_scene`):
   the statements are sorted into the blocks of the IDL template
   (mcm331cv.pro / mcm_converter.py easy_temp):

       INCLUDES → WALLDEP → SOURCE → BGREAC → DIL → MEASURED → FILES/FACS
       → INITIALS → TRACER DILUTION → CHAMBER SOURCES → INJECTIONS → EXTRA

   `write_scene` emits the canonical IDL-style text with the usual banners
   and comments.  Everything a recogniser does not understand goes to
   `SceneSpec.extra` and is re-emitted verbatim, so the scene level is also
   lossless at the statement level.

Equivalence checks (`active_statements`, `statements_equivalent`):
   two EASY files are equivalent for FACSIMILE/JUICE when they contain the
   same multiset of active statements (comments, blank lines and whitespace
   ignored) and, for every definition that occurs more than once, the same
   last occurrence ("last wins").

CLI:
    python easy_scene.py scene.easy            # round-trip report
    python easy_scene.py scene.easy out.easy   # write canonical regeneration
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ─────────────────────────────────────────────────────────────────────────────
# Level 1 — document model
# ─────────────────────────────────────────────────────────────────────────────

_RX_INCLUDE  = re.compile(r"^#include\s+(\S+)\s*$", re.I)
_RX_FILES    = re.compile(r"^FILES\s*\[\s*(\w+)\s*\]\s*=\s*(.+?)\s*$", re.I)
_RX_FACS     = re.compile(r"^FACS\s*\[\s*(\w+)\s*\]\s*=\s*(.+?)\s*$", re.I)
_RX_INITIAL  = re.compile(r"^INITIAL\s*\[\s*([\w*]+)\s*\]\s*=\s*(.+?)\s*$", re.I)
_RX_RATE     = re.compile(r"^k\s*\[(.+?)\]\s*=\s*(.+?)\s*$", re.I)
_RX_ASSIGN   = re.compile(r"^([A-Za-z_][\w]*)\s*=\s*(.+?)\s*$")
_RX_REACTION = re.compile(r"^(.*?)(-/->|-->)(.*)$")


@dataclass
class Line:
    """One physical line of an EASY file."""
    raw: str                      # verbatim, without the trailing newline
    kind: str                     # blank|comment|banner|include|assign|reaction|rate|files|facs|initial|other
    code: str = ""                # statement text without comment / surrounding whitespace
    comment: str = ""             # text after ';' (without the ';'), '' if none
    fields: Tuple[str, ...] = ()  # kind-specific parsed fields (see _classify)


def _split_comment(raw: str) -> Tuple[str, str]:
    """Split a raw line into (code, comment).  EASY comments start with ';'."""
    if ";" in raw:
        code, comment = raw.split(";", 1)
        return code.strip(), comment.rstrip("\r\n")
    return raw.strip(), ""


def _classify(code: str) -> Tuple[str, Tuple[str, ...]]:
    if not code:
        return "blank", ()
    if set(code) <= {"*"}:
        return "banner", ()
    m = _RX_INCLUDE.match(code)
    if m:
        return "include", (m.group(1),)
    m = _RX_FILES.match(code)
    if m:
        return "files", (m.group(1).upper(), m.group(2))
    m = _RX_FACS.match(code)
    if m:
        return "facs", (m.group(1).upper(), m.group(2))
    m = _RX_INITIAL.match(code)
    if m:
        return "initial", (m.group(1), m.group(2))
    m = _RX_RATE.match(code)
    if m:
        return "rate", (m.group(1).strip(), m.group(2))
    m = _RX_REACTION.match(code)
    if m and "=" not in code.split(m.group(2))[0]:
        return "reaction", (m.group(1).strip(), m.group(3).strip(), m.group(2))
    m = _RX_ASSIGN.match(code)
    if m:
        return "assign", (m.group(1), m.group(2))
    return "other", (code,)


def parse_line(raw: str) -> Line:
    raw = raw.rstrip("\r\n")
    stripped = raw.strip()
    if stripped.startswith("﻿"):          # BOM on first line (scene_08)
        stripped = stripped.lstrip("﻿")
    if not stripped:
        return Line(raw, "blank")
    if stripped.startswith(";"):
        return Line(raw, "comment", comment=stripped[1:])
    code, comment = _split_comment(stripped)
    kind, fields = _classify(code)
    return Line(raw, kind, code, comment, fields)


@dataclass
class EasyDocument:
    path: str
    lines: List[Line]
    newline: str = "\n"
    bom: bool = False
    encoding: str = "latin-1"
    trailing_newline: bool = True

    @property
    def statements(self) -> List[Line]:
        return [l for l in self.lines if l.kind not in ("blank", "comment", "banner")]


def parse_document(path: str) -> EasyDocument:
    with open(path, "rb") as fh:
        data = fh.read()
    bom = data.startswith(b"\xef\xbb\xbf")
    if bom:
        encoding, text = "utf-8", data[3:].decode("utf-8")
    else:
        try:
            encoding, text = "utf-8", data.decode("utf-8")
        except UnicodeDecodeError:
            encoding, text = "latin-1", data.decode("latin-1")
    newline = "\r\n" if "\r\n" in text else "\n"
    trailing = text.endswith(newline)
    lines = [parse_line(l) for l in text.split(newline)]
    if lines and lines[-1].raw == "" and trailing:
        lines.pop()                           # trailing newline, not an extra blank line
    return EasyDocument(path, lines, newline, bom, encoding, trailing)


def render_document(doc: EasyDocument) -> bytes:
    """Byte-identical reproduction of a parsed document."""
    text = doc.newline.join(l.raw for l in doc.lines) + (doc.newline if doc.trailing_newline else "")
    return (b"\xef\xbb\xbf" if doc.bom else b"") + text.encode(doc.encoding)


# ─────────────────────────────────────────────────────────────────────────────
# Statement normalisation / equivalence
# ─────────────────────────────────────────────────────────────────────────────

def normalise(code: str) -> str:
    """Canonical form for comparing statements: no whitespace, upper-case."""
    return re.sub(r"\s+", "", code).upper()


def _definition_key(line: Line) -> Optional[str]:
    """Key under which FACSIMILE/JUICE apply 'last definition wins'."""
    if line.kind == "assign":
        return "A:" + line.fields[0].upper()
    if line.kind == "rate":
        return "K:" + normalise(line.fields[0])
    if line.kind == "initial":
        return "I:" + line.fields[0].upper()
    if line.kind in ("files", "facs"):
        return line.kind.upper() + ":" + line.fields[0]
    return None


def _is_noise(l: Line) -> bool:
    """Statements without ASCII content (e.g. a UTF-8 BOM read as latin-1 by an
    upstream tool) are ignored by FACSIMILE/JUICE and are not compared."""
    return not l.code or ord(l.code[0]) >= 0x7f or not re.search(r"[A-Za-z0-9#\[]", l.code)


def active_statements(path_or_doc) -> List[str]:
    doc = path_or_doc if isinstance(path_or_doc, EasyDocument) else parse_document(path_or_doc)
    return [normalise(l.code) for l in doc.statements if not _is_noise(l)]


def statements_equivalent(a, b) -> Tuple[bool, List[str]]:
    """Compare two EASY files/documents at the active-statement level.

    Returns (equivalent, list_of_differences).  Equivalent means: identical
    multiset of normalised statements and identical last definition for every
    repeated definition key.
    """
    da = a if isinstance(a, EasyDocument) else parse_document(a)
    db = b if isinstance(b, EasyDocument) else parse_document(b)
    diffs: List[str] = []
    from collections import Counter
    # Definitions (assign / rate / initial / files / facs) may be repeated without
    # effect ("last wins"), so they are compared as sets; reactions, includes and
    # other statements are compared with multiplicity.
    def _parts(doc):
        defs, rest = set(), Counter()
        for l in doc.statements:
            if _is_noise(l):
                continue
            if _definition_key(l):
                defs.add(normalise(l.code))
            else:
                rest[normalise(l.code)] += 1
        return defs, rest
    sa, ra = _parts(da); sb, rb = _parts(db)
    for s in sorted(sa ^ sb):
        diffs.append(f"definition only in {'A' if s in sa else 'B'}: {s[:100]}")
    for s in sorted(set(ra) | set(rb)):
        if ra[s] != rb[s]:
            diffs.append(f"count {ra[s]} vs {rb[s]}: {s[:100]}")
    last_a = {_definition_key(l): normalise(l.code) for l in da.statements if _definition_key(l)}
    last_b = {_definition_key(l): normalise(l.code) for l in db.statements if _definition_key(l)}
    for k in sorted(set(last_a) | set(last_b)):
        if last_a.get(k) != last_b.get(k):
            diffs.append(f"last-wins differs for {k}: {last_a.get(k)} vs {last_b.get(k)}")
    return (not diffs), diffs


# ─────────────────────────────────────────────────────────────────────────────
# Level 2 — scene model
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Stmt:
    """A statement kept verbatim (code + optional trailing comment)."""
    code: str
    comment: str = ""

    @classmethod
    def from_line(cls, l: Line) -> "Stmt":
        return cls(l.code, l.comment)

    def render(self) -> str:
        return self.code + (f"   ;{self.comment}" if self.comment else "")


@dataclass
class ReactionDef:
    """A reaction line plus (optionally) its k[...] rate line."""
    lhs: str
    rhs: str
    rate: Optional[str] = None          # expression inside CONST( ) or full RHS of k[]=
    arrow: str = "-->"
    comment: str = ""
    rate_comment: str = ""
    rate_raw: Optional[str] = None      # full text after '=' in the k[] line (verbatim)
    text_raw: Optional[str] = None      # reaction text exactly as written (spacing kept)

    @property
    def text(self) -> str:
        if self.text_raw is not None:
            return self.text_raw
        return f"{self.lhs} {self.arrow} {self.rhs}".strip()

    def render(self, rxn_text: Optional[str] = None) -> List[str]:
        out = []
        rt = rxn_text if rxn_text is not None else self.text
        out.append(rt + (f"   ;{self.comment}" if self.comment else ""))
        if self.rate_raw is not None:
            out.append(f"k[{rt}]={self.rate_raw}" + (f"   ;{self.rate_comment}" if self.rate_comment else ""))
        return out


@dataclass
class Setting:
    """NAME = CONST(expr) | INPUT(col) | FUNCTION[...] | CALL(...)  (kept verbatim)."""
    name: str
    rhs: str                             # verbatim right-hand side
    comment: str = ""

    def render(self) -> str:
        return f"{self.name} = {self.rhs}" + (f"   ;{self.comment}" if self.comment else "")


@dataclass
class Injection:
    """One ENZ-switched source: `switch = input(COL)` + channel reactions."""
    switch: str                          # e.g. qNO, Q_ethanol, qVOC
    enz_column: str                      # e.g. QNO, qVCP2
    channels: List[ReactionDef] = field(default_factory=list)   # switch --> products
    chain: List[ReactionDef] = field(default_factory=list)      # follow-up reactions (VOCin1 --> ...)
    comment: str = ""                    # comment on the input() line
    pre_comments: List[str] = field(default_factory=list)       # comment lines directly above the block


@dataclass
class SceneSpec:
    title: str = ""
    header_comments: List[str] = field(default_factory=list)
    # INCLUDES
    includes_pre: List[Stmt] = field(default_factory=list)      # e.g. default_minimum.incl
    anorg_includes: List[Stmt] = field(default_factory=list)
    mech_stem: str = ""                                          # e.g. ethanol.fac
    mech_includes: List[Stmt] = field(default_factory=list)     # <stem>.saphir.incl/.ro2.incl/.mech/.rate
    j_clear_sky: Optional[Stmt] = None
    includes_other: List[Stmt] = field(default_factory=list)
    # WALLDEP
    walldep_include: Optional[Stmt] = None
    walldep: Optional[Setting] = None
    # SOURCE
    wvc2rh: Optional[Setting] = None
    rh: Optional[Setting] = None
    hono_source: Optional[ReactionDef] = None
    hono_const: Optional[Setting] = None
    hcho_source: Optional[ReactionDef] = None
    hcho_const: Optional[Setting] = None
    # BGREAC
    y_reaction: Optional[ReactionDef] = None
    y_const: Optional[Setting] = None
    # DIL
    dil_includes: List[Stmt] = field(default_factory=list)
    dilute: Optional[Setting] = None
    dil_flow: Optional[Setting] = None                           # FLOW = CONST(6.) inside the DIL block
    dil: Optional[Setting] = None
    # MEASURED
    jmeasured_include: Optional[Stmt] = None
    j_inputs: List[Setting] = field(default_factory=list)       # jno2 = INPUT(jno2) / PHTFKT(...)
    j_overrides: List[ReactionDef] = field(default_factory=list) # k[NO2 + HV --> ...] = CONST(jno2)  (rate-only)
    settings: List[Setting] = field(default_factory=list)       # T_in, T, P, SZA, SZ, HV_in, HV, FLOW, H2O_in, H2O
    # FILES / FACS
    files: List[Tuple[str, str]] = field(default_factory=list)
    facs: List[Tuple[str, str]] = field(default_factory=list)
    # INITIALS
    initials: List[Tuple[str, str, str]] = field(default_factory=list)   # (species, expr, comment)
    # TRACER DILUTION  (XNOy+DIL-->, NA+DIL-->)
    tracer_dilution: List[ReactionDef] = field(default_factory=list)
    # CHAMBER SOURCES: dark HONO, SNO2 light source, mechanism overrides
    dark_hono: Optional[ReactionDef] = None
    light_sources: List[ReactionDef] = field(default_factory=list)       # SNO2 --> NO2 (+XNOy)
    # INJECTIONS
    injections: List[Injection] = field(default_factory=list)
    # everything else, verbatim, in original order
    extra: List[Stmt] = field(default_factory=list)
    extra_reactions: List[ReactionDef] = field(default_factory=list)
    # Encounter-ordered views used by the writer so that the regenerated file
    # keeps the original statement order where it matters for downstream tools:
    #   measured_items: Setting | ReactionDef (rate-only J override) of the MEASURED block
    #   tail:           Injection | Stmt | ReactionDef | ('include', Stmt) after the fixed blocks
    measured_items: List[object] = field(default_factory=list)
    tail: List[object] = field(default_factory=list)


_J_INPUT_NAMES = {"JNO2", "JO1D", "JHCHO_M", "JHCHO_R", "JH2O2", "JHONO", "JHNO3", "JNO3"}
_SETTING_NAMES = {"T_IN", "T", "P", "P_IN", "SZA", "SZ", "HV_IN", "HV", "FLOW", "H2O_IN", "H2O",
                  "TEMP", "PRESS"}


def _rate_parts(rhs: str) -> str:
    return rhs.strip()


def extract_scene(doc: EasyDocument) -> SceneSpec:
    """Classify the statements of a document into SceneSpec blocks."""
    spec = SceneSpec()
    lines = doc.lines

    # ── header: TITLE and leading comments ──────────────────────────────────
    for l in lines:
        if l.kind == "comment":
            c = l.comment
            if c.upper().startswith("TITLE=") and not spec.title:
                spec.title = c[len("TITLE="):]
            elif c.strip().upper().startswith("INCLUDES"):
                continue                                  # written by write_scene itself
            else:
                spec.header_comments.append(c)
        elif l.kind in ("blank", "banner"):
            continue
        else:
            break

    # ── pair reactions with their k[] lines ─────────────────────────────────
    stmts = doc.statements
    # k[] lines per reaction key, in file order; each reaction occurrence takes
    # the next unconsumed k[] line with the same key (duplicated blocks stay paired).
    rate_for: Dict[str, List[Line]] = {}
    for l in stmts:
        if l.kind == "rate":
            rate_for.setdefault(normalise(l.fields[0]), []).append(l)
    consumed_rates = set()

    def take_rate(rxn_text: str) -> Tuple[Optional[str], str]:
        for rl in rate_for.get(normalise(rxn_text), []):
            if id(rl) not in consumed_rates:
                consumed_rates.add(id(rl))
                return rl.fields[1], rl.comment
        return None, ""

    def mk_reaction(l: Line) -> ReactionDef:
        lhs, rhs, arrow = l.fields
        rate_raw, rc = take_rate(l.code)
        rd = ReactionDef(lhs=lhs, rhs=rhs, arrow=arrow, comment=l.comment,
                         rate_raw=rate_raw, rate_comment=rc, text_raw=l.code)
        rd.rate = rate_raw
        return rd

    # pending comments directly above a statement (used for injection blocks)
    pending_comments: List[str] = []

    # switch names (X = input(COL)) — collected in a first pass
    switches: Dict[str, Tuple[str, str]] = {}      # upper name → (orig name, column)
    for l in stmts:
        if l.kind == "assign":
            m = re.match(r"INPUT\s*\(\s*(\w+)\s*\)\s*$", l.fields[1], re.I)
            if m and l.fields[0].upper() not in _SETTING_NAMES and l.fields[0].upper() not in _J_INPUT_NAMES \
                    and l.fields[0].upper() != "DIL":
                switches[l.fields[0].upper()] = (l.fields[0], m.group(1))
    injections: Dict[str, Injection] = {}
    chain_species: set = set()

    seen_mech_stem = ""
    for l in lines:
        if l.kind == "comment":
            pending_comments.append(l.comment)
            continue
        if l.kind in ("blank", "banner"):
            pending_comments = []
            continue
        if l.kind == "include":
            name = l.fields[0]; low = name.lower(); st = Stmt.from_line(l)
            if low.startswith("anorg.") and not low.endswith((".dil.mech", ".dil.rate")):
                spec.anorg_includes.append(st)
            elif low == "j_clear_sky.incl":
                spec.j_clear_sky = st
            elif low == "saphir.walldep.incl":
                spec.walldep_include = st
            elif low.startswith("saphir.jmeasured"):
                spec.jmeasured_include = st
            elif low.endswith((".dil.mech", ".dil.rate")):
                spec.dil_includes.append(st)
            elif re.search(r"\.fac\.(saphir\.incl|ro2\.incl|mech|rate)$", low):
                stem = re.sub(r"\.(saphir\.incl|ro2\.incl|mech|rate)$", "", name, flags=re.I)
                seen_mech_stem = seen_mech_stem or stem
                spec.mech_stem = seen_mech_stem
                spec.mech_includes.append(st)
            elif not spec.anorg_includes and not spec.mech_includes:
                spec.includes_pre.append(st)
            elif spec.files or spec.facs or spec.tail:
                # include after the fixed blocks (e.g. VOC_Quellen_CHANEL2024.incl
                # right before the injections): keep its position in the tail
                spec.tail.append(("include", st))
            else:
                spec.includes_other.append(st)
            pending_comments = []
            continue
        if l.kind == "files":
            spec.files.append((l.fields[0], l.fields[1])); pending_comments = []; continue
        if l.kind == "facs":
            spec.facs.append((l.fields[0], l.fields[1])); pending_comments = []; continue
        if l.kind == "initial":
            spec.initials.append((l.fields[0], l.fields[1], l.comment)); pending_comments = []; continue
        if l.kind == "rate":
            if id(l) in consumed_rates:
                continue
            # rate-only statements: J overrides on mechanism reactions (k[NO2 + HV --> ...] = CONST(jno2))
            rxn = l.fields[0]
            rd = ReactionDef(lhs=rxn.split("-->")[0].strip(), rhs=rxn.split("-->", 1)[1].strip() if "-->" in rxn else "",
                             rate_raw=l.fields[1], rate_comment=l.comment)
            rd.rate = l.fields[1]
            rd._rate_only_text = rxn            # type: ignore[attr-defined]
            if "HV" in rxn.upper() and re.match(r"CONST\s*\(\s*J\w*\s*\)", l.fields[1], re.I):
                spec.j_overrides.append(rd); spec.measured_items.append(rd)
            else:
                spec.extra_reactions.append(rd); spec.tail.append(rd)
            consumed_rates.add(id(l)); pending_comments = []
            continue
        if l.kind == "assign":
            name, rhs = l.fields; up = name.upper(); s = Setting(name, rhs, l.comment)
            if up == "WALLDEP":
                spec.walldep = s
            elif up == "WVC2RH":
                spec.wvc2rh = s
            elif up == "RH":
                spec.rh = s
            elif up == "SHONO":
                spec.hono_const = s
            elif up in ("SHCHO", "QSHCHO"):
                spec.hcho_const = s
            elif up in ("Y", "QY"):
                spec.y_const = s
            elif up == "DILUTE":
                spec.dilute = s
            elif up == "DIL":
                spec.dil = s
            elif up == "FLOW" and spec.dilute is not None and not spec.jmeasured_include and spec.dil is None:
                spec.dil_flow = s
            elif up in _J_INPUT_NAMES or re.match(r"(INPUT|PHTFKT|PHTFKT_CS)\s*\(", rhs, re.I) and up.startswith("J"):
                spec.j_inputs.append(s); spec.measured_items.append(s)
            elif up in _SETTING_NAMES:
                spec.settings.append(s); spec.measured_items.append(s)
            elif up in switches:
                inj = injections.get(up)
                if inj is None:
                    inj = Injection(switch=name, enz_column=switches[up][1], comment=l.comment,
                                    pre_comments=list(pending_comments))
                    injections[up] = inj; spec.injections.append(inj); spec.tail.append(inj)
                else:
                    inj.comment = inj.comment or l.comment
            else:
                st = Stmt.from_line(l)
                spec.extra.append(st); spec.tail.append(st)
            pending_comments = []
            continue
        if l.kind == "reaction":
            rd = mk_reaction(l)
            lhs_up = normalise(rd.lhs); rhs_up = normalise(rd.rhs)
            lhs_sp = [p for p in re.split(r"\+", lhs_up) if p]
            if rd.arrow == "-/->":
                spec.extra_reactions.append(rd); spec.tail.append(rd)
            elif "SHONO" in lhs_sp:
                spec.hono_source = rd
            elif "SHCHO" in lhs_sp or "QSHCHO" in lhs_sp:
                spec.hcho_source = rd
            elif ("Y" in lhs_sp or "QY" in lhs_sp) and "OH" in lhs_sp:
                spec.y_reaction = rd
            elif "DIL" in lhs_sp and not rhs_up:
                spec.tracer_dilution.append(rd)
            elif lhs_up == "H2O" and rhs_up.startswith("HONO"):
                spec.dark_hono = rd
            elif lhs_up in ("SNO2", "SNO", "SOH", "SHO2"):
                spec.light_sources.append(rd)
            elif lhs_up in switches:
                inj = injections.get(lhs_up)
                if inj is None:
                    inj = Injection(switch=switches[lhs_up][0], enz_column=switches[lhs_up][1],
                                    pre_comments=list(pending_comments))
                    injections[lhs_up] = inj; spec.injections.append(inj); spec.tail.append(inj)
                inj.channels.append(rd)
                for term in re.split(r"\+", rhs_up):
                    sp = re.sub(r"^[\d.E+\-]*\*", "", term)
                    chain_species.add(sp)
            elif lhs_up in chain_species and spec.injections:
                spec.injections[-1].chain.append(rd)
                for term in re.split(r"\+", rhs_up):
                    chain_species.add(re.sub(r"^[\d.E+\-]*\*", "", term))
            else:
                spec.extra_reactions.append(rd); spec.tail.append(rd)
            pending_comments = []
            continue
        # other
        st = Stmt.from_line(l)
        spec.extra.append(st); spec.tail.append(st)
        pending_comments = []
    return spec


# ─────────────────────────────────────────────────────────────────────────────
# Canonical writer (IDL template style)
# ─────────────────────────────────────────────────────────────────────────────

_BANNER = "*********************************************************************"


def write_scene(spec: SceneSpec) -> str:
    out: List[str] = []
    P = out.append

    if spec.title:
        P(f";TITLE={spec.title}")
    for c in spec.header_comments:
        P(f";{c}")
    P("   ")
    P("; INCLUDES: ")
    for s in spec.includes_pre:
        P(s.render())
    for s in spec.anorg_includes:
        P(s.render())
    for s in spec.mech_includes:
        P(s.render())
    if spec.j_clear_sky:
        P(spec.j_clear_sky.render())
    for s in spec.includes_other:
        P(s.render())
    P(_BANNER)
    P(" ")

    if spec.walldep_include or spec.walldep:
        P(";/WALLDEP, Wall deposition in SAPHIR chamber")
        if spec.walldep_include:
            P(spec.walldep_include.render())
        if spec.walldep:
            P(spec.walldep.render())
        P(" ")
        P(_BANNER)
        P(" ")

    if any((spec.wvc2rh, spec.rh, spec.hono_source, spec.hcho_source)):
        P(";/SOURCE, HONO and HCHO source of SAPHIR chamber")
        P("")
        if spec.wvc2rh:
            P("; Funktion die die Wasserkonzentration in die relative")
            P("; Feuchtigkeit umrechnet")
            P(spec.wvc2rh.render())
            P("")
        if spec.rh:
            P("; Relative Feuchtigkeit")
            P(spec.rh.render())
            P("")
        if spec.hono_source:
            P("; HONO Quelle aus Paper Rohrer et al.")
            for s in spec.hono_source.render():
                P(s)
            if spec.hono_const:
                P(spec.hono_const.render())
            P("")
        if spec.hcho_source:
            P("; HCHO Quelle (lichtabhängig, photochemie unabhängig)")
            for s in spec.hcho_source.render():
                P(s)
            if spec.hcho_const:
                P(spec.hcho_const.render())
        P(" ")
        P(_BANNER)
        P(" ")

    if spec.y_reaction or spec.y_const:
        P(";/BGREAC, Background reactivity (see saphir.bgreac.incl) ")
        if spec.y_reaction:
            for s in spec.y_reaction.render():
                P(s)
        if spec.y_const:
            P(spec.y_const.render())
        P(" ")
        P(_BANNER)
        P(" ")

    if spec.dil_includes or spec.dilute or spec.dil:
        P("; DILUTION, /DIL")
        for s in spec.dil_includes:
            P(s.render())
        if spec.dilute:
            P(spec.dilute.render())
        if spec.dil_flow:
            P(spec.dil_flow.render())
        if spec.dil:
            P(spec.dil.render())
        P(" ")
        P(_BANNER)
        P(" ")

    P(";/MEASURED, use measured photolysis frequencies (given in ENZ-file)")
    P(";           and other measured parameters]")
    if spec.jmeasured_include:
        P(spec.jmeasured_include.render())
    # J inputs, their k[] overrides and the T/P/SZ/HV/FLOW/H2O settings in the
    # original encounter order (matters for the elog statement sequence).
    _in_settings = False
    for item in spec.measured_items:
        if isinstance(item, Setting):
            if item.name.upper() in _SETTING_NAMES and not _in_settings:
                P(""); _in_settings = True
            P(item.render())
        else:
            P(f"k[{item._rate_only_text}] = {item.rate_raw}" + (f"   ;{item.rate_comment}" if item.rate_comment else ""))
    P(_BANNER)
    P(" ")

    if spec.files:
        P("; FILES: ")
        for k, v in spec.files:
            P(f"FILES[{k}]   = {v}")
        P(" ")
    if spec.facs:
        P("; FACSIMILE PARAMETERS")
        for k, v in spec.facs:
            P(f"FACS[{k}]   = {v}")
        P(" ")

    if spec.initials:
        P("; INITIALS:")
        for sp, expr, c in spec.initials:
            P(f"INITIAL[{sp}] = {expr}" + (f"   ;{c}" if c else ""))
        P("")

    if spec.tracer_dilution:
        P("; Tracer dilution")
        for rd in spec.tracer_dilution:
            for s in rd.render():
                P(s)
            P("")

    if spec.dark_hono or spec.light_sources:
        P("; Chamber sources")
        if spec.dark_hono:
            for s in spec.dark_hono.render():
                P(s)
            P("")
        for rd in spec.light_sources:
            for s in rd.render():
                P(s)
            P("")

    if spec.tail:
        P(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;")
        P("; INJECTIONS (ENZ switch columns) and additional definitions")
        for item in spec.tail:
            if isinstance(item, Injection):
                for c in item.pre_comments:
                    P(f";{c}")
                P(f"{item.switch}=input({item.enz_column})" + (f"   ;{item.comment}" if item.comment else ""))
                for rd in item.channels:
                    for s in rd.render():
                        P(s)
                for rd in item.chain:
                    for s in rd.render():
                        P(s)
                P("")
            elif isinstance(item, tuple):               # ('include', Stmt)
                P(item[1].render())
            elif isinstance(item, Stmt):
                P(item.render())
            elif hasattr(item, "_rate_only_text"):
                P(f"k[{item._rate_only_text}] = {item.rate_raw}" + (f"   ;{item.rate_comment}" if item.rate_comment else ""))
            else:
                for s in item.render():
                    P(s)
        P(";;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;")
    return "\n".join(out) + "\n"


# ─────────────────────────────────────────────────────────────────────────────
# Round-trip report
# ─────────────────────────────────────────────────────────────────────────────

def roundtrip_report(path: str, out_path: Optional[str] = None) -> dict:
    doc = parse_document(path)
    byte_identical = render_document(doc) == open(path, "rb").read()
    spec = extract_scene(doc)
    text = write_scene(spec)
    tmp = out_path or (os.path.splitext(path)[0] + ".regen.easy")
    with open(tmp, "w", encoding=doc.encoding, newline="\n") as fh:
        fh.write(text)
    ok, diffs = statements_equivalent(doc, parse_document(tmp))
    n_extra = len(spec.extra) + len(spec.extra_reactions)
    return {"file": path, "level1_byte_identical": byte_identical,
            "statements": len(doc.statements), "unstructured": n_extra,
            "level2_equivalent": ok, "diffs": diffs, "regenerated": tmp,
            "injections": [(i.switch, i.enz_column, len(i.channels), len(i.chain)) for i in spec.injections],
            "extra": [s.code for s in spec.extra] + [r.lhs + r.arrow + r.rhs for r in spec.extra_reactions]}


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    rep = roundtrip_report(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else None)
    for k, v in rep.items():
        print(f"{k:22s}: {v}")
