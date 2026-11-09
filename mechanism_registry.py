"""
mechanism_registry.py — which chemistry, which version, from where.

Two kinds of chemistry, kept apart:

  mechanisms/base/<ID>/manifest.json      an MCM download exactly as published (never edited)
      files as downloaded (+ checksums), version, source URL and selection, download date,
      citation, status.  Today: MCM_v3.3.1.
      kind "gecko": a GECKO-A generator output converted once by gecko_import.py into EASY
      includes (organic chemistry, GECKO-A inorganic core, aliases, dilution pair, partitioning,
      species table) with the generator files' SHA-256 and the import options in the manifest.
      An experiment with such a base chooses the inorganic core and the partitioning include.

  mechanisms/addons/<ID>/manifest.json    chemistry added on top, as EASY include files
      version, status (draft | reviewed), compatible base IDs, the files, the species it adds
      (name + formula, so the GUI can resolve injections to them), the base species it needs
      (seed for the subset reachability), the base reactions it overrides, references.

Rules the registry enforces or reports:
  * a base file whose checksum differs from the manifest is reported (modified download);
  * an add-on is offered only for compatible bases;
  * add-ons never edit the base: they add reactions, override a rate with k[...] = CONST(...)
    after the base (last definition wins) or switch a reaction off with a zero rate.  Reactions
    are never deduplicated, so an add-on must not re-declare a base reaction with a bare line.

`provenance_lines()` gives the comment lines written into every generated EASY file, so a run
always states its chemistry.
"""
from __future__ import annotations

import hashlib
import json
import os
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional

HERE = os.path.dirname(os.path.abspath(__file__))
MECH_DIR = os.path.join(HERE, "mechanisms")
BASE_DIR = os.path.join(MECH_DIR, "base")
ADDON_DIR = os.path.join(MECH_DIR, "addons")
DEFAULT_BASE = "MCM_v3.3.1"


def sha256_of(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


@dataclass
class BaseMechanism:
    id: str
    name: str
    version: str
    folder: str
    files: Dict[str, str]                      # role → file name (fac, species_tsv, species_export_tsv, kpp, …)
    checksums: Dict[str, str] = field(default_factory=dict)   # file name → sha256
    source: Dict[str, str] = field(default_factory=dict)      # url, selection, downloaded, notes
    citation: str = ""
    status: str = "reviewed"
    notes: str = ""
    kind: str = "base"                         # "base" (an MCM download) | "gecko" (a converted GECKO-A generator output)
    gecko: Dict = field(default_factory=dict)  # kind gecko: precursor, primary_code, n_species, n_reactions_organic, ...

    @property
    def is_gecko(self) -> bool:
        return self.kind == "gecko"

    def path(self, role: str) -> str:
        return os.path.join(self.folder, self.files[role]) if role in self.files else ""

    @property
    def fac(self) -> str:
        return self.path("fac")

    def verify(self) -> List[str]:
        """Problems with the download: missing files or checksums that no longer match."""
        out = []
        for role, fn in self.files.items():
            p = os.path.join(self.folder, fn)
            if not os.path.isfile(p):
                out.append(f"{self.id}: {role} file missing ({fn})"); continue
            want = self.checksums.get(fn)
            if want and sha256_of(p) != want:
                out.append(f"{self.id}: {fn} differs from the recorded download (sha256 mismatch) — the base was modified")
        return out


@dataclass
class AddOn:
    id: str
    name: str
    version: str
    folder: str
    files: List[str]                           # EASY include files, in load order
    compatible_with: List[str] = field(default_factory=list)
    status: str = "draft"                      # draft | reviewed
    description: str = ""
    adds_species: Dict[str, str] = field(default_factory=dict)     # NAME → formula
    needs_base_species: List[str] = field(default_factory=list)    # must be in the subset (reachability seed)
    precursors: List[str] = field(default_factory=list)            # species of the add-on that can be injected
    overrides: List[str] = field(default_factory=list)             # base reactions it re-rates or switches off
    references: List[str] = field(default_factory=list)
    notes: str = ""

    def paths(self) -> List[str]:
        return [os.path.join(self.folder, f) for f in self.files]

    def compatible(self, base_id: str) -> bool:
        return not self.compatible_with or base_id in self.compatible_with


def _read_manifest(folder: str) -> Optional[dict]:
    p = os.path.join(folder, "manifest.json")
    if not os.path.isfile(p):
        return None
    with open(p, encoding="utf-8") as fh:
        return json.load(fh)


def load_bases(base_dir: str = "") -> Dict[str, BaseMechanism]:
    out: Dict[str, BaseMechanism] = {}
    root = base_dir or BASE_DIR
    if os.path.isdir(root):
        for d in sorted(os.listdir(root)):
            folder = os.path.join(root, d)
            m = _read_manifest(folder)
            if m and m.get("kind", "base") in ("base", "gecko"):
                keys = {k for k in BaseMechanism.__dataclass_fields__} - {"folder"}
                out[m["id"]] = BaseMechanism(folder=folder, **{k: v for k, v in m.items() if k in keys})
    return out


def load_addons() -> Dict[str, AddOn]:
    out: Dict[str, AddOn] = {}
    if os.path.isdir(ADDON_DIR):
        for d in sorted(os.listdir(ADDON_DIR)):
            folder = os.path.join(ADDON_DIR, d)
            m = _read_manifest(folder)
            if m and m.get("kind") == "addon":
                keys = {k for k in AddOn.__dataclass_fields__} - {"folder"}
                out[m["id"]] = AddOn(folder=folder, **{k: v for k, v in m.items() if k in keys})
    return out


def base(base_id: str = "") -> BaseMechanism:
    bases = load_bases()
    bid = base_id or DEFAULT_BASE
    if bid not in bases:
        if bases:
            raise KeyError(f"base mechanism '{bid}' not installed; available: {', '.join(bases)}")
        raise KeyError(f"no base mechanism installed under {BASE_DIR}")
    return bases[bid]


def is_gecko(base_id: str) -> bool:
    try:
        return base(base_id).is_gecko
    except KeyError:
        return False


def species_names(base_id: str) -> Dict[str, str]:
    """NAME -> formula of a GECKO-A base's species table ({} for an MCM base, whose species come from mcm_species_db)."""
    try:
        b = base(base_id)
    except KeyError:
        return {}
    p = b.path("species_tsv") if b.is_gecko else ""
    out: Dict[str, str] = {}
    if p and os.path.isfile(p):
        with open(p, encoding="utf-8", errors="replace") as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_name = header.index("Name") if "Name" in header else 0
            i_formula = header.index("Formula") if "Formula" in header else 4
            for line in fh:
                cols = line.rstrip("\n").split("\t")
                if len(cols) > max(i_name, i_formula) and cols[i_name]:
                    out[cols[i_name].upper()] = cols[i_formula]
    return out


def addons_for(base_id: str, ids: Optional[List[str]] = None) -> List[AddOn]:
    """Add-ons compatible with `base_id`; if `ids` is given, only those (in that order)."""
    all_ = load_addons()
    if ids is None:
        return [a for a in all_.values() if a.compatible(base_id)]
    out = []
    for i in ids:
        a = all_.get(i)
        if a is None:
            raise KeyError(f"add-on '{i}' not installed; available: {', '.join(all_) or 'none'}")
        if not a.compatible(base_id):
            raise ValueError(f"add-on '{i}' ({a.name}) is not compatible with {base_id}")
        out.append(a)
    return out


def addon_species(ids: List[str]) -> Dict[str, str]:
    """NAME → formula of every species the chosen add-ons add (for completer / resolver)."""
    out: Dict[str, str] = {}
    all_ = load_addons()
    for i in ids or []:
        if i in all_:
            out.update(all_[i].adds_species)
    return out


def provenance_lines(base_id: str, addon_ids: List[str]) -> List[str]:
    """Comment lines for the generated EASY file: what chemistry this run used."""
    b = base(base_id)
    if b.is_gecko:
        src = b.source.get("generator_files_sha256", {}) or {}
        g = b.gecko or {}
        lines = [f"; CHEMISTRY: base {b.id} ({b.name}; {b.version}) - GECKO-A generator output {b.source.get('folder', '?')}",
                 "; CHEMISTRY:   generator files sha256[:12] " + ", ".join(f"{k} {v[:12]}" for k, v in src.items()),
                 f"; CHEMISTRY:   {g.get('n_species', '?')} species, {g.get('n_reactions_inorganic', '?')} inorganic + "
                 f"{g.get('n_reactions_organic', '?')} organic reactions, {g.get('n_partitioning', '?')} partitioning species; "
                 f"primary {g.get('primary_code', '?')} = {g.get('precursor', '?')}; options {b.source.get('options', {})}"]
        return [l.encode("ascii", "replace").decode() for l in lines]
    sha = b.checksums.get(b.files.get("fac", ""), "")[:12]
    lines = [f"; CHEMISTRY: base {b.id} ({b.name}, {b.version}; fac sha256 {sha or 'n/a'}; "
             f"source {b.source.get('url', '?')}, downloaded {b.source.get('downloaded', '?')})"]
    for a in addons_for(base_id, addon_ids):
        lines.append(f"; CHEMISTRY: add-on {a.id} v{a.version} [{a.status}] - {a.name}: {', '.join(a.files)}")
    return [l.encode("ascii", "replace").decode() for l in lines]      # EASY files stay ASCII


def write_manifest(folder: str, data: dict) -> str:
    os.makedirs(folder, exist_ok=True)
    p = os.path.join(folder, "manifest.json")
    with open(p, "w", encoding="utf-8") as fh:
        json.dump(data, fh, indent=1, ensure_ascii=False)
    return p


if __name__ == "__main__":
    print("base mechanisms:")
    for b in load_bases().values():
        probs = b.verify()
        print(f"  {b.id:14s} {b.name} v{b.version} [{b.status}]  fac={os.path.basename(b.fac)}  " + ("OK" if not probs else "; ".join(probs)))
    print("add-ons:")
    for a in load_addons().values():
        print(f"  {a.id:20s} v{a.version} [{a.status}] for {', '.join(a.compatible_with) or 'any'}: {a.name} — files {a.files}")
