import re
import os
from mechanism import Mechanism, Reaction


# ---------------------------------------------------------------------------
# FUNCTION macro expander
# ---------------------------------------------------------------------------

def _split_top_level_args(s: str) -> list:
    """Split s by top-level commas (not inside parentheses)."""
    args, current, depth = [], [], 0
    for ch in s:
        if ch == ',' and depth == 0:
            args.append(''.join(current).strip())
            current = []
        else:
            if ch == '(': depth += 1
            elif ch == ')': depth -= 1
            current.append(ch)
    if current:
        args.append(''.join(current).strip())
    return [a for a in args if a]


_INCLUDE_NAME_RE = re.compile(r'[ \t]+["\']?([^"\';\r\n]+?)["\']?[ \t]*(?=;|\r?$)')   # the file name after '#include'


def expand_macros(expr: str, functions: dict, max_depth: int = 10) -> str:
    """Expand user-defined FUNCTION macros inline.

    For each NAME=FUNCTION[body] definition, replaces NAME(a,b,c) with the
    body where %1→a, %2→b, %3→c.  Applied recursively until stable.
    """
    if not functions:
        return expr
    for _ in range(max_depth):
        changed = False
        for fname, template in functions.items():
            pat = re.compile(rf'\b{re.escape(fname)}\s*\(', re.IGNORECASE)
            new_expr = ''
            pos = 0
            for m in pat.finditer(expr):
                new_expr += expr[pos:m.start()]
                arg_start = m.end()
                depth = 1
                p = arg_start
                while p < len(expr) and depth > 0:
                    if expr[p] == '(':   depth += 1
                    elif expr[p] == ')': depth -= 1
                    p += 1
                args = _split_top_level_args(expr[arg_start:p - 1])
                expanded = template
                for i, arg in enumerate(args, 1):
                    expanded = re.sub(rf'%{i}\b', f'({arg.strip()})', expanded)
                new_expr += f'({expanded})'
                pos = p
                changed = True
            new_expr += expr[pos:]
            expr = new_expr
        if not changed:
            break
    return expr

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir
        self.rates_map = {}
        self.functions = {}   # NAME → FUNCTION[body] template (from EASY files)
        # FUNCTION macros are expanded where they are USED (EASY semantics): a redefinition only
        # applies to the statements after it.  The first definition keeps its name; a different
        # later body becomes NAME__2, NAME__3, …, and raw definitions stored for the evaluator
        # are rewritten to call the version in effect at their line.
        self._fn_current = {}   # NAME → versioned name in effect (NAME, NAME__2, …)
        self._fn_bodies = {}    # versioned name → body
        # Canon keys of reactions explicitly disabled by -/-> across ALL parse() calls.
        # scenario_runner uses this to retroactively drop already-materialised reactions.
        self.disabled_rxn_keys = set()
        # Fallback rates: when -/-> disables A+B→C, a replacement bare reaction
        # A+B→D with no explicit k[] definition inherits the disabled rate.
        # Maps frozenset(reactant_names_upper) → rate_expression.
        self._reactant_fallback_rates = {}

    @staticmethod
    def _canon_rxn_key(lhs: str, rhs: str, aliases: dict = None) -> str:
        """Return a fully normalised reaction key (homogeniser).

        Applied in order:
          1. Alias substitution  — HO2NO2 → HNO4, HNO2 → HONO, etc.
          2. Coefficient merging — 0.5*NO + 0.5*NO → NO
          3. Alphabetical sort  — O3+NO → NO+O3
          4. Canonical notation — NA+NA → 2*NA
        All four steps must happen together so that every call site
        (pre-scans, Pass 1 rate defs, Pass 2 bare reactions, comparison
        script) produces the same canonical string regardless of order.
        """
        def _expand(side: str) -> list:
            coeff_map: dict = {}
            for term in side.upper().replace(' ', '').split('+'):
                if not term:
                    continue
                m = re.match(r'^([\d.]+)\*(.+)$', term)
                if m:
                    try:
                        coeff = float(m.group(1))
                    except ValueError:
                        coeff = 1.0
                    species = m.group(2)
                else:
                    coeff = 1.0
                    species = term
                # Step 1: alias substitution (HO2NO2 → HNO4 etc.)
                if aliases:
                    species = aliases.get(species, species)
                coeff_map[species] = coeff_map.get(species, 0.0) + coeff
            tokens = []
            for sp in sorted(coeff_map.keys()):
                c = coeff_map[sp]
                c_int = int(round(c))
                if abs(c - c_int) < 1e-9:
                    if c_int == 1:
                        tokens.append(sp)
                    else:
                        tokens.append(f'{c_int}*{sp}')
                else:
                    tokens.append(f'{c:.6g}*{sp}')
            return tokens

        lhs_tokens = _expand(lhs)
        rhs_tokens = _expand(rhs)
        return '+'.join(lhs_tokens) + '-->' + '+'.join(rhs_tokens)

    def _find_file(self, filename: str) -> str:
        if os.path.isabs(filename) and os.path.exists(filename): return filename
        local_path = os.path.join(self.base_path, filename)
        if os.path.exists(local_path): return local_path
        if self.include_dir and os.path.exists(self.include_dir):
            for root, _, files in os.walk(self.include_dir):
                if filename in files: return os.path.join(root, filename)
        raise FileNotFoundError(f"Include file '{filename}' not found.")

    def preprocess(self, filename: str) -> str:
        try:
            path = self._find_file(filename)
            with open(path, 'r', encoding='latin-1') as f:
                content = f.read()
        except FileNotFoundError:
            print(f"[WARNING] Include file '{filename}' not found. Skipping.")
            return ""
        # content = re.sub(r';.*', '', content) # DEACTIVATED: Keep comments for line-by-line parsing
        # Expand #include directives that are NOT commented out (no ';' before them on the line).  The name
        # runs to the end of the line or a ';' comment and may contain spaces.  Line-based instead of a
        # whole-file regex with a greedy ^[^;]* prefix: that prefix backtracked over every character of
        # every line, quadratic in the line length - the chunked COA/PERO sums of the GECKO-A includes are
        # 10-kB lines and the full toluene mechanism took minutes here (F. Rohrer, 85057db).  Semantics
        # unchanged: the text from the line start up to the include name is replaced by the included file,
        # a trailing comment is kept.
        if '#include' in content:
            out = []
            for line in content.split('\n'):
                idx = line.find('#include')
                if idx >= 0 and ';' not in line[:idx]:
                    m = _INCLUDE_NAME_RE.match(line, idx + len('#include'))
                    if m:
                        out.append(self.preprocess(m.group(1).strip()) + line[m.end():])
                        continue
                out.append(line)
            content = '\n'.join(out)
        content = re.sub(r'\$\s*\n', '', content) 
        return content

    def parse(self, filename: str):
        # Create a new mechanism object for each file, as the runner accumulates them.
        self.mech = Mechanism()
        # Use an ordered dict to deduplicate reactions by normalised key (last definition wins,
        # matching FACSIMILE override semantics — duplicate k[A-->B] lines in .mech files are
        # treated as overrides, not additions).
        self._rxn_map = {}       # normalised_key -> (reacs, prods, rate_expr)
        self._k_raw_defs = {}    # normalised_key -> original definition string (Phase 1)
        _species_seen = set(self.mech.species)   # membership set mirroring the ordered species list (process_side)
        full_content = self.preprocess(filename)
        lines = full_content.split('\n')

        # RESERVED: names excluded from parameter capture.
        # Always reserved: math functions and FACSIMILE keywords.
        # Physical env variables (T, P, M, H2O, HV, RH, DIL) are reserved when
        # parsing individual mechanism files — prevents include-file defaults from
        # overriding EASY-master or ENZ values.
        # In elog mode the flat file already has correct order (EASY master comes
        # last), so physical variables are NOT reserved — "last wins" handles it.
        # Example: HV=INPUT(hv) then HV=CONST(0) at the end → HV=0 wins.
        _is_elog = filename.lower().endswith('.elog')
        _PHYS = {'T', 'TEMP', 'P', 'PRESS', 'M', 'H2O', 'HV', 'RH', 'DIL'}
        RESERVED = {'EXP', 'LOG', 'LOG10',
                    'CONST', 'SUM', 'VARIABLE', 'FUNCTION', 'INPUT',
                    'PHTFKT', 'TROE', 'STOP', 'PHOTMMSC'}
        if not _is_elog:
            RESERVED = RESERVED | _PHYS   # non-elog: keep physical vars reserved

        if 'JSCALE' not in self.mech.parameters:
             self.mech.parameters['JSCALE'] = 'JNO2 / 1.165E-02'

        # ── TWO-PASS PARSE ───────────────────────────────────────────────────────
        # EASY files often place bare reaction lines (A --> B) *before* the
        # corresponding k[A-->B]=CONST(...) rate definition in the same file
        # (e.g. saphir.walldep.incl: reactions on lines 11-16, rates on 18-23).
        # A single-pass parser picks up the wrong rate (or "0.0") for those.
        #
        # Pass 1: collect parameters + ALL rate definitions into rates_map.
        # Pass 2: register bare reaction lines using the fully-populated rates_map.
        # ─────────────────────────────────────────────────────────────────────────

        # Pre-scan 0: collect ALL FUNCTION definitions before Pass 1 so that
        # parameters defined before their function macro (e.g. KFPAN=COMPRC(...)
        # appearing in an include before COMPRC=FUNCTION[...]) are still expanded.
        # Only the FIRST definition of a name is taken here (a use before any definition binds to
        # the first one); redefinitions take effect in Pass 1 at their own line.
        _func_re = re.compile(r'^(\w+)\s*=\s*FUNCTION\s*\[(.+)\]\s*$', re.I)
        for _l in lines:
            _s = _l.strip()
            if _s.startswith((';', '*', '#')):
                continue
            _fm = _func_re.match(_s.split(';')[0].strip())
            if _fm:
                _fname = _fm.group(1).strip().upper()
                _fbody = _fm.group(2).strip()
                if _fname not in self.functions:
                    self.functions[_fname] = _fbody
                    self.mech.functions[_fname] = _fbody

        def _versioned(text: str) -> str:
            """Rewrite calls NAME( → NAME__n( for macros redefined before this line."""
            for _n, _v in self._fn_current.items():
                if _v != _n:
                    text = re.sub(rf'\b{_n}\s*\(', f'{_v}(', text, flags=re.I)
            return text

        # Pre-scan 1b: collect ALIAS directives so species substitution is available
        # before Pass 1 processes reactions and parameters.
        _alias_re = re.compile(r'^ALIAS\s*\[\s*(\w+)\s*\]\s*=\s*(.+)$', re.I)
        for _l in lines:
            _s = _l.strip()
            if _s.startswith((';', '*', '#')):
                continue
            _am = _alias_re.match(_s.split(';')[0].strip())
            if _am:
                _canonical = _am.group(1).strip().upper()
                for _alt in _am.group(2).split():
                    _alt_up = _alt.strip().upper()
                    if _alt_up and _alt_up != _canonical:
                        self.mech.aliases[_alt_up] = _canonical

        # Pre-scan 1c: collect -/-> directives and build disabled_rxn_keys BEFORE
        # Pass 1/2 touch the mechanism.  A disabled reaction never enters _rxn_map
        # and therefore never reaches mechanism.reactions, kv, or pd output.
        # Must run after ALIAS pre-scan (canonical names needed for canon key).
        # Must also collect bare-reaction LHS tokens for replacement detection.
        _bare_lhs_prescan = []   # (line_idx, frozenset(lhs_tokens))
        for _idx2, _l2 in enumerate(lines):
            _s2 = _l2.strip()
            if not _s2 or _s2.startswith((';', '*', '#')):
                continue
            _c2 = _s2.split(';')[0].strip()
            if '-->' in _c2 and '=' not in _c2 and '-/->' not in _c2:
                _lhs2 = _c2.split('-->', 1)[0]
                _tok2 = frozenset(
                    self.mech.aliases.get(t.upper(), t.upper())
                    for t in _lhs2.replace(' ', '').split('+') if t)
                _bare_lhs_prescan.append((_idx2, _tok2))

        def _has_replacement(line_idx, lhs_tokens):
            for _bi, _blhs in _bare_lhs_prescan:
                if _bi > line_idx and _blhs == lhs_tokens:
                    return True
            return False

        for _idx3, _l3 in enumerate(lines):
            _s3 = _l3.strip()
            if not _s3 or _s3.startswith((';', '*', '#')):
                continue
            _c3 = _s3.split(';')[0].strip()
            if '-/->' not in _c3 or '=' in _c3:
                continue
            _ld, _rd = _c3.split('->', 1)
            _ld = _ld.rstrip('-/')
            _rd = _rd.lstrip('>')
            _lhs_tok = frozenset(
                self.mech.aliases.get(t.upper(), t.upper())
                for t in _ld.replace(' ', '').split('+') if t)
            _canon3 = self._canon_rxn_key(_ld, _rd, self.mech.aliases)
            if not _has_replacement(_idx3, _lhs_tok):
                # Lone -/->: disable without replacement
                self.disabled_rxn_keys.add(_canon3)
            else:
                # Has replacement with different products: disable original
                # (same-canon replacements overwrite _rxn_map naturally)
                _repl_same = False
                for _bi3, _blhs3 in _bare_lhs_prescan:
                    if _bi3 > _idx3 and _blhs3 == _lhs_tok:
                        _repl_l3 = lines[_bi3].strip().split(';')[0].strip()
                        if '-->' in _repl_l3:
                            _rl3, _rr3 = _repl_l3.split('-->', 1)
                            if self._canon_rxn_key(_rl3, _rr3, self.mech.aliases) == _canon3:
                                _repl_same = True
                        break
                if not _repl_same:
                    self.disabled_rxn_keys.add(_canon3)

        # Pre-scan: collect (line_index, lhs_tokens) for every BARE `-->` reaction
        # in this file. Used by the -/-> handler to honour FAC's rule: a `-/->`
        # line only disables the reaction if a bare `-->` replacement with
        # matching LHS appears LATER in the same file. A lone `-/->` with no
        # subsequent replacement is ignored (FAC keeps the original rate).
        # Example that required this (scene_03, Apr 2026):
        #   `HV + O3 -/-> O` with no follow-up → FAC kept k[O3+hv-->O] active.
        bare_rxn_lhs = []  # list of (line_idx, frozenset(lhs_tokens_upper))
        for _idx, _ln in enumerate(lines):
            _l = _ln.strip()
            if not _l or _l.startswith(('*', ';', '#')):
                continue
            _lr = _l.split(';')[0].strip()
            if '-/->' in _lr or '=' in _lr or _lr.startswith('%'):
                continue
            if '-->' in _lr:
                _left = _lr.split('-->', 1)[0]
                _tokens = frozenset(t.upper() for t in
                                    _left.replace(' ', '').split('+') if t)
                bare_rxn_lhs.append((_idx, _tokens))

        def _has_subsequent_bare_replacement(line_idx, lhs_tokens):
            for _bi, _blhs in bare_rxn_lhs:
                if _bi > line_idx and _blhs == lhs_tokens:
                    return True
            return False

        # ── Builtin function expansion helpers ──────────────────────────────────
        # These are defined at parse() scope so they are available throughout
        # _parse_line_pass1 regardless of execution order.
        # _BUILTIN_FNS, _resolve_builtin, _check_unexpanded_builtins are the
        # SINGLE control point for all FACSIMILE built-in function calls.

        def _phtfkt_builtin(args_str):
            al = [a.strip() for a in args_str.split(',')]
            if len(al) < 3:
                return None
            a, b, c = al[0], al[1], al[2]
            out = (f"({a}) * AMAX(COSX,0.01)**({b}) * EXP(({c})/AMAX(COSX,0.01))"
                   f" * AMIN(JNO2 / (LNO2 * AMAX(COSX,0.01)**MNO2"
                   f" * EXP(NNO2/AMAX(COSX,0.01))), 1.0)")
            return out.upper()

        def _photmmsc_builtin(args_str):
            al = [a.strip() for a in args_str.split(',')]
            if len(al) < 7:
                return None
            a, b, c, d, e, f, g = al[:7]
            out = (f"({a}) * AMAX(COSX,0.01)**({b}) * EXP(({c})/AMAX(COSX,0.01))"
                   f" * AMIN(({g}) / (({d}) * AMAX(COSX,0.01)**({e})"
                   f" * EXP(({f})/AMAX(COSX,0.01))), 1.0)")
            return out.upper()

        _BUILTIN_FNS = {'PHTFKT': _phtfkt_builtin, 'PHOTMMSC': _photmmsc_builtin}

        def _resolve_builtin(func_name, args_str, context_line=''):
            """Single control point for FACSIMILE built-in function calls.
            EASY script FUNCTION wins; ELOG mode + missing def → STOP."""
            fn_up = func_name.upper()
            if fn_up in self.functions:
                return expand_macros(f"{fn_up}({args_str})", self.functions).upper()
            if _is_elog:
                raise RuntimeError(
                    f"\n[PARSER ERROR] ELOG '{filename}':\n"
                    f"  '{fn_up}(...)' is used but '{fn_up} = FUNCTION[...]'\n"
                    f"  is not defined in the EASY script.\n"
                    f"  All functions must be defined by the script itself.\n"
                    + (f"  Line: {context_line}" if context_line else "")
                )
            if fn_up in _BUILTIN_FNS:
                return _BUILTIN_FNS[fn_up](args_str)
            raise RuntimeError(
                f"\n[PARSER ERROR] '{filename}':\n"
                f"  '{fn_up}(...)' is used but has no FUNCTION definition\n"
                f"  and no builtin expansion exists for it.\n"
                + (f"  Line: {context_line}" if context_line else "")
            )

        def _check_unexpanded_builtins(expr, context):
            """After expand_macros: verify no FACSIMILE built-in calls remain."""
            for fn in tuple(_BUILTIN_FNS) + ('TROE', 'COMPRC'):
                if re.search(rf'\b{fn}\s*\(', expr, re.I):
                    if _is_elog:
                        raise RuntimeError(
                            f"\n[PARSER ERROR] ELOG '{filename}':\n"
                            f"  {context}: '{fn}(...)' remains unexpanded after\n"
                            f"  macro expansion — '{fn} = FUNCTION[...]' is missing\n"
                            f"  or defined after first use."
                        )
                    print(f"  [WARNING] {context}: '{fn}(...)' unexpanded "
                          f"('{fn}' not defined in '{filename}').")

        def _parse_line_pass1(line_idx, line):
            """Pass 1: handle parameters, rate defs, photolysis defs, initials."""
            original_line = line.strip()
            if not original_line:
                return
            is_commented = original_line.startswith(('*', ';', '#'))
            clean_line = re.sub(r'^\s*[*;#]+\s*', '', original_line).strip()
            line_for_rxn = original_line.split(';')[0].strip()
            if not clean_line:
                return

            # FACSIMILE "-/->" notation: disable-plus-replace marker. Only applies
            # when a bare `-->` with matching LHS follows later in this file.
            if '-/->' in line_for_rxn and '=' not in line_for_rxn and not is_commented:
                left_d, right_d = line_for_rxn.split('->',  1)   # splits on first '->'
                left_d = left_d.rstrip('-/')                      # strip trailing -/
                right_d = right_d.lstrip('>')                     # strip leading >
                lhs_d_tokens = frozenset(
                    t.upper() for t in left_d.replace(' ', '').split('+') if t)
                if not _has_subsequent_bare_replacement(line_idx, lhs_d_tokens):
                    # Lone -/-> with no replacement — DISABLE the reaction.
                    # FACSIMILE removes the reaction pathway; no bare replacement
                    # means it is simply deactivated, not replaced.
                    canon_d = self._canon_rxn_key(left_d, right_d, self.mech.aliases)
                    original_rate = self.rates_map.get(canon_d)
                    if original_rate and original_rate != '0.0':
                        self._reactant_fallback_rates[lhs_d_tokens] = original_rate
                    self.rates_map[canon_d] = '0.0'
                    self.disabled_rxn_keys.add(canon_d)
                    return
                # Normalize to canonical form — single key throughout.
                canon_d = self._canon_rxn_key(left_d, right_d, self.mech.aliases)
                # Capture original rate BEFORE zeroing — replacement reactions with the
                # same reactant pair (different products) can inherit it (FACSIMILE
                # rate-inheritance convention).
                original_rate = self.rates_map.get(canon_d)
                if original_rate and original_rate != '0.0':
                    self._reactant_fallback_rates[lhs_d_tokens] = original_rate
                # Check whether the replacement has the same canon key as the original
                # (e.g. reactants in different order: NO+O3-/->NO2 / O3+NO-->NO2).
                # In that case the replacement naturally overwrites _rxn_map —
                # no disable needed (and adding the key would block the replacement).
                _replacement_same_canon = False
                for _bi, _blhs in bare_rxn_lhs:
                    if _bi > line_idx and _blhs == lhs_d_tokens:
                        _repl_raw = lines[_bi].strip().split(';')[0].strip()
                        if '-->' in _repl_raw and '-/->' not in _repl_raw:
                            _repl_l, _repl_r = _repl_raw.split('-->', 1)
                            if self._canon_rxn_key(_repl_l, _repl_r, self.mech.aliases) == canon_d:
                                _replacement_same_canon = True
                        break
                if not _replacement_same_canon:
                    # Zero the canonical rate so the replacement (different products)
                    # doesn't inherit the old rate via fallback lookup.
                    self.rates_map[canon_d] = "0.0"
                    # Track across all parse() calls so scenario_runner can
                    # retroactively drop reactions materialised in earlier files.
                    self.disabled_rxn_keys.add(canon_d)
                return

            # ALIAS directive: ALIAS[canonical] = alt1 [alt2 ...]
            # Maps alternative species names to a single canonical name so that
            # HO2NO2 and HNO4, HONO and HNO2, etc. are treated identically.
            alias_m = re.match(r'ALIAS\s*\[\s*(\w+)\s*\]\s*=\s*(.+)$', clean_line, re.I)
            if alias_m and not is_commented:
                canonical = alias_m.group(1).strip().upper()
                for alt in alias_m.group(2).split():
                    alt_up = alt.strip().upper()
                    if alt_up and alt_up != canonical:
                        self.mech.aliases[alt_up] = canonical
                return

            # FUNCTION macro definition: NAME = FUNCTION[body]   (an inline "; comment" after the
            # closing bracket is allowed, as in the pre-scan above — e.g. scene_01's
            # DENSITY=FUNCTION[(%1)/((%2)*1.379E-19)]  ; P/(T*1.379E-19))
            func_m = re.match(r'(\w+)\s*=\s*FUNCTION\s*\[(.+)\]\s*$', clean_line.split(';')[0].strip(), re.I)
            if func_m:
                if not is_commented:
                    fname = func_m.group(1).strip().upper()
                    body  = _versioned(func_m.group(2).strip())   # nested macros: the versions in effect now
                    cur = self._fn_current.get(fname)
                    if cur is None:
                        vname = fname                              # first definition keeps the name
                    elif self._fn_bodies.get(cur) == body:
                        vname = cur                                # identical redefinition: nothing changes
                    else:
                        n = 2
                        while f"{fname}__{n}" in self._fn_bodies:
                            n += 1
                        vname = f"{fname}__{n}"                    # a new body from this line on
                    self._fn_current[fname] = vname
                    self._fn_bodies[vname] = body
                    # inline expansion (mech.parameters, rates) uses the definition in effect now
                    self.functions[fname] = body
                    self.mech.functions[fname] = body
                    from mechanism import FunctionEntry
                    # the evaluator resolves calls by name: every version stays available
                    self.mech.function_entries = [
                        e for e in self.mech.function_entries if e.name != vname]
                    self.mech.function_entries.append(FunctionEntry(name=vname, definition=body))
                return

            # Parameter definitions
            if ('=' in clean_line
                    and not clean_line.lower().startswith(('k[', 'initial', 'facs', 'files', 'alias'))
                    and '-->' not in clean_line
                    and not clean_line.startswith('%')):
                name = clean_line.split('=')[0].strip().upper()
                if name and name not in RESERVED and re.match(r'^[A-Z][A-Z0-9_]*$', name):
                    # Use line_for_rxn (inline comment stripped) for the expression
                    expr = (line_for_rxn.split('=', 1)[1].strip()
                            if '=' in line_for_rxn else clean_line.split('=', 1)[1].strip())
                    expr = re.sub(r'^(CONST|INPUT)\(\s*(.*?)\s*\)\s*$', r'\2', expr, flags=re.I)
                    # Expand any user-defined FUNCTION macros in the expression
                    expr = expand_macros(expr, self.functions)
                    # Validate: no FACSIMILE built-in function calls may remain
                    _check_unexpanded_builtins(expr, f"parameter '{name}'")
                    target_dict = self.mech.commented_parameters if is_commented else self.mech.parameters
                    target_dict[name] = expr.upper()
                    if not is_commented:
                        from mechanism import SEntry
                        raw_def = (line_for_rxn.split('=', 1)[1].strip()
                                   if '=' in line_for_rxn else expr)
                        self.mech.s_entries = [
                            e for e in self.mech.s_entries if e.name != name]
                        self.mech.s_entries.append(SEntry(name=name, definition=_versioned(raw_def)))
                return

            if is_commented:
                return

            # Photolysis with J<N>
            # Full FACSIMILE PHTFKT(a,b,c) expansion:
            #   HV × (a · COSX^b · exp(c/COSX)) × min(JNO2 / (LNO2·COSX^MNO2·exp(NNO2/COSX)), 1)
            # The clear-sky spectral factor (COSX^b · exp(c/COSX)) is specific to each J_i and
            # must NOT be dropped — the old `JSCALE * a * HV` collapse (valid only if every J_i
            # shared JNO2's spectrum) caused systematic 1.5–2× bias in scenarios driven by
            # clear-sky photolysis. Full PHOTMMSC(a,b,c,d,e,f,g) takes its own reference tuple
            # (d,e,f,g) instead of the global (LNO2,MNO2,NNO2,JNO2).
            photo_match = re.search(r'k\[(.*?)\](?:\{J<(\d+)>)\}\s*=\s*(PHTFKT|PHOTMMSC)\((.*?)\)', line_for_rxn, re.I)
            if photo_match:
                rxn_str, j_num, func, args = photo_match.groups()
                rate_param_name = f"J{j_num}"
                new_rate = _resolve_builtin(func, args)
                if new_rate is None:
                    new_rate = f"JSCALE * {args.split(',')[0].strip()}"
                self.mech.parameters[rate_param_name] = new_rate
                raw_rxn_key = rxn_str.replace(' ', '').upper()
                if '-->' in raw_rxn_key:
                    lhs_p, rhs_p = raw_rxn_key.split('-->', 1)
                    norm_key = self._canon_rxn_key(lhs_p, rhs_p, self.mech.aliases)
                else:
                    norm_key = raw_rxn_key
                self.rates_map[norm_key] = rate_param_name
                self._k_raw_defs[norm_key] = _versioned(f"{func.upper()}({args})")
                # Phase 1: register in s_entries if not already there (from a separate param line).
                # Mark synthetic so it is excluded from cs export (not in EASY script explicitly).
                from mechanism import SEntry
                if not any(e.name == rate_param_name for e in self.mech.s_entries):
                    self.mech.s_entries.append(
                        SEntry(name=rate_param_name, definition=_versioned(f"{func.upper()}({args})"),
                               is_synthetic=True))
                return

            # PHOTMMSC without J<N>
            phot_no_j_match = re.search(r'k\[(.*?)\]\s*=\s*PHOTMMSC\((.*?)\)', line_for_rxn, re.I)
            if phot_no_j_match:
                rxn_str, args = phot_no_j_match.groups()
                # Build a unique parameter name from the full reaction key so that
                # different PHOTMMSC branches (e.g. NO3+hv→NO and NO3+hv→NO2+O3)
                # each get their own J parameter and don't overwrite each other.
                raw_rxn_key = rxn_str.replace(' ', '').upper()
                if '-->' in raw_rxn_key:
                    lhs_p, rhs_p = raw_rxn_key.split('-->', 1)
                    norm_key = self._canon_rxn_key(lhs_p, rhs_p, self.mech.aliases)
                else:
                    norm_key = raw_rxn_key
                param_suffix = norm_key.replace('-->', '_').replace('+', '_').replace('*', 'x')
                rate_param_name = f"J_{param_suffix}"
                new_rate = _resolve_builtin('PHOTMMSC', args)
                if new_rate is None:
                    new_rate = f"JSCALE * {args.split(',')[0].strip()}"
                self.mech.parameters[rate_param_name] = new_rate
                self.rates_map[norm_key] = rate_param_name
                self._k_raw_defs[norm_key] = _versioned(f"PHOTMMSC({args})")
                # Phase 1: register in s_entries so the evaluator can use it.
                # Mark synthetic: J_HV_* names never appear in the EASY script explicitly.
                from mechanism import SEntry
                self.mech.s_entries = [e for e in self.mech.s_entries
                                       if e.name != rate_param_name]
                self.mech.s_entries.append(
                    SEntry(name=rate_param_name, definition=_versioned(f"PHOTMMSC({args})"),
                           is_synthetic=True))
                return

            # INITIAL conditions
            init_match = re.search(r'INITIAL\s*\[(.*?)\]\s*=\s*CONST\((.*)\)', line_for_rxn, re.I)
            if init_match:
                spec, expr = init_match.groups()
                self.mech.initial_conditions[spec.strip().upper()] = expr.strip().upper()
                return

            # Rate definitions  k[...] = CONST/INPUT/FUNCTION(...)
            rate_def = re.match(r'k\[(.*?)\](?:\{[^}]*\})?\s*=\s*(?:CONST|INPUT|FUNCTION)?\s*\((.*)\)', line_for_rxn, re.I)
            if rate_def:
                raw_key = rate_def.group(1).replace(' ', '').upper()
                rate_val = expand_macros(rate_def.group(2).strip(), self.functions).upper()
                # HV must NOT be injected into the rate constant — it is already a
                # reactant in the reaction equation (reactants_env) and gates the
                # reaction there at runtime.  Injecting HV here would cause the kv
                # export to show time-averaged J-values instead of instantaneous ones,
                # making them incomparable with FACSIMILE kv output.
                # Normalize to canonical form BEFORE storing — single key, no dual-key.
                if '-->' in raw_key:
                    lhs_r, rhs_r = raw_key.split('-->', 1)
                    norm_key = self._canon_rxn_key(lhs_r, rhs_r, self.mech.aliases)
                else:
                    norm_key = raw_key
                self.rates_map[norm_key] = rate_val
                raw_def = line_for_rxn.split('=', 1)[1].strip() if '=' in line_for_rxn else rate_val
                self._k_raw_defs[norm_key] = _versioned(raw_def)

        def _parse_line_pass2(line):
            """Pass 2: register bare reaction lines using the complete rates_map."""
            original_line = line.strip()
            if not original_line or original_line.startswith(('*', ';', '#')):
                return
            line_for_rxn = original_line.split(';')[0].strip()

            def process_side(side_str):
                parts = [s.strip() for s in side_str.split('+') if s.strip()]
                normalised = []
                for p in parts:
                    coeff_part, _, sp_part = p.rpartition('*')
                    name = sp_part.strip().upper()
                    # Apply alias substitution: HO2NO2 → HNO4, HNO2 → HONO, etc.
                    name = self.mech.aliases.get(name, name)
                    norm_p = (coeff_part + '*' + name) if coeff_part else name
                    normalised.append(norm_p)
                    if name not in RESERVED and name not in self.mech.parameters and len(name) > 0:
                        if name not in _species_seen:          # set lookup; `name in self.mech.species` was a
                            _species_seen.add(name)            # linear list scan -> O(species x terms), the
                            self.mech.species.append(name)     # dominant setup cost for 1e5-species mechanisms
                return normalised

            if '-->' in line_for_rxn and '=' not in line_for_rxn:
                left, right = line_for_rxn.split('-->', 1)
                reacs = process_side(left)
                prods = process_side(right)
                # Normalize to canonical form first — single key lookup, no fallback.
                norm_key = self._canon_rxn_key(left, right, self.mech.aliases)
                # Skip reactions disabled by pre-scan — they must not enter _rxn_map.
                if norm_key in self.disabled_rxn_keys:
                    return
                rate = self.rates_map.get(norm_key)
                if rate is None:
                    # FACSIMILE rate-inheritance: bare reaction with no explicit k[]
                    # inherits the rate of a -/->-disabled reaction with same reactants.
                    reac_key = frozenset(
                        t.upper() for t in left.replace(' ', '').split('+') if t)
                    rate = self._reactant_fallback_rates.get(reac_key, "0.0")
                self._rxn_map[norm_key] = (reacs, prods, rate)

            elif line_for_rxn.startswith('%'):
                parts = line_for_rxn.split(':', 1)
                if len(parts) == 2:
                    rate_expr = parts[0].replace('%', '').strip().upper()
                    rxn_part = parts[1]
                    if '=' in rxn_part:
                        left, right = rxn_part.split('=', 1)
                        reacs = process_side(left)
                        prods = process_side(right)
                        key_pct = rxn_part.replace(' ', '').upper()
                        self._rxn_map[key_pct] = (reacs, prods, rate_expr)

        for _i, line in enumerate(lines):
            _parse_line_pass1(_i, line)
        for line in lines:
            _parse_line_pass2(line)

        # Materialise reactions from dedup map (preserves insertion order of last occurrence).
        # Skip reactions disabled by -/-> within this same file (or any previously
        # parsed file — disabled_rxn_keys persists across calls).
        for canon_key, (reacs, prods, rate_expr) in self._rxn_map.items():
            if canon_key not in self.disabled_rxn_keys:
                self.mech.reactions.append(Reaction(reacs, prods, rate_expr))

        # ── Phase 1: populate k_entries from _rxn_map + _k_raw_defs ─────────
        from mechanism import KEntry
        for canon_key, (reacs, prods, rate_expr) in self._rxn_map.items():
            if canon_key in self.disabled_rxn_keys:
                continue
            raw_def = self._k_raw_defs.get(canon_key)
            if raw_def is None:
                # Bare reaction with no explicit k[]: definition from rates_map value
                raw_def = f"CONST({rate_expr})" if rate_expr != '0.0' else 'CONST(0)'
            self.mech.k_entries.append(KEntry(
                name=canon_key,
                definition=raw_def,
                reactants=reacs,
                products=prods,
            ))

        return self.mech