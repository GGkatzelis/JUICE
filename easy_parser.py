import re
import os
from mechanism import Mechanism, Reaction

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir
        self.rates_map = {}
        # Canon keys of reactions explicitly disabled by -/-> across ALL parse() calls.
        # scenario_runner uses this to retroactively drop already-materialised reactions.
        self.disabled_rxn_keys = set()
        # Fallback rates: when -/-> disables A+B→C, a replacement bare reaction
        # A+B→D with no explicit k[] definition inherits the disabled rate.
        # Maps frozenset(reactant_names_upper) → rate_expression.
        self._reactant_fallback_rates = {}

    @staticmethod
    def _canon_rxn_key(lhs: str, rhs: str) -> str:
        """Return a stoichiometry-normalised reaction key.

        Expands integer stoichiometric coefficients so that e.g.
        '2*NA' and 'NA+NA' both produce the same sorted token list.
        This lets k[N2O5-->2*NA] override the bare 'N2O5 --> NA + NA' line
        even though the literal strings differ.
        """
        def _expand(side: str) -> list:
            tokens = []
            for term in side.upper().replace(' ', '').split('+'):
                m = re.match(r'^(\d+)\*(.+)$', term)
                if m:
                    tokens.extend([m.group(2)] * int(m.group(1)))
                elif term:
                    tokens.append(term)
            return sorted(tokens)

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
        def include_replacer(match):
            return self.preprocess(match.group(1).strip())
        # Only match #include that is NOT preceded by a semicolon (i.e. not commented out).
        content = re.sub(r'(?m)^[^;]*#include\s+["\']?([^"\'\s]+)["\']?', include_replacer, content)
        content = re.sub(r'\$\s*\n', '', content) 
        return content

    def parse(self, filename: str):
        # Create a new mechanism object for each file, as the runner accumulates them.
        self.mech = Mechanism()
        # Use an ordered dict to deduplicate reactions by normalised key (last definition wins,
        # matching FACSIMILE override semantics — duplicate k[A-->B] lines in .mech files are
        # treated as overrides, not additions).
        self._rxn_map = {}   # normalised_key -> (reacs, prods, rate_expr)
        full_content = self.preprocess(filename)
        lines = full_content.split('\n')

        RESERVED = {'EXP', 'LOG', 'LOG10', 'TEMP', 'M', 'P', 'DIL', 'HV', 'RH', 'T', 'H2O', 'O2', 'N2',
                    'CONST', 'SUM', 'VARIABLE', 'FUNCTION', 'INPUT', 'PHTFKT', 'TROE', 'STOP', 'PHOTMMSC'}

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
                    # Lone -/-> with no replacement — FAC keeps original rate.
                    return
                canon_d = self._canon_rxn_key(left_d, right_d)
                raw_d = (left_d.replace(' ', '').upper()
                         + '-->' + right_d.replace(' ', '').upper())
                # Capture original rate BEFORE zeroing — replacement reactions with the
                # same reactant pair (different products) can inherit it (FACSIMILE
                # rate-inheritance convention).
                original_rate = (self.rates_map.get(raw_d)
                                 or self.rates_map.get(canon_d))
                if original_rate and original_rate != '0.0':
                    self._reactant_fallback_rates[lhs_d_tokens] = original_rate
                # Zero both the canon key (sorted) and the natural-order raw key so
                # that the pass-2 lookup (which checks raw key first) also returns 0.0.
                self.rates_map[canon_d] = "0.0"
                self.rates_map[raw_d]   = "0.0"
                # Track across all parse() calls so scenario_runner can retroactively
                # drop reactions that were materialised in an earlier file parse.
                self.disabled_rxn_keys.add(canon_d)
                return

            # ALIAS directive: ALIAS new_name = 'canonical_name'  or  ALIAS new = canonical
            # Maps one species name to another for display/output purposes.
            alias_m = re.match(r'ALIAS\s+(\w+)\s*=\s*[\'"]?(\w+)[\'"]?', clean_line, re.I)
            if alias_m and not is_commented:
                self.mech.aliases[alias_m.group(1).upper()] = alias_m.group(2).upper()
                return

            # Parameter definitions
            if ('=' in clean_line
                    and not clean_line.lower().startswith(('k[', 'initial', 'facs', 'files', 'alias'))
                    and '-->' not in clean_line
                    and not clean_line.startswith('%')):
                name = clean_line.split('=')[0].strip().upper()
                if name and name not in RESERVED and re.match(r'^[A-Z][A-Z0-9_]*$', name):
                    expr = clean_line.split('=')[1].strip()
                    expr = re.sub(r'^(CONST|INPUT|FUNCTION)\(\s*(.*?)\s*\)\s*$', r'\2', expr, flags=re.I)
                    target_dict = self.mech.commented_parameters if is_commented else self.mech.parameters
                    target_dict[name] = expr.upper()
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
            def _phtfkt_expansion(args_str):
                al = [a.strip() for a in args_str.split(',')]
                if len(al) < 3:
                    return None
                a, b, c = al[0], al[1], al[2]
                out = (f"HV * (({a}) * AMAX(COSX,0.01)**({b}) * EXP(({c})/AMAX(COSX,0.01)))"
                       f" * AMIN(JNO2 / (LNO2 * AMAX(COSX,0.01)**MNO2 * EXP(NNO2/AMAX(COSX,0.01))), 1.0)")
                # ctx keys are uppercase (JNO2, LNO2, COSX, HV, ...); args may reference
                # lowercase names from the easy file (e.g. `jno2`). Match the convention
                # used elsewhere in this parser where stored expressions are uppercased.
                return out.upper()
            def _photmmsc_expansion(args_str):
                al = [a.strip() for a in args_str.split(',')]
                if len(al) < 7:
                    return None
                a, b, c, d, e, f, g = al[:7]
                out = (f"HV * (({a}) * AMAX(COSX,0.01)**({b}) * EXP(({c})/AMAX(COSX,0.01)))"
                       f" * AMIN(({g}) / (({d}) * AMAX(COSX,0.01)**({e}) * EXP(({f})/AMAX(COSX,0.01))), 1.0)")
                return out.upper()

            photo_match = re.search(r'k\[(.*?)\](?:\{J<(\d+)>)\}\s*=\s*(PHTFKT|PHOTMMSC)\((.*?)\)', line_for_rxn, re.I)
            if photo_match:
                rxn_str, j_num, func, args = photo_match.groups()
                rate_param_name = f"J{j_num}"
                new_rate = (_photmmsc_expansion(args) if func.upper() == 'PHOTMMSC'
                            else _phtfkt_expansion(args))
                if new_rate is None:
                    # Fallback to legacy collapse if arg count is unexpected
                    new_rate = f"JSCALE * {args.split(',')[0].strip()} * HV"
                self.mech.parameters[rate_param_name] = new_rate
                raw_rxn_key = rxn_str.replace(' ', '').upper()
                self.rates_map[raw_rxn_key] = rate_param_name
                # Also store canon key so -/-> inheritance lookup can find it.
                if '-->' in raw_rxn_key:
                    lhs_p, rhs_p = raw_rxn_key.split('-->', 1)
                    canon_p = self._canon_rxn_key(lhs_p, rhs_p)
                    if canon_p != raw_rxn_key:
                        self.rates_map[canon_p] = rate_param_name
                return

            # PHOTMMSC without J<N>
            phot_no_j_match = re.search(r'k\[(.*?)\]\s*=\s*PHOTMMSC\((.*?)\)', line_for_rxn, re.I)
            if phot_no_j_match:
                rxn_str, args = phot_no_j_match.groups()
                # Build a unique parameter name from the full reaction key so that
                # different PHOTMMSC branches (e.g. NO3+hv→NO and NO3+hv→NO2+O3)
                # each get their own J parameter and don't overwrite each other.
                raw_rxn_key = rxn_str.replace(' ', '').upper()
                # Replace arrow and '+' to get a valid Python identifier
                param_suffix = raw_rxn_key.replace('-->', '_').replace('+', '_')
                rate_param_name = f"J_{param_suffix}"
                new_rate = _photmmsc_expansion(args)
                if new_rate is None:
                    new_rate = f"JSCALE * {args.split(',')[0].strip()} * HV"
                self.mech.parameters[rate_param_name] = new_rate
                self.rates_map[raw_rxn_key] = rate_param_name
                # Also store canon key.
                if '-->' in raw_rxn_key:
                    lhs_p, rhs_p = raw_rxn_key.split('-->', 1)
                    canon_p = self._canon_rxn_key(lhs_p, rhs_p)
                    if canon_p != raw_rxn_key:
                        self.rates_map[canon_p] = rate_param_name
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
                rate_val = rate_def.group(2).strip().upper()
                # Auto HV-gate: FACSIMILE treats `hv` on the reactant side as a
                # pseudo-species (0 in dark, 1 in light), so `k[X+hv-->Y]=<rate>`
                # is implicitly gated. Python strips HV from reactants (reserved),
                # so an explicit `* HV` must be injected when absent.
                if '-->' in raw_key:
                    _lhs_r = raw_key.split('-->', 1)[0]
                    _lhs_parts = {p.strip() for p in _lhs_r.split('+') if p.strip()}
                    if 'HV' in _lhs_parts and not re.search(r'\bHV\b', rate_val):
                        rate_val = f"HV*({rate_val})"
                self.rates_map[raw_key] = rate_val
                # Also store under stoichiometry-normalised key so that
                # k[N2O5-->2*NA] matches "N2O5 --> NA + NA" bare reaction.
                if '-->' in raw_key:
                    lhs_r, rhs_r = raw_key.split('-->', 1)
                    canon = self._canon_rxn_key(lhs_r, rhs_r)
                    if canon != raw_key:
                        self.rates_map[canon] = rate_val

        def _parse_line_pass2(line):
            """Pass 2: register bare reaction lines using the complete rates_map."""
            original_line = line.strip()
            if not original_line or original_line.startswith(('*', ';', '#')):
                return
            line_for_rxn = original_line.split(';')[0].strip()

            def process_side(side_str):
                parts = [s.strip() for s in side_str.split('+') if s.strip()]
                for p in parts:
                    name = p.split('*')[-1].strip().upper()
                    if name not in RESERVED and name not in self.mech.parameters and len(name) > 0:
                        if name not in self.mech.species:
                            self.mech.species.append(name)
                return parts

            if '-->' in line_for_rxn and '=' not in line_for_rxn:
                left, right = line_for_rxn.split('-->', 1)
                reacs = process_side(left)
                prods = process_side(right)
                key = line_for_rxn.replace(' ', '').upper()
                canon_key = self._canon_rxn_key(left, right)
                rate = self.rates_map.get(key) or self.rates_map.get(canon_key)
                if rate is None:
                    # FACSIMILE rate-inheritance: bare reaction with no explicit k[]
                    # inherits the rate of a -/->-disabled reaction with same reactants.
                    reac_key = frozenset(
                        t.upper() for t in left.replace(' ', '').split('+') if t)
                    rate = self._reactant_fallback_rates.get(reac_key, "0.0")
                self._rxn_map[canon_key] = (reacs, prods, rate)

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

        return self.mech