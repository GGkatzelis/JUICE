import re
import os
from mechanism import Mechanism, Reaction

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir
        self.rates_map = {}

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

        def _parse_line_pass1(line):
            """Pass 1: handle parameters, rate defs, photolysis defs, initials."""
            original_line = line.strip()
            if not original_line:
                return
            is_commented = original_line.startswith(('*', ';', '#'))
            clean_line = re.sub(r'^\s*[*;#]+\s*', '', original_line).strip()
            line_for_rxn = original_line.split(';')[0].strip()
            if not clean_line:
                return

            # Parameter definitions
            if ('=' in clean_line
                    and not clean_line.lower().startswith(('k[', 'initial', 'facs', 'files'))
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
            photo_match = re.search(r'k\[(.*?)\](?:\{J<(\d+)>)\}\s*=\s*(PHTFKT|PHOTMMSC)\((.*?)\)', line_for_rxn, re.I)
            if photo_match:
                rxn_str, j_num, func, args = photo_match.groups()
                rate_param_name = f"J{j_num}"
                j_arg = args.split(',')[0].strip()
                self.mech.parameters[rate_param_name] = f"JSCALE * {j_arg} * HV"
                self.rates_map[rxn_str.replace(' ', '').upper()] = rate_param_name
                return

            # PHOTMMSC without J<N>
            phot_no_j_match = re.search(r'k\[(.*?)\]\s*=\s*PHOTMMSC\((.*?)\)', line_for_rxn, re.I)
            if phot_no_j_match:
                rxn_str, args = phot_no_j_match.groups()
                rate_param_name = f"J_{rxn_str.split('+')[0].strip()}"
                j_arg = args.split(',')[0].strip()
                self.mech.parameters[rate_param_name] = f"JSCALE * {j_arg} * HV"
                self.rates_map[rxn_str.replace(' ', '').upper()] = rate_param_name
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
                rate = self.rates_map.get(key, self.rates_map.get(canon_key, "0.0"))
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

        for line in lines:
            _parse_line_pass1(line)
        for line in lines:
            _parse_line_pass2(line)

        # Materialise reactions from dedup map (preserves insertion order of last occurrence)
        for reacs, prods, rate_expr in self._rxn_map.values():
            self.mech.reactions.append(Reaction(reacs, prods, rate_expr))

        return self.mech