import re
import os
from mechanism import Mechanism, Reaction

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir

    def _find_file(self, filename: str) -> str:
        if os.path.isabs(filename) and os.path.exists(filename): return filename
        local_path = os.path.join(self.base_path, filename)
        if os.path.exists(local_path): return local_path
        if self.include_dir and os.path.exists(self.include_dir):
            for root, _, files in os.walk(self.include_dir):
                if filename in files: return os.path.join(root, filename)
        raise FileNotFoundError(f"Include file '{filename}' not found.")

    def preprocess(self, filename: str) -> str:
        path = self._find_file(filename)
        with open(path, 'r', encoding='latin-1') as f:
            content = f.read()
        content = re.sub(r';.*', '', content) 
        def include_replacer(match):
            return self.preprocess(match.group(1).strip())
        content = re.sub(r'#include\s+["\']?([^"\'\s]+)["\']?', include_replacer, content)
        content = re.sub(r'\$\s*\n', '', content) 
        return content

    def parse(self, filename: str):
        full_content = self.preprocess(filename)
        lines = full_content.split('\n')
        rates_map = {}
        
        # Reserved math/env keywords (Must be Upper)
        RESERVED = {'EXP', 'LOG', 'LOG10', 'TEMP', 'M', 'P', 'DIL', 'HV', 'RH', 'T', 'H2O', 'O2', 'N2',
                    'CONST', 'SUM', 'VARIABLE', 'FUNCTION', 'INPUT', 'PHTFKT', 'TROE', 'STOP'}

        # First pass: collect all parameters
        for line in lines:
            line = line.strip()
            if '=' in line and not line.lower().startswith(('k[', 'initial', 'facs', 'files')):
                name = line.split('=')[0].strip().upper()
                if name not in RESERVED:
                    expr = line.split('=')[1].strip()
                    expr = re.sub(r'^(CONST|INPUT|FUNCTION)\((.*)\)$', r'\2', expr, flags=re.I)
                    self.mech.parameters[name] = expr.upper()

        # Second pass: collect reactions and initial conditions
        for line in lines:
            line = line.strip()
            if not line or line.startswith(('*', ';')): continue

            # 1. Capture INITIAL
            init_match = re.search(r'INITIAL\s*\[(.*?)\]\s*=\s*CONST\((.*)\)', line, re.I)
            if init_match:
                spec, expr = init_match.groups()
                self.mech.initial_conditions[spec.strip().upper()] = expr.strip().upper()
                continue

            # 2. Capture Rate Definitions
            rate_def = re.match(r'k\[(.*?)\]\s*=\s*(?:CONST|INPUT|FUNCTION)?\s*\((.*)\)', line, re.I)
            if rate_def:
                rates_map[rate_def.group(1).replace(' ', '').upper()] = rate_def.group(2).strip().upper()
                continue

            # 3. Capture Reactions
            if '-->' in line and '=' not in line:
                clean_line = re.split(r'[;{]', line)[0].strip()
                left, right = clean_line.split('-->', 1)
                
                def process_side(side_str):
                    parts = [s.strip() for s in side_str.split('+') if s.strip()]
                    for p in parts:
                        name = p.split('*')[-1].strip().upper()
                        # Only add to species if it's not a parameter or reserved word
                        if name not in RESERVED and name not in self.mech.parameters and len(name) > 0:
                            if name not in self.mech.species: self.mech.species.append(name)
                    return parts

                reacs = process_side(left)
                prods = process_side(right)
                key = clean_line.replace(' ', '').upper()
                self.mech.reactions.append(Reaction(reacs, prods, rates_map.get(key, "0.0")))
        
        return self.mech