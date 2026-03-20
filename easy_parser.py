import re
import os
from mechanism import Mechanism, Reaction

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir

    def _find_file(self, filename: str) -> str:
        local_path = os.path.join(self.base_path, filename)
        if os.path.exists(local_path): return local_path
        if self.include_dir and os.path.exists(self.include_dir):
            for root, _, files in os.walk(self.include_dir):
                if filename in files: return os.path.join(root, filename)
        raise FileNotFoundError(f"Include file '{filename}' not found.")

    def preprocess(self, filename: str) -> str:
        if os.path.isabs(filename) and os.path.exists(filename): path = filename
        else: path = self._find_file(filename)

        with open(path, 'r', encoding='latin-1') as f:
            content = f.read()

        content = re.sub(r';.*', '', content) # Strip comments first
        def include_replacer(match):
            return self.preprocess(match.group(1).strip())
        content = re.sub(r'#include\s+["\']?([^"\'\s]+)["\']?', include_replacer, content)
        content = re.sub(r'\$\s*\n', '', content) # Handle continuations
        return content

    def parse(self, filename: str):
        full_content = self.preprocess(filename)

        # 1. Extract Aliases
        alias_matches = re.findall(r'ALIAS\[(.*?)\]\s*=\s*(.*)', full_content)
        for original, alias in alias_matches:
            for a in alias.split(): self.mech.aliases[a.strip()] = original.strip()

        # 2. Extract Functions (Macros)
        func_matches = re.findall(r'(\w+)\s*=\s*FUNCTION\[(.*?)\]', full_content, re.DOTALL)
        for name, body in func_matches: self.mech.functions[name] = body.strip()

        # 3. Extract Initial Conditions & Parameters
        lines = full_content.split('\n')
        for line in lines:
            line = line.strip()
            if not line: continue
            
            # Capture INITIAL[Species] = Value
            init_match = re.match(r'INITIAL\[(.*?)\]\s*=\s*(.*)', line)
            if init_match:
                spec, val = init_match.groups()
                self.mech.initial_conditions[spec.strip()] = val.strip()
                continue

            # Capture Parameters (T, P, etc.)
            param_match = re.match(r'^([a-zA-Z0-9_]+)\s*=\s*(.+)', line)
            if param_match and "FUNCTION" not in line and "-->" not in line and not line.startswith('k['):
                self.mech.parameters[param_match.group(1)] = param_match.group(2)
                continue

            # 4. Extract Reactions
            if '-->' in line and not line.startswith('k['):
                left, right = line.split('-->', 1)
                reactants = [r.strip() for r in left.split('+') if r.strip()]
                products = [p.strip() for p in right.split('+') if p.strip()]
                
                for s in reactants + products:
                    clean_s = s.split('*', 1)[1].strip() if '*' in s else s.strip()
                    resolved_s = self.mech.aliases.get(clean_s, clean_s)
                    if resolved_s not in self.mech.species and resolved_s.upper() not in ['HV', 'M', 'O2', 'N2', 'H2O']:
                        self.mech.species.append(resolved_s)
                self.mech.reactions.append(Reaction(reactants, products, ""))
            
        return self.mech