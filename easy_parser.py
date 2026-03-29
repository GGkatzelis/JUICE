import re
import os
from mechanism import Mechanism, Reaction

class EasyParser:
    def __init__(self, base_path=".", include_dir=None):
        self.mech = Mechanism()
        self.base_path = base_path
        self.include_dir = include_dir
        self.rates_map = {}

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
        content = re.sub(r'#include\s+["\']?([^"\'\s]+)["\']?', include_replacer, content)
        content = re.sub(r'\$\s*\n', '', content) 
        return content

    def parse(self, filename: str):
        # Create a new mechanism object for each file, as the runner accumulates them.
        self.mech = Mechanism()
        full_content = self.preprocess(filename)
        lines = full_content.split('\n')
        
        RESERVED = {'EXP', 'LOG', 'LOG10', 'TEMP', 'M', 'P', 'DIL', 'HV', 'RH', 'T', 'H2O', 'O2', 'N2',
                    'CONST', 'SUM', 'VARIABLE', 'FUNCTION', 'INPUT', 'PHTFKT', 'TROE', 'STOP', 'PHOTMMSC'}

        if 'JSCALE' not in self.mech.parameters:
             self.mech.parameters['JSCALE'] = 'JNO2 / 1.165E-02'
        
        for line in lines:
            original_line = line.strip()
            if not original_line:
                continue

            is_commented = original_line.startswith(('*', ';', '#'))
            
            # Clean the line for parameter detection
            clean_line = re.sub(r'^\s*[*;#]+\s*', '', original_line).strip()
            
            # Remove trailing comments for reaction parsing
            line_for_rxn = original_line.split(';')[0].strip()

            if not clean_line:
                continue

            # Handle parameter definitions (check in the cleaned line)
            if '=' in clean_line and not clean_line.lower().startswith(('k[', 'initial', 'facs', 'files')) and '-->' not in clean_line and not clean_line.startswith('%'):
                name = clean_line.split('=')[0].strip().upper()
                if name and name not in RESERVED and re.match(r'^[A-Z][A-Z0-9_]*$', name):
                    expr = clean_line.split('=')[1].strip()
                    expr = re.sub(r'^(CONST|INPUT|FUNCTION)\(\s*(.*?)\s*\)\s*$', r'\2', expr, flags=re.I)
                    
                    target_dict = self.mech.commented_parameters if is_commented else self.mech.parameters
                    target_dict[name] = expr.upper()
                continue
            
            # From here, we only care about non-commented lines for reactions
            if is_commented:
                continue

            # New photolysis handler
            photo_match = re.search(r'k\[(.*?)\](?:\{J<(\d+)>)\}\s*=\s*(PHTFKT|PHOTMMSC)\((.*?)\)', line_for_rxn, re.I)
            if photo_match:
                rxn_str, j_num, func, args = photo_match.groups()
                rate_param_name = f"J{j_num}"
                j_arg = args.split(',')[0].strip()
                self.mech.parameters[rate_param_name] = f"JSCALE * {j_arg} * HV"
                self.rates_map[rxn_str.replace(' ', '').upper()] = rate_param_name
                continue

            # Handle PHOTMMSC calls that don't have a J<> number
            phot_no_j_match = re.search(r'k\[(.*?)\]\s*=\s*PHOTMMSC\((.*?)\)', line_for_rxn, re.I)
            if phot_no_j_match:
                rxn_str, args = phot_no_j_match.groups()
                rate_param_name = f"J_{rxn_str.split('+')[0].strip()}"
                j_arg = args.split(',')[0].strip()
                self.mech.parameters[rate_param_name] = f"JSCALE * {j_arg} * HV"
                self.rates_map[rxn_str.replace(' ', '').upper()] = rate_param_name
                continue

            init_match = re.search(r'INITIAL\s*\[(.*?)\]\s*=\s*CONST\((.*)\)', line_for_rxn, re.I)
            if init_match:
                spec, expr = init_match.groups()
                self.mech.initial_conditions[spec.strip().upper()] = expr.strip().upper()
                continue

            rate_def = re.match(r'k\[(.*?)\]\s*=\s*(?:CONST|INPUT|FUNCTION)?\s*\((.*)\)', line_for_rxn, re.I)
            if rate_def:
                self.rates_map[rate_def.group(1).replace(' ', '').upper()] = rate_def.group(2).strip().upper()
                continue

            def process_side(side_str):
                parts = [s.strip() for s in side_str.split('+') if s.strip()]
                for p in parts:
                    name = p.split('*')[-1].strip().upper()
                    if name not in RESERVED and name not in self.mech.parameters and len(name) > 0:
                        if name not in self.mech.species: self.mech.species.append(name)
                return parts

            if '-->' in line_for_rxn and '=' not in line_for_rxn:
                left, right = line_for_rxn.split('-->', 1)
                reacs = process_side(left)
                prods = process_side(right)
                key = line_for_rxn.replace(' ', '').upper()
                self.mech.reactions.append(Reaction(reacs, prods, self.rates_map.get(key, "0.0")))

            elif line_for_rxn.startswith('%'):
                parts = line_for_rxn.split(':', 1)
                if len(parts) == 2:
                    rate_expr = parts[0].replace('%', '').strip().upper()
                    rxn_part = parts[1]
                    if '=' in rxn_part:
                        left, right = rxn_part.split('=', 1)
                        reacs = process_side(left)
                        prods = process_side(right)
                        self.mech.reactions.append(Reaction(reacs, prods, rate_expr))
        
        return self.mech