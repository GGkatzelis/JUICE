import re
from mechanism import Mechanism, Reaction

class EasyParser:
    """Parses an .easy file and populates a Mechanism object."""

    def __init__(self, file_path: str, include_dirs: list = None):
        self.file_path = file_path
        import os
        # Add the directory of the root file to the include paths
        self.base_dir = os.path.dirname(os.path.abspath(file_path))
        self.include_dirs = include_dirs or []
        if self.base_dir not in self.include_dirs:
            self.include_dirs.insert(0, self.base_dir)

        self.mechanism = Mechanism()
        self.lines = self._read_lines_with_includes(self.file_path)

    def _read_lines_with_includes(self, current_file: str, processed_files: set = None) -> list[str]:
        """Reads lines from the file and recursively processes #include directives."""
        import os
        
        if processed_files is None:
            processed_files = set()

        abs_current_file = os.path.abspath(current_file)
        if abs_current_file in processed_files:
            return []
        
        processed_files.add(abs_current_file)

        lines = []
        
        try:
            with open(abs_current_file, 'r') as f:
                raw_lines = f.readlines()
        except FileNotFoundError:
            found = False
            for inc_dir in self.include_dirs:
                path_to_try = os.path.join(inc_dir, os.path.basename(current_file))
                if os.path.exists(path_to_try):
                    return self._read_lines_with_includes(path_to_try, processed_files)
            
            # If not found directly in include dirs, walk them
            for inc_dir in self.include_dirs:
                for root, _, files in os.walk(inc_dir):
                    if os.path.basename(current_file) in files:
                        found_path = os.path.join(root, os.path.basename(current_file))
                        return self._read_lines_with_includes(found_path, processed_files)

            print(f"Warning: Could not find include file '{current_file}'.")
            return lines

        current_dir = os.path.dirname(abs_current_file)

        for line in raw_lines:
            stripped = line.strip()
            if stripped.startswith('#include'):
                parts = stripped.split()
                if len(parts) > 1:
                    inc_filename = parts[1]
                    
                    # 1. Try relative path first
                    path_to_try = os.path.join(current_dir, inc_filename)
                    if os.path.exists(path_to_try):
                        lines.extend(self._read_lines_with_includes(path_to_try, processed_files))
                        continue

                    # 2. If not, search in include directories (including subdirectories)
                    found_in_includes = False
                    for d in self.include_dirs:
                        for root, _, files in os.walk(d):
                            if inc_filename in files:
                                found_path = os.path.join(root, inc_filename)
                                lines.extend(self._read_lines_with_includes(found_path, processed_files))
                                found_in_includes = True
                                break
                        if found_in_includes:
                            break
                    
                    if not found_in_includes:
                        print(f"Warning: Included file '{inc_filename}' not found.")
            else:
                lines.append(line)
                
        return lines

    def parse(self) -> Mechanism:
        """Parses the .easy file."""
        self._parse_sections()
        self._parse_reactions_and_rates()
        self.mechanism._extract_species() # Recalculate species
        return self.mechanism

    def _parse_sections(self):
        """Parses sections like FILES, FACS, INITIALS, FUNCTIONS."""
        facs_regex = re.compile(r"FACS\[(\w+)\]\s*=\s*(.*)")
        files_regex = re.compile(r"FILES\[(\w+)\]\s*=\s*(.*)")
        initial_regex = re.compile(r"INITIAL\[(.*?)\]\s*=\s*(.*)")
        function_regex = re.compile(r"(\w+)\s*=\s*FUNCTION\[(.*)") # Modified to be non-greedy
        input_regex = re.compile(r"(\w+)\s*=\s*INPUT\((.*?)\)")
        const_regex = re.compile(r"(\w+)\s*=\s*CONST\((.*?)\)")
        var_regex = re.compile(r"(\w+)\s*=\s*(\w+)\((.*)\)")

        i = 0
        while i < len(self.lines):
            line = self.lines[i].strip()
            i += 1

            if not line or line.startswith(';'):
                continue

            if match := facs_regex.match(line):
                self.mechanism.facsimile_params[match.group(1)] = match.group(2).strip()
            elif match := files_regex.match(line):
                self.mechanism.enz_files[match.group(1)] = match.group(2).strip()
            elif match := initial_regex.match(line):
                self.mechanism.initial_conditions[match.group(1)] = match.group(2).strip()
            elif match := function_regex.match(line):
                name = match.group(1)
                body = match.group(2)
                open_brackets = body.count('[') - body.count(']')
                
                # Handle multi-line functions
                while open_brackets > 0 and i < len(self.lines):
                    next_line = self.lines[i].strip()
                    body += next_line
                    open_brackets += next_line.count('[') - next_line.count(']')
                    i += 1
                
                # Clean up the body
                body = body.replace('$', '').strip()
                if body.endswith(']'):
                    body = body[:-1]

                self.mechanism.functions[name] = f"FUNCTION[{body}]"

            elif match := input_regex.match(line):
                self.mechanism.variables[match.group(1)] = f"INPUT({match.group(2).strip()})"
            elif match := const_regex.match(line):
                self.mechanism.variables[match.group(1)] = f"CONST({match.group(2).strip()})"
            elif match := var_regex.match(line):
                 self.mechanism.variables[match.group(1)] = f"{match.group(2)}({match.group(3).strip()})"


    def _parse_reactions_and_rates(self):
        """Parses reactions and their rate constants."""
        # Regex to find reaction (e.g., A+B-->C+D) and rate (e.g., k[A+B-->C+D]=rate)
        reaction_regex = re.compile(r"(.+?)\s*-->\s*(.+)")
        rate_regex = re.compile(r"k\[(.+?)\s*-->\s*(.+)\]\s*=\s*(.+)")

        parsed_reactions = []

        def parse_species_list(s_list_str):
            species_list = []
            for part in s_list_str.split('+'):
                part = part.strip()
                if '*' in part:
                    coeff_str, name = part.split('*', 1)
                    try:
                        # We don't use the coefficient for now, just the species name
                        species_list.append(name.strip())
                    except ValueError:
                        species_list.append(part) # Not a valid coefficient, treat as part of the name
                else:
                    species_list.append(part)
            return tuple(species_list)

        for line in self.lines:
            line = line.split(';')[0].strip()
            if not line:
                continue

            rate_match = rate_regex.match(line)
            if rate_match:
                reactants_str = rate_match.group(1).strip()
                products_str = rate_match.group(2).strip()
                rate = rate_match.group(3).strip()

                reactants = parse_species_list(reactants_str)
                products = parse_species_list(products_str)
                
                reaction_key = (reactants, products)
                
                # Find the first matching reaction without a rate assigned yet
                assigned = False
                for r in parsed_reactions:
                    if (r['reactants'], r['products']) == reaction_key and r.get('rate') is None:
                        r['rate'] = rate
                        assigned = True
                        break
                
                if not assigned:
                    parsed_reactions.append({'reactants': reactants, 'products': products, 'rate': rate})
                continue

            reaction_match = reaction_regex.match(line)
            if reaction_match:
                # Check if it's not a rate line
                if line.startswith('k['):
                    continue

                reactants_str = reaction_match.group(1).strip()
                products_str = reaction_match.group(2).strip()
                
                reactants = parse_species_list(reactants_str)
                products = parse_species_list(products_str)

                parsed_reactions.append({'reactants': reactants, 'products': products, 'rate': None})
        
        for data in parsed_reactions:
            rate = data.get('rate')
            if rate is None:
                rate = '0.0'
            self.mechanism.reactions.append(Reaction(data['reactants'], data['products'], rate))

if __name__ == '__main__':
    parser = EasyParser('NOx_O3.easy')
    mechanism = parser.parse()
    
    print("--- Parsed Mechanism ---")
    print("\nENZ Files:")
    print(mechanism.enz_files)
    print("\nFACSIMILE Params:")
    print(mechanism.facsimile_params)
    print("\nInitial Conditions:")
    print(mechanism.initial_conditions)
    print("\nFunctions:")
    print(mechanism.functions)
    print("\nVariables:")
    print(mechanism.variables)
    print("\nReactions:")
    for reaction in mechanism.reactions:
        print(f"  {'+'.join(reaction.reactants)} --> {'+'.join(reaction.products)}  (k={reaction.rate_constant})")
    print("\nSpecies:")
    for species in mechanism.species:
        print(f"  {species.name}")
    print("\n------------------------")
