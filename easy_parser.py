import re
from mechanism import Mechanism, Reaction

class EasyParser:
    """Parses an .easy file and populates a Mechanism object."""

    def __init__(self, file_path: str, include_dirs: list = None):
        self.file_path = file_path
        import os
        self.include_dirs = include_dirs or [os.path.dirname(os.path.abspath(file_path))]
        self.mechanism = Mechanism()
        self.lines = self._read_lines_with_includes(self.file_path)

    def _read_lines_with_includes(self, current_file: str) -> list[str]:
        """Reads lines from the file and recursively processes #include directives."""
        import os
        lines = []
        try:
            with open(current_file, 'r') as f:
                raw_lines = f.readlines()
        except FileNotFoundError:
            print(f"Warning: Could not read file '{current_file}'.")
            return lines
            
        for line in raw_lines:
            stripped = line.strip()
            if stripped.startswith('#include'):
                parts = stripped.split()
                if len(parts) > 1:
                    inc_filename = parts[1]
                    # Find the included file in the provided directories recursively
                    found_path = None
                    for d in self.include_dirs:
                        for root, _, files in os.walk(d):
                            if inc_filename in files:
                                found_path = os.path.join(root, inc_filename)
                                break
                        if found_path:
                            break
                            
                    if found_path:
                        lines.extend(self._read_lines_with_includes(found_path))
                    else:
                        print(f"Warning: Included file '{inc_filename}' not found in include directories.")
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
        section_regex = re.compile(r"FACS\[(\w+)\]\s*=\s*(.*)|"
                                   r"FILES\[(\w+)\]\s*=\s*(.*)|"
                                   r"INITIAL\[(.*?)\]\s*=\s*(.*)|"
                                   r"(\w+)\s*=\s*FUNCTION\[(.*?)\]|"
                                   r"(\w+)\s*=\s*INPUT\((\w+)\)|"
                                   r"(\w+)\s*=\s*(\w+)\((.*)\)")

        for line in self.lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue

            match = section_regex.match(line)
            if not match:
                continue

            groups = match.groups()
            if groups[0] is not None:
                self.mechanism.facsimile_params[groups[0]] = groups[1].strip()
            elif groups[2] is not None:
                self.mechanism.enz_files[groups[2]] = groups[3].strip()
            elif groups[4] is not None:
                self.mechanism.initial_conditions[groups[4]] = groups[5].strip()
            elif groups[6] is not None:
                self.mechanism.functions[groups[6]] = f"FUNCTION[{groups[7].strip()}]"
            elif groups[8] is not None:
                self.mechanism.variables[groups[8]] = f"INPUT({groups[9].strip()})"
            elif groups[10] is not None and groups[11] is not None:
                self.mechanism.variables[groups[10]] = f"{groups[11]}({groups[12].strip()})"


    def _parse_reactions_and_rates(self):
        """Parses reactions and their rate constants."""
        # Regex to find reaction (e.g., A+B-->C+D) and rate (e.g., k[A+B-->C+D]=rate)
        reaction_regex = re.compile(r"(.+?)\s*-->\s*(.+)")
        rate_regex = re.compile(r"k\[(.+?)\s*-->\s*(.+)\]\s*=\s*(.+)")

        parsed_reactions = []

        for line in self.lines:
            line = line.split(';')[0].strip()
            if not line:
                continue

            rate_match = rate_regex.match(line)
            if rate_match:
                reactants_str = rate_match.group(1).strip()
                products_str = rate_match.group(2).strip()
                rate = rate_match.group(3).strip()

                reactants = tuple(r.strip() for r in reactants_str.split('+'))
                products = tuple(p.strip() for p in products_str.split('+'))
                
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
                
                reactants = tuple(r.strip() for r in reactants_str.split('+'))
                products = tuple(p.strip() for p in products_str.split('+'))

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
