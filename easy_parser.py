import re
from mechanism import Mechanism, Reaction

class EasyParser:
    """Parses an .easy file and populates a Mechanism object."""

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.mechanism = Mechanism()
        self.lines = self._read_lines()

    def _read_lines(self) -> list[str]:
        """Reads lines from the file."""
        with open(self.file_path, 'r') as f:
            return f.readlines()

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

        reactions = {}

        for line in self.lines:
            line = line.strip()
            if not line or line.startswith(';'):
                continue

            rate_match = rate_regex.match(line)
            if rate_match:
                reactants_str = rate_match.group(1).strip()
                products_str = rate_match.group(2).strip()
                rate = rate_match.group(3).split(';')[0].strip() # Remove comments

                reactants = tuple(r.strip() for r in reactants_str.split('+'))
                products = tuple(p.strip() for p in products_str.split('+'))
                
                # Store rate with a canonical representation of the reaction
                reaction_key = (reactants, products)
                if reaction_key in reactions:
                    reactions[reaction_key]['rate'] = rate
                else:
                    reactions[reaction_key] = {'rate': rate}
                continue # Go to next line after finding a rate

            reaction_match = reaction_regex.match(line)
            if reaction_match:
                # Check if it's not a rate line
                if line.startswith('k['):
                    continue

                reactants_str = reaction_match.group(1).strip()
                products_str = reaction_match.group(2).strip()
                
                reactants = tuple(r.strip() for r in reactants_str.split('+'))
                products = tuple(p.strip() for p in products_str.split('+'))

                reaction_key = (reactants, products)
                if reaction_key not in reactions:
                    reactions[reaction_key] = {}
        
        for (reactants, products), data in reactions.items():
            rate = data.get('rate', '0.0') # Default rate if not found
            self.mechanism.reactions.append(Reaction(reactants, products, rate))

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
