from Bio.PDB import *
from Bio.PDB.NACCESS import NACCESS_atomic
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from collections import defaultdict
import numpy as np
from pathlib import Path
import json
import subprocess
import math
from Bio.PDB import PDBParser

class InterfaceAnalyzer:
    def __init__(self, pdb_paths, output_dir="output", cutoff_distance=10.0):
        """Initialize the interface analyzer with PDB file path(s).
        
        Args:
            pdb_paths: Either a single PDB path (string) for dimer analysis
                      or a tuple/list of two PDB paths for separate chains
            output_dir: Directory to store analysis results
            cutoff_distance: Distance cutoff for interface detection in Angstroms
        """
        self.parser = PDBParser(QUIET=True)
        self.cutoff_distance = cutoff_distance
        self.interface_residues = defaultdict(list)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Handle single PDB (dimer) or dual PDB (separate chains) cases
        if isinstance(pdb_paths, (str, Path)):
            self.mode = "dimer"
            self.structure = self.parser.get_structure("complex", pdb_paths)
        else:
            self.mode = "separate"
            self.structure_A = self.parser.get_structure("chainA", pdb_paths[0])
            self.structure_B = self.parser.get_structure("chainB", pdb_paths[1])
            # Combine structures for analysis
            self.structure = self.merge_structures()
    
    def merge_structures(self):
        """Merge two separate PDB structures into one for analysis."""
        # Create a new empty structure
        merged = Structure("merged")
        model = Model(0)
        
        # Add chains from structure A
        for chain in self.structure_A[0]:
            chain.id = "A"  # Force chain ID to A
            model.add(chain)
            
        # Add chains from structure B
        for chain in self.structure_B[0]:
            chain.id = "B"  # Force chain ID to B
            model.add(chain)
            
        return merged
        
    def calculate_distance(self, atom1, atom2):
        """Calculate distance between two atoms."""
        return np.linalg.norm(atom1.coord - atom2.coord)
    
    def is_interface_residue(self, residue1, chain2):
        """Check if a residue is part of the interface with another chain."""
        for atom1 in residue1:
            if atom1.name == "CA":  # Only check alpha carbons for efficiency
                for residue2 in chain2.get_residues():
                    for atom2 in residue2:
                        if atom2.name == "CA":
                            if self.calculate_distance(atom1, atom2) <= self.cutoff_distance:
                                return True
        return False
    
    def is_perimeter_residue(self, residue, chain, interface_residues):
        """Check if a residue is on the perimeter of the interface.
        
        A residue is considered to be on the perimeter if it has fewer than 
        4 interface residue neighbors within 8 Angstroms.
        """
        neighbor_count = 0
        res_id = residue.get_id()[1]
        
        for other_res in chain.get_residues():
            if other_res.get_id()[1] == res_id:
                continue
            
            other_tuple = (other_res.get_resname(), chain.id, other_res.get_id()[1])
            if other_tuple in interface_residues:
                # Check if they're neighbors (using CA atoms)
                ca1 = residue["CA"]
                ca2 = other_res["CA"]
                if self.calculate_distance(ca1, ca2) <= 8.0:  # 8A neighbor cutoff
                    neighbor_count += 1
                
        return neighbor_count < 4  # Perimeter residues have fewer neighbors
    
    def find_interface_residues(self):
        """Identify interface residues between chains."""
        chains = list(self.structure[0].get_chains())
        
        # First pass - find all interface residues
        for i, chain1 in enumerate(chains):
            for chain2 in chains[i+1:]:
                for residue in chain1.get_residues():
                    if residue.get_resname() == "HOH":  # Skip water molecules
                        continue
                    if self.is_interface_residue(residue, chain2):
                        self.interface_residues[chain1.id].append(
                            (residue.get_resname(), 
                             chain1.id,
                             residue.get_id()[1])
                        )
                        
                for residue in chain2.get_residues():
                    if residue.get_resname() == "HOH":
                        continue
                    if self.is_interface_residue(residue, chain1):
                        self.interface_residues[chain2.id].append(
                            (residue.get_resname(), 
                             chain2.id,
                             residue.get_id()[1])
                        )
        
        # Second pass - filter for perimeter residues
        perimeter_residues = defaultdict(list)
        for chain_id, residues in self.interface_residues.items():
            chain = next(chain for chain in self.structure[0].get_chains() if chain.id == chain_id)
            for res_name, chain_id, res_num in residues:
                residue = next(res for res in chain.get_residues() if res.get_id()[1] == res_num)
                if self.is_perimeter_residue(residue, chain, set(residues)):
                    perimeter_residues[chain_id].append((res_name, chain_id, res_num))
        
        self.interface_residues = perimeter_residues
    
    def get_interface_residues(self):
        """Return the interface residues in a formatted way."""
        if not self.interface_residues:
            self.find_interface_residues()
            
        formatted_residues = []
        for chain_id in sorted(self.interface_residues.keys()):
            residues = sorted(self.interface_residues[chain_id], key=lambda x: x[2])
            for res_name, chain, res_num in residues:
                formatted_residues.append(f"{res_name} {chain}{res_num}")
        return formatted_residues
    
    def save_results(self):
        """Save analysis results to JSON file."""
        # Format residues as strings
        perimeter_residues = []
        for chain_id in sorted(self.interface_residues.keys()):
            for res_name, chain, res_num in sorted(self.interface_residues[chain_id], key=lambda x: (x[1], x[2])):
                perimeter_residues.append(f"{res_name} {chain}{res_num}")
        
        results = {
            "perimeter_residues": perimeter_residues
        }
        
        output_file = self.output_dir / "interface_analysis.json"
        with open(output_file, "w") as f:
            json.dump(results, f, indent=2)
        
        return output_file

def analyze_interface_lining(pdb_path: str, output_dir: str = "analyses") -> dict:
    """Analyze protein interface from PDB file."""
    analyzer = InterfaceAnalyzer(pdb_path, output_dir)
    analyzer.find_interface_residues()  # Make sure interface residues are found
    output_file = analyzer.save_results()
    
    # Format residues as strings in the desired format
    perimeter_residues = []
    for chain_id in sorted(analyzer.interface_residues.keys()):
        for res_name, chain, res_num in sorted(analyzer.interface_residues[chain_id], key=lambda x: (x[1], x[2])):
            perimeter_residues.append(f"{res_name} {chain}{res_num}")
    
    return {
        "perimeter_residues": perimeter_residues,
        "output_file": str(output_file)
    }

def analyze_pockets(pdb_path: str, output_dir: str = "analyses") -> dict:
    """Analyze pockets in both chains using fpocket."""
    # Setup paths
    pdb_path = Path(pdb_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    base_name = pdb_path.stem
    fpocket_dir = pdb_path.parent / f'{base_name}_out'
    info_path = fpocket_dir / f'{base_name}_info.txt'
    pocket_pdb = fpocket_dir / f'{base_name}_out.pdb'
    
    # Run fpocket and capture output
    try:
        subprocess.run(['fpocket', '-f', str(pdb_path)], check=True, capture_output=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running fpocket: {e}")
        return {}
    
    # Read and parse pocket info
    try:
        with open(info_path, 'r') as f:
            pocket_info = parse_pocket_info(f.read())
    except FileNotFoundError:
        print(f"Could not find fpocket output at {info_path}")
        return {}
    
    # Get pocket residues
    pocket_residues = find_all_pocket_residues(str(pocket_pdb))
    
    # Format results in the desired structure
    formatted_results = {}
    for pocket_num in sorted(pocket_info.keys()):
        formatted_residues = [
            f"{res[0]} {res[2]}{res[1]}" 
            for res in pocket_residues.get(pocket_num, [])
        ]
        
        formatted_results[str(pocket_num)] = {
            "info": pocket_info[pocket_num]['info'],
            "residues": formatted_residues
        }
    
    # Save results
    output_file = output_dir / "pocket_analysis.json"
    with open(output_file, "w") as f:
        json.dump(formatted_results, f, indent=2)
        
    return formatted_results

def parse_pdb_line(line):
    """Parse a PDB file line and return relevant atomic information."""
    if line.startswith(("ATOM", "HETATM")):
        return {
            'record_type': line[0:6].strip(),
            'atom_num': int(line[6:11]),
            'atom_name': line[12:16].strip(),
            'res_name': line[17:20].strip(),
            'chain_id': line[21],
            'res_num': int(line[22:26]),
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54])
        }
    return None

def distance(p1, p2):
    """Calculate Euclidean distance between two 3D points."""
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(p1, p2)))

def parse_pocket_info(info_text):
    """Parse the fpocket info file and return pocket information."""
    pockets = {}
    current_pocket = None
    
    for line in info_text.split('\n'):
        line = line.strip()
        if line.startswith('Pocket'):
            if ':' in line:
                current_pocket = line.split(':')[0].strip()
                pocket_num = int(current_pocket.split()[-1])
                pockets[pocket_num] = {'header': current_pocket, 'info': []}
        elif line and current_pocket and '\t' in line:
            pockets[int(current_pocket.split()[-1])]['info'].append(line)
            
    return pockets

def find_all_pocket_residues(pdb_file, cutoff=4.0):
    """Find residues lining all pockets within a cutoff distance."""
    pocket_points = defaultdict(list)
    protein_atoms = []
    pocket_residues = defaultdict(set)
    
    with open(pdb_file, 'r') as f:
        for line in f:
            entry = parse_pdb_line(line)
            if not entry:
                continue
                
            if entry['record_type'] == "HETATM":
                if ("POL STP" in line or "APOL STP" in line):
                    pocket_num = entry['res_num']
                    pocket_points[pocket_num].append((entry['x'], entry['y'], entry['z']))
            
            elif entry['record_type'] == "ATOM":
                protein_atoms.append({
                    'res_name': entry['res_name'],
                    'res_num': entry['res_num'],
                    'chain_id': entry['chain_id'],
                    'coords': (entry['x'], entry['y'], entry['z'])
                })

    for pocket_num, points in pocket_points.items():
        for atom in protein_atoms:
            for point in points:
                if distance(atom['coords'], point) <= cutoff:
                    residue = (atom['res_name'], atom['res_num'], atom['chain_id'])
                    pocket_residues[pocket_num].add(residue)
                    break
                
    return {k: sorted(v, key=lambda x: (x[2], x[1])) for k, v in pocket_residues.items()}

def get_missing_residues(pdb_path: str, output_dir: str = "analyses") -> dict:
    """Detect missing residues by finding gaps in residue numbering."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    missing = {}
    
    for chain in structure[0]:
        chain_id = chain.id
        residue_nums = []
        
        # Get all residue numbers in this chain
        for residue in chain:
            if residue.id[0] == ' ':  # Only standard amino acids
                residue_nums.append(residue.id[1])
        
        # Sort residue numbers
        residue_nums.sort()
        
        # Find gaps in numbering
        missing_nums = []
        for i in range(len(residue_nums) - 1):
            current = residue_nums[i]
            next_res = residue_nums[i + 1]
            # If there's a gap larger than 1, those residues are missing
            if next_res - current > 1:
                missing_nums.extend(range(current + 1, next_res))
        
        if missing_nums:
            missing[chain_id] = missing_nums
    
    # Save to separate JSON file
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    output_file = output_dir / "missing_residues.json"
    
    with open(output_file, "w") as f:
        json.dump(missing, f, indent=2)
            
    return {
        "missing_residues": missing,
        "output_file": str(output_file)
    }

def combine_analyses(interface_analysis: dict, pocket_analysis: dict, missing_residues: dict, output_dir: str = "analyses") -> dict:
    """Combine all analyses into a single JSON file."""
    combined = {
        "interface_analysis": interface_analysis.get("perimeter_residues", []),
        "pocket_analysis": pocket_analysis,
        "missing_residues": missing_residues.get("missing_residues", {}),
        "pdb_path": interface_analysis.get("output_file", "").replace("interface_analysis.json", "").rstrip("/")
    }
    
    output_file = Path(output_dir) / "combined_analysis.json"
    with open(output_file, "w") as f:
        json.dump(combined, f, indent=2)
    
    return combined

if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1:
        pdb_file = sys.argv[1]
        interface_analysis = analyze_interface_lining(pdb_file)
        pocket_analysis = analyze_pockets(pdb_file)
        missing_residues = get_missing_residues(pdb_file)
        combine_analyses(interface_analysis, pocket_analysis, missing_residues)
