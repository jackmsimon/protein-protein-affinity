import requests
import os
import json
from pathlib import Path
import sys

key = '' # NVIDIA API key would need to be added here

def get_reduced_pdb():
    pdb = Path("/Users/jacksimon/Bioengineering DSL/files/FGFR2_dimer.pdb")
    lines = []
    first_res = {'A': None, 'B': None}
    
    for line in pdb.read_text().split("\n"):
        if line.startswith("ATOM"):
            chain = line[21]
            res_num = int(line[22:26])
            if first_res[chain] is None:
                first_res[chain] = res_num
            lines.append(line)
    
    print(f"First residue numbers: {first_res}")  # Debug info
    return "\n".join(lines), first_res

pdb_lines, first_residues = get_reduced_pdb()

# Adjust contig string based on first residue numbers
r = requests.post(
    url='https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate',
    headers={'Authorization': f'Bearer {key}'},
    json={
        'input_pdb': pdb_lines,
        'contigs': 'A16-144 A149-357 B16-144 B149-292/0 50-80',
        'hotspot_res': ['A175', 'A198', 'A207', 'B198', 'B223'],
        'diffusion_steps': 20
    }
)

# Create outputs directory if it doesn't exist
outputs_dir = Path(__file__).parent
outputs_dir.mkdir(parents=True, exist_ok=True)

# Write the output PDB with error handling
try:
    response_data = json.loads(r.text)
    if r.status_code != 200:
        print(f"Error: API request failed with status {r.status_code}")
        print(f"Response: {response_data}")
        sys.exit(1)
    
    if 'output_pdb' not in response_data:
        print(f"Error: No output_pdb in response. Full response: {response_data}")
        sys.exit(1)
        
    (outputs_dir / "design_output.pdb").write_text(response_data["output_pdb"])
    print("Design successfully generated and saved")
    
except Exception as e:
    print(f"Error processing API response: {e}")
    print(f"Response text: {r.text}")
    sys.exit(1)
