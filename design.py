import requests
import os
import json
from pathlib import Path
import subprocess
import anthropic
from dotenv import load_dotenv
import sys

# Load environment variables from .env file
load_dotenv()

# Create client with API key from environment
client = anthropic.Anthropic(
    api_key=os.getenv('ANTHROPIC_API_KEY')
)

def get_design_params(analysis_json: dict, instructions_path: Path, pdb_path: str) -> tuple[str, str]:
    """Use LLM to determine contig and hotspot parameters based on analysis and instructions."""
    
    # First get the PDB residue numbering
    def get_residue_range():
        pdb = Path(pdb_path)
        first_res = {'A': None, 'B': None}
        last_res = {'A': None, 'B': None}
        
        for line in pdb.read_text().split("\n"):
            if line.startswith("ATOM"):
                chain = line[21]
                res_num = int(line[22:26])
                if first_res[chain] is None:
                    first_res[chain] = res_num
                last_res[chain] = res_num
        return first_res, last_res

    first_res, last_res = get_residue_range()
    
    # Read instructions
    with open(instructions_path, 'r') as f:
        instructions = f.read()
    
    # Create prompt combining analysis and instructions
    prompt = f"""Based on this protein analysis:
{json.dumps(analysis_json, indent=2)}

And following these instructions:
{instructions}

IMPORTANT PDB INFORMATION:
Chain A starts at residue {first_res['A']}
Chain B starts at residue {first_res['B']}

Determine the optimal contig and hotspot residues for RFdiffusion design.
Return ONLY a JSON object in this exact format:
{{
    "contig": "string defining the contig regions (e.g. 'B16-357/0 50-80')",
    "hotspots": ["A50","A51","A52"]
}}

IMPORTANT: 
1. The hotspots must be a list of strings
2. Use exact PDB residue numbers (chain A and B both start at 16)
3. Return only the JSON object, no other text
"""

    # Get LLM response
    message = client.messages.create(
        model="claude-3-7-sonnet-20250219",
        max_tokens=2000,
        temperature=0,
        system="You are an expert in protein design and RFdiffusion.",
        messages=[{"role": "user", "content": prompt}]
    )

    # Parse response
    try:
        text = message.content[0].text.strip()
        # Extract JSON part (between curly braces)
        json_start = text.find('{')
        json_end = text.rfind('}') + 1
        if json_start == -1 or json_end == 0:
            raise ValueError("No JSON object found in response")
            
        json_text = text[json_start:json_end]
        response = json.loads(json_text)
        
        # Ensure hotspots is a comma-separated string
        if isinstance(response["hotspots"], list):
            response["hotspots"] = ",".join(response["hotspots"])
            
        return response["contig"], response["hotspots"]
        
    except Exception as e:
        print(f"Error parsing LLM response: {e}")
        print(f"Raw response:\n{text}")
        raise

def get_rfdiffusion_script(contig: str, hotspots: str, pdb_path: str) -> str:
    """Generate RFdiffusion script with specified parameters."""
    # Convert hotspots string to list
    hotspot_list = [h.strip() for h in hotspots.split(',')]
    
    return f"""import requests
import os
import json
from pathlib import Path
import sys

key = '' # Add NVIDIA API key here

def get_reduced_pdb():
    pdb = Path("{pdb_path}")
    lines = []
    first_res = {{'A': None, 'B': None}}
    
    for line in pdb.read_text().split("\\n"):
        if line.startswith("ATOM"):
            chain = line[21]
            res_num = int(line[22:26])
            if first_res[chain] is None:
                first_res[chain] = res_num
            lines.append(line)
    
    print(f"First residue numbers: {{first_res}}")  # Debug info
    return "\\n".join(lines), first_res

pdb_lines, first_residues = get_reduced_pdb()

# Adjust contig string based on first residue numbers
r = requests.post(
    url='https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate',
    headers={{'Authorization': f'Bearer {{key}}'}},
    json={{
        'input_pdb': pdb_lines,
        'contigs': '{contig}',
        'hotspot_res': {hotspot_list},
        'diffusion_steps': 20
    }}
)

# Create outputs directory if it doesn't exist
outputs_dir = Path(__file__).parent
outputs_dir.mkdir(parents=True, exist_ok=True)

# Write the output PDB with error handling
try:
    response_data = json.loads(r.text)
    if r.status_code != 200:
        print(f"Error: API request failed with status {{r.status_code}}")
        print(f"Response: {{response_data}}")
        sys.exit(1)
    
    if 'output_pdb' not in response_data:
        print(f"Error: No output_pdb in response. Full response: {{response_data}}")
        sys.exit(1)
        
    (outputs_dir / "design_output.pdb").write_text(response_data["output_pdb"])
    print("Design successfully generated and saved")
    
except Exception as e:
    print(f"Error processing API response: {{e}}")
    print(f"Response text: {{r.text}}")
    sys.exit(1)
"""

def main():
    if len(sys.argv) != 2:
        print("Usage: python design.py <pdb_path>")
        sys.exit(1)
        
    pdb_path = sys.argv[1]
    
    # Setup paths
    base_dir = Path(__file__).parent
    designs_dir = base_dir / "designs"
    designs_dir.mkdir(parents=True, exist_ok=True)
    
    # Load analysis results
    analysis_file = base_dir / "analyses" / "combined_analysis.json"
    with open(analysis_file, 'r') as f:
        analysis = json.load(f)
    
    # Get instructions file
    instructions_file = base_dir / "RFdiffusion_Context.txt"
    
    # Get design parameters from LLM
    contig, hotspots = get_design_params(analysis, instructions_file, pdb_path)
    
    # Generate design script
    script = get_rfdiffusion_script(
        contig=contig,
        hotspots=hotspots,
        pdb_path=pdb_path
    )
    
    # Save script
    script_path = designs_dir / "rfdiffusion_design.py"
    with open(script_path, "w") as f:
        f.write(script)
    
    # Run the script
    result = subprocess.run(
        ['python', str(script_path)],
        capture_output=True,
        text=True
    )
    print("\nDesign Results:")
    print("STDOUT:", result.stdout)
    print("STDERR:", result.stderr)

if __name__ == "__main__":
    main()
