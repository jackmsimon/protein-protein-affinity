import os
import sys
from pathlib import Path
import json
import subprocess

# Import from local analysis module
from analysis import analyze_interface_lining, analyze_pockets, get_missing_residues

def run_pipeline(interaction_pdb: str) -> None:
    """Run binding affinity increase pipeline.
    
    Args:
        interaction_pdb: Path to the PDB file to analyze
    """
    try:
        # Convert to absolute path if relative
        pdb_path = str(Path(interaction_pdb).resolve())
        
        # Setup output directory
        output_dir = Path(__file__).parent / "analyses"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Run analyses
        print("Running interface lining analysis...")
        interface_results = analyze_interface_lining(pdb_path, output_dir)
        
        print("Running pocket analysis...")
        pocket_results = analyze_pockets(pdb_path, output_dir)
        
        print("Running missing residues analysis...")
        missing_residues = get_missing_residues(pdb_path, output_dir)
        
        # Combine results
        combined_results = {
            "interface_analysis": interface_results,
            "pocket_analysis": pocket_results,
            "missing_residues": missing_residues,
            "pdb_path": pdb_path
        }
        
        # Save combined results
        combined_output = output_dir / "combined_analysis.json"
        with open(combined_output, "w") as f:
            json.dump(combined_results, f, indent=2)
            
        print(f"Analysis complete. Results saved to: {output_dir}")
        
        # Run design with the same PDB path
        design_script = Path(__file__).parent / "design.py"
        subprocess.run(['python', str(design_script), pdb_path], check=True)
        
    except Exception as e:
        print(f"Error in pipeline: {e}")
        raise

# Only run if called directly, not when imported
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python main.py <pdb_path>")
        sys.exit(1)
    run_pipeline(sys.argv[1]) 