# Increasing Protein-Protein Affinity

Autonomous system for enhancing protein-protein binding affinity using LLM-generated RFdiffusion inference scripts.

## What it does
Analyzes protein-protein interfaces and automatically generates optimized binders with enhanced binding affinity. Takes a PDB as input, identifies hotspots and druggable pockets, then uses an LLM to determine optimal RFdiffusion parameters for binder generation.

## Quick Start
```bash
python main.py path/to/protein_complex.pdb
```

## Requirements
- Python 3.8+
- BioPython
- fpocket (installed and in PATH)
- Anthropic API key (set `ANTHROPIC_API_KEY` in .env)
- NVIDIA API key (add to generated script)

## Dependencies
```bash
uv install biopython anthropic python-dotenv requests
```

## Workflow
1. **Analysis**: Interface analysis, pocket detection, missing residue identification
2. **Design**: AI-powered parameter selection for RFdiffusion
3. **Generation**: Automated binder creation via NVIDIA RFdiffusion API

## Output
- `analyses/`: Structural analysis results
- `designs/rfdiffusion_design.py`: Generated design script
- `designs/design_output.pdb`: Optimized binder structure

## Warning
Highly experimental. Not recommended for production protein engineering applications. 