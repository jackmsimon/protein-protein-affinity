# Protein-Protein Affinity Enhancement

An autonomous protein-protein interaction enhancement system that combines structural analysis with LLM-orchestrated RFdiffusion design to discover optimized binding partners with increased affinity. This system embraces computational approaches to dynamically allocate design resources, pushing towards more agentic approaches in protein engineering.

The system integrates BioPython for structural analysis, fpocket for druggable pocket identification, Anthropic's Claude for scientific reasoning and parameter optimization, and NVIDIA's RFdiffusion API for protein structure generation. The pipeline automatically analyzes protein-protein interfaces, identifies critical binding hotspots, and generates complete RFdiffusion scripts for each design hypothesis.

The workflow requires no manual intervention once initiated, automatically progressing from structural analysis through hypothesis generation to final binder design output. Please note that this work is highly exploratory, and should be altered considerably for production use.

## System Architecture

The autonomous pipeline consists of three integrated components:

1. **Structural Analysis Engine**: Comprehensive interface analysis using BioPython, including binding hotspot identification, missing residue detection, and druggable pocket characterization via fpocket integration.

2. **AI-Driven Design Orchestrator**: LLM-powered parameter optimization that analyzes structural features and generates hypothesis-driven RFdiffusion configurations tailored to specific binding enhancement strategies.

3. **Automated Binder Generation**: Direct integration with NVIDIA's RFdiffusion API for autonomous protein structure generation based on optimized design parameters.

## Quick Start

```bash
python main.py path/to/protein_complex.pdb
```

## Installation

```bash
# Install Python dependencies
uv install biopython anthropic python-dotenv requests

# Install fpocket (required for pocket analysis)
# macOS:
brew install fpocket
# Linux: download from https://fpocket.sourceforge.net/
```

## Configuration

Set up your environment variables:

```bash
# .env file
ANTHROPIC_API_KEY=your_anthropic_api_key_here
```

Note: NVIDIA API key must be added to the generated RFdiffusion script.

## Workflow

1. **Interface Analysis**: Automated identification of protein-protein binding interfaces, calculation of interaction energies, and characterization of binding hotspots.

2. **Pocket Detection**: Integration with fpocket to identify druggable binding pockets and potential allosteric sites for enhanced binding.

3. **Missing Residue Analysis**: Detection and analysis of missing residues that could impact binding interface completeness.

4. **AI-Powered Design Strategy**: LLM analysis of structural features to determine optimal RFdiffusion parameters, including contiguous length, guidance scale, and binding constraints.

5. **Automated Script Generation**: Creation of complete, executable RFdiffusion inference scripts with optimized parameters.

6. **Binder Generation**: Execution via NVIDIA's RFdiffusion API to generate optimized binding partners.

## Output Structure

```
analyses/
├── interface_analysis.json      # Binding interface characterization
├── pocket_analysis.json         # Druggable pocket identification
├── missing_residues.json       # Missing residue analysis
└── combined_analysis.json       # Integrated structural analysis

designs/
├── rfdiffusion_design.py       # Generated RFdiffusion script
└── design_output.pdb           # Final optimized binder structure
```

## Technical Requirements

- Python 3.8+
- BioPython for structural analysis
- fpocket (installed and accessible in PATH)
- Anthropic API access for Claude integration
- NVIDIA API access for RFdiffusion execution

## Research Context

This system builds on advances in protein design methodology, particularly the RFdiffusion framework for protein structure generation. The approach integrates automated structural analysis with AI-driven parameter optimization to enable autonomous discovery of enhanced protein-protein binding partners.

## Limitations and Warnings

This work is highly experimental and represents an exploratory approach to autonomous protein design. The system is not validated for production protein engineering applications and should be used primarily for research and development purposes. Generated designs require extensive experimental validation before any practical application. 