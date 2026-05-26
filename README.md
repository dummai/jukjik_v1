# jukjik_v1

Anti-DENV drug combination efficacy and safety evaluation.

## Overview

This project provides a complete pipeline for evaluating drug combinations against Dengue virus (DENV). It calculates percent inhibition from raw experimental data, converts to SynergyFinder format, and performs drug synergy analysis using multiple mathematical models.

## Installation

```bash
pip install -r requirements.txt
```

Dependencies:
- openpyxl >= 3.1.0
- numpy >= 1.24.0
- pandas >= 2.0.0
- scipy >= 1.10.0
- matplotlib >= 3.7.0
- seaborn >= 0.12.0

## Usage

### Step 1: Generate Experiment Template (Optional)

```bash
python example_usage.py
```

This reads from `GenExcel.xlsx` and generates a template for data entry.

### Step 2: Calculate Percent Inhibition

```bash
python calculate_inhibition.py <input.xlsx> [output.xlsx]
```

- `input.xlsx`: Raw experimental data with replicates
- `output.xlsx`: Output file with % inhibition (default: PercInhibition.xlsx)

### Step 3: Convert to SynergyFinder Format

```bash
python synergyfinder_input.py <PercInhibition.xlsx> [output.csv]
```

Options:
- `--drug-pairs DrugA:DrugB DrugA:DrugC`: Specify drug pairs (optional, auto-detected if omitted)

Output:
- `SynergyFinder_input.csv`: Main input file
- `SynergyFinder_BlockID.csv`: Block summary

### Step 4: Calculate Synergy Scores

```bash
python synergy_analysis.py <input.csv> [output.csv]
```

Calculates synergy using four models:
- **ZIP** (Zero Interaction Potency)
- **Bliss Independence**
- **Loewe Additivity**
- **HSA** (Highest Single Agent)

Output:
- `synergy_score_table.csv`: Detailed scores per concentration
- `synergy_summary.csv`: Summary statistics

### Step 5: Generate Heatmap

```bash
python heatmap_generator.py
```

Creates publication-quality synergy heatmaps.

## Input/Output Formats

### RawInput.xlsx
| Single drugs | Drugs | Rep1 | Rep2 | Rep3 |
|--------------|-------|------|------|------|
| Control | | 100 | 105 | 98 |
| DrugA | 0.78 | 85 | 82 | 88 |
| ... | ... | ... | ... | ... |

### SynergyFinder Input CSV
| block_id | drug1 | drug2 | conc1 | conc2 | response | conc_unit |
|----------|-------|-------|-------|-------|----------|-----------|
| 1 | DrugA | DrugB | 0 | 0 | 0 | uM |
| 1 | DrugA | DrugB | 1.56 | 0 | 12.5 | uM |
| ... | ... | ... | ... | ... | ... | ... |

### Synergy Output CSV
| block_id | conc1 | conc2 | ZIP_synergy | Bliss_synergy | Loewe_synergy | HSA_synergy |
|----------|-------|-------|-------------|---------------|---------------|-------------|
| 1 | 1.56 | 3.12 | 5.2 | 3.8 | 4.1 | 2.5 |
| ... | ... | ... | ... | ... | ... | ... |

## Synergy Interpretation

- **Synergy score > 5**: Synergistic effect
- **Synergy score < -5**: Antagonistic effect
- **-5 <= score <= 5**: Additive effect

## References

- Yadav B, et al. (2015) Searching for Drug Synergy in Complex Dose-Response Landscape Using an Interaction Potency Model. Computational and Structural Biotechnology Journal.
- Ritz C, et al. (2015) Dose-Response Analysis Using R. PLoS ONE.

## License

MIT License
