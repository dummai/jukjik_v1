"""
Generate SynergyFinder Plus Input Format

Converts PercInhibition.xlsx data to SynergyFinder Plus long table format.

Input: PercInhibition.xlsx (from calculate_inhibition.py)
Output: SynergyFinder_input.csv

Required columns for SynergyFinder Plus:
- block_id: Identifier for drug combination blocks
- drug1: Name of first drug
- drug2: Name of second drug
- conc1: Concentration of first drug
- conc2: Concentration of second drug
- response: % inhibition value
- conc_unit: Unit of concentration (uM)
"""

import re
from openpyxl import load_workbook
import csv


def parse_drug_conc(value):
    """
    Parse drug_concentration format (e.g., "DrugA_0.78").
    
    Returns:
        tuple: (drug_name, concentration) or (None, None) if not parseable
    """
    if value is None:
        return None, None
    
    value_str = str(value).strip()
    
    # Match pattern: DrugName_Concentration
    match = re.match(r'^(Drug[A-Z])[_\s]+(\d+\.?\d*)$', value_str)
    if match:
        drug_name = match.group(1)
        concentration = float(match.group(2))
        return drug_name, concentration
    
    return None, None


def generate_synergyfinder_input(
    input_path: str,
    output_path: str = "SynergyFinder_input.csv"
):
    """
    Convert PercInhibition.xlsx to SynergyFinder Plus format.
    
    Args:
        input_path: Path to PercInhibition.xlsx
        output_path: Output CSV file path
        
    Output format:
        block_id,drug1,drug2,conc1,conc2,response,conc_unit
    """
    
    # Load input workbook
    wb = load_workbook(input_path)
    
    # Collect all data rows
    all_rows = []
    
    # Track block_id
    block_id = 0
    
    # Define drug pairs in order
    drug_pairs = [
        ("DrugA", "DrugB"),
        ("DrugA", "DrugC"),
        ("DrugB", "DrugC")
    ]
    
    # Track block summary for BlockID.csv
    block_summary_data = []
    
    for sheet_name in wb.sheetnames:
        ws = wb[sheet_name]
        print(f"Processing: {sheet_name}")
        
        # Read all rows
        rows = list(ws.iter_rows(values_only=True))
        
        # Get header to find Avg_% column index
        header = rows[0]
        avg_col_idx = None
        for i, col in enumerate(header):
            if col == "Avg_%":
                avg_col_idx = i
                break
        
        if avg_col_idx is None:
            print(f"  Warning: Avg_% column not found, skipping sheet")
            continue
        
        # Track block_ids for this experiment
        exp_block_ids = {}
        
        # Initialize block_ids for each drug pair in this experiment
        for pair in drug_pairs:
            block_id += 1
            exp_block_ids[pair] = block_id
            # Add to block summary
            block_summary_data.append({
                "block_id": block_id,
                "drug1": pair[0],
                "drug2": pair[1],
                "experiment": sheet_name
            })
        
        # Track state for combination rows
        current_drug1 = None
        current_conc1 = None
        
        # Process each data row
        for row_idx, row in enumerate(rows[1:], start=2):
            col_a = row[0]
            col_b = row[1]
            response = row[avg_col_idx]
            
            # Skip control row or None response
            if col_a == "Control" or response is None:
                continue
            
            if col_a is not None:
                # New entry - check if it's combination format
                drug1_check, conc1_check = parse_drug_conc(col_a)
                
                if drug1_check and conc1_check:
                    # Combination row
                    current_drug1 = drug1_check
                    current_conc1 = conc1_check
                    
                    drug2, conc2 = parse_drug_conc(col_b)
                    
                    if drug2 and conc2 is not None:
                        # Find matching drug pair
                        pair_key = None
                        for pair in drug_pairs:
                            if (pair[0] == current_drug1 and pair[1] == drug2) or \
                               (pair[1] == current_drug1 and pair[0] == drug2):
                                pair_key = pair
                                break
                        
                        if pair_key and pair_key in exp_block_ids:
                            all_rows.append({
                                "block_id": exp_block_ids[pair_key],
                                "drug1": current_drug1,
                                "drug2": drug2,
                                "conc1": current_conc1,
                                "conc2": conc2,
                                "response": response,
                                "conc_unit": "uM",
                                "experiment": sheet_name
                            })
            
            elif col_a is None and col_b is not None and current_drug1:
                # Continuation row - same drug1, different conc1/conc2
                drug2, conc2 = parse_drug_conc(col_b)
                
                if drug2 and conc2 and current_drug1 and current_conc1:
                    pair_key = None
                    for pair in drug_pairs:
                        if (pair[0] == current_drug1 and pair[1] == drug2) or \
                           (pair[1] == current_drug1 and pair[0] == drug2):
                            pair_key = pair
                            break
                    
                    if pair_key and pair_key in exp_block_ids:
                        all_rows.append({
                            "block_id": exp_block_ids[pair_key],
                            "drug1": current_drug1,
                            "drug2": drug2,
                            "conc1": current_conc1,
                            "conc2": conc2,
                            "response": response,
                            "conc_unit": "uM",
                            "experiment": sheet_name
                        })
        
        print(f"  Processed rows for {sheet_name}")
    
    # Remove duplicates (keep unique combinations)
    seen = set()
    unique_rows = []
    for row in all_rows:
        key = (row["block_id"], row["drug1"], row["drug2"], row["conc1"], row["conc2"])
        if key not in seen:
            seen.add(key)
            unique_rows.append(row)
    
    # Write to CSV
    with open(output_path, 'w', newline='') as csvfile:
        fieldnames = ["block_id", "drug1", "drug2", "conc1", "conc2", "response", "conc_unit"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in unique_rows:
            writer.writerow({
                "block_id": row["block_id"],
                "drug1": row["drug1"],
                "drug2": row["drug2"],
                "conc1": row["conc1"],
                "conc2": row["conc2"],
                "response": row["response"],
                "conc_unit": row["conc_unit"]
            })
    
    print(f"\nTotal unique rows: {len(unique_rows)}")
    print(f"Saved to: {output_path}")
    
    # Print summary
    print("\nBlock ID Summary:")
    block_summary = {}
    for row in unique_rows:
        key = (row["block_id"], row["drug1"], row["drug2"], row["experiment"])
        if key not in block_summary:
            block_summary[key] = 0
        block_summary[key] += 1
    
    for key, count in sorted(block_summary.items()):
        print(f"  Block {key[0]}: {key[1]}-{key[2]} ({key[3]}) - {count} combinations")
    
    # Write BlockID.csv
    block_id_path = output_path.replace(".csv", "_BlockID.csv").replace("_input", "")
    with open(block_id_path, 'w', newline='') as csvfile:
        fieldnames = ["block_id", "drug1", "drug2", "experiment"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        
        writer.writeheader()
        for row in block_summary_data:
            writer.writerow({
                "block_id": row["block_id"],
                "drug1": row["drug1"],
                "drug2": row["drug2"],
                "experiment": row["experiment"]
            })
    
    print(f"\nBlock ID summary saved to: {block_id_path}")
    
    return output_path


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
        output_file = sys.argv[2] if len(sys.argv) > 2 else "SynergyFinder_input.csv"
    else:
        input_file = "/Users/kwan/Documents/dengue/DENV_AI/PercInhibition.xlsx"
        output_file = "/Users/kwan/Documents/dengue/DENV_AI/jukjik_v1/SynergyFinder_input.csv"
    
    generate_synergyfinder_input(input_file, output_file)
