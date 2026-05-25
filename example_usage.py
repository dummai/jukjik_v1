"""
Example usage of drug_combination_template.py

This script demonstrates how to use the drug combination template generator.
"""

from drug_combination_template import generate_drug_combination_excel


def main():
    # Define your experiments
    experiments = [
        {
            "name": "DrugA_DrugB",
            "drug_a": "DrugA",
            "drug_b": "DrugB",
            "conc_a": [10, 20, 50],      # μM concentrations
            "conc_b": [5, 10, 20],       # μM concentrations
            "replicates": 3
        },
        {
            "name": "DrugA_DrugC",
            "drug_a": "DrugA",
            "drug_b": "DrugC",
            "conc_a": [10, 20, 50],
            "conc_b": [1, 5, 10],
            "replicates": 3
        },
        {
            "name": "DrugB_DrugC",
            "drug_a": "DrugB",
            "drug_b": "DrugC",
            "conc_a": [5, 10, 20],
            "conc_b": [1, 5, 10],
            "replicates": 3
        }
    ]
    
    # Generate the Excel file
    output_path = "drug_combination_results.xlsx"
    generate_drug_combination_excel(
        experiments=experiments,
        output_path=output_path
    )
    
    print(f"\nGenerated: {output_path}")
    print(f"Worksheets: {len(experiments)}")
    for exp in experiments:
        print(f"  - {exp['name']}: {exp['drug_a']} x {exp['drug_b']}")


if __name__ == "__main__":
    main()
