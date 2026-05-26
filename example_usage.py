"""
Example usage of drug_combination_template.py

This script demonstrates how to use the drug combination template generator.
"""

from drug_combination_template import generate_from_excel


def main():
    print("=== Generating Drug Combination Template ===")
    generate_from_excel(
        input_path="GenExcel.xlsx",
        output_path="drug_combination_results.xlsx"
    )


if __name__ == "__main__":
    main()
