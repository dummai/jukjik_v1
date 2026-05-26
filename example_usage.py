"""
Example usage of drug_combination_template.py

This script demonstrates how to use the drug combination template generator.
"""

from drug_combination_template import generate_from_excel


def main():
    # Method 1: Read from Excel input file
    print("=== Method 1: Read from GenExcel.xlsx ===")
    print("Usage: generate_from_excel('GenExcel.xlsx', 'output.xlsx')")
 


def main():
    # Generate from input Excel file
    generate_from_excel(
        input_path="GenExcel.xlsx",
        output_path="drug_combination_results.xlsx"
    )


if __name__ == "__main__":
    main()
