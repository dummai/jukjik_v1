"""
Drug Combination Excel Template Generator

Generates an Excel file with worksheets for drug combination experiments.
Each worksheet contains a template for entering drug combination results.
"""

from openpyxl import Workbook
from openpyxl.styles import Font, Alignment, Border, Side, PatternFill
from openpyxl.utils import get_column_letter


def generate_drug_combination_excel(
    experiments: list,
    output_path: str = "drug_combination_results.xlsx"
):
    """
    Generate Excel file with drug combination template.
    
    Args:
        experiments: List of experiment dictionaries
        output_path: Output Excel file path
    
    Each experiment dict:
        {
            "name": "Experiment Name",      # Worksheet name
            "drug_a": "DrugA",              # First drug name
            "drug_b": "DrugB",              # Second drug name
            "conc_a": [10, 20, 50],         # Concentrations for Drug A
            "conc_b": [5, 10, 20],          # Concentrations for Drug B
            "replicates": 3                 # Number of replicates
        }
    
    Structure per worksheet:
        - Row 1: No Drug control
        - Rows 2-N: Drug A alone (all concentrations)
        - Rows N+1-M: Drug B alone (all concentrations)
        - Remaining rows: All combinations (Drug A x Drug B)
        
        Columns:
        - Col A: Drug A (drug_concentration format)
        - Col B: Drug B (drug_concentration format or "No Drug")
        - Col C onwards: Replicate results (empty for manual entry)
    """
    
    wb = Workbook()
    wb.remove(wb.active)  # Remove default sheet
    
    for exp in experiments:
        ws = wb.create_sheet(title=exp["name"])
        
        drug_a = exp["drug_a"]
        drug_b = exp["drug_b"]
        conc_a = sorted(exp["conc_a"])
        conc_b = sorted(exp["conc_b"])
        replicates = exp["replicates"]
        
        # Style definitions
        header_font = Font(bold=True)
        header_fill = PatternFill(start_color="D9E1F2", end_color="D9E1F2", fill_type="solid")
        no_drug_fill = PatternFill(start_color="FCE4D6", end_color="FCE4D6", fill_type="solid")
        thin_border = Border(
            left=Side(style='thin'),
            right=Side(style='thin'),
            top=Side(style='thin'),
            bottom=Side(style='thin')
        )
        center_align = Alignment(horizontal='center', vertical='center')
        
        # Write headers
        headers = ["Drug A", "Drug B"] + [f"Rep{i+1}" for i in range(replicates)]
        for col_idx, header in enumerate(headers, start=1):
            cell = ws.cell(row=1, column=col_idx, value=header)
            cell.font = header_font
            cell.fill = header_fill
            cell.border = thin_border
            cell.alignment = center_align
        
        # Track current row
        current_row = 2
        
        # Row 1: No Drug control
        cell_a = ws.cell(row=current_row, column=1, value="No Drug")
        cell_a.fill = no_drug_fill
        cell_a.border = thin_border
        cell_a.alignment = center_align
        
        cell_b = ws.cell(row=current_row, column=2, value="No Drug")
        cell_b.fill = no_drug_fill
        cell_b.border = thin_border
        cell_b.alignment = center_align
        
        # Empty replicate cells
        for col_idx in range(3, 3 + replicates):
            cell = ws.cell(row=current_row, column=col_idx, value="")
            cell.border = thin_border
        
        current_row += 1
        
        # Drug A alone (Drug A concentrations with Drug B = No Drug)
        for conc in conc_a:
            cell_a = ws.cell(row=current_row, column=1, value=f"{drug_a}_{conc}")
            cell_a.border = thin_border
            cell_a.alignment = center_align
            
            cell_b = ws.cell(row=current_row, column=2, value="No Drug")
            cell_b.fill = no_drug_fill
            cell_b.border = thin_border
            cell_b.alignment = center_align
            
            for col_idx in range(3, 3 + replicates):
                cell = ws.cell(row=current_row, column=col_idx, value="")
                cell.border = thin_border
            
            current_row += 1
        
        # Drug B alone (Drug B concentrations in Column A with Drug B column = No Drug)
        for conc in conc_b:
            cell_a = ws.cell(row=current_row, column=1, value=f"{drug_b}_{conc}")
            cell_a.border = thin_border
            cell_a.alignment = center_align
            
            cell_b = ws.cell(row=current_row, column=2, value="No Drug")
            cell_b.fill = no_drug_fill
            cell_b.border = thin_border
            cell_b.alignment = center_align
            
            for col_idx in range(3, 3 + replicates):
                cell = ws.cell(row=current_row, column=col_idx, value="")
                cell.border = thin_border
            
            current_row += 1
        
        # All combinations (Drug A x Drug B)
        for conc_a_val in conc_a:
            for conc_b_val in conc_b:
                cell_a = ws.cell(row=current_row, column=1, value=f"{drug_a}_{conc_a_val}")
                cell_a.border = thin_border
                cell_a.alignment = center_align
                
                cell_b = ws.cell(row=current_row, column=2, value=f"{drug_b}_{conc_b_val}")
                cell_b.border = thin_border
                cell_b.alignment = center_align
                
                for col_idx in range(3, 3 + replicates):
                    cell = ws.cell(row=current_row, column=col_idx, value="")
                    cell.border = thin_border
                
                current_row += 1
        
        # Adjust column widths
        ws.column_dimensions['A'].width = 15
        ws.column_dimensions['B'].width = 15
        for col_idx in range(3, 3 + replicates):
            ws.column_dimensions[get_column_letter(col_idx)].width = 10
    
    # Save workbook
    wb.save(output_path)
    print(f"Excel file saved to: {output_path}")
    return output_path


if __name__ == "__main__":
    # Example usage
    experiments = [
        {
            "name": "DrugA_DrugB",
            "drug_a": "DrugA",
            "drug_b": "DrugB",
            "conc_a": [10, 20, 50],
            "conc_b": [5, 10, 20],
            "replicates": 3
        },
        {
            "name": "DrugA_DrugC",
            "drug_a": "DrugA",
            "drug_b": "DrugC",
            "conc_a": [10, 20, 50],
            "conc_b": [1, 5, 10],
            "replicates": 3
        }
    ]
    
    generate_drug_combination_excel(
        experiments=experiments,
        output_path="drug_combination_results.xlsx"
    )
