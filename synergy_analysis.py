#!/usr/bin/env python3
"""
Drug Synergy Analysis - Matching SynergyFinder Plus Output Format

Calculates drug synergy scores using ZIP, Bliss, Loewe, and HSA models.
Implements the same algorithms as R's SynergyFinder package.

Input: CSV with columns: block_id, drug1, drug2, conc1, conc2, response, conc_unit
Output: CSV with columns: block_id, conc1, conc2, ZIP_fit, ZIP_ref, ZIP_synergy,
        HSA_ref, HSA_synergy, Loewe_ref, Loewe_synergy, Loewe_ci, Bliss_ref, Bliss_synergy

References:
- Yadav B, Wennerberg K, Aittokallio T, Tang J. Searching for Drug Synergy in 
  Complex Dose-Response Landscape Using an Interaction Potency Model.
  Computational and Structural Biotechnology Journal 2015; 13: 504-513.
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve, curve_fit
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')


def hill_equation_4pl(dose: np.ndarray, hill: float, e0: float, emax: float, ec50: float) -> np.ndarray:
    """
    4-parameter logistic (4PL) equation - standard form used in SynergyFinder.
    
    E(d) = e0 + (emax - e0) / (1 + (dose/ec50)^hill)
    
    Parameters are ordered as: hill (slope), e0 (min), emax (max), ec50
    """
    dose = np.atleast_1d(dose).astype(float)
    result = np.zeros_like(dose)
    
    zero_mask = dose == 0
    nonzero_mask = dose > 0
    
    result[zero_mask] = e0
    
    if np.any(nonzero_mask):
        d = dose[nonzero_mask]
        result[nonzero_mask] = e0 + (emax - e0) / (1.0 + np.power(d / ec50, hill))
    
    return result


def fit_4pl_curve(doses: np.ndarray, responses: np.ndarray, 
                  fixed_lower: float = None, fixed_upper: float = 100.0) -> Optional[Dict]:
    """
    Fit 4PL curve to dose-response data matching R's drc::drm with L.4.
    
    Parameters:
    -----------
    doses : array of concentrations (>0)
    responses : array of responses
    fixed_lower : if provided, fix the lower bound (E0) to this value
    fixed_upper : if provided, fix the upper bound (Emax) to this value (default 100)
    
    Returns:
    --------
    dict with parameters: Hill, E0, Emax, EC50
    """
    valid_mask = ~(np.isnan(doses) | np.isnan(responses)) & (doses > 0)
    doses = doses[valid_mask]
    responses = responses[valid_mask]
    
    if len(doses) < 2:
        return None
    
    responses = np.maximum(responses, 0)
    
    if fixed_lower is not None:
        e0_fixed = fixed_lower
    else:
        e0_fixed = np.min(responses) if len(responses) > 0 else 0
    
    emax_fixed = fixed_upper
    
    ec50_init = np.median(doses)
    hill_init = 1.0
    
    def fit_func(d, hill, ec50):
        return hill_equation_4pl(d, hill, e0_fixed, emax_fixed, ec50)
    
    try:
        popt, _ = curve_fit(
            fit_func, doses, responses,
            p0=[hill_init, ec50_init],
            bounds=([0.01, min(doses)/10], [10.0, max(doses)*10]),
            maxfev=10000
        )
        
        return {
            'Hill': popt[0],
            'E0': e0_fixed,
            'Emax': emax_fixed,
            'EC50': popt[1]
        }
    except Exception:
        return None


def predict_4pl(dose: float, params: Dict) -> float:
    """Predict response at given dose using 4PL parameters."""
    if params is None or dose <= 0:
        return params['E0'] if params else 0.0
    result = hill_equation_4pl(np.array([dose]), params['Hill'], params['E0'], params['Emax'], params['EC50'])
    return float(result[0])


def iso_effective_dose(E: float, params: Dict) -> float:
    """
    Calculate the dose that produces effect E (inverse of 4PL).
    
    D = EC50 * ((E - E0) / (Emax - E))^(1/Hill)
    """
    if params is None:
        return float('inf')
    
    e0 = params['E0']
    emax = params['Emax']
    ec50 = params['EC50']
    hill = params['Hill']
    
    E = np.clip(E, e0 + 0.001, emax - 0.001)
    
    try:
        ratio = (E - e0) / (emax - E)
        if ratio <= 0:
            return float('inf')
        dose = ec50 * np.power(ratio, 1.0 / hill)
        return dose
    except Exception:
        return float('inf')


def create_response_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[float], List[float]]:
    """
    Convert input CSV to dose-response matrix format.
    
    Returns:
    --------
    matrix : 2D numpy array (rows = conc1, cols = conc2)
    conc1_vals : sorted list of unique conc1 values
    conc2_vals : sorted list of unique conc2 values
    """
    conc1_vals = sorted(df['conc1'].unique())
    conc2_vals = sorted(df['conc2'].unique())
    
    matrix = np.full((len(conc1_vals), len(conc2_vals)), np.nan)
    
    for _, row in df.iterrows():
        i = conc1_vals.index(row['conc1'])
        j = conc2_vals.index(row['conc2'])
        matrix[i, j] = row['response']
    
    return matrix, conc1_vals, conc2_vals


def calculate_zip(response_matrix: np.ndarray, conc1_vals: List[float], conc2_vals: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate ZIP synergy scores matching R implementation.
    
    The R ZIP algorithm:
    1. Fit single-drug curves for Drug1 (row 0, cols 1:n) and Drug2 (col 0, rows 1:n)
       - Fix lower bound to control (response_matrix[0,0]) and upper bound to 100
    2. Create updated_single matrix with FITTED single-drug values
    3. For each column j>0: fit curve through rows, fixing lower bound to Drug2 response at that conc
    4. For each row i>0: fit curve through cols, fixing lower bound to Drug1 response at that conc
    5. ZIP_fit = average of column-fitted and row-fitted matrices
    6. ZIP_ref = Bliss independence: updated_single[i,0] + updated_single[0,j] - product/100
    7. ZIP_synergy = ZIP_fit - ZIP_ref
    
    Key: R uses drm/drc package with LL.4 (4PL) where params are: Hill, E0, Emax, EC50
    """
    n_rows, n_cols = response_matrix.shape
    
    zip_fit_matrix = np.zeros_like(response_matrix)
    zip_ref = np.zeros_like(response_matrix)
    
    control = response_matrix[0, 0]
    
    drug1_doses = np.array(conc1_vals[1:])
    drug1_responses = response_matrix[1:, 0]
    drug1_params = fit_4pl_curve(drug1_doses, drug1_responses, fixed_lower=control)
    
    drug2_doses = np.array(conc2_vals[1:])
    drug2_responses = response_matrix[0, 1:]
    drug2_params = fit_4pl_curve(drug2_doses, drug2_responses, fixed_lower=control)
    
    updated_single = np.zeros_like(response_matrix)
    updated_single[0, 0] = 0
    updated_single[1:, 0] = drug1_responses
    updated_single[0, 1:] = drug2_responses
    
    if drug1_params:
        for i, d in enumerate(drug1_doses):
            updated_single[i+1, 0] = predict_4pl(d, drug1_params)
    
    if drug2_params:
        for j, d in enumerate(drug2_doses):
            updated_single[0, j+1] = predict_4pl(d, drug2_params)
    
    updated_col = np.zeros_like(response_matrix)
    updated_col[0, :] = updated_single[0, :]
    updated_col[:, 0] = updated_single[:, 0]
    
    for j in range(1, n_cols):
        col_responses = response_matrix[1:, j]
        col_doses = drug1_doses
        
        lower_bound = updated_single[0, j]
        params = fit_4pl_curve(col_doses, col_responses, fixed_lower=lower_bound)
        
        if params:
            for i, d in enumerate(col_doses):
                updated_col[i+1, j] = predict_4pl(d, params)
        else:
            updated_col[1:, j] = col_responses
    
    updated_row = np.zeros_like(response_matrix)
    updated_row[0, :] = updated_single[0, :]
    updated_row[:, 0] = updated_single[:, 0]
    
    for i in range(1, n_rows):
        row_responses = response_matrix[i, 1:]
        row_doses = drug2_doses
        
        lower_bound = updated_single[i, 0]
        params = fit_4pl_curve(row_doses, row_responses, fixed_lower=lower_bound)
        
        if params:
            for j, d in enumerate(row_doses):
                updated_row[i, j+1] = predict_4pl(d, params)
        else:
            updated_row[i, 1:] = row_responses
    
    zip_fit_matrix = (updated_col + updated_row) / 2.0
    zip_fit_matrix = np.clip(zip_fit_matrix, 0, 100)
    zip_fit_matrix[0, 0] = 0
    
    for i in range(n_rows):
        zip_fit_matrix[i, 0] = response_matrix[i, 0]
    for j in range(n_cols):
        zip_fit_matrix[0, j] = response_matrix[0, j]
    
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            e1 = updated_single[i, 0]
            e2 = updated_single[0, j]
            zip_ref[i, j] = e1 + e2 - (e1 * e2) / 100.0
    
    zip_synergy = zip_fit_matrix - zip_ref
    
    return zip_fit_matrix, zip_ref, zip_synergy, updated_single


def calculate_loewe(response_matrix: np.ndarray, conc1_vals: List[float], 
                    conc2_vals: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate Loewe synergy scores matching R implementation.
    
    Loewe equation: d1/D1(E) + d2/D2(E) = 1
    where D1(E) and D2(E) are doses that produce effect E if used alone.
    
    This is solved numerically for each combination to find the expected
    additive effect (Loewe_ref). The synergy is obs - Loewe_ref.
    """
    n_rows, n_cols = response_matrix.shape
    
    loewe_ref = np.zeros_like(response_matrix)
    loewe_synergy = np.zeros_like(response_matrix)
    loewe_ci = np.full_like(response_matrix, np.nan)
    
    drug1_doses = np.array(conc1_vals[1:])
    drug1_responses = response_matrix[1:, 0]
    drug1_params = fit_4pl_curve(drug1_doses, drug1_responses, fixed_lower=response_matrix[0, 0])
    
    drug2_doses = np.array(conc2_vals[1:])
    drug2_responses = response_matrix[0, 1:]
    drug2_params = fit_4pl_curve(drug2_doses, drug2_responses, fixed_lower=response_matrix[0, 0])
    
    if drug1_params is None or drug2_params is None:
        return loewe_ref, loewe_synergy, loewe_ci
    
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            d1 = conc2_vals[j]  # Drug2 (column) concentration
            d2 = conc1_vals[i]  # Drug1 (row) concentration
            obs_response = response_matrix[i, j]
            
            def loewe_equation(E):
                if E <= drug1_params['E0'] or E >= drug1_params['Emax']:
                    return 1e10
                if E <= drug2_params['E0'] or E >= drug2_params['Emax']:
                    return 1e10
                
                D1 = iso_effective_dose(E, drug1_params)
                D2 = iso_effective_dose(E, drug2_params)
                
                if D1 <= 0 or D2 <= 0 or np.isinf(D1) or np.isinf(D2):
                    return 1e10
                
                return d1 / D1 + d2 / D2 - 1.0
            
            E_guess = min(drug1_params['Emax'], drug2_params['Emax']) - 5
            
            try:
                E_solution = fsolve(loewe_equation, E_guess, full_output=True)
                E_sol = float(E_solution[0])
                
                if E_sol < 0 or E_sol > 100 or np.isnan(E_sol):
                    E_sol = max(drug1_params['E0'], drug2_params['E0']) + 10
                
                D1_sol = iso_effective_dose(E_sol, drug1_params)
                D2_sol = iso_effective_dose(E_sol, drug2_params)
                
                if D1_sol > 0 and D2_sol > 0:
                    ci = d1 / D1_sol + d2 / D2_sol
                else:
                    ci = np.nan
                
            except Exception:
                E_sol = 0
                ci = np.nan
            
            loewe_ref[i, j] = E_sol
            loewe_synergy[i, j] = obs_response - E_sol
            loewe_ci[i, j] = ci
    
    return loewe_ref, loewe_synergy, loewe_ci


def calculate_bliss(response_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate Bliss synergy scores.
    
    Bliss_ref = E1 + E2 - E1*E2/100 (Bliss independence)
    Bliss_synergy = observed - Bliss_ref
    """
    n_rows, n_cols = response_matrix.shape
    
    bliss_ref = np.zeros_like(response_matrix)
    bliss_synergy = np.zeros_like(response_matrix)
    
    drug1_responses = response_matrix[:, 0]
    drug2_responses = response_matrix[0, :]
    
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            e1 = drug1_responses[i]
            e2 = drug2_responses[j]
            bliss_ref[i, j] = e1 + e2 - (e1 * e2) / 100.0
            bliss_synergy[i, j] = response_matrix[i, j] - bliss_ref[i, j]
    
    return bliss_ref, bliss_synergy


def calculate_hsa(response_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Calculate HSA (Highest Single Agent) synergy scores.
    
    HSA_ref = max(E1, E2)
    HSA_synergy = observed - HSA_ref
    """
    n_rows, n_cols = response_matrix.shape
    
    hsa_ref = np.zeros_like(response_matrix)
    hsa_synergy = np.zeros_like(response_matrix)
    
    drug1_responses = response_matrix[:, 0]
    drug2_responses = response_matrix[0, :]
    
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            e1 = drug1_responses[i]
            e2 = drug2_responses[j]
            hsa_ref[i, j] = max(e1, e2)
            hsa_synergy[i, j] = response_matrix[i, j] - hsa_ref[i, j]
    
    return hsa_ref, hsa_synergy


def generate_output_order(conc1_vals: List[float], conc2_vals: List[float]) -> List[Tuple[float, float]]:
    """
    Generate output order matching SynergyFinder format.
    
    Order:
    1. Control (0, 0)
    2. Drug1 alone (all conc1 > 0 with conc2=0)
    3. For each conc2 > 0:
       - Drug2 alone (0, conc2)
       - All combinations (conc1 > 0, conc2)
    """
    grid = []
    
    grid.append((0.0, 0.0))
    
    for c1 in conc1_vals:
        if c1 > 0:
            grid.append((float(c1), 0.0))
    
    for c2 in conc2_vals:
        if c2 > 0:
            grid.append((0.0, float(c2)))
            for c1 in conc1_vals:
                if c1 > 0:
                    grid.append((float(c1), float(c2)))
    
    return grid


def calculate_synergy(input_file: str, output_file: str):
    """
    Main function to calculate drug synergy from input file.
    
    Parameters:
    -----------
    input_file : str
        Path to input CSV file
    output_file : str
        Path to output CSV file
    """
    print(f"Loading input data from: {input_file}")
    df = pd.read_csv(input_file)
    print(f"  Total rows: {len(df)}")
    print(f"  Columns: {list(df.columns)}")
    
    block_id = int(df['block_id'].iloc[0])
    
    response_matrix, conc1_vals, conc2_vals = create_response_matrix(df)
    
    input_lookup = {}
    for _, row in df.iterrows():
        key = (float(row['conc1']), float(row['conc2']))
        input_lookup[key] = float(row['response'])
    
    print(f"\nResponse matrix shape: {response_matrix.shape}")
    print(f"  Conc1 values: {conc1_vals}")
    print(f"  Conc2 values: {conc2_vals}")
    
    print(f"\nCalculating ZIP synergy...")
    zip_fit, zip_ref, zip_synergy, updated_single = calculate_zip(
        response_matrix, conc1_vals, conc2_vals)
    
    print(f"Calculating Loewe synergy...")
    loewe_ref, loewe_synergy, loewe_ci = calculate_loewe(
        response_matrix, conc1_vals, conc2_vals)
    
    print(f"Calculating Bliss synergy...")
    bliss_ref, bliss_synergy = calculate_bliss(response_matrix)
    
    print(f"Calculating HSA synergy...")
    hsa_ref, hsa_synergy = calculate_hsa(response_matrix)
    
    output_order = generate_output_order(conc1_vals, conc2_vals)
    
    conc1_to_idx = {v: i for i, v in enumerate(conc1_vals)}
    conc2_to_idx = {v: i for i, v in enumerate(conc2_vals)}
    
    results = []
    for (c1, c2) in output_order:
        i = conc1_to_idx[c1]
        j = conc2_to_idx[c2]
        
        is_single = (c1 == 0 and c2 == 0) or (c1 > 0 and c2 == 0) or (c1 == 0 and c2 > 0)
        obs_response = input_lookup.get((c1, c2), response_matrix[i, j])
        
        if is_single:
            row = {
                'block_id': block_id,
                'conc1': c1,
                'conc2': c2,
                'ZIP_fit': round(obs_response, 3),
                'ZIP_ref': round(obs_response, 3),
                'ZIP_synergy': 0,
                'HSA_ref': round(obs_response, 3),
                'HSA_synergy': 0,
                'Loewe_ref': round(obs_response, 3),
                'Loewe_synergy': 0,
                'Loewe_ci': '',
                'Bliss_ref': round(obs_response, 3),
                'Bliss_synergy': 0
            }
        else:
            ci_val = round(loewe_ci[i, j], 3) if not np.isnan(loewe_ci[i, j]) else ''
            row = {
                'block_id': block_id,
                'conc1': c1,
                'conc2': c2,
                'ZIP_fit': round(zip_fit[i, j], 3),
                'ZIP_ref': round(zip_ref[i, j], 3),
                'ZIP_synergy': round(zip_synergy[i, j], 3),
                'HSA_ref': round(hsa_ref[i, j], 3),
                'HSA_synergy': round(hsa_synergy[i, j], 3),
                'Loewe_ref': round(loewe_ref[i, j], 3),
                'Loewe_synergy': round(loewe_synergy[i, j], 3),
                'Loewe_ci': ci_val,
                'Bliss_ref': round(bliss_ref[i, j], 3),
                'Bliss_synergy': round(bliss_synergy[i, j], 3)
            }
        results.append(row)
    
    cols = ['block_id', 'conc1', 'conc2', 'ZIP_fit', 'ZIP_ref', 'ZIP_synergy',
            'HSA_ref', 'HSA_synergy', 'Loewe_ref', 'Loewe_synergy', 'Loewe_ci',
            'Bliss_ref', 'Bliss_synergy']
    output_df = pd.DataFrame(results)[cols]
    
    output_df.to_csv(output_file, index=False)
    print(f"\nSaved synergy results to: {output_file}")
    
    combo_results = [r for r in results if r['conc1'] > 0 and r['conc2'] > 0]
    
    print(f"\n{'='*60}")
    print("Synergy Summary:")
    print(f"{'='*60}")
    
    for model in ['ZIP', 'Bliss', 'HSA', 'Loewe']:
        syn_key = f'{model}_synergy'
        scores = [r[syn_key] for r in combo_results]
        mean_syn = np.mean(scores)
        min_syn = np.min(scores)
        max_syn = np.max(scores)
        print(f"  {model}: Mean={mean_syn:.2f}, Range=[{min_syn:.2f}, {max_syn:.2f}]")
    
    print(f"{'='*60}")
    
    return output_df


def main():
    """Main entry point."""
    import argparse
    
    parser = argparse.ArgumentParser(description='Calculate drug synergy scores')
    parser.add_argument('input', nargs='?', 
                        default='/Users/kwan/Documents/dengue/DENV_AI/SynergyFinder_input_block6.csv',
                        help='Input CSV file path')
    parser.add_argument('output', nargs='?',
                        default='/Users/kwan/Documents/dengue/DENV_AI/synergy_score_table.csv',
                        help='Output CSV file path')
    
    args = parser.parse_args()
    
    calculate_synergy(args.input, args.output)


if __name__ == '__main__':
    main()
