#!/usr/bin/env python3
"""
Drug Synergy Analysis - Matching SynergyFinder Plus Output Format

Calculates drug synergy scores using ZIP, Bliss, Loewe, and HSA models.
Implements the same algorithms as R's SynergyFinder package using drc::LL.4.

Input: CSV with columns: block_id, drug1, drug2, conc1, conc2, response, conc_unit
Output: CSV with columns: block_id, conc1, conc2, ZIP_fit, ZIP_ref, ZIP_synergy,
        HSA_ref, HSA_synergy, Loewe_ref, Loewe_synergy, Loewe_ci, Bliss_ref, Bliss_synergy

References:
- Yadav B, Wennerberg K, Aittokallio T, Tang J. Searching for Drug Synergy in 
  Complex Dose-Response Landscape Using an Interaction Potency Model.
  Computational and Structural Biotechnology Journal 2015; 13: 504-513.
- Ritz C, Baty F, Streibig JC, Gerhard D. Dose-Response Analysis Using R. PLoS ONE 2015.
"""

import numpy as np
import pandas as pd
from scipy.optimize import fsolve, curve_fit, minimize
from typing import Dict, List, Tuple, Optional
import warnings
warnings.filterwarnings('ignore')


def ll4_model(dose: np.ndarray, b: float, c: float, d: float, e: float) -> np.ndarray:
    """
    4-parameter logistic model (LL.4) from R's drc package.
    
    f(x) = c + (d - c) / (1 + exp(b * (log(x) - log(e))))
    
    Parameters:
    -----------
    b : slope parameter (negative for increasing dose-response)
    c : lower limit (minimum response)
    d : upper limit (maximum response)  
    e : EC50 (dose producing 50% of maximum effect)
    
    Note: In R's drc, parameter order is (b, c, d, e)
    """
    dose = np.atleast_1d(dose).astype(float)
    result = np.zeros_like(dose)
    
    zero_mask = (dose == 0) | np.isnan(dose)
    nonzero_mask = ~zero_mask
    
    if np.any(zero_mask):
        result[zero_mask] = c
    
    if np.any(nonzero_mask):
        log_dose = np.log(dose[nonzero_mask])
        result[nonzero_mask] = c + (d - c) / (1.0 + np.exp(b * (log_dose - np.log(e))))
    
    return result


def fit_ll4_curve(doses: np.ndarray, responses: np.ndarray,
                  fixed_lower: float = None, fixed_upper: float = 100.0) -> Optional[Dict]:
    """
    Fit LL.4 curve matching R's drc::drm with fct = LL.4(fixed = c(NA, c, d, NA)).
    
    Parameters:
    -----------
    doses : array of concentrations (>0)
    responses : array of responses
    fixed_lower : if provided, fix parameter c (lower limit) to this value
    fixed_upper : if provided, fix parameter d (upper limit) to this value (default 100)
    
    Returns:
    --------
    dict with parameters: b (slope), c (lower), d (upper), e (EC50)
    """
    valid_mask = ~(np.isnan(doses) | np.isnan(responses)) & (doses > 0)
    doses = doses[valid_mask]
    responses = responses[valid_mask]
    
    if len(doses) < 2:
        return None
    
    responses = np.clip(responses, 0, 150)
    
    if fixed_lower is not None:
        c_fixed = fixed_lower
    else:
        c_fixed = np.min(responses) if len(responses) > 0 else 0
    
    d_fixed = fixed_upper
    
    e_init = np.median(doses)
    b_init = -2.0
    
    def fit_func(dose, b, e):
        return ll4_model(dose, b, c_fixed, d_fixed, e)
    
    try:
        popt, _ = curve_fit(
            fit_func, doses, responses,
            p0=[b_init, e_init],
            bounds=([-20.0, min(doses)/100], [0.01, max(doses)*100]),
            maxfev=10000,
            method='trf'
        )
        
        return {
            'b': popt[0],
            'c': c_fixed,
            'd': d_fixed,
            'e': popt[1]
        }
    except Exception:
        try:
            popt, _ = curve_fit(
                fit_func, doses, responses,
                p0=[-5.0, np.mean(doses)],
                maxfev=10000
            )
            return {
                'b': popt[0],
                'c': c_fixed,
                'd': d_fixed,
                'e': popt[1]
            }
        except Exception:
            return None


def predict_ll4(dose: float, params: Dict) -> float:
    """Predict response at given dose using LL.4 parameters."""
    if params is None:
        return 0.0
    if dose <= 0:
        return params['c']
    result = ll4_model(np.array([dose]), params['b'], params['c'], params['d'], params['e'])
    return float(np.clip(result[0], 0, 100))


def predict_hill(dose: float, params: Dict) -> float:
    """
    Predict response using standard Hill equation.
    
    y = E0 + (Emax - E0) * (dose/EC50)^hill / (1 + (dose/EC50)^hill)
    
    This is used for Loewe calculations in R's SynergyFinder.
    """
    if params is None:
        return 0.0
    
    e0 = params['c']
    emax = params['d']
    ec50 = params['e']
    hill = -params['b']  # R's LL.4 b is negative, Hill equation uses positive
    
    if dose <= 0:
        return e0
    
    ratio = (dose / ec50) ** hill
    result = e0 + (emax - e0) * ratio / (1.0 + ratio)
    return float(np.clip(result, 0, 100))


def iso_effective_dose_ll4(E: float, params: Dict) -> float:
    """
    Calculate the dose that produces effect E using LL.4 inverse.
    
    From LL.4: E = c + (d-c)/(1 + exp(b*(log(D) - log(e))))
    Solving for D: D = e * exp((log((d-c)/(E-c) - 1)) / b)
    """
    if params is None:
        return float('inf')
    
    c = params['c']
    d = params['d']
    e = params['e']
    b = params['b']
    
    if E <= c:
        return float('inf')
    if E >= d:
        return 0.0
    
    try:
        ratio = (d - c) / (E - c) - 1.0
        if ratio <= 0:
            return float('inf')
        dose = e * np.power(ratio, 1.0 / b)
        return max(0.0, dose)
    except Exception:
        return float('inf')


def iso_effective_dose_hill(E: float, params: Dict) -> float:
    """
    Calculate the dose that produces effect E using Hill equation inverse.
    
    From Hill: E = E0 + (Emax-E0) * (D/EC50)^h / (1 + (D/EC50)^h)
    Solving for D: D = EC50 * ((E - E0)/(Emax - E))^(1/hill)
    """
    if params is None:
        return float('inf')
    
    e0 = params['c']
    emax = params['d']
    ec50 = params['e']
    hill = -params['b']  # Positive slope for Hill equation
    
    if E <= e0:
        return float('inf')
    if E >= emax:
        return 0.0
    
    try:
        ratio = (E - e0) / (emax - E)
        if ratio <= 0:
            return float('inf')
        dose = ec50 * np.power(ratio, 1.0 / hill)
        return max(0.0, dose)
    except Exception:
        return float('inf')


def create_response_matrix(df: pd.DataFrame) -> Tuple[np.ndarray, List[float], List[float]]:
    """Convert input CSV to dose-response matrix format."""
    conc1_vals = sorted(df['conc1'].unique())
    conc2_vals = sorted(df['conc2'].unique())
    
    matrix = np.full((len(conc1_vals), len(conc2_vals)), np.nan)
    
    for _, row in df.iterrows():
        i = conc1_vals.index(row['conc1'])
        j = conc2_vals.index(row['conc2'])
        matrix[i, j] = row['response']
    
    return matrix, conc1_vals, conc2_vals


def calculate_zip(response_matrix: np.ndarray, conc1_vals: List[float], 
                  conc2_vals: List[float]) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Calculate ZIP synergy scores matching R's SynergyFinder implementation.
    
    Algorithm from R's ZIP() function:
    1. Fit single-drug curves with drm(..., fct = LL.4(fixed = c(NA, 0, 100, NA)))
    2. Create updated_single matrix with fitted single-drug values
    3. For each column (fixed conc2), fit LL.4 through rows with lower bound = Drug2 response
    4. For each row (fixed conc1), fit LL.4 through columns with lower bound = Drug1 response
    5. ZIP_fit = average of column-fitted and row-fitted matrices
    6. ZIP_ref = Bliss independence: updated_single[i,0] + updated_single[0,j] - product/100
    7. ZIP_synergy = ZIP_fit - ZIP_ref
    """
    n_rows, n_cols = response_matrix.shape
    
    zip_fit_matrix = np.zeros_like(response_matrix)
    zip_ref = np.zeros_like(response_matrix)
    
    control = response_matrix[0, 0]
    
    drug1_doses = np.array(conc1_vals[1:])
    drug1_responses = response_matrix[1:, 0]
    drug1_params = fit_ll4_curve(drug1_doses, drug1_responses, fixed_lower=control)
    
    drug2_doses = np.array(conc2_vals[1:])
    drug2_responses = response_matrix[0, 1:]
    drug2_params = fit_ll4_curve(drug2_doses, drug2_responses, fixed_lower=control)
    
    updated_single = np.zeros_like(response_matrix)
    updated_single[0, 0] = 0
    updated_single[1:, 0] = drug1_responses
    updated_single[0, 1:] = drug2_responses
    
    if drug1_params:
        for i, d in enumerate(drug1_doses):
            updated_single[i+1, 0] = predict_ll4(d, drug1_params)
    
    if drug2_params:
        for j, d in enumerate(drug2_doses):
            updated_single[0, j+1] = predict_ll4(d, drug2_params)
    
    updated_col = np.zeros_like(response_matrix)
    updated_col[0, :] = updated_single[0, :]
    updated_col[:, 0] = updated_single[:, 0]
    
    for j in range(1, n_cols):
        col_responses = response_matrix[1:, j]
        lower_bound = updated_single[0, j]
        params = fit_ll4_curve(drug1_doses, col_responses, fixed_lower=lower_bound)
        
        if params:
            for i, d in enumerate(drug1_doses):
                updated_col[i+1, j] = predict_ll4(d, params)
        else:
            updated_col[1:, j] = col_responses
    
    updated_row = np.zeros_like(response_matrix)
    updated_row[0, :] = updated_single[0, :]
    updated_row[:, 0] = updated_single[:, 0]
    
    for i in range(1, n_rows):
        row_responses = response_matrix[i, 1:]
        lower_bound = updated_single[i, 0]
        params = fit_ll4_curve(drug2_doses, row_responses, fixed_lower=lower_bound)
        
        if params:
            for j, d in enumerate(drug2_doses):
                updated_row[i, j+1] = predict_ll4(d, params)
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
    Calculate Loewe synergy scores matching R's SynergyFinder implementation.
    
    Loewe equation: d1/D1(E) + d2/D2(E) = 1
    where D1(E) and D2(E) are doses that produce effect E if used alone.
    
    Uses standard Hill equation for inverse calculations (not LL.4).
    """
    n_rows, n_cols = response_matrix.shape
    
    loewe_ref = np.zeros_like(response_matrix)
    loewe_synergy = np.zeros_like(response_matrix)
    loewe_ci = np.full_like(response_matrix, np.nan)
    
    control = response_matrix[0, 0]
    
    drug1_doses = np.array(conc1_vals[1:])
    drug1_responses = response_matrix[1:, 0]
    drug1_params = fit_ll4_curve(drug1_doses, drug1_responses, fixed_lower=control)
    
    drug2_doses = np.array(conc2_vals[1:])
    drug2_responses = response_matrix[0, 1:]
    drug2_params = fit_ll4_curve(drug2_doses, drug2_responses, fixed_lower=control)
    
    if drug1_params is None or drug2_params is None:
        for i in range(1, n_rows):
            for j in range(1, n_cols):
                max_single = max(response_matrix[i, 0], response_matrix[0, j])
                loewe_ref[i, j] = max_single
                loewe_synergy[i, j] = response_matrix[i, j] - max_single
        return loewe_ref, loewe_synergy, loewe_ci
    
    for i in range(1, n_rows):
        for j in range(1, n_cols):
            d1 = conc2_vals[j]
            d2 = conc1_vals[i]
            obs_response = response_matrix[i, j]
            
            def ci_at_E_Hill(E):
                """Calculate CI using Hill equation."""
                if E <= drug1_params['c'] or E >= drug1_params['d']:
                    return float('inf')
                if E <= drug2_params['c'] or E >= drug2_params['d']:
                    return float('inf')
                
                D1 = iso_effective_dose_hill(E, drug1_params)
                D2 = iso_effective_dose_hill(E, drug2_params)
                
                if D1 <= 0 or D2 <= 0 or np.isinf(D1) or np.isinf(D2):
                    return float('inf')
                
                return d1 / D1 + d2 / D2
            
            def loewe_equation(E):
                return ci_at_E_Hill(E) - 1.0
            
            E_guess = 50.0
            try:
                E_solution = fsolve(loewe_equation, E_guess, full_output=True)
                E_sol = float(E_solution[0][0])
                ci = ci_at_E_Hill(E_sol)
                
                if E_sol < 0 or E_sol > 100 or np.isnan(E_sol) or np.isinf(ci):
                    pred_d1 = predict_hill(d1 + d2, drug1_params)
                    pred_d2 = predict_hill(d1 + d2, drug2_params)
                    E_sol = max(pred_d1, pred_d2)
                    ci = np.nan
            
            except Exception:
                pred_d1 = predict_hill(d1 + d2, drug1_params)
                pred_d2 = predict_hill(d1 + d2, drug2_params)
                E_sol = max(pred_d1, pred_d2)
                ci = np.nan
            
            loewe_ref[i, j] = E_sol
            loewe_synergy[i, j] = obs_response - E_sol
            loewe_ci[i, j] = ci
    
    return loewe_ref, loewe_synergy, loewe_ci


def calculate_bliss(response_matrix: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Calculate Bliss synergy scores."""
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
    """Calculate HSA (Highest Single Agent) synergy scores."""
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
    """Generate output order matching SynergyFinder format."""
    grid = [(0.0, 0.0)]
    
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


def generate_synergy_summary(block_id: int, results: list, loewe_ci: np.ndarray,
                             conc1_vals: List[float], conc2_vals: List[float],
                             output_file: str):
    """
    Generate summary CSV with max/min synergy information for each model.
    
    Output columns:
    - block_id: Block identifier
    - model: Synergy model name
    - max_synergy: Maximum synergy score
    - max_conc1: Drug 1 concentration at max synergy
    - max_conc2: Drug 2 concentration at max synergy
    - min_synergy: Minimum synergy score
    - min_conc1: Drug 1 concentration at min synergy
    - min_conc2: Drug 2 concentration at min synergy
    - mean_synergy: Mean synergy across all combinations
    - synergy_category: "Synergy", "Antagonism", or "Additive"
    - ci_at_max_synergy: CI value at max synergy (Loewe only)
    """
    combo_results = [r for r in results if r['conc1'] > 0 and r['conc2'] > 0]
    
    if not combo_results:
        return
    
    conc1_to_idx = {v: i for i, v in enumerate(conc1_vals)}
    conc2_to_idx = {v: i for i, v in enumerate(conc2_vals)}
    
    summary_rows = []
    
    for model in ['ZIP', 'Bliss', 'HSA', 'Loewe']:
        syn_key = f'{model}_synergy'
        scores = [r[syn_key] for r in combo_results]
        
        mean_syn = np.mean(scores)
        max_syn = np.max(scores)
        min_syn = np.min(scores)
        
        max_idx = scores.index(max_syn)
        min_idx = scores.index(min_syn)
        
        max_conc1 = combo_results[max_idx]['conc1']
        max_conc2 = combo_results[max_idx]['conc2']
        min_conc1 = combo_results[min_idx]['conc1']
        min_conc2 = combo_results[min_idx]['conc2']
        
        if mean_syn > 5:
            category = "Synergy"
        elif mean_syn < -5:
            category = "Antagonism"
        else:
            category = "Additive"
        
        ci_at_max = ''
        if model == 'Loewe':
            i = conc1_to_idx.get(max_conc1)
            j = conc2_to_idx.get(max_conc2)
            if i is not None and j is not None:
                ci_val = loewe_ci[i, j]
                if not np.isnan(ci_val):
                    ci_at_max = round(ci_val, 3)
        
        summary_rows.append({
            'block_id': block_id,
            'model': model,
            'max_synergy': round(max_syn, 3),
            'max_conc1': max_conc1,
            'max_conc2': max_conc2,
            'min_synergy': round(min_syn, 3),
            'min_conc1': min_conc1,
            'min_conc2': min_conc2,
            'mean_synergy': round(mean_syn, 3),
            'synergy_category': category,
            'ci_at_max_synergy': ci_at_max
        })
    
    summary_df = pd.DataFrame(summary_rows)
    
    cols = ['block_id', 'model', 'max_synergy', 'max_conc1', 'max_conc2',
            'min_synergy', 'min_conc1', 'min_conc2', 'mean_synergy',
            'synergy_category', 'ci_at_max_synergy']
    summary_df = summary_df[cols]
    
    summary_df.to_csv(output_file, index=False)
    print(f"Saved synergy summary to: {output_file}")
    
    return summary_df


def calculate_synergy(input_file: str, output_file: str):
    """Main function to calculate drug synergy from input file."""
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
    
    import os
    summary_file = os.path.join(os.path.dirname(output_file), 'synergy_summary.csv')
    generate_synergy_summary(block_id, results, loewe_ci, conc1_vals, conc2_vals, summary_file)
    
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
