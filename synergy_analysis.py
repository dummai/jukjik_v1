#!/usr/bin/env python3
"""
Drug Synergy Analysis - Matching SynergyFinder Plus Output Format

Calculates drug synergy scores using ZIP, Bliss, Loewe, and HSA models.
Generates output in the exact format of SynergyFinder Plus.

Input: CSV with columns: block_id, drug1, drug2, conc1, conc2, response, conc_unit
Output: CSV with columns: block_id, conc1, conc2, ZIP_fit, ZIP_ref, ZIP_synergy,
        HSA_ref, HSA_synergy, Loewe_ref, Loewe_synergy, Loewe_ci, Bliss_ref, Bliss_synergy

Author: Drug Combination Analysis Pipeline
"""

import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from typing import Dict, List, Tuple, Optional


def hill_equation(d: np.ndarray, E0: float, Emax: float, EC50: float, h: float) -> np.ndarray:
    """
    4-parameter logistic (Hill) equation for dose-response curve.
    
    E(d) = E0 + (Emax - E0) / (1 + (EC50/d)^h)
    
    Parameters:
    -----------
    d : array-like
        Drug concentrations
    E0 : float
        Baseline effect (minimum response)
    Emax : float
        Maximum effect (maximum response)
    EC50 : float
        Concentration producing 50% of maximum effect
    h : float
        Hill slope (steepness)
    
    Returns:
    --------
    array : Response values
    """
    # Handle d=0 case
    d = np.atleast_1d(d)
    result = np.zeros_like(d, dtype=float)
    nonzero_mask = d > 0
    zero_mask = d == 0
    
    if np.any(nonzero_mask):
        dn = d[nonzero_mask]
        result[nonzero_mask] = E0 + (Emax - E0) / (1.0 + np.power(EC50 / dn, h))
    
    if np.any(zero_mask):
        result[zero_mask] = E0
    
    return result


def fit_hill_curve(concentrations: np.ndarray, responses: np.ndarray) -> Optional[Dict]:
    """
    Fit Hill equation to single-drug dose-response data.
    
    Parameters:
    -----------
    concentrations : array
        Drug concentrations
    responses : array
        Observed responses (% inhibition)
    
    Returns:
    --------
    dict : Hill parameters {E0, Emax, EC50, Hill} or None if fitting fails
    """
    # Remove NaN values
    valid_mask = ~(np.isnan(concentrations) | np.isnan(responses))
    concentrations = concentrations[valid_mask]
    responses = responses[valid_mask]
    
    if len(concentrations) < 3:
        return None
    
    # Initial parameter guesses
    E0_init = np.min(responses)
    Emax_init = np.max(responses)
    EC50_init = np.median(concentrations[concentrations > 0])
    h_init = 1.0
    
    p0 = [E0_init, Emax_init, EC50_init, h_init]
    
    # Parameter bounds
    bounds_low = [0.0, 50.0, 0.001, 0.1]
    bounds_high = [50.0, 150.0, 100.0, 20.0]
    
    try:
        popt, _ = curve_fit(
            hill_equation,
            concentrations,
            responses,
            p0=p0,
            bounds=(bounds_low, bounds_high),
            maxfev=10000,
            method='trf'
        )
        
        return {
            'E0': popt[0],
            'Emax': popt[1],
            'EC50': popt[2],
            'Hill': popt[3]
        }
    except Exception as e:
        print(f"    Warning: Hill curve fitting failed: {e}")
        return None


def predict_hill(dose: float, params: Dict) -> float:
    """Predict response using Hill equation parameters."""
    if params is None:
        return 0.0
    return hill_equation(np.array([dose]), params['E0'], params['Emax'], params['EC50'], params['Hill'])[0]


def find_isoeffective_dose(E_target: float, params: Dict, max_dose: float = 100.0) -> float:
    """
    Find dose that gives target effect using Hill equation.
    
    Solve: E_target = E0 + (Emax - E0) / (1 + (EC50/d)^h)
    Rearranged: d = EC50 * ((Emax - E0) / (E_target - E0) - 1)^(1/h)
    
    Returns:
    --------
    float : Dose giving E_target, or max_dose if not achievable
    """
    if params is None:
        return max_dose
    
    E0 = params['E0']
    Emax = params['Emax']
    EC50 = params['EC50']
    h = params['Hill']
    
    # Check if target is achievable
    if E_target <= E0:
        return max_dose  # Would need infinite dose
    if E_target >= Emax:
        return EC50 * 0.01  # Already at max effect
    
    # Calculate dose
    try:
        ratio = (Emax - E0) / (E_target - E0) - 1.0
        if ratio <= 0:
            return max_dose
        dose = EC50 * np.power(ratio, 1.0 / h)
        return min(dose, max_dose)
    except:
        return max_dose


def calculate_loewe_ref(d1: float, d2: float, 
                        hill1: Dict, hill2: Dict,
                        E_comb: float) -> Tuple[float, float]:
    """
    Calculate Loewe reference response and synergy.
    
    Loewe additivity: d1/D1 + d2/D2 = 1
    where D1 and D2 are doses that would give effect E if used alone.
    
    Parameters:
    -----------
    d1, d2 : float
        Drug concentrations in combination
    hill1, hill2 : dict
        Hill parameters for drugs 1 and 2
    E_comb : float
        Observed combination effect
    
    Returns:
    --------
    tuple : (Loewe_ref, Loewe_synergy)
    """
    if d1 == 0 or d2 == 0:
        return E_comb, 0.0
    
    if hill1 is None or hill2 is None:
        return 0.0, E_comb
    
    # Find isoeffective doses
    D1 = find_isoeffective_dose(E_comb, hill1)
    D2 = find_isoeffective_dose(E_comb, hill2)
    
    # Calculate Loewe combination index
    CI = d1 / D1 + d2 / D2
    
    # Loewe synergy
    Loewe_synergy = E_comb * (1.0 - 1.0 / CI) if CI > 0 else 0.0
    
    # Loewe_ref is the effect on the isobologram
    Loewe_ref = E_comb / CI if CI > 0 else 0.0
    
    return Loewe_ref, Loewe_synergy


def fit_2d_surface_zip(conc1_data: np.ndarray, conc2_data: np.ndarray, 
                        response_data: np.ndarray,
                        hill1: Dict, hill2: Dict) -> Dict:
    """
    Fit 2D dose-response surface for ZIP model.
    
    ZIP uses a linear model on transformed dose scale:
    E(d1, d2) = a + b1*log(d1/conc1_EC50) + b2*log(d2/conc2_EC50) 
                + interaction terms
    
    For single drugs, ZIP_fit = observed response.
    For combinations, ZIP_fit is predicted from the surface model.
    
    Returns model parameters.
    """
    # Predict single-drug responses using Hill curves
    E1_pred = np.array([predict_hill(c, hill1) for c in conc1_data])
    E2_pred = np.array([predict_hill(c, hill2) for c in conc2_data])
    
    # ZIP model: assumes multiplicative action on effect scale
    # For each point, predict from 2D surface
    
    # Create feature matrix for linear regression
    # Features: E1_pred, E2_pred, interaction
    # E_model = a*E1 + b*E2 + c*E1*E2
    
    # Filter to valid points
    valid_mask = ~np.isnan(response_data)
    E1_valid = E1_pred[valid_mask]
    E2_valid = E2_pred[valid_mask]
    resp_valid = response_data[valid_mask]
    
    # Build design matrix
    # Features: [1, E1, E2, E1*E2]
    X = np.column_stack([
        np.ones(len(resp_valid)),
        E1_valid,
        E2_valid,
        E1_valid * E2_valid / 100.0  # Interaction term
    ])
    
    # Solve least squares: X @ params = response
    try:
        params, _, _, _ = np.linalg.lstsq(X, resp_valid, rcond=None)
    except:
        # Fallback parameters
        params = [0, 0.5, 0.5, 1.0]
    
    return {
        'params': params,
        'E1_pred': E1_pred,
        'E2_pred': E2_pred,
        'conc1_data': conc1_data,
        'conc2_data': conc2_data,
        'response_data': response_data
    }


def get_zip_fit(conc1: float, conc2: float, 
                surface: Dict,
                input_lookup: Dict,
                hill1: Dict,
                hill2: Dict,
                is_single: bool) -> float:
    """
    Get ZIP_fit value at (conc1, conc2).
    
    For single drugs: return observed response (or Hill-predicted)
    For combinations: use 2D ZIP surface model prediction
    """
    # Single drug or control: return observed value
    key = (conc1, conc2)
    if is_single:
        if key in input_lookup:
            return input_lookup[key]
        # Fallback to Hill curve prediction
        if conc1 > 0:
            return predict_hill(conc1, hill1)
        elif conc2 > 0:
            return predict_hill(conc2, hill2)
        return 0.0
    
    # Combination: use ZIP surface model
    # Predict E1 and E2 at these concentrations
    E1 = predict_hill(conc1, hill1)
    E2 = predict_hill(conc2, hill2)
    
    # Apply ZIP model: E_fit = params[0] + params[1]*E1 + params[2]*E2 + params[3]*E1*E2/100
    params = surface['params']
    zip_fit = params[0] + params[1]*E1 + params[2]*E2 + params[3]*E1*E2/100.0
    
    # Clamp to [0, 100] range
    zip_fit = max(0.0, min(100.0, zip_fit))
    
    return zip_fit


def generate_dose_grid(data: pd.DataFrame) -> List[Tuple[float, float]]:
    """
    Generate complete dose-response grid in SynergyFinder order.
    
    Order: Control → Drug1 alone → Grid (for each conc2, all conc1)
    
    Returns:
    --------
    list : [(conc1, conc2), ...] grid points
    """
    # Extract unique concentrations
    conc1_values = sorted(data['conc1'].unique())
    conc2_values = sorted(data['conc2'].unique())
    
    grid = []
    
    # Control (0, 0)
    grid.append((0.0, 0.0))
    
    # Drug1 alone (all conc1 with conc2=0)
    for c1 in conc1_values:
        if c1 > 0:
            grid.append((float(c1), 0.0))
    
    # Grid: for each conc2 > 0, all conc1 including 0 (Drug2 alone included in grid)
    for c2 in conc2_values:
        if c2 > 0:
            # Add Drug2 alone first (0, c2)
            grid.append((0.0, float(c2)))
            # Then combinations
            for c1 in conc1_values:
                if c1 > 0:
                    grid.append((float(c1), float(c2)))
    
    return grid


def calculate_synergy_row(conc1: float, conc2: float,
                          input_lookup: Dict,
                          hill1: Dict,
                          hill2: Dict,
                          surface: Dict,
                          block_id: int) -> Dict:
    """
    Calculate all synergy metrics for a single (conc1, conc2) point.
    
    Returns:
    --------
    dict : All columns for output CSV row
    """
    # Get observed response (exists in input)
    key = (conc1, conc2)
    response = input_lookup.get(key, None)
    
    # Predict responses from single-drug curves
    E1 = predict_hill(conc1, hill1) if conc1 > 0 else 0.0
    E2 = predict_hill(conc2, hill2) if conc2 > 0 else 0.0
    
    # Determine if single drug or combination
    is_single = (conc1 == 0 and conc2 == 0) or (conc1 > 0 and conc2 == 0) or (conc1 == 0 and conc2 > 0)
    
    # Get ZIP_fit
    ZIP_fit = get_zip_fit(conc1, conc2, surface, input_lookup, hill1, hill2, is_single)
    
    if is_single:
        # Single drug: synergy = 0, reference = response
        ZIP_ref = ZIP_fit
        ZIP_synergy = 0.0
        
        HSA_ref = ZIP_fit
        HSA_synergy = 0.0
        
        Loewe_ref = ZIP_fit
        Loewe_synergy = 0.0
        Loewe_ci = None
        
        Bliss_ref = ZIP_fit
        Bliss_synergy = 0.0
    else:
        # Combination: calculate all synergies
        
        # For Bliss and HSA: use RAW response (not ZIP_fit)
        # For ZIP: use ZIP_fit
        
        # Use response if available, otherwise use ZIP_fit (for interpolated points)
        if response is None:
            response = ZIP_fit
        
        # ZIP_ref: expected from independent action (Bliss-like)
        ZIP_ref = E1 + E2 - (E1 * E2) / 100.0
        ZIP_synergy = ZIP_fit - ZIP_ref
        
        # HSA: uses raw response
        HSA_ref = max(E1, E2)
        HSA_synergy = response - HSA_ref
        
        # Loewe: uses raw response
        Loewe_ref, Loewe_synergy = calculate_loewe_ref(conc1, conc2, hill1, hill2, response)
        Loewe_ci = None  # NA
        
        # Bliss: uses raw response
        Bliss_ref = E1 + E2 - (E1 * E2) / 100.0
        Bliss_synergy = response - Bliss_ref
    
    return {
        'block_id': block_id,
        'conc1': conc1,
        'conc2': conc2,
        'ZIP_fit': ZIP_fit,
        'ZIP_ref': ZIP_ref,
        'ZIP_synergy': ZIP_synergy,
        'HSA_ref': HSA_ref,
        'HSA_synergy': HSA_synergy,
        'Loewe_ref': Loewe_ref,
        'Loewe_synergy': Loewe_synergy,
        'Loewe_ci': Loewe_ci,
        'Bliss_ref': Bliss_ref,
        'Bliss_synergy': Bliss_synergy
    }


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
    data = pd.read_csv(input_file)
    print(f"  Total rows: {len(data)}")
    print(f"  Columns: {list(data.columns)}")
    
    # Get block_id (assume single block per file)
    block_id = int(data['block_id'].iloc[0])
    
    # Create lookup dictionary for input responses
    input_lookup = {}
    for _, row in data.iterrows():
        key = (float(row['conc1']), float(row['conc2']))
        input_lookup[key] = float(row['response'])
    
    # Extract single-drug data
    drug1_alone = data[(data['conc1'] > 0) & (data['conc2'] == 0)]
    drug2_alone = data[(data['conc1'] == 0) & (data['conc2'] > 0)]
    
    print(f"\nFitting Hill curves...")
    print(f"  Drug1 ({drug1_alone['drug1'].iloc[0] if len(drug1_alone) > 0 else 'N/A'})")
    
    # Fit Hill curves
    hill1 = fit_hill_curve(drug1_alone['conc1'].values, drug1_alone['response'].values)
    if hill1:
        print(f"    E0={hill1['E0']:.2f}, Emax={hill1['Emax']:.2f}, EC50={hill1['EC50']:.4f}, Hill={hill1['Hill']:.2f}")
    else:
        print(f"    Warning: Could not fit Hill curve for Drug1")
        hill1 = {'E0': 0, 'Emax': 100, 'EC50': 1.0, 'Hill': 1.0}
    
    print(f"  Drug2 ({drug2_alone['drug2'].iloc[0] if len(drug2_alone) > 0 else 'N/A'})")
    hill2 = fit_hill_curve(drug2_alone['conc2'].values, drug2_alone['response'].values)
    if hill2:
        print(f"    E0={hill2['E0']:.2f}, Emax={hill2['Emax']:.2f}, EC50={hill2['EC50']:.4f}, Hill={hill2['Hill']:.2f}")
    else:
        print(f"    Warning: Could not fit Hill curve for Drug2")
        hill2 = {'E0': 0, 'Emax': 100, 'EC50': 1.0, 'Hill': 1.0}
    
    # Fit 2D surface for ZIP
    print(f"\nFitting 2D dose-response surface for ZIP model...")
    surface = fit_2d_surface_zip(data['conc1'].values, data['conc2'].values, 
                                  data['response'].values, hill1, hill2)
    print(f"  ZIP surface parameters: intercept={surface['params'][0]:.2f}, "
          f"E1_coef={surface['params'][1]:.3f}, E2_coef={surface['params'][2]:.3f}, "
          f"interaction={surface['params'][3]:.3f}")
    
    # Generate complete dose grid
    print(f"\nGenerating dose grid...")
    grid = generate_dose_grid(data)
    print(f"  Grid points: {len(grid)}")
    
    # Calculate synergy for each grid point
    print(f"\nCalculating synergy scores...")
    results = []
    for (conc1, conc2) in grid:
        row = calculate_synergy_row(conc1, conc2, input_lookup, hill1, hill2, surface, block_id)
        results.append(row)
    
    # Create output DataFrame
    output_df = pd.DataFrame(results)
    
    # Round numeric columns
    numeric_cols = ['ZIP_fit', 'ZIP_ref', 'ZIP_synergy', 'HSA_ref', 'HSA_synergy', 
                    'Loewe_ref', 'Loewe_synergy', 'Bliss_ref', 'Bliss_synergy']
    for col in numeric_cols:
        output_df[col] = output_df[col].round(2)
    
    # Save to CSV
    output_df.to_csv(output_file, index=False)
    print(f"\nSaved synergy results to: {output_file}")
    
    # Print summary
    print(f"\n{'='*60}")
    print("Synergy Summary:")
    print(f"{'='*60}")
    
    combo_results = [r for r in results if r['conc1'] > 0 and r['conc2'] > 0]
    
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
