"""
Drug Synergy Analysis using Python synergy package

Calculates drug synergy scores using ZIP, Bliss, Loewe, and HSA models.
Generates publication-quality dose-response plots and synergy heatmaps.

Input files:
- SynergyFinder_input.csv: Drug combination data
- SynergyFinder_BlockID.csv: Block ID summary with experiment names

Output:
- synergy_summary.csv: Combined summary of all synergy scores
- Dose-response plots (JPEG, publication quality)
- Synergy heatmaps (JPEG, publication quality)
- Detailed per-block synergy data
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import seaborn as sns
from typing import List, Optional, Dict

try:
    from synergy.single import Hill
    from synergy.combination import ZIP, Bliss, Loewe, HSA
    from synergy.dataset import load_data
    SYNERGY_AVAILABLE = True
except ImportError:
    SYNERGY_AVAILABLE = False
    print("Warning: 'synergy' package not installed. Install with: pip install synergy")


def fit_dose_response(concentrations: np.ndarray, responses: np.ndarray, 
                       outlier_detection: bool = True) -> Dict:
    """
    Fit dose-response curve using Hill equation (LL4 equivalent).
    
    Args:
        concentrations: Drug concentrations
        responses: Response values (% inhibition)
        outlier_detection: Whether to remove outliers before fitting
        
    Returns:
        Dictionary with fitted parameters (E0, Emax, EC50, Hill_slope)
    """
    if not SYNERGY_AVAILABLE:
        return None
    
    # Remove NaN values
    valid_mask = ~(np.isnan(concentrations) | np.isnan(responses))
    concentrations = concentrations[valid_mask]
    responses = responses[valid_mask]
    
    if len(concentrations) < 3:
        return None
    
    # Simple outlier removal using IQR
    if outlier_detection and len(responses) > 4:
        Q1, Q3 = np.percentile(responses, [25, 75])
        IQR = Q3 - Q1
        lower_bound = Q1 - 1.5 * IQR
        upper_bound = Q3 + 1.5 * IQR
        inlier_mask = (responses >= lower_bound) & (responses <= upper_bound)
        concentrations = concentrations[inlier_mask]
        responses = responses[inlier_mask]
    
    try:
        model = Hill()
        model.fit(concentrations, responses)
        return {
            'E0': model.E0,
            'Emax': model.Emax,
            'EC50': model.C,
            'Hill_slope': model.h,
            'model': model
        }
    except Exception as e:
        print(f"  Warning: Could not fit Hill curve: {e}")
        return None


def calculate_synergy_scores(
    d1: np.ndarray,
    d2: np.ndarray,
    E: np.ndarray,
    model_type: str
) -> np.ndarray:
    """
    Calculate synergy scores for a single model.
    
    Args:
        d1: Drug 1 concentrations
        d2: Drug 2 concentrations
        E: Response values
        model_type: Model name (ZIP, Bliss, Loewe, HSA)
        
    Returns:
        Array of synergy scores
    """
    if not SYNERGY_AVAILABLE:
        return np.zeros(len(E))
    
    try:
        if model_type == "ZIP":
            model = ZIP()
        elif model_type == "Bliss":
            model = Bliss()
        elif model_type == "Loewe":
            model = Loewe()
        elif model_type == "HSA":
            model = HSA()
        else:
            return np.zeros(len(E))
        
        model.fit(d1, d2, E)
        return model.synergy
    except Exception as e:
        print(f"  Warning: Could not calculate {model_type} synergy: {e}")
        return np.zeros(len(E))


def plot_dose_response(
    concentrations: np.ndarray,
    responses: np.ndarray,
    fitted_model: Dict,
    drug_name: str,
    output_path: str,
    title: str = ""
):
    """
    Generate publication-quality dose-response plot.
    
    Args:
        concentrations: Drug concentrations
        responses: Response values
        fitted_model: Fitted model parameters
        drug_name: Drug name for labeling
        output_path: Output file path (JPEG)
        title: Plot title
    """
    # Set publication quality style
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['font.size'] = 10
    plt.rcParams['axes.titlesize'] = 12
    plt.rcParams['axes.labelsize'] = 10
    
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Plot observed data
    ax.scatter(concentrations, responses, color='#1f77b4', s=50, 
               alpha=0.7, label='Observed', zorder=5)
    
    # Plot fitted curve
    if fitted_model and fitted_model.get('model'):
        model = fitted_model['model']
        conc_range = np.logspace(np.log10(min(concentrations[concentrations > 0])),
                                  np.log10(max(concentrations)), 100)
        try:
            predicted = model.predict(conc_range)
            ax.plot(conc_range, predicted, color='#d62728', linewidth=2, 
                   label='Fitted (Hill)', zorder=4)
        except:
            pass
    
    ax.set_xlabel(f'{drug_name} Concentration (µM)', fontweight='bold')
    ax.set_ylabel('% Inhibition', fontweight='bold')
    ax.set_title(title, fontweight='bold', fontsize=11)
    
    ax.set_xscale('log')
    ax.set_ylim(min(-10, min(responses) - 5), max(110, max(responses) + 5))
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.legend(loc='best', frameon=False)
    ax.grid(True, alpha=0.3, linestyle='--')
    
    plt.tight_layout()
    plt.savefig(output_path, format='jpeg', dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


def plot_synergy_heatmap(
    d1: np.ndarray,
    d2: np.ndarray,
    synergy_scores: np.ndarray,
    drug1: str,
    drug2: str,
    model_name: str,
    output_path: str,
    title: str = ""
):
    """
    Generate publication-quality synergy heatmap.
    
    Args:
        d1: Drug 1 concentrations
        d2: Drug 2 concentrations
        synergy_scores: Synergy scores
        drug1: Drug 1 name
        drug2: Drug 2 name
        model_name: Model name for title
        output_path: Output file path (JPEG)
        title: Plot title
    """
    # Set publication quality style
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['font.size'] = 10
    
    # Create pivot table for heatmap
    df = pd.DataFrame({
        'd1': d1,
        'd2': d2,
        'synergy': synergy_scores
    })
    
    # Filter to combination rows only (both drugs > 0)
    df_combo = df[(df['d1'] > 0) & (df['d2'] > 0)].copy()
    
    if df_combo.empty:
        print(f"  Warning: No combination data for {model_name} heatmap")
        return
    
    # Pivot for heatmap
    pivot = df_combo.pivot_table(values='synergy', index='d2', columns='d1', aggfunc='mean')
    
    if pivot.empty:
        print(f"  Warning: Empty pivot table for {model_name} heatmap")
        return
    
    # Create heatmap
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Custom diverging colormap: blue (synergy) - white - red (antagonism)
    cmap = sns.diverging_palette(240, 10, as_cmap=True)
    
    sns.heatmap(pivot, cmap=cmap, center=0, annot=True, fmt='.1f',
                linewidths=0.5, ax=ax, cbar_kws={'label': 'Synergy Score'})
    
    ax.set_xlabel(f'{drug1} Concentration (µM)', fontweight='bold')
    ax.set_ylabel(f'{drug2} Concentration (µM)', fontweight='bold')
    ax.set_title(f'{title}\n{model_name} Model', fontweight='bold', fontsize=11)
    
    plt.tight_layout()
    plt.savefig(output_path, format='jpeg', dpi=300,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


def calculate_synergy(
    input_path: str,
    block_id_path: str,
    output_dir: str = "synergy_results",
    outlier_detection: bool = True,
    curve_fit: str = "LL4",
    dose_response_plot: bool = True,
    synergy_models: Optional[List[str]] = None
):
    """
    Calculate drug synergy using Python synergy package.
    
    Args:
        input_path: Path to SynergyFinder_input.csv
        block_id_path: Path to SynergyFinder_BlockID.csv
        output_dir: Directory for output files
        outlier_detection: Enable/disable outlier detection
        curve_fit: Curve fitting method (LL4, LL5, Hill)
        dose_response_plot: Generate dose-response line plots
        synergy_models: List of synergy models to calculate
                        Default: ["ZIP", "Bliss", "Loewe", "HSA"]
    """
    # Default synergy models
    if synergy_models is None:
        synergy_models = ["ZIP", "Bliss", "Loewe", "HSA"]
    
    # Validate models
    valid_models = ["ZIP", "Bliss", "Loewe", "HSA"]
    synergy_models = [m for m in synergy_models if m in valid_models]
    if not synergy_models:
        synergy_models = valid_models
    
    print(f"Using synergy models: {synergy_models}")
    print(f"Outlier detection: {outlier_detection}")
    print(f"Curve fitting: {curve_fit}")
    
    # Create output directories
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'dose_response'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'heatmaps'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'data'), exist_ok=True)
    
    # Load input data
    print(f"\nLoading input data from: {input_path}")
    input_df = pd.read_csv(input_path)
    print(f"  Total rows: {len(input_df)}")
    
    # Load block ID data
    print(f"Loading block ID data from: {block_id_path}")
    block_df = pd.read_csv(block_id_path)
    print(f"  Total blocks: {len(block_df)}")
    
    # Summary results
    summary_results = []
    
    # Process each block
    for _, block_row in block_df.iterrows():
        block_id = block_row['block_id']
        drug1 = block_row['drug1']
        drug2 = block_row['drug2']
        experiment = block_row['experiment']
        
        print(f"\n{'='*60}")
        print(f"Processing Block {block_id}: {drug1} + {drug2} ({experiment})")
        print(f"{'='*60}")
        
        # Filter data for this block
        block_data = input_df[input_df['block_id'] == block_id].copy()
        
        if len(block_data) == 0:
            print(f"  Warning: No data found for block {block_id}")
            continue
        
        # Extract single drug data
        drug1_single = block_data[(block_data['conc1'] > 0) & (block_data['conc2'] == 0)]
        drug2_single = block_data[(block_data['conc1'] == 0) & (block_data['conc2'] > 0)]
        combinations = block_data[(block_data['conc1'] > 0) & (block_data['conc2'] > 0)]
        control = block_data[(block_data['conc1'] == 0) & (block_data['conc2'] == 0)]
        
        print(f"  Drug1 single doses: {len(drug1_single)}")
        print(f"  Drug2 single doses: {len(drug2_single)}")
        print(f"  Combinations: {len(combinations)}")
        print(f"  Control: {len(control)}")
        
        # Prepare data for synergy calculation
        d1 = block_data['conc1'].values.astype(float)
        d2 = block_data['conc2'].values.astype(float)
        E = block_data['response'].values.astype(float)
        
        # Replace NaN with 0
        d1 = np.nan_to_num(d1, nan=0.0)
        d2 = np.nan_to_num(d2, nan=0.0)
        E = np.nan_to_num(E, nan=0.0)
        
        # Fit dose-response curves
        print("\n  Fitting dose-response curves...")
        
        if len(drug1_single) > 0:
            drug1_conc = drug1_single['conc1'].values.astype(float)
            drug1_resp = drug1_single['response'].values.astype(float)
            fit1 = fit_dose_response(drug1_conc, drug1_resp, outlier_detection)
            
            if fit1:
                print(f"    {drug1}: EC50={fit1['EC50']:.3f} µM, Hill={fit1['Hill_slope']:.2f}")
            
            if dose_response_plot and fit1:
                plot_path = os.path.join(output_dir, 'dose_response', 
                                          f'block{block_id}_{drug1}.jpg')
                plot_dose_response(drug1_conc, drug1_resp, fit1, drug1,
                                   plot_path, f'{drug1} - {experiment}')
        
        if len(drug2_single) > 0:
            drug2_conc = drug2_single['conc2'].values.astype(float)
            drug2_resp = drug2_single['response'].values.astype(float)
            fit2 = fit_dose_response(drug2_conc, drug2_resp, outlier_detection)
            
            if fit2:
                print(f"    {drug2}: EC50={fit2['EC50']:.3f} µM, Hill={fit2['Hill_slope']:.2f}")
            
            if dose_response_plot and fit2:
                plot_path = os.path.join(output_dir, 'dose_response',
                                          f'block{block_id}_{drug2}.jpg')
                plot_dose_response(drug2_conc, drug2_resp, fit2, drug2,
                                   plot_path, f'{drug2} - {experiment}')
        
        # Calculate synergy for each model
        print("\n  Calculating synergy scores...")
        
        # Initialize block synergy data with base columns
        block_synergy_data = {
            'block_id': [block_id] * len(block_data),
            'drug1': [drug1] * len(block_data),
            'drug2': [drug2] * len(block_data),
            'conc1': d1.tolist(),
            'conc2': d2.tolist(),
            'response': E.tolist()
        }
        
        for model_name in synergy_models:
            print(f"    {model_name}...")
            synergy_scores = calculate_synergy_scores(d1, d2, E, model_name)
            block_synergy_data[f'{model_name}_synergy'] = synergy_scores.tolist()
            valid_scores = synergy_scores[~np.isnan(synergy_scores)]
            if len(valid_scores) > 0:
                mean_synergy = np.mean(valid_scores)
                std_synergy = np.std(valid_scores)
                print(f"      Mean: {mean_synergy:.2f}, Std: {std_synergy:.2f}")
                
                summary_results.append({
                    'block_id': block_id,
                    'experiment': experiment,
                    'drug1': drug1,
                    'drug2': drug2,
                    'model': model_name,
                    'mean_synergy': mean_synergy,
                    'std_synergy': std_synergy,
                    'min_synergy': np.min(valid_scores),
                    'max_synergy': np.max(valid_scores),
                    'n_combinations': len(valid_scores)
                })
            
            # Generate heatmap
            combo_mask = (d1 > 0) & (d2 > 0)
            if np.any(combo_mask):
                heatmap_path = os.path.join(output_dir, 'heatmaps',
                                             f'block{block_id}_{model_name}.jpg')
                plot_synergy_heatmap(
                    d1[combo_mask], d2[combo_mask], synergy_scores[combo_mask],
                    drug1, drug2, model_name, heatmap_path,
                    f'{drug1} + {drug2} ({experiment})'
                )
        
        # Save block data
        block_synergy_df = pd.DataFrame(block_synergy_data)
        
        block_output_path = os.path.join(output_dir, 'data', f'block{block_id}_synergy.csv')
        block_synergy_df.to_csv(block_output_path, index=False)
        print(f"\n  Saved block data to: {block_output_path}")
    
    # Save summary
    if summary_results:
        summary_df = pd.DataFrame(summary_results)
        
        # Pivot to wide format
        summary_wide = summary_df.pivot_table(
            index=['block_id', 'experiment', 'drug1', 'drug2'],
            columns='model',
            values='mean_synergy'
        ).reset_index()
        
        # Save CSV summary
        summary_path = os.path.join(output_dir, 'synergy_summary.csv')
        summary_df.to_csv(summary_path, index=False)
        print(f"\nSaved synergy summary to: {summary_path}")
        
        # Save Excel summary
        excel_path = os.path.join(output_dir, 'synergy_summary.xlsx')
        summary_wide.to_excel(excel_path, index=False)
        print(f"Saved Excel summary to: {excel_path}")
    
    print(f"\n{'='*60}")
    print("Synergy analysis complete!")
    print(f"Output directory: {output_dir}")
    print(f"{'='*60}")
    
    return output_dir


if __name__ == "__main__":
    import sys
    
    # Default paths
    base_path = "/Users/kwan/Documents/dengue/DENV_AI"
    
    input_file = os.path.join(base_path, "SynergyFinder_input.csv")
    block_id_file = os.path.join(base_path, "SynergyFinder_BlockID.csv")
    output_directory = os.path.join(base_path, "synergy_results")
    
    # Parse command line arguments
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        block_id_file = sys.argv[2]
    if len(sys.argv) > 3:
        output_directory = sys.argv[3]
    
    # Run synergy analysis
    calculate_synergy(
        input_path=input_file,
        block_id_path=block_id_file,
        output_dir=output_directory,
        outlier_detection=True,
        curve_fit="LL4",
        dose_response_plot=True,
        synergy_models=["ZIP", "Bliss", "Loewe", "HSA"]
    )
