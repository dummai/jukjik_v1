#!/usr/bin/env python3
"""
Generate synergy heatmap from summary and block CSV files.

Usage:
    python heatmap_generator.py \n        --summary /path/to/synergy_summary.csv \n        --blocks /path/to/SynergyFinder_BlockID.csv \n        --output /desired/path/heatmap.jpeg
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os


def load_files(summary_path: str, block_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Read the two CSV files and return their DataFrames."""
    df_sum = pd.read_csv(summary_path)
    df_block = pd.read_csv(block_path)
    return df_sum, df_block


def build_matrix(df_sum: pd.DataFrame, df_block: pd.DataFrame) -> pd.DataFrame:
    """Create the heatmap matrix.

    Parameters
    ----------
    df_sum : DataFrame
        Contains ``block_id``, ``model`` and ``max_synergy``.
    df_block : DataFrame
        Contains ``block_id``, ``drug1``, ``drug2`` and ``experiment``.

    Returns
    -------
    DataFrame
        Pivoted table: rows = models, columns = drug1_drug2_experiment.
    """
    # Create composite label
    df_block = df_block.copy()
    df_block["label"] = (
        df_block["drug1"] + "_" + df_block["drug2"] + "_" + df_block["experiment"]
    )

    # Merge with summary on block_id
    df = df_sum.merge(df_block[["block_id", "label"]], on="block_id", how="left")

    # Warn if any label missing
    if df["label"].isna().any():
        missing = df[df["label"].isna()]["block_id"].unique()
        print(f"Warning: missing block labels for IDs: {missing}")

    # Pivot to matrix
    matrix = df.pivot(index="model", columns="label", values="max_synergy")
    return matrix


def plot_and_save(matrix: pd.DataFrame, out_path: str, dpi: int = 300) -> None:
    """Render a heatmap and save it as a publication‑quality JPEG."""
    plt.figure(figsize=(max(8, len(matrix.columns) * 0.4), max(6, len(matrix.index) * 0.4)))
    # Diverging green/red palette
    cmap = sns.diverging_palette(150, 10, as_cmap=True)
    ax = sns.heatmap(
        matrix,
        cmap=cmap,
        center=0,
        linewidths=0.5,
        cbar_kws={"label": "max_synergy"},
        square=False,
        annot=True,
        fmt=".2f",
        annot_kws={"size": 7},
    )
    ax.set_xlabel("Drug1_Drug2_Experiment")
    ax.set_ylabel("Model")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    # Ensure parent dir exists
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    plt.savefig(out_path, dpi=dpi, format="jpeg", bbox_inches="tight")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Generate synergy heatmap")
    parser.add_argument("--summary", default="../synergy_summary.csv", help="Path to synergy_summary.csv")
    parser.add_argument("--blocks", default="../SynergyFinder_BlockID.csv", help="Path to SynergyFinder_BlockID.csv")
    parser.add_argument("--output", default="../heatmap.jpeg", help="Output JPEG file path")
    parser.add_argument("--dpi", type=int, default=300, help="Output image resolution")
    args = parser.parse_args()

    df_sum, df_block = load_files(args.summary, args.blocks)
    matrix = build_matrix(df_sum, df_block)
    plot_and_save(matrix, args.output, dpi=args.dpi)
    print(f"Heatmap written to {args.output}")


if __name__ == "__main__":
    main()
