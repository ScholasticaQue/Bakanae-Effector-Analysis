"""
Pathogen Effector Analysis: Genomic Distribution Visualization
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
This script identifies and compares the genomic localization of classical 
Small Cysteine-Rich Proteins (SCRPs) between the Kenyan isolates 5B and 5C. 
It parses GFF3 files and cross-references them with the refined secretome 
to demonstrate conserved expansion across specific scaffolds.
"""
import pandas as pd
from Bio import SeqIO
from pathlib import Path
import matplotlib.pyplot as plt

# --- CONFIGURATION (Relative Paths for Reproducibility) ---
# Ensure these files are placed in the 'data/' or 'results/' directories
EFFECTOR_FASTA = Path("results/refined_classical_effectors.fasta")

GFF_FILES = {
    "5B": Path("data/genomes/Fusarium_fujikuroi_5B.gff3"),
    "5C": Path("data/genomes/Fusarium_fujikuroi_5C.gff3")
}
def get_effector_distribution(isolate, gff_path):
    target_ids = set()
    for record in SeqIO.parse(effector_fasta, "fasta"):
        if record.id.startswith(isolate):
            clean_id = record.id.split('|')[-1]
            target_ids.add(clean_id)
    
    results = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip(): continue
            parts = line.split('\t')
            if parts[2] == 'gene':
                attributes = parts[8]
                if 'ID=' in attributes:
                    gff_id = attributes.split('ID=')[1].split(';')[0]
                    if gff_id in target_ids:
                        results.append(parts[0]) # Scaffold name
    
    return pd.Series(results).value_counts()

def main():
    print("🚀 Comparing 5B/5C and generating distribution plot...")
    
    dist_5B = get_effector_distribution("5B", gff_files["5B"])
    dist_5C = get_effector_distribution("5C", gff_files["5C"])
    
    # Merge data for comparison
    comparison = pd.DataFrame({'5B': dist_5B, '5C': dist_5C}).fillna(0)
    comparison.index.name = 'Scaffold'
    comparison = comparison.sort_values(by='5B', ascending=False)

    # 1. Print the Table to Console
    print("\n" + "="*45)
    print(f"{'Scaffold':<15} | {'5B Count':<10} | {'5C Count':<10}")
    print("-" * 45)
    for index, row in comparison.head(15).iterrows():
        print(f"{index:<15} | {int(row['5B']):<10} | {int(row['5C']):<10}")
    print("="*45)

    # 2. Generate the Plot
    # We will plot the top 12 scaffolds for clarity
    plot_df = comparison.head(12)
    
    ax = plot_df.plot(kind='bar', figsize=(12, 7), color=['#1f77b4', '#ff7f0e'], width=0.8)
    
    plt.title('Conserved Multi-Loci Effector Expansion (5B vs 5C)', fontsize=15, fontweight='bold')
    plt.ylabel('Number of Classical Effectors (SCRPs)', fontsize=12)
    plt.xlabel('Genomic Scaffold', fontsize=12)
    plt.xticks(rotation=45)
    plt.legend(title='Kenyan Isolate', frameon=True)
    plt.grid(axis='y', linestyle='--', alpha=0.6)
    
    # Add value labels on top of bars
    for p in ax.patches:
        ax.annotate(str(int(p.get_height())), (p.get_x() + 0.05, p.get_height() + 0.3), fontsize=9)

    plt.tight_layout()
    output_fig = 'effector_distribution_comparison.png'
    plt.savefig(output_fig, dpi=300)
    print(f"✅ Comparison complete. Figure saved as: {output_fig}")

if __name__ == "__main__":
    main()

