"""
Pathogen Effector Analysis: Presence/Absence Heatmap Generation
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
This script parses CD-HIT clustering results to visualize the conservation 
of effector families across 13 global F. fujikuroi isolates. It categorizes 
clusters into 'Core', 'Accessory', and 'Unique' suites to highlight 
regional effector expansions.
"""
#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

# ==============================
# USER INPUT & PATHS
# ==============================

INPUT_FILE = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/mapped_effector_orthogroups.tsv"
EFFECTOR_FILE = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/refined_effector_ids.txt"
OUTPUT_MATRIX = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/effector_presence_absence_matrix.csv"
HEATMAP_OUT = "/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results/effector_heatmap.png"

ISOLATES = ["5B", "5C", "Baka", "B20", "C1995", "CFF1", "CFF2",
            "E282", "FSU48", "IMI58289", "KSU3368", "M567",
            "MRC2276", "NCIM1100"]

# ==============================
# HELPER FUNCTIONS
# ==============================

def normalize_id(gene_id):
    """
    Cleans gene IDs by removing whitespace, quotes, 
    and prefixes like '5B|'.
    """
    gene_id = str(gene_id).strip().replace('"', '')
    
    # Remove prefix if present (e.g., '5B|5B_005636' -> '5B_005636')
    if "|" in gene_id:
        gene_id = gene_id.split("|")[-1]
    
    return gene_id

# ==============================
# 1. LOAD EFFECTOR LIST
# ==============================

print("--- Step 1: Loading Effector IDs ---")
if not os.path.exists(EFFECTOR_FILE):
    print(f"❌ Error: {EFFECTOR_FILE} not found.")
    exit()

with open(EFFECTOR_FILE) as f:
    # Normalize every ID in the text file
    effector_ids = {normalize_id(line) for line in f if line.strip()}

print(f"Total unique effector IDs loaded: {len(effector_ids)}")

# ==============================
# 2. LOAD ORTHOGROUP DATA
# ==============================

print("\n--- Step 2: Loading Orthogroups ---")
df = pd.read_csv(INPUT_FILE, sep="\t", dtype=str).fillna("")
print(f"Total Orthogroups loaded: {df.shape[0]}")

# ==============================
# 3. BUILD PRESENCE/ABSENCE MATRIX
# ==============================

print("\n--- Step 3: Building Presence/Absence Matrix ---")
matrix = []

for idx, row in df.iterrows():
    og = row["Orthogroup"]
    result = {"Orthogroup": og}

    for isolate in ISOLATES:
        cell_content = row.get(isolate, "")
        
        if cell_content:
            # Split comma-separated genes and normalize them
            genes_in_cell = [normalize_id(g) for g in str(cell_content).split(",") if g.strip()]
            
            # Check if any gene in this orthogroup/isolate combo is in our effector list
            is_present = any(gene in effector_ids for gene in genes_in_cell)
            result[isolate] = 1 if is_present else 0
        else:
            result[isolate] = 0

    matrix.append(result)

matrix_df = pd.DataFrame(matrix)

# Add a helper column for statistics and sorting
matrix_df["Total_Isolates"] = matrix_df[ISOLATES].sum(axis=1)

# Filter out orthogroups that don't contain any effectors at all
effector_matrix = matrix_df[matrix_df["Total_Isolates"] > 0].copy()

# ==============================
# 4. SUMMARY STATISTICS
# ==============================

print("\n--- Step 4: Summary Statistics ---")
core = effector_matrix[effector_matrix["Total_Isolates"] == len(ISOLATES)]
accessory = effector_matrix[(effector_matrix["Total_Isolates"] > 1) & 
                            (effector_matrix["Total_Isolates"] < len(ISOLATES))]
singleton = effector_matrix[effector_matrix["Total_Isolates"] == 1]

print(f"Core Effector OGs (in all {len(ISOLATES)}): {len(core)}")
print(f"Accessory Effector OGs: {len(accessory)}")
print(f"Singleton Effector OGs: {len(singleton)}")

# Save CSV
effector_matrix.to_csv(OUTPUT_MATRIX, index=False)
print(f"Matrix saved to: {OUTPUT_MATRIX}")

# ==============================
# 5. GENERATE HEATMAP
# ==============================

print("\n--- Step 5: Generating Heatmap ---")
if not effector_matrix.empty:
    # Sort by Total_Isolates so Core is at the top, Singletons at the bottom
    heatmap_data = effector_matrix.sort_values(by="Total_Isolates", ascending=False)
    
    # Set index for the plot and drop the helper column
    plot_df = heatmap_data.set_index("Orthogroup").drop(columns=["Total_Isolates"])

    plt.figure(figsize=(10, 12))
    sns.heatmap(
        plot_df, 
        cmap=["#f8f9fa", "#1a5276"],  # Off-white for 0, Deep Blue for 1
        cbar_kws={'label': 'Presence (1) / Absence (0)', 'ticks': [0, 1]},
        linewidths=0.05, 
        linecolor='#dee2e6'
    )

    plt.title(f"Effector Presence/Absence Matrix\n(n = {len(plot_df)} Orthogroups)", fontsize=14)
    plt.xlabel("Fusarium Isolate", fontsize=12)
    plt.ylabel("Orthogroup ID", fontsize=12)
    
    plt.tight_layout()
    plt.savefig(HEATMAP_OUT, dpi=300)
    print(f"✅ Heatmap saved to: {HEATMAP_OUT}")
else:
    print("⚠️ No effector orthogroups found to plot.")

print("\nAnalysis Complete.")
