import pandas as pd
from pathlib import Path

# --- CONFIGURATION ---
results_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
output_matrix = results_dir / "effector_presence_absence_matrix.csv"

ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

# 1. Load all unique IDs
master_list_path = results_dir / "all_unique_effector_ids.txt"
with open(master_list_path, 'r') as f:
    all_ids = [line.strip() for line in f if line.strip()]

# 2. Build the matrix
matrix_data = []

for eid in sorted(all_ids):
    row = {"Effector_ID": eid}
    # Check each isolate file for this specific ID
    for iso in ISOLATES:
        file_path = results_dir / f"{iso}_effectorp_results.txt"
        found = 0
        if file_path.exists():
            with open(file_path, 'r') as f:
                if eid in f.read():
                    found = 1
        row[iso] = found
    matrix_data.append(row)

# 3. Save to CSV
df = pd.DataFrame(matrix_data)
df.to_csv(output_matrix, index=False)

print(f" Matrix created with {len(all_ids)} effectors.")
print(f" Saved to: {output_matrix}")

# Quick Summary Statistics
df['Total_Count'] = df[ISOLATES].sum(axis=1)
core_effectors = len(df[df['Total_Count'] == len(ISOLATES)])
singleton_effectors = len(df[df['Total_Count'] == 1])

print(f"\n--- Preliminary Insights ---")
print(f"Core Effectors (in all 13): {core_effectors}")
print(f"Accessory/Unique Effectors: {len(all_ids) - core_effectors}")
print(f"Singletons (only in 1 isolate): {singleton_effectors}")
