import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import re

# --- CONFIGURATION ---
results_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
clstr_file = results_dir / "effectors_clustered.clstr"
matrix_output = results_dir / "clustered_effector_matrix.csv"
heatmap_output = results_dir / "clustered_effector_heatmap.png"

ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

# 1. Parse CD-HIT .clstr file
clusters = {}
current_cluster = None

with open(clstr_file, 'r') as f:
    for line in f:
        if line.startswith(">Cluster"):
            current_cluster = line.strip().replace(">", "")
            clusters[current_cluster] = set()
        else:
            # Extract isolate name from the tagged ID (e.g., 5B|5B_000012)
            match = re.search(r'>(\w+)\|', line)
            if match:
                isolate = match.group(1)
                clusters[current_cluster].add(isolate)

# 2. Create the Presence/Absence Matrix
matrix_data = []
for cluster, isolate_set in clusters.items():
    row = {"Cluster_ID": cluster}
    for iso in ISOLATES:
        row[iso] = 1 if iso in isolate_set else 0
    matrix_data.append(row)

df = pd.DataFrame(matrix_data)
df['Conservation'] = df[ISOLATES].sum(axis=1)

# 3. Categorize
def get_cat(count):
    if count == len(ISOLATES): return "Core"
    if count == 1: return "Unique"
    return "Accessory"

df['Category'] = df['Conservation'].apply(get_cat)
df.sort_values(by=['Conservation', 'Cluster_ID'], ascending=[False, True], inplace=True)
df.to_csv(matrix_output, index=False)

# 4. Generate the Summary
print("\n--- Final Effector Conservation Summary ---")
print(df['Category'].value_counts())

# 5. Generate Heatmap
plt.figure(figsize=(12, 10))
sns.set_theme(style="white")
plot_data = df[ISOLATES]

sns.heatmap(plot_data, cmap="YlGnBu", cbar=False, yticklabels=False)
plt.title(f"Clustered Effector Presence/Absence\n({len(df)} Unique Effector Families)", fontsize=16)
plt.xlabel("Isolate", fontsize=12)
plt.ylabel("Effector Clusters (Sorted by Conservation)", fontsize=12)

plt.tight_layout()
plt.savefig(heatmap_output, dpi=300)
print(f"\n Heatmap saved to: {heatmap_output}")
