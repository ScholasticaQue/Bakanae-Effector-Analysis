import os
import subprocess
from pathlib import Path

# -------- CONFIGURATION --------
# The directory where our previous script saved the secretome files
input_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_input")

# Path to the EffectorP 3.0 script
effectorp_script = "/mnt/e/Fusarium_fujikuroi_analysis/scripts/EffectorP-3.0/EffectorP.py"

# Where to save results
output_dir = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
output_dir.mkdir(parents=True, exist_ok=True)

# Final summary list path
effector_gene_list_path = output_dir / "all_unique_effector_ids.txt"

ISOLATES = ["5B", "5C", "B20", "C1995", "CFF1", "CFF2", "E282", "FSU48", 
            "IMI58289", "KSU3368", "M567", "MRC2276", "NCIM1100"]

all_effector_ids = set()


print(f"{'Isolate':<12} | {'Status':<12} | {'Apoplastic':<11} | {'Cytoplasmic':<11} | {'Dual/Other'}")
print("-" * 75)

for isolate in ISOLATES:
    output_path = output_dir / f"{isolate}_effectorp_results.txt"
    
    # Counts for this specific isolate
    apo, cyto, dual = 0, 0, 0
    
    if not output_path.exists():
        print(f"{isolate:<12} |  Missing")
        continue

    with open(output_path, 'r') as f:
        for line in f:
            line = line.strip()
            # Skip empty lines, comments, or headers
            if not line or line.startswith("#") or "Identifier" in line:
                continue
            
            # Use whitespace splitting to handle the space/tab inconsistency
            parts = line.split()
            prediction = " ".join(parts[1:]).lower() # Joins the rest of the line to search for keywords
            
            # Logic: If it's an effector and NOT a non-effector
            if "effector" in prediction and "non-effector" not in prediction:
                gene_id = parts[0]
                all_effector_ids.add(gene_id)
                
                # Sort into categories for your report
                if "apoplastic/cytoplasmic" in prediction or "cytoplasmic/apoplastic" in prediction:
                    dual += 1
                elif "apoplastic" in prediction:
                    apo += 1
                elif "cytoplasmic" in prediction:
                    cyto += 1
        
    print(f"{isolate:<12} |  Parsed  | {apo:<11} | {cyto:<11} | {dual}")

# Write the Master ID list
with open(effector_gene_list_path, "w") as f:
    for eid in sorted(all_effector_ids):
        f.write(f"{eid}\n")

print("-" * 75)
print(f" Total unique effector IDs across all isolates: {len(all_effector_ids)}")
print(f" Master list: {effector_gene_list_path}")
