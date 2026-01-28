"""
Pathogen Effector Analysis: TargetP 2.0 Subcellular Localization Wrapper
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
This script automates TargetP 2.0 to predict subcellular localization for 
13 Fusarium fujikuroi isolates. It specifically filters for the presence of 
Signal Peptides (SP) while excluding proteins targeted to other organelles 
(e.g., Mitochondria), ensuring a high-confidence secretome.
"""

import subprocess
from pathlib import Path
from datetime import datetime

# ----------- CONFIGURATION (Relative Paths) -----------
INPUT_DIR = Path("data/protein_input")
OUTPUT_DIR = Path("results/targetp_outputs")

# Ensure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Find all .faa files (representing 13 global isolates)
faa_files = sorted(list(INPUT_DIR.glob("*.faa")))
print(f"TargetP 2.0: Found {len(faa_files)} files. Starting at {datetime.now().strftime('%H:%M:%S')}")

for faa_file in faa_files:
    base = faa_file.stem
    output_subdir = OUTPUT_DIR / base
    
    # Skip if already processed (Efficiency for large genomic datasets)
    if output_subdir.exists() and any(output_subdir.iterdir()):
        print(f"[-] Skipping {base}: Already processed.")
        continue

    output_subdir.mkdir(parents=True, exist_ok=True)
    print(f"\n[+] Processing: {faa_file.name}")

    # TargetP 2.0 Command:
    # -org euk: Configured for fungi (Eukaryote)
    # -format short: Generates the summary table for downstream parsing
    cmd = [
        "targetp", 
        "-fasta", str(faa_file),
        "-org", "euk",
        "-format", "short",
        "-out_dir", str(output_subdir)
    ]

    try:
        # Execute TargetP via subprocess
        subprocess.run(cmd, check=True)
        print(f"    Finished {base} successfully.")
    
    except subprocess.CalledProcessError as e:
        print(f"    !!! Error running TargetP on {base} !!!")
        continue

print(f"\nTargetP runs completed at {datetime.now().strftime('%H:%M:%S')}")
