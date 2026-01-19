import os
import subprocess
from pathlib import Path
from datetime import datetime

# Paths - Ensure these match your environment
INPUT_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")
OUTPUT_DIR = INPUT_DIR / "targetp_results"

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Find all .faa files
faa_files = sorted(list(INPUT_DIR.glob("*.faa")))
print(f"TargetP 2.0: Found {len(faa_files)} files. Starting at {datetime.now().strftime('%H:%M:%S')}")

for faa_file in faa_files:
    base = faa_file.stem
    output_subdir = OUTPUT_DIR / base
    
    # Check if this isolate was already processed
    if output_subdir.exists() and any(output_subdir.iterdir()):
        print(f"[-] Skipping {base}: Already processed.")
        continue

    output_subdir.mkdir(parents=True, exist_ok=True)
    print(f"\n[+] Processing: {faa_file.name}")

    # TargetP 2.0 Command
    # -org euk: for fungi/eukaryotes
    # -format short: produces the summary .txt table we need
    cmd = [
        "targetp", 
        "-fasta", str(faa_file),
        "-org", "euk",
        "-format", "short",
        "-out_dir", str(output_subdir)
    ]

    try:
        # TargetP usually outputs to stdout, but we use out_dir for the files
        subprocess.run(cmd, check=True)
        print(f"    Finished {base} successfully.")
    
    except subprocess.CalledProcessError as e:
        print(f"    !!! Error running TargetP on {base} !!!")
        continue

print(f"\nTargetP runs completed at {datetime.now().strftime('%H:%M:%S')}")

