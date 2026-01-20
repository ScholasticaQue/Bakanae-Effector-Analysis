import os
import subprocess
from pathlib import Path

# Paths to your sanitized data
INPUT_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")
OUTPUT_DIR = INPUT_DIR / "deeptmhmm_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Find all sanitized files
fasta_files = sorted(list(INPUT_DIR.glob("*.faa")))



for faa_file in fasta_files:
    print(f"[+] Running DeepTMHMM via BioLib: {faa_file.name}")
    
    # The correct command for the BioLib.com CLI
    cmd = [
        "biolib", "run", "DTU/DeepTMHMM",
        "--fasta", str(faa_file.absolute())
    ]
    
    try:
        # We use check=True so the script stops if there's a connection/login error
        subprocess.run(cmd, check=True)
        print(f"    Finished {base} successfully.")
    
    except subprocess.CalledProcessError as e:
        print(f"    !!! Error running DeepTMHMM on {base} !!!")
        print(f"    Reason: {e}")
        continue

print(f"\nAll DeepTMHMM runs completed at {datetime.now().strftime('%H:%M:%S')}")

