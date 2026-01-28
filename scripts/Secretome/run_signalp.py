"""
Pathogen Effector Analysis: SignalP 6.0 Batch Processor
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
Automates the prediction of signal peptides across 13 global Fusarium fujikuroi 
isolates using SignalP 6.0 in 'fast' mode. This is Step 1 of the secretome 
filtering pipeline.
"""
import subprocess
from pathlib import Path

# ----------- CONFIGURATION (Relative Paths) -----------
INPUT_DIR = Path("data/raw_proteins")  
OUTPUT_DIR = Path("results/signalp_outputs")

 # Ensure SignalP directory is in the Python path
sys.path.append(str(SIGNALP_DIR))

# Create output directory
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Find all .faa files and sort them
faa_files = sorted(list(INPUT_DIR.glob("*.faa")))
print(f"Found {len(faa_files)} .faa files to process.")

for faa_file in faa_files:
    base = faa_file.stem
    output_subdir = OUTPUT_DIR / base
    
    # Check if this file has already been processed to save time
    if output_subdir.exists() and any(output_subdir.iterdir()):
        print(f"Skipping {base}: Output already exists.")
        continue

    output_subdir.mkdir(parents=True, exist_ok=True)

    print(f">>> Running SignalP on: {faa_file.name}")

    # Build and run command
    # We point '--fastafile' directly to the original file path
    cmd = [
        "python3", "-m", "signalp",
        "--fastafile", str(faa_file),
        "--output_dir", str(output_subdir),
        "--format", "txt",
        "--organism", "eukarya",
        "--mode", "fast"
    ]

    try:
        # Run the command and capture output
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"Success: {base} completed.")
    
    except subprocess.CalledProcessError as e:
        print(f"!!! Error running SignalP on {base} !!!")
        print("STDERR:", e.stderr)
        # We continue to the next file even if one fails
        continue

print("\nProcessing complete.")

