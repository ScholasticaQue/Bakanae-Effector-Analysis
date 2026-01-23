"""
Pathogen Effector Analysis: DeepTMHMM Transmembrane Prediction Wrapper
Author: Scholastica Quaicoe 
Institution: Africa Rice Center

Description:
This script automates transmembrane domain prediction for 13 F. fujikuroi isolates. 
It implements a 'split-run-merge' strategy:
1. Splits large proteomes into 1000-sequence chunks to optimize BioLib performance.
2. Executes DeepTMHMM (DTU/DeepTMHMM) via BioLib.
3. Merges chunked results back into a single comprehensive output per isolate.
"""

import os
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO
from datetime import datetime

# ===================== CONFIG (Relative Paths) =====================
INPUT_DIR = Path("data/protein_input")
FINAL_OUTPUT_DIR = Path("results/deeptmhmm_final")
WORK_DIR = Path("temp/biolib_runs")

# Ensure directories exist
FINAL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
WORK_DIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 1000
# ===================================================================

def split_fasta(file_path):
    """Split a FASTA into chunks to prevent BioLib memory errors."""
    records = list(SeqIO.parse(file_path, "fasta"))
    base_name = file_path.stem
    chunk_files = []

    for i in range(0, len(records), CHUNK_SIZE):
        chunk = records[i:i + CHUNK_SIZE]
        chunk_num = (i // CHUNK_SIZE) + 1
        output_name = INPUT_DIR / f"{base_name}_part{chunk_num}.faa"

        if not output_name.exists():
            SeqIO.write(chunk, output_name, "fasta")
        chunk_files.append(output_name)

    return chunk_files

def merge_results(isolate_name, result_dirs):
    """Merges prediction files from chunks into one isolate file."""
    merged_file = FINAL_OUTPUT_DIR / f"{isolate_name}_all_predictions.3line"
    with open(merged_file, "w") as outfile:
        for folder in result_dirs:
            # DeepTMHMM usually outputs a .3line or T3SE file
            prediction_file = folder / "predicted_topologies.3line"
            if prediction_file.exists():
                with open(prediction_file, "r") as infile:
                    outfile.write(infile.read())
    print(f"  [Merge] Created consolidated file: {merged_file.name}")

# ===================== MAIN PIPELINE =====================
print(f" DeepTMHMM Pipeline started: {datetime.now()}")

original_files = sorted([f for f in INPUT_DIR.glob("*.faa") if "_part" not in f.stem])

for faa_file in original_files:
    isolate_name = faa_file.stem
    print(f"\n[+] Processing: {isolate_name}")

    parts = split_fasta(faa_file)
    processed_dirs = []

    for part in parts:
        part_base = part.stem
        chunk_dest = FINAL_OUTPUT_DIR / "chunks" / part_base
        chunk_dest.mkdir(parents=True, exist_ok=True)

        # BioLib execution
        before_dirs = set(os.listdir(WORK_DIR))
        cmd = ["biolib", "run", "DTU/DeepTMHMM", "--fasta", str(part.resolve())]

        try:
            subprocess.run(cmd, cwd=WORK_DIR, check=True, capture_output=True, text=True)
            
            after_dirs = set(os.listdir(WORK_DIR))
            new_dirs = list(after_dirs - before_dirs)
            
            if new_dirs:
                run_dir = WORK_DIR / new_dirs[0]
                shutil.move(run_dir, chunk_dest)
                processed_dirs.append(chunk_dest)
                print(f"+ Chunk {part_base} complete")

        except subprocess.CalledProcessError as e:
            print(f"- Failed: {part_base}")

    # Merge all chunks for this isolate
    if processed_dirs:
        merge_results(isolate_name, processed_dirs)

print(f"\n Pipeline finished: {datetime.now()}")
