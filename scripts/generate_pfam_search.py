"""
Pathogen Effector Annotation: Pfam Domain Search Generator
Author: Scholastica Quaicoe
Institution: Africa Rice Center

Description:
This script automates the generation of a high-throughput HMMER search 
against the Pfam-A database for candidate effector proteins identified 
in 13 global Fusarium fujikuroi isolates.
"""
from pathlib import Path

# ----------- CONFIGURATION (Relative Paths for Portability) -----------
# Ensure the Pfam-A.hmm file and input directory are correctly mapped
pfam_db = Path("databases/Pfam-A.hmm")
input_dir = Path("results/effector_candidates")
output_dir = input_dir / "pfam_annotations"
output_dir.mkdir(parents=True, exist_ok=True)

# ----------- BUILD COMMANDS -----------
commands = []

# Search for all candidate FASTA files
faa_files = list(input_dir.glob("*_filtered_candidates.faa"))

for faa_file in faa_files:
    sample = faa_file.stem.replace("_filtered_candidates", "")
    domtblout = output_dir / f"{sample}_pfam.domtblout"
    rawout = output_dir / f"{sample}_pfam_full.txt"
    
    cmd = f"hmmsearch --domtblout {domtblout} {pfam_db} {faa_file} > {rawout}"
    commands.append(cmd)

# ----------- WRITE BASH SCRIPT -----------
bash_script_path = output_dir / "run_all_pfam.sh"
with open(bash_script_path, "w") as f:
    f.write("#!/bin/bash\n\n")
    for cmd in commands:
        f.write(cmd + "\n")

print(f" Bash script written: {bash_script_path}")
