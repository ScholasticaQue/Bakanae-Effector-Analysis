import os
import subprocess
from pathlib import Path

# --- CONFIGURATION ---
base_path = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input/effectorp_results")
input_dir = base_path / "dnds_input_groups"
output_dir = base_path / "dnds_alignments"
os.makedirs(output_dir, exist_ok=True)

fasta_files = list(input_dir.glob("*.fasta"))

print(f" Starting alignments for {len(fasta_files)} clusters...")

for fasta in fasta_files:
    output_file = output_dir / f"{fasta.stem}_aligned.faa"
    # Using --auto to let MAFFT choose the best strategy
    cmd = f"mafft --auto {fasta} > {output_file}"
    
    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f" Aligned: {fasta.name}")
    except subprocess.CalledProcessError:
        print(f" Error aligning {fasta.name}")

print(f"\n All alignments saved to: {output_dir}")
