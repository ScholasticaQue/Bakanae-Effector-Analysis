import os
import subprocess
from pathlib import Path
from Bio import SeqIO  # Essential import
from datetime import datetime

# Paths
INPUT_DIR = Path("/mnt/e/Fusarium_fujikuroi_analysis/cds_data/protein_input")
OUTPUT_DIR = INPUT_DIR / "deeptmhmm_results"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 5000 

def split_fasta(file_path):
    """Splits a large FASTA into smaller chunks for BioLib."""
    records = list(SeqIO.parse(file_path, "fasta"))
    base_name = file_path.stem
    chunk_files = []
    
    for i in range(0, len(records), CHUNK_SIZE):
        chunk = records[i:i + CHUNK_SIZE]
        chunk_num = (i // CHUNK_SIZE) + 1
        output_name = INPUT_DIR / f"{base_name}_part{chunk_num}.faa"
        
        # Only write if it doesn't exist
        if not output_name.exists():
            SeqIO.write(chunk, output_name, "fasta")
            print(f"  [Split] Created {output_name.name} ({len(chunk)} seqs)")
        
        chunk_files.append(output_name)
    return chunk_files

# 1. First, find all original sanitized files
original_files = sorted(list(INPUT_DIR.glob("*.faa")))

print(f"Starting DeepTMHMM Pipeline at {datetime.now().strftime('%H:%M:%S')}")

for faa_file in original_files:
    print(f"\n[+] Processing Isolate: {faa_file.stem}")
    
    # 2. Split the file into chunks
    parts = split_fasta(faa_file)
    
    # 3. Run DeepTMHMM on each chunk
    for part in parts:
        part_base = part.stem
        part_output = OUTPUT_DIR / part_base
        
        if part_output.exists() and any(part_output.iterdir()):
            print(f"  [-] Skipping {part_base}: Already processed.")
            continue

        print(f"  [+] BioLib Run: {part.name}...")
        
        cmd = [
            "biolib", "run", "DTU/DeepTMHMM", 
            "--fasta", str(part.absolute()), 
            "--output_dir", str(part_output)
        ]

        try:
            subprocess.run(cmd, check=True)
            print(f"  [Success] Finished {part_base}")
        except subprocess.CalledProcessError:
            print(f"  [!] Failed to run BioLib on {part_base}")

print(f"\nPipeline finished at {datetime.now().strftime('%H:%M:%S')}")


