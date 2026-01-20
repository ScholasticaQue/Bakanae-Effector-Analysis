import os
import subprocess
import shutil
from pathlib import Path
from Bio import SeqIO
from datetime import datetime

# ===================== CONFIG =====================
INPUT_DIR = Path("/protein_input")
FINAL_OUTPUT_DIR = INPUT_DIR / "deeptmhmm_results"
FINAL_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CHUNK_SIZE = 1000
WORK_DIR = INPUT_DIR / "biolib_runs"
WORK_DIR.mkdir(exist_ok=True)

# ==================================================

def split_fasta(file_path):
    """Split a FASTA into chunks of CHUNK_SIZE."""
    records = list(SeqIO.parse(file_path, "fasta"))
    base_name = file_path.stem
    chunk_files = []

    for i in range(0, len(records), CHUNK_SIZE):
        chunk = records[i:i + CHUNK_SIZE]
        chunk_num = (i // CHUNK_SIZE) + 1
        output_name = INPUT_DIR / f"{base_name}_part{chunk_num}.faa"

        if not output_name.exists():
            SeqIO.write(chunk, output_name, "fasta")
            print(f"  [Split] {output_name.name} ({len(chunk)} seqs)")

        chunk_files.append(output_name)

    return chunk_files


# ===================== MAIN =====================
print(f"\n DeepTMHMM Pipeline started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

# Only ORIGINAL FASTA files (exclude _part files)
original_files = sorted([
    f for f in INPUT_DIR.glob("*.faa")
    if "_part" not in f.stem
])

for faa_file in original_files:
    print(f"\n[+] Processing isolate: {faa_file.stem}")

    parts = split_fasta(faa_file)

    for part in parts:
        part_base = part.stem
        final_dest = FINAL_OUTPUT_DIR / part_base

        if final_dest.exists():
            print(f"  [-] Skipping {part_base} (already processed)")
            continue

        print(f"  [+] Running DeepTMHMM on {part.name}")

        # Capture directories before run
        before_dirs = set(os.listdir(WORK_DIR))

        cmd = [
            "biolib", "run", "DTU/DeepTMHMM",
            "--fasta", str(part.resolve())
        ]

        try:
            subprocess.run(
                cmd,
                cwd=WORK_DIR,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # Capture directories after run
            after_dirs = set(os.listdir(WORK_DIR))
            new_dirs = list(after_dirs - before_dirs)

            if not new_dirs:
                print(f"  [!] No output directory detected for {part_base}")
                continue

            # BioLib creates exactly ONE new directory per run
            run_dir = WORK_DIR / new_dirs[0]
            shutil.move(run_dir, final_dest)

            print(f"  [✓] Results saved → {final_dest}")

        except subprocess.CalledProcessError as e:
            print(f"  [✗] DeepTMHMM failed for {part_base}")
            print(e.stderr)

print(f"\n Pipeline finished at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
