import subprocess
from pathlib import Path

# Paths - adjust as needed
SIGNALP_DIR = Path("/Fusarium_fujikuroi_analysis/signalp6_fast/signalp-6-package")
INPUT_DIR = Path("/Fusarium_fujikuroi_analysis/processed_data/proteins")
OUTPUT_DIR = SIGNALP_DIR / "signalp_results"

# Make sure output directory exists
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Glob all .faa files under Country_*/Isolate_*/ recursively
faa_files = list(INPUT_DIR.glob("Country_*/Isolate_*/*.faa"))

print(f"Found {len(faa_files)} .faa files to process.")

for faa_file in faa_files:
    base = faa_file.stem
    output_subdir = OUTPUT_DIR / base
    output_subdir.mkdir(parents=True, exist_ok=True)
    
    print(f"Running SignalP on {base}...")
    
    # Build the command as a list
    cmd = [
        "python", "-m", "signalp",
        "--fastafile", str(faa_file),
        "--output_dir", str(output_subdir),
        "--format", "txt",
        "--organism", "eukarya",
        "--mode", "fast"
    ]
    
    # Run the command and wait for it to finish
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print(f"Finished: {base}")
        print(f"Output stored in: {output_subdir}")
    else:
        print(f"Error running SignalP on {base}")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)

print("All done.")
