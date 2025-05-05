#!/usr/bin/env python

"""
Wrapper to launch PDC file downloads from a merged manifest.

This script performs the following steps:

1. Defines paths and parameters for download.
2. Ensures the output directory exists.
3. Submits the merged manifest to the download engine in parallel,
   showing progress with a tqdm bar.
4. Reports completion status (test or full run).
"""

import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm  # Progress bar

# ------------------------------------------------------------------
# Configuration
# ------------------------------------------------------------------
output_dir  = "data/raw/TMT11_downloads"
script_path = os.path.join(os.getcwd(), "scr/python/01_downloads/02_download_pdc_data.py")
merged_manifest = "data/raw/manifest_files/TMT11_manifests"

# Set to True to run a small test batch
TEST_MODE = False

# Number of worker threads for manifest processing
max_workers = 5

# Ensure the download directory exists
os.makedirs(output_dir, exist_ok=True)

def process_manifest(manifest, progress_bar):
    """
    Invoke the detailed download script on the merged manifest.
    """
    cmd = [
        "python", script_path,
        os.path.abspath(manifest),
        "1"  # Option 1 = Download only
    ]
    subprocess.run(cmd, check=True)
    progress_bar.update(1)

if __name__ == "__main__":
    # Process the single merged manifest with a progress bar
    with tqdm(total=1, desc="Processing Merged Manifest", unit="file") as pbar:
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            future = executor.submit(process_manifest, merged_manifest, pbar)
            future.result()

    print("Test run completed." if TEST_MODE else "Download process completed.")
