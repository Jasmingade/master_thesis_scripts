#!/usr/bin/env python

"""
Download and organize PDC files based on a manifest.

This script performs the following steps:

1. Reads the manifest (CSV/TSV) and filters by allowed categories/types.
2. Builds download tasks, sanitizing folder names.
3. Downloads files in batches with a thread pool, showing per-batch progress.
4. Logs any failed downloads to a CSV report.
"""

import sys
import os
import csv
import shutil
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd

# Allowed data filters
ALLOWED_CATEGORIES = {"Peptide Spectral Matches", "Other Metadata", "Protein Assembly"}
ALLOWED_FILE_TYPES = {"Text", "Document"}

# Download parameters
BATCH_SIZE  = 50
max_workers = 5

# Collect failures
failed_downloads = []

def sanitize_folder_name(name):
    """Replace disallowed characters and spaces for filesystem safety."""
    return "".join(c if c.isalnum() or c in (" ", "_", "-") else "_" for c in name).replace(" ", "_")

def download_file(task, progress_bar):
    """
    Download a single file; skip if exists, handle errors, and update progress.
    """
    fname, url, dest_folder, study_id, study_name, file_id, manifest_folder = task
    dest_path = os.path.join(dest_folder, fname)

    if os.path.isfile(dest_path):
        progress_bar.update(1)
        return f"Skipped (exists): {fname}"

    os.makedirs(dest_folder, exist_ok=True)
    try:
        resp = requests.get(url, stream=True, timeout=130, headers={'Cache-Control': 'no-cache'})
        resp.raise_for_status()
        with open(dest_path, 'wb') as f:
            shutil.copyfileobj(resp.raw, f)
        progress_bar.update(1)
        return f"Downloaded: {fname}"
    except Exception as e:
        progress_bar.update(1)
        failed_downloads.append({
            "PDC Study ID": study_id,
            "Study Name":   study_name,
            "File ID":      file_id,
            "Manifest":     manifest_folder,
            "File Name":    fname,
            "URL":          url,
            "Error":        str(e)
        })
        return f"Failed: {fname} ({e})"

def download_organize(manifest_file, delimiter, manifest_folder):
    """
    Parse manifest, filter entries, and download in batches.
    """
    tasks = []
    with open(manifest_file, newline='') as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            if row['Data Category'] not in ALLOWED_CATEGORIES or row['File Type'] not in ALLOWED_FILE_TYPES:
                continue

            study_folder = f"{sanitize_folder_name(row['PDC Study ID'])}_{sanitize_folder_name(row['Study Name'])}"
            base_dir = os.path.join(os.getcwd(), "data/raw/TMT11_downloads", study_folder)
            subfolder = row['Run Metadata ID'] if row['Run Metadata ID'] not in ("", "null") else ""
            dest = os.path.join(base_dir, row['PDC Study Version'], row['Data Category'], subfolder, row['File Type'])
            tasks.append((
                row['File Name'], row['File Download Link'],
                dest, row['PDC Study ID'], row['Study Name'],
                row['File ID'], manifest_folder
            ))

    # Process in batches
    for i in range(0, len(tasks), BATCH_SIZE):
        batch = tasks[i:i+BATCH_SIZE]
        with tqdm(total=len(batch), desc=f"Batch {i//BATCH_SIZE+1}", unit="file") as pbar:
            with ThreadPoolExecutor(max_workers=max_workers) as executor:
                futures = [executor.submit(download_file, task, pbar) for task in batch]
                for fut in as_completed(futures):
                    print(fut.result())
        print(f"✅ Batch {i//BATCH_SIZE+1} complete.")

    # Write failures, if any
    if failed_downloads:
        df = pd.DataFrame(failed_downloads)
        log_path = os.path.join(os.getcwd(), "reports/logs/failed_downloads.csv")
        df.to_csv(log_path, index=False)
        print(f"⚠ Failed downloads logged to: {log_path}")

def main():
    """Parse arguments and dispatch to download_organize."""
    if len(sys.argv) != 3:
        print("Usage: python 02_download_pdc_data.py <manifest_file> <option: 1>")
        sys.exit(1)

    file_name, option = sys.argv[1], sys.argv[2]
    if not os.path.isfile(file_name):
        sys.exit(f"Error: Manifest '{file_name}' not found.")

    delim = ',' if file_name.lower().endswith('.csv') else '\t' if file_name.lower().endswith('.tsv') else None
    if delim is None:
        sys.exit("Error: Manifest must be .csv or .tsv")

    if option == "1":
        manifest_folder = os.path.basename(os.path.dirname(file_name))
        download_organize(file_name, delim, manifest_folder)
    else:
        sys.exit("Error: Invalid option. Use '1' to download.")

if __name__ == "__main__":
    main()
