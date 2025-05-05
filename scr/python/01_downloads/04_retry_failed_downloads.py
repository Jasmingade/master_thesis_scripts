#!/usr/bin/env python
"""
Retry failed PDC downloads.

This script retries downloading files listed in `failed_downloads.csv`,
applying configurable retry attempts with exponential backoff
and reporting progress via a tqdm progress bar.

Steps:
1. Configure retry parameters (max attempts, backoff factor, timeout).
2. Read the CSV of previously failed downloads.
3. Build download tasks, sanitizing study folder names.
4. Attempt each download up to MAX_RETRIES, waiting BACKOFF_FACTOR**attempt between tries.
5. Report success or ultimate failure for each file.
"""
import os
import shutil
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm  # Progress bar

# ------------------------------------------------------------------
# 1. Retry configuration
# ------------------------------------------------------------------
MAX_RETRIES    = 3   # Number of attempts per file
BACKOFF_FACTOR = 2   # Exponential backoff multiplier (seconds)
TIMEOUT        = 60  # Request timeout in seconds

def sanitize_folder_name(name):
    """Sanitize folder names by replacing spaces and invalid characters."""
    return "".join(
        c if c.isalnum() or c in (" ", "_", "-") else "_"
        for c in name
    ).replace(" ", "_")

def download_file(file_info, progress_bar):
    """
    Attempt to download a single file with retry logic.

    file_info: tuple of
        (filename, url, dest_folder, study_id, study_name, file_id, manifest_folder)
    """
    fname, url_link, folder_name, study_id, study_name, file_id, manifest_folder = file_info
    file_path = os.path.join(folder_name, fname)

    # Skip if already present
    if os.path.isfile(file_path):
        progress_bar.update(1)
        return f"{fname} already exists. Skipping."

    os.makedirs(folder_name, exist_ok=True)

    # Retry loop
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.get(url_link, stream=True, timeout=TIMEOUT)
            response.raise_for_status()

            with open(file_path, 'wb') as f:
                shutil.copyfileobj(response.raw, f)

            progress_bar.update(1)
            return f"Downloaded: {fname}"

        except requests.exceptions.RequestException as e:
            print(f"Attempt {attempt}/{MAX_RETRIES} failed for {fname}: {e}")
            if attempt < MAX_RETRIES:
                wait = BACKOFF_FACTOR ** attempt
                print(f"Retrying in {wait} seconds...")
                time.sleep(wait)
            else:
                progress_bar.update(1)
                return f"❌ Failed to download {fname} after {MAX_RETRIES} attempts."

def retry_failed_downloads():
    """Load failed_downloads.csv and retry each entry in parallel."""
    failed_csv = "reports/logs/failed_downloads.csv"
    if not os.path.exists(failed_csv):
        print("⚠️ No failed downloads found.")
        return

    df = pd.read_csv(failed_csv)
    tasks = []
    for _, row in df.iterrows():
        study_folder = f"{sanitize_folder_name(row['PDC Study ID'])}_{sanitize_folder_name(row['Study Name'])}"
        dest
