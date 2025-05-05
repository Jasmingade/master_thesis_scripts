#!/usr/bin/env python
import sys
import os
import csv
import shutil
import requests
from glob import glob
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
import pandas as pd

# ────────────────────────────────────────────────────────────────────────
#  CONFIG
# ────────────────────────────────────────────────────────────────────────
ALLOWED_CATEGORIES = {
    "Peptide Spectral Matches",
    "Other Metadata",
    "Protein Assembly"
}
ALLOWED_FILE_TYPES = {"Text", "Document"}
BATCH_SIZE = 50
failed_downloads = []

# ────────────────────────────────────────────────────────────────────────
def sanitize_folder_name(name):
    return "".join(c if c.isalnum() or c in (" ", "_", "-") else "_" for c in name) \
               .replace(" ", "_")

def download_file(file_info, progress_bar):
    fname, url_link, folder_name, study_id, study_name, file_id = file_info
    os.makedirs(folder_name, exist_ok=True)
    path = os.path.join(folder_name, fname)
    if os.path.isfile(path):
        progress_bar.update(1)
        return f"✔ {fname} exists"
    try:
        resp = requests.get(url_link, stream=True, timeout=130)
        resp.raise_for_status()
        with open(path, "wb") as f:
            shutil.copyfileobj(resp.raw, f)
        progress_bar.update(1)
        return f"Downloaded {fname}"
    except Exception as e:
        progress_bar.update(1)
        failed_downloads.append({
            "PDC Study ID": study_id,
            "Study Name": study_name,
            "File ID": file_id,
            "File Name": fname,
            "Download Link": url_link,
            "Error": str(e)
        })
        return f"✖ {fname} failed: {e}"

def downloadOrganize(file_name, delimiter, manifest_folder):
    # infer plex from the manifest folder name, e.g. "TMT11_manifests" → "TMT11"
    plex = manifest_folder.split("_")[0]
    file_list = []

    with open(file_name) as f:
        reader = csv.DictReader(f, delimiter=delimiter)
        for row in reader:
            if (row["Data Category"] not in ALLOWED_CATEGORIES
                or row["File Type"] not in ALLOWED_FILE_TYPES):
                continue

            fname = row["File Name"]
            url = row["File Download Link"]
            study_id   = row["PDC Study ID"]
            study_name = row["Study Name"]
            version    = row["PDC Study Version"]
            folder = row.get("Run Metadata ID","") or version

            study_folder = sanitize_folder_name(f"{study_id}_{study_name}")
            base_dir     = os.path.join(os.getcwd(),
                                        f"data/raw/{plex}_downloads",
                                        study_folder,
                                        version,
                                        row["Data Category"],
                                        row["File Type"])
            file_list.append((fname, url, base_dir, study_id, study_name, row["File ID"]))

    for i in range(0, len(file_list), BATCH_SIZE):
        batch = file_list[i:i+BATCH_SIZE]
        with tqdm(batch, desc=f"{plex} {manifest_folder} batch {i//BATCH_SIZE+1}",
                  unit="file") as bar:
            with ThreadPoolExecutor(max_workers=5) as ex:
                for res in ex.map(lambda fi: download_file(fi, bar), batch):
                    print(res)

    if failed_downloads:
        df = pd.DataFrame(failed_downloads)
        df.to_csv("failed_downloads.csv", index=False)
        print("⚠ Failed → failed_downloads.csv")

def main():
    if len(sys.argv) != 3:
        print("Usage: python 00_downloadPDCData.py <manifest> <1>")
        sys.exit(1)
    manifest = sys.argv[1]
    opt      = sys.argv[2]
    if not os.path.isfile(manifest):
        sys.exit("Manifest not found.")
    delim = "," if manifest.lower().endswith(".csv") else "\t"
    if opt == "1":
        folder = os.path.basename(os.path.dirname(manifest))
        downloadOrganize(manifest, delim, folder)
    else:
        sys.exit("Invalid option.")

if __name__ == "__main__":
    main()
