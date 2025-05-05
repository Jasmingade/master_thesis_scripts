#!/usr/bin/env python
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm

# ────────────────────────────────────────────────────────
#  Adjust these if you move your manifest folders
# ────────────────────────────────────────────────────────
PLEXES = ["TMT10", "TMT11"]
MANIFEST_BASE = "data/raw/manifest_files"   # contains e.g. TMT10_manifests, TMT11_manifests
SCRIPT       = os.path.join(os.getcwd(), "scr/python/00_downloadPDCData.py")

def process_manifest(manifest):
    cmd = ["python", SCRIPT, os.path.abspath(manifest), "1"]
    subprocess.run(cmd, check=True)

def main():
    tasks = []
    for plex in PLEXES:
        folder = f"{MANIFEST_BASE}/{plex}_manifests"
        files  = [os.path.join(folder, f)
                  for f in os.listdir(folder)
                  if f.endswith((".csv",".tsv"))]
        for mf in files:
            tasks.append(mf)

    with tqdm(total=len(tasks), desc="All Manifests", unit="manifest") as pbar:
        with ThreadPoolExecutor(max_workers=4) as ex:
            for _ in ex.map(lambda m: (process_manifest(m), pbar.update()), tasks):
                pass

    print("✅ All downloads triggered.")

if __name__ == "__main__":
    main()
