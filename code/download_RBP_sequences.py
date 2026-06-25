#!/usr/bin/env python3
"""
Download human RBP sequences (cDNA and protein) based on ENSG gene IDs.
Uses Ensembl REST API with batch queries.
"""

import requests
import time
import os
import sys
import re

# --- Config ---
ENSG_FILE = "./data/raw/RBP_geneid_list.txt"
OUTPUT_CDNA = "./data/raw/human_RBP_cdna.fasta"
OUTPUT_PROTEIN = "./data/raw/human_RBP_protein.fasta"
OUTPUT_LOG = "./data/raw/RBP_download_log.txt"
BATCH_SIZE = 50  # Ensembl max batch size
RATE_LIMIT_SLEEP = 1.0  # seconds between batches (conservative)

BASE_URL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "text/x-fasta"}


def read_ensg_ids(filepath):
    """Read ENSG IDs from file. Format: each line has 'index\\tENSG...'"""
    ids = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split()
            for p in parts:
                if re.match(r'^ENSG\d{11}(\.\d+)?$', p):
                    ids.append(p)
                    break
    return ids


def fetch_batch(ids, seq_type):
    """Fetch sequences for a batch of IDs from Ensembl REST API."""
    url = f"{BASE_URL}/sequence/id?type={seq_type}"
    data = {"ids": ids}

    try:
        r = requests.post(url, headers=HEADERS, json=data, timeout=30)
        r.raise_for_status()
        return r.text
    except requests.exceptions.RequestException as e:
        return None


def write_fasta_entries(fasta_text, f_out, log_handle):
    """Parse and write FASTA entries, skipping empty sequences."""
    if not fasta_text or fasta_text.strip() == "":
        return 0

    entries = re.findall(r'(>.*?)(?=\n>|\Z)', fasta_text, re.DOTALL)
    count = 0
    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:]).replace(' ', '').replace('\n', '')
        if len(seq) > 0:
            # Reformat to 60-char lines
            f_out.write(header + '\n')
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + '\n')
            count += 1
        else:
            log_handle.write(f"WARNING: Empty sequence for {header}\n")
    return count


def main():
    print("=" * 60)
    print("Human RBP Sequence Downloader")
    print("=" * 60)

    # Read ENSG IDs
    ensg_ids = read_ensg_ids(ENSG_FILE)
    print(f"\nRead {len(ensg_ids)} ENSG IDs from {ENSG_FILE}")

    if not ensg_ids:
        print("ERROR: No valid ENSG IDs found!")
        sys.exit(1)

    # Open output files
    total_cdna_seqs = 0
    total_protein_seqs = 0
    failed_batches = 0
    total_batches = (len(ensg_ids) + BATCH_SIZE - 1) // BATCH_SIZE

    with open(OUTPUT_CDNA, 'w') as f_cdna, \
         open(OUTPUT_PROTEIN, 'w') as f_prot, \
         open(OUTPUT_LOG, 'w') as f_log:

        f_log.write(f"RBP sequence download log\n")
        f_log.write(f"Total ENSG IDs: {len(ensg_ids)}\n")
        f_log.write(f"Batch size: {BATCH_SIZE}, Sleep: {RATE_LIMIT_SLEEP}s\n\n")

        for i in range(0, len(ensg_ids), BATCH_SIZE):
            batch = ensg_ids[i:i+BATCH_SIZE]
            batch_num = i // BATCH_SIZE + 1
            print(f"\rBatch {batch_num}/{total_batches} (IDs {i+1}-{i+len(batch)})...", end="", flush=True)

            # --- Download cDNA ---
            cdna_raw = fetch_batch(batch, "cdna")
            if cdna_raw is None:
                f_log.write(f"Batch {batch_num} (cDNA): FAILED - HTTP error\n")
                failed_batches += 1
            else:
                n = write_fasta_entries(cdna_raw, f_cdna, f_log)
                total_cdna_seqs += n

            # --- Download Protein ---
            prot_raw = fetch_batch(batch, "protein")
            if prot_raw is None:
                f_log.write(f"Batch {batch_num} (protein): FAILED - HTTP error\n")
                failed_batches += 1
            else:
                n = write_fasta_entries(prot_raw, f_prot, f_log)
                total_protein_seqs += n

            # Rate limiting
            time.sleep(RATE_LIMIT_SLEEP)

        f_log.write(f"\n--- Summary ---\n")
        f_log.write(f"Total cDNA sequences: {total_cdna_seqs}\n")
        f_log.write(f"Total protein sequences: {total_protein_seqs}\n")
        f_log.write(f"Failed batches: {failed_batches}\n")

    print(f"\n\nDownload complete!")
    print(f"  cDNA:     {OUTPUT_CDNA}  ({total_cdna_seqs} sequences)")
    print(f"  Protein:  {OUTPUT_PROTEIN}  ({total_protein_seqs} sequences)")
    print(f"  Log:      {OUTPUT_LOG}")


if __name__ == "__main__":
    main()
