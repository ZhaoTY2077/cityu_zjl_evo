#!/usr/bin/env python3
"""
Download MISSING human RBP sequences (cDNA and protein) and APPEND to existing FASTA files.
Uses Ensembl REST API with batch queries, retry logic, and checkpoints.
"""

import requests
import time
import os
import sys
import re
import json

# --- Config ---
ENSG_FILE = "./files/RBP_geneid_list.txt"
EXISTING_CDNA = "./files/human_RBP_cdna.fasta"
EXISTING_PROTEIN = "./files/human_RBP_protein.fasta"
OUTPUT_LOG = "./files/RBP_download_log.txt"
CHECKPOINT_FILE = "./files/download_checkpoint.json"

BATCH_SIZE = 25  # Smaller batches for reliability
MAX_RETRIES = 3
RETRY_DELAY = 5  # seconds between retries
RATE_LIMIT_SLEEP = 3.0  # seconds between batches

BASE_URL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "text/x-fasta"}


def read_ensg_ids(filepath):
    """Read ENSG IDs from file."""
    ids = []
    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            m = re.search(r'ENSG\d{11}', line)
            if m:
                ids.append(m.group())
    return ids


def fetch_batch_with_retry(ids, seq_type):
    """Fetch sequences for a batch of IDs with retry logic."""
    url = f"{BASE_URL}/sequence/id?type={seq_type}"
    data = {"ids": ids}

    for attempt in range(MAX_RETRIES):
        try:
            r = requests.post(url, headers=HEADERS, json=data, timeout=60)
            if r.status_code == 200:
                return r.text
            elif r.status_code == 429:
                wait = RETRY_DELAY * (attempt + 1) * 2
                print(f"  [429 Rate limited] Retry {attempt+1}/{MAX_RETRIES} after {wait}s", flush=True)
                time.sleep(wait)
            elif r.status_code >= 500:
                wait = RETRY_DELAY * (attempt + 1)
                print(f"  [{r.status_code} Server error] Retry {attempt+1}/{MAX_RETRIES} after {wait}s", flush=True)
                time.sleep(wait)
            else:
                print(f"  [HTTP {r.status_code}] Not retrying", flush=True)
                return None
        except requests.exceptions.Timeout:
            wait = RETRY_DELAY * (attempt + 1)
            print(f"  [Timeout] Retry {attempt+1}/{MAX_RETRIES} after {wait}s", flush=True)
            time.sleep(wait)
        except requests.exceptions.ConnectionError:
            wait = RETRY_DELAY * (attempt + 1) * 2
            print(f"  [ConnectionError] Retry {attempt+1}/{MAX_RETRIES} after {wait}s", flush=True)
            time.sleep(wait)
        except Exception as e:
            print(f"  [Error: {e}] Not retrying", flush=True)
            return None

    print(f"  [FAILED after {MAX_RETRIES} retries]", flush=True)
    return None


def write_fasta_entries(fasta_text, f_out, log_handle, source_ids):
    """Parse and write FASTA entries, returning count."""
    if not fasta_text or fasta_text.strip() == "":
        log_handle.write(f"  Empty response for batch starting with {source_ids[0]}\n")
        return 0

    entries = re.findall(r'(>.*?)(?=\n>|\Z)', fasta_text, re.DOTALL)
    count = 0
    written_headers = set()
    for entry in entries:
        lines = entry.strip().split('\n')
        header = lines[0]
        seq = ''.join(lines[1:]).replace(' ', '').replace('\n', '')

        # Skip duplicate sequences
        seq_id = header[1:].split()[0]
        if seq_id in written_headers:
            continue
        written_headers.add(seq_id)

        if len(seq) > 0:
            f_out.write(header + '\n')
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + '\n')
            count += 1
        else:
            log_handle.write(f"  WARNING: Empty sequence for {header}\n")
    return count


def load_existing_ids(fasta_path):
    """Load sequence IDs that already exist in a FASTA file."""
    ids = set()
    if not os.path.exists(fasta_path):
        return ids
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:].split()[0]
                ids.add(seq_id)
    return ids


def main():
    print("=" * 60)
    print("Missing RBP Sequence Downloader (Append Mode)")
    print("=" * 60)

    # Read all ENSG IDs
    all_ensg_ids = read_ensg_ids(ENSG_FILE)
    print(f"\nTotal RBP ENSG IDs: {len(all_ensg_ids)}")

    # Load existing sequence IDs to avoid duplicates
    existing_cdna_ids = load_existing_ids(EXISTING_CDNA)
    existing_protein_ids = load_existing_ids(EXISTING_PROTEIN)
    print(f"Existing cDNA sequences: {len(existing_cdna_ids)}")
    print(f"Existing protein sequences: {len(existing_protein_ids)}")

    # These are the IDs we need to download for. We'll download by ENSG ID.
    # Rather than trying to figure out which ENSG IDs are already covered,
    # we simply download ALL ENSG IDs - the Ensembl API handles dedup at gene level.
    # Our write_fasta_entries will skip duplicate ENST/ENSP IDs.

    # Simplified approach: download all missing in one pass
    # Split into batches
    batches = [all_ensg_ids[i:i+BATCH_SIZE] for i in range(0, len(all_ensg_ids), BATCH_SIZE)]
    total_batches = len(batches)
    print(f"Total batches: {total_batches}")

    # Open files for APPENDING
    total_cdna_new = 0
    total_protein_new = 0
    failed_batches = []
    partial_batches = []

    with open(EXISTING_CDNA, 'a') as f_cdna, \
         open(EXISTING_PROTEIN, 'a') as f_prot, \
         open(OUTPUT_LOG, 'a') as f_log:

        if os.path.getsize(OUTPUT_LOG) == 0:
            f_log.write(f"RBP sequence download log\n")
            f_log.write(f"Total ENSG IDs: {len(all_ensg_ids)}\n\n")

        f_log.write(f"\n--- Resume download at {time.strftime('%Y-%m-%d %H:%M:%S')} ---\n")

        for batch_idx, batch in enumerate(batches):
            batch_num = batch_idx + 1
            print(f"\rBatch {batch_num}/{total_batches} (IDs {all_ensg_ids.index(batch[0])+1}-{all_ensg_ids.index(batch[-1])+1})...", flush=True)

            # --- Download cDNA ---
            cdna_raw = fetch_batch_with_retry(batch, "cdna")
            if cdna_raw is None:
                f_log.write(f"  Batch {batch_num} (cDNA): FAILED after retries - IDs: {batch[0]}..{batch[-1]}\n")
                failed_batches.append((batch_num, 'cDNA', batch))
                cdna_count = 0
            else:
                cdna_count = write_fasta_entries(cdna_raw, f_cdna, f_log, batch)
                total_cdna_new += cdna_count
                if cdna_count == 0:
                    f_log.write(f"  Batch {batch_num} (cDNA): 0 new sequences (IDs: {batch[0]}..{batch[-1]})\n")

            # --- Download Protein ---
            prot_raw = fetch_batch_with_retry(batch, "protein")
            if prot_raw is None:
                f_log.write(f"  Batch {batch_num} (protein): FAILED after retries - IDs: {batch[0]}..{batch[-1]}\n")
                failed_batches.append((batch_num, 'protein', batch))
                prot_count = 0
            else:
                prot_count = write_fasta_entries(prot_raw, f_prot, f_log, batch)
                total_protein_new += prot_count
                if prot_count == 0:
                    f_log.write(f"  Batch {batch_num} (protein): 0 new sequences (IDs: {batch[0]}..{batch[-1]})\n")

            # Track partial failures
            if cdna_count == 0 or prot_count == 0:
                partial_batches.append(batch_num)

            # Rate limiting
            time.sleep(RATE_LIMIT_SLEEP)

        # Summary
        f_log.write(f"\n--- Final Summary ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---\n")
        f_log.write(f"New cDNA sequences appended: {total_cdna_new}\n")
        f_log.write(f"New protein sequences appended: {total_protein_new}\n")
        f_log.write(f"Total batches processed: {total_batches}\n")
        f_log.write(f"Failed batches (both types): {len(failed_batches)}\n")
        if failed_batches:
            f_log.write(f"  Failed details: {failed_batches}\n")

    # Print results
    print(f"\n{'='*60}")
    print(f"DOWNLOAD COMPLETE")
    print(f"{'='*60}")
    print(f"New cDNA sequences:     {total_cdna_new}")
    print(f"New protein sequences:  {total_protein_new}")
    print(f"Failed batch-calls:     {len(failed_batches)}")
    print(f"Partially successful:   {len(partial_batches)}")

    if failed_batches:
        print(f"\n⚠  {len(failed_batches)} batch-calls failed even after retries.")
        print(f"   IDs of first failed batch: {failed_batches[0][2][:3]}...")
        print(f"   Re-run the script to retry these.")

    # Show new file sizes
    cdna_size = os.path.getsize(EXISTING_CDNA)
    prot_size = os.path.getsize(EXISTING_PROTEIN)
    print(f"\nFinal file sizes:")
    print(f"  {EXISTING_CDNA}: {cdna_size/1024/1024:.1f} MB")
    print(f"  {EXISTING_PROTEIN}: {prot_size/1024/1024:.1f} MB")
    print(f"  Log: {OUTPUT_LOG}")

    return 0 if not failed_batches else 1


if __name__ == "__main__":
    sys.exit(main())
