#!/usr/bin/env python3
"""
Download remaining missing RBP sequences (cDNA and protein) and APPEND to existing FASTA.
Small batches + conservative rate limiting for reliability.
"""
import requests, re, time, sys, os

ENSG_FILE = "./data/raw/missing_cdna_ensg.txt"  # IDs needing cDNA
EXISTING_CDNA = "./data/raw/human_RBP_cdna.fasta"
EXISTING_PROTEIN = "./data/raw/human_RBP_protein.fasta"
LOG_FILE = "./data/raw/RBP_download_log.txt"

BATCH_SIZE = 10        # small batches to avoid timeouts
MAX_RETRIES = 5
RETRY_DELAY = 10
RATE_LIMIT_SLEEP = 5.0  # 5s between batches

BASE_URL = "https://rest.ensembl.org"
HEADERS = {"Content-Type": "application/json", "Accept": "text/x-fasta"}


def read_ensg_ids(filepath):
    with open(filepath) as f:
        return [line.strip() for line in f if line.strip()]


def fetch_with_retry(ids, seq_type):
    url = f"{BASE_URL}/sequence/id?type={seq_type}"
    data = {"ids": ids}
    for attempt in range(1, MAX_RETRIES + 1):
        try:
            r = requests.post(url, headers=HEADERS, json=data, timeout=120)
            if r.status_code == 200:
                return r.text
            elif r.status_code == 429:
                print(f"  [429] retry {attempt}/{MAX_RETRIES}", flush=True)
                time.sleep(RETRY_DELAY * attempt * 2)
            elif r.status_code >= 500:
                print(f"  [{r.status_code}] retry {attempt}/{MAX_RETRIES}", flush=True)
                time.sleep(RETRY_DELAY * attempt)
            else:
                print(f"  [HTTP {r.status_code}] abort", flush=True)
                return None
        except requests.exceptions.Timeout:
            print(f"  [Timeout] retry {attempt}/{MAX_RETRIES}", flush=True)
            time.sleep(RETRY_DELAY * attempt)
        except requests.exceptions.ConnectionError:
            print(f"  [ConnErr] retry {attempt}/{MAX_RETRIES}", flush=True)
            time.sleep(RETRY_DELAY * attempt * 2)
        except Exception as e:
            print(f"  [{e}] abort", flush=True)
            return None
    return None


def append_fasta(fasta_text, f_out, seen_ids):
    """Write new entries, skipping duplicates."""
    if not fasta_text or not fasta_text.strip():
        return 0
    count = 0
    for entry in re.findall(r'(>.*?)(?=\n>|\Z)', fasta_text, re.DOTALL):
        lines = entry.strip().split('\n')
        header = lines[0]
        seq_id = header[1:].split()[0]
        if seq_id in seen_ids:
            continue
        seq = ''.join(lines[1:]).replace(' ', '').replace('\n', '')
        if len(seq) > 0:
            f_out.write(header + '\n')
            for i in range(0, len(seq), 60):
                f_out.write(seq[i:i+60] + '\n')
            seen_ids.add(seq_id)
            count += 1
    return count


def load_seen_ids(fasta_path):
    if not os.path.exists(fasta_path):
        return set()
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                ids.add(line[1:].split()[0])
    return ids


def main():
    print("=" * 60)
    print("Remaining RBP Downloader (Append, Conservative)")
    print("=" * 60)

    ensg_ids = read_ensg_ids(ENSG_FILE)
    print(f"\nENSG IDs to download: {len(ensg_ids)}")

    # Also read IDs needing protein
    prot_file = "./data/raw/missing_prot_ensg.txt"
    needs_protein = set(read_ensg_ids(prot_file))
    needs_cdna = set(ensg_ids)

    # Pre-load existing IDs to skip duplicates
    seen_cdna = load_seen_ids(EXISTING_CDNA)
    seen_prot = load_seen_ids(EXISTING_PROTEIN)
    print(f"Existing: cDNA={len(seen_cdna)}, protein={len(seen_prot)}")

    batches = [ensg_ids[i:i+BATCH_SIZE] for i in range(0, len(ensg_ids), BATCH_SIZE)]
    print(f"Batches: {len(batches)} (size {BATCH_SIZE}, {RATE_LIMIT_SLEEP}s gap)\n")

    total_cdna, total_prot = 0, 0
    failed = []

    with open(EXISTING_CDNA, 'a') as f_cdna, \
         open(EXISTING_PROTEIN, 'a') as f_prot, \
         open(LOG_FILE, 'a') as f_log:

        f_log.write(f"\n--- Remaining download {time.strftime('%Y-%m-%d %H:%M:%S')} ---\n")

        for idx, batch in enumerate(batches, 1):
            print(f"Batch {idx}/{len(batches)} ({batch[0]:>15}...) ", end="", flush=True)

            # cDNA
            if any(g in needs_cdna for g in batch):
                cdna_raw = fetch_with_retry(batch, "cdna")
                if cdna_raw is not None:
                    n = append_fasta(cdna_raw, f_cdna, seen_cdna)
                    total_cdna += n
                    print(f"cDNA+{n}", end=" ", flush=True)
                else:
                    print("cDNA:FAIL", end=" ", flush=True)
                    failed.append((idx, "cDNA", batch[0]))
            else:
                print("cDNA:skip", end=" ", flush=True)

            # Protein
            if any(g in needs_protein for g in batch):
                prot_raw = fetch_with_retry(batch, "protein")
                if prot_raw is not None:
                    n = append_fasta(prot_raw, f_prot, seen_prot)
                    total_prot += n
                    print(f"prot+{n}", end=" ", flush=True)
                else:
                    print("prot:FAIL", end=" ", flush=True)
                    failed.append((idx, "protein", batch[0]))
            else:
                print("prot:skip", end=" ", flush=True)

            print(flush=True)
            time.sleep(RATE_LIMIT_SLEEP)

        f_log.write(f"New: cDNA={total_cdna}, protein={total_prot}\n")
        if failed:
            f_log.write(f"Failed: {failed}\n")

    print(f"\n{'='*60}")
    print(f"Done! Added: {total_cdna} cDNA, {total_prot} protein sequences")
    if failed:
        print(f"⚠ {len(failed)} batch failures — re-run to retry")
    else:
        print(f"✓ All batches succeeded!")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
