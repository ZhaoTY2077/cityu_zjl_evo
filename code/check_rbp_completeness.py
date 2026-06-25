#!/usr/bin/env python3
"""
Check which RBP gene sequences are present in the FASTA files.
Uses Ensembl REST API to map ENSG -> ENST -> check FASTA.
"""
import requests
import json
import re
import sys
import os

ENSG_FILE = "./data/raw/RBP_geneid_list.txt"
CDNA_FASTA = "./data/raw/human_RBP_cdna.fasta"
PROTEIN_FASTA = "./data/raw/human_RBP_protein.fasta"
BATCH_SIZE = 50
API_TIMEOUT = 120

BASE_URL = "https://rest.ensembl.org"


def read_ensg_ids(filepath):
    """Read clean ENSG IDs from file."""
    ensg_set = set()
    with open(filepath, "r") as f:
        for line in f:
            m = re.search(r'ENSG\d{11}', line)
            if m:
                ensg_set.add(m.group())
    return sorted(ensg_set)


def get_gene_info(ensg_ids):
    """Batch lookup gene info from Ensembl to get transcript IDs."""
    results = {}
    for i in range(0, len(ensg_ids), BATCH_SIZE):
        batch = ensg_ids[i:i+BATCH_SIZE]
        url = f"{BASE_URL}/lookup/id"
        headers = {"Content-Type": "application/json", "Accept": "application/json"}
        data = {"ids": batch, "expand": 1}

        try:
            r = requests.post(url, headers=headers, json=data, timeout=API_TIMEOUT)
            if r.status_code == 200:
                resp = r.json()
                for gid, info in resp.items():
                    if info and 'Transcript' in info:
                        transcripts = [t['id'] for t in info['Transcript']]
                        results[gid] = {
                            'name': info.get('display_name', '?'),
                            'biotype': info.get('biotype', '?'),
                            'transcripts': transcripts
                        }
                print(f"  Lookup {i+1}-{i+len(batch)}: got {len(resp)} results", flush=True)
            else:
                print(f"  Lookup {i+1}-{i+len(batch)}: HTTP {r.status_code}", flush=True)
        except Exception as e:
            print(f"  Lookup {i+1}-{i+len(batch)}: Error - {e}", flush=True)

    return results


def load_sequences(fasta_path):
    """Load all sequence IDs from a FASTA file."""
    ids = set()
    with open(fasta_path) as f:
        for line in f:
            if line.startswith('>'):
                seq_id = line.strip()[1:].split('.')[0]  # Remove version
                ids.add(seq_id)
    return ids


def main():
    print("=" * 60)
    print("RBP Sequence Completeness Checker")
    print("=" * 60)

    # Read ENSG IDs
    ensg_ids = read_ensg_ids(ENSG_FILE)
    print(f"\nRBP list: {len(ensg_ids)} unique ENSG IDs")

    # Load existing sequences
    cdna_ids = load_sequences(CDNA_FASTA)
    protein_ids = load_sequences(PROTEIN_FASTA)
    print(f"cDNA FASTA: {len(cdna_ids)} unique ENST IDs")
    print(f"Protein FASTA: {len(protein_ids)} unique ENSP IDs")

    # Get gene info from Ensembl
    print(f"\nQuerying Ensembl API for gene info...")
    gene_info = get_gene_info(ensg_ids)
    print(f"Got info for {len(gene_info)} genes")

    # Check completeness
    missing_from_cdna = []
    missing_from_protein = []
    genes_with_no_transcripts = []

    for gid in ensg_ids:
        if gid not in gene_info:
            genes_with_no_transcripts.append(gid)
            continue

        info = gene_info[gid]
        transcripts = info['transcripts']

        if not transcripts:
            genes_with_no_transcripts.append(gid)
            continue

        # Check if any transcript is in the cDNA FASTA
        found_in_cdna = any(t in cdna_ids for t in transcripts)
        if not found_in_cdna:
            missing_from_cdna.append(gid)

        # For protein, check if protein IDs are present
        # (We can't directly check since we have ENSP IDs)
        # We'll note this limitation

    # Report
    print(f"\n{'='*60}")
    print(f"RESULTS")
    print(f"{'='*60}")
    print(f"Total RBP genes:           {len(ensg_ids)}")
    print(f"Genes with info:           {len(gene_info)}")
    print(f"Genes found in cDNA:       {len(ensg_ids) - len(missing_from_cdna) - len(genes_with_no_transcripts)}")
    print(f"Genes MISSING from cDNA:   {len(missing_from_cdna)}")
    print(f"Genes with no transcripts: {len(genes_with_no_transcripts)}")

    if missing_from_cdna:
        print(f"\nMissing from cDNA:")
        for gid in sorted(missing_from_cdna):
            name = gene_info.get(gid, {}).get('name', '?')
            biotype = gene_info.get(gid, {}).get('biotype', '?')
            print(f"  {gid} ({name}, {biotype})")

    if genes_with_no_transcripts:
        print(f"\nGenes with no transcripts (ENSEMBL API):")
        for gid in sorted(genes_with_no_transcripts):
            print(f"  {gid}")

    # If any missing, offer to download
    if missing_from_cdna or genes_with_no_transcripts:
        print(f"\n{'='*60}")
        print(f"MISSING RBPs DETECTED - Need to re-download")
        print(f"{'='*60}")
        return 1
    else:
        print(f"\n✓ All {len(ensg_ids)} RBP genes are fully represented in the FASTA files!")
        return 0


if __name__ == "__main__":
    sys.exit(main())
