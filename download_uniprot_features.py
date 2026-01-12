#!/usr/bin/env python3
import time
import csv
import requests

# --- 20 PROTEIN ---
ACCESSIONS = [
    "Q45614","P35164","Q99027","O32193","O31516","P94414","O34638","O05250","O06979","Q795K2",
    "O31433","P16497","Q08430","O31671","O31661","P70954","O34989","O32198","O34757","P42422"
]

UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{}"

KEEP_TYPES = {
    "Domain",
    "Region",
    "Repeat",
    "Motif",
    "Coiled coil",
    "Topological domain",
    "Transmembrane"
}

session = requests.Session()
session.headers.update({
    "Accept": "application/json",
    "User-Agent": "slow-fast-codon-project/1.0"
})

def fetch_uniprot(acc):
    r = session.get(UNIPROT_URL.format(acc), timeout=30)
    r.raise_for_status()
    return r.json()

def extract_features(acc, data):
    feats = []
    for f in data.get("features", []):
        ftype = f.get("type")
        if ftype not in KEEP_TYPES:
            continue

        loc = f.get("location", {})
        start = (loc.get("start") or {}).get("value")
        end = (loc.get("end") or {}).get("value")
        if not start or not end:
            continue

        desc = f.get("description") or f.get("featureId") or ""
        feats.append({
            "accession": acc,
            "feature_type": ftype,
            "start_aa": int(start),
            "end_aa": int(end),
            "description": desc
        })
    return feats

def main():
    all_features = []

    for acc in ACCESSIONS:
        try:
            data = fetch_uniprot(acc)
            feats = extract_features(acc, data)
            all_features.extend(feats)
            print(f"[OK] {acc}: {len(feats)} features")
            time.sleep(0.2)
        except Exception as e:
            print(f"[ERR] {acc}: {e}")

    with open("uniprot_domain_features.tsv", "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=["accession","feature_type","start_aa","end_aa","description"],
            delimiter="\t"
        )
        writer.writeheader()
        for row in all_features:
            writer.writerow(row)

    print("Wrote: uniprot_domain_features.tsv")

if __name__ == "__main__":
    main()
