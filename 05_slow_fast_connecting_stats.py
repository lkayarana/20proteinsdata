#!/usr/bin/env python3
import csv
import math
from collections import defaultdict, Counter
from typing import Dict, List, Tuple
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Data import CodonTable

# ---- CONFIG ----
ACCESSIONS = [
    "Q45614","P35164","Q99027","O32193","O31516","P94414","O34638","O05250","O06979","Q795K2",
    "O31433","P16497","Q08430","O31671","O31661","P70954","O34989","O32198","O34757","P42422"
]

CODON_TABLE_ID = 11      # bacterial code
MARGIN_AA = 5            # ±5 aa connecting margin
SLOW_Q = 0.25            # bottom 25%
FAST_Q = 0.25            # top 25%

CDS_FASTA = "cds_from_ncbi.fna"
FEATURES_TSV = "uniprot_domain_features.tsv"

KEEP_FEATURE_TYPES = {
    "Domain","Region","Topological domain","Transmembrane",
    "Coiled coil","Repeat","Motif"
}

# ---- HELPERS ----
def read_cds_fasta(path: str) -> Dict[str, str]:
    cds = {}
    for rec in SeqIO.parse(path, "fasta"):
        acc = rec.id.split("|")[0]
        cds[acc] = str(rec.seq).upper()
    return cds

def read_features(path: str) -> Dict[str, List[Tuple[int,int]]]:
    feats = defaultdict(list)
    with open(path, newline="", encoding="utf-8") as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            acc = row["accession"]
            ftype = row["feature_type"]
            if ftype not in KEEP_FEATURE_TYPES:
                continue
            s, e = row["start_aa"], row["end_aa"]
            if s and e:
                feats[acc].append((int(s), int(e)))
    for acc in feats:
        feats[acc].sort()
    return feats

def build_connecting_regions(domains: List[Tuple[int,int]], margin: int) -> List[Tuple[int,int]]:
    regs = []
    for i in range(len(domains)-1):
        cs = max(1, domains[i][1] + 1 - margin)
        ce = domains[i+1][0] - 1 + margin
        if cs <= ce:
            regs.append((cs, ce))
    return regs

def aa_in_regions(pos: int, regs: List[Tuple[int,int]]) -> bool:
    return any(s <= pos <= e for s,e in regs)

def split_codons(cds: str) -> List[str]:
    return [cds[i:i+3] for i in range(0, len(cds)-2, 3)]

def codon_to_aa_map():
    table = CodonTable.unambiguous_dna_by_id[CODON_TABLE_ID]
    return table.forward_table

CODON2AA = codon_to_aa_map()

# ---- MAIN ----
def main():
    cds_by_acc = read_cds_fasta(CDS_FASTA)
    feats_by_acc = read_features(FEATURES_TSV)

    # connecting regions
    connecting = {
        acc: build_connecting_regions(feats_by_acc.get(acc, []), MARGIN_AA)
        for acc in ACCESSIONS
    }

    # CAI-like codon weights (relative within dataset)
    counts = defaultdict(Counter)
    for cds in cds_by_acc.values():
        for codon in split_codons(cds):
            if codon in CODON2AA:
                counts[CODON2AA[codon]][codon] += 1

    w = {}
    for aa, cnts in counts.items():
        maxf = max(cnts.values())
        for codon, f in cnts.items():
            w[codon] = f / maxf

    # thresholds
    ws = sorted(w.values())
    slow_thr = ws[int(len(ws) * SLOW_Q)]
    fast_thr = ws[int(len(ws) * (1 - FAST_Q))]

    rows = []
    for acc, cds in cds_by_acc.items():
        regs = connecting.get(acc, [])
        codons = split_codons(cds)
        for i, codon in enumerate(codons, start=1):
            if codon not in CODON2AA:
                continue
            ww = w.get(codon, 0.0)
            region = "connecting" if aa_in_regions(i, regs) else "nonconnecting"
            rows.append([
                acc, i, codon, CODON2AA[codon], ww, region,
                int(ww <= slow_thr), int(ww >= fast_thr)
            ])

    with open("per_codon_table.tsv", "w", newline="", encoding="utf-8") as f:
        wr = csv.writer(f, delimiter="\t")
        wr.writerow(["accession","aa_pos","codon","aa","w","region","is_slow","is_fast"])
        wr.writerows(rows)

    a=b=c=d=0
    for r in rows:
        if r[5]=="connecting":
            a+=r[6]; b+=1-r[6]
        else:
            c+=r[6]; d+=1-r[6]

    OR = (a*d)/(b*c) if b*c>0 else float("inf")

    with open("stats_summary.txt","w") as f:
        f.write(f"Slow threshold w <= {slow_thr}\n")
        f.write(f"Connecting slow={a}, not_slow={b}\n")
        f.write(f"Nonconnecting slow={c}, not_slow={d}\n")
        f.write(f"Odds Ratio={OR}\n")

    print("DONE")
    print(f"a={a}, b={b}, c={c}, d={d}")
    print(f"OR={OR}")

if __name__ == "__main__":
    main()
