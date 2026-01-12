#!/usr/bin/env python3
import re
import time
import requests
import xml.etree.ElementTree as ET
from typing import Dict, Any, List, Optional, Tuple
from Bio.Seq import Seq

# --- SENİN 20 ACCESSION ---
ACCESSIONS = [
    "Q45614","P35164","Q99027","O32193","O31516","P94414","O34638","O05250","O06979","Q795K2",
    "O31433","P16497","Q08430","O31671","O31661","P70954","O34989","O32198","O34757","P42422"
]

# --- AYARLAR ---
NCBI_EMAIL = "rana.kaya@live.acibadem.edu.tr"  # <- burayı değiştirmen iyi olur
NCBI_API_KEY = None  # varsa yaz (rate limit artar)
SLEEP_S = 0.34       # NCBI'ye nazik istek aralığı

UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{}"
NCBI_EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"

session = requests.Session()
session.headers.update({
    "Accept": "application/json",
    "User-Agent": "slow-fast-codon-project/1.0"
})

def fetch_uniprot_json(acc: str) -> Dict[str, Any]:
    r = session.get(UNIPROT_ENTRY.format(acc), headers={"Accept": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()

def pick_refseq_protein_xref(uniprot_entry: Dict[str, Any]) -> Optional[str]:
    """
    UniProt JSON'da cross-reference yapısı zamanla değişebilir.
    En pratik: crossReferences içinde RefSeq id'lerini ara.
    """
    candidates = []

    for key in ["uniProtKBCrossReferences", "crossReferences", "dbReferences"]:
        xrefs = uniprot_entry.get(key)
        if isinstance(xrefs, list):
            for x in xrefs:
                db = (x.get("database") or x.get("db") or "").lower()
                if db == "refseq":
                    xid = x.get("id") or x.get("primaryId")
                    if xid:
                        candidates.append(xid)

    # tercih sırası: NP_/YP_/WP_
    pref = ["NP_", "YP_", "WP_"]
    for p in pref:
        for c in candidates:
            if c.startswith(p):
                return c
    return candidates[0] if candidates else None

def fetch_ncbi_protein_genpept(protein_id: str) -> str:
    """
    protein_id hem NP_... gibi accession olabilir, hem de NCBI UID olabilir.
    """
    params = {
        "db": "protein",
        "id": protein_id,
        "rettype": "gp",
        "retmode": "text",
        "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = requests.get(NCBI_EFETCH, params=params, timeout=60)
    r.raise_for_status()
    return r.text

def esearch_protein_ids(query: str, retmax: int = 20) -> List[str]:
    """
    coded_by olmayan durumlarda NCBI'de alternatif protein kayıtlarını aramak için.
    """
    params = {
        "db": "protein",
        "term": query,
        "retmode": "xml",
        "retmax": str(retmax),
        "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = requests.get(NCBI_ESEARCH, params=params, timeout=60)
    r.raise_for_status()
    root = ET.fromstring(r.text)
    return [e.text for e in root.findall(".//IdList/Id") if e.text]

def parse_coded_by(genpept_text: str) -> Tuple[str, List[Tuple[int, int]], int]:
    """
    GenPept CDS feature altında coded_by satırı genelde şöyle olur:
      coded_by="NC_000964.3:12345..13000"
      coded_by="complement(NC_000964.3:12345..13000)"
      coded_by="join(NC_000964.3:123..200,NC_000964.3:300..400)"
    """
    m = re.search(r'coded_by="([^"]+)"', genpept_text, flags=re.DOTALL)
    if not m:
        raise ValueError("coded_by not found in GenPept record.")

    coded = m.group(1).strip()

    strand = +1
    if coded.startswith("complement(") and coded.endswith(")"):
        strand = -1
        coded = coded[len("complement("):-1].strip()

    if coded.startswith("join(") and coded.endswith(")"):
        coded = coded[len("join("):-1].strip()
        parts = [p.strip() for p in coded.split(",")]
    else:
        parts = [coded]

    nuc_acc = None
    segs: List[Tuple[int, int]] = []

    for p in parts:
        mm = re.match(r"^([A-Za-z0-9_\.]+):<?(\d+)\.\.>?(\d+)$", p)
        if not mm:
            raise ValueError(f"Unrecognized coded_by segment: {p}")
        a, s, e = mm.group(1), int(mm.group(2)), int(mm.group(3))
        nuc_acc = nuc_acc or a
        if a != nuc_acc:
            raise ValueError("coded_by join across different nucleotide accessions not supported in this script.")
        segs.append((s, e))

    return nuc_acc, segs, strand

def fetch_nuccore_segment_fasta(nuc_acc: str, start: int, end: int, strand: int) -> str:
    params = {
        "db": "nuccore",
        "id": nuc_acc,
        "rettype": "fasta",
        "retmode": "text",
        "seq_start": str(start),
        "seq_stop": str(end),
        "strand": "2" if strand == -1 else "1",
        "email": NCBI_EMAIL,
    }
    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY
    r = requests.get(NCBI_EFETCH, params=params, timeout=60)
    r.raise_for_status()
    return r.text

def fasta_to_seq(fasta_text: str) -> str:
    lines = [ln.strip() for ln in fasta_text.splitlines() if ln.strip()]
    seq = "".join([ln for ln in lines if not ln.startswith(">")])
    return seq.upper()

def translate_cds_bacteria(cds_seq: str) -> str:
    """
    Bacterial translation:
    - table=11
    - force first aa to M (alt start codons)
    - stop handling
    """
    prot_full = str(Seq(cds_seq).translate(table=11))
    prot_full = prot_full.rstrip("*")
    if prot_full:
        prot_full = "M" + prot_full[1:]
    # first stop'a kadar (bazı kayıtlarda internal stop yoktur ama güvenli)
    return prot_full.split("*")[0] if prot_full else ""

def get_coded_by_with_fallback(acc: str, refseq_prot: str) -> Tuple[str, List[Tuple[int, int]], int, str]:
    """
    Önce RefSeq protein accession'dan coded_by dene.
    Bulamazsa NCBI protein'de arama yapıp coded_by içeren ilk kaydı kullan.
    """
    gp = fetch_ncbi_protein_genpept(refseq_prot)
    try:
        nuc_acc, segs, strand = parse_coded_by(gp)
        return nuc_acc, segs, strand, refseq_prot
    except ValueError as e:
        if "coded_by not found" not in str(e):
            raise

    # fallback search:
    # Organism kısıtlaması çok işe yarıyor
    org = "Bacillus subtilis subsp. subtilis str. 168"
    query = f'"{org}"[Organism] AND {acc}[All Fields]'
    ids = esearch_protein_ids(query, retmax=30)

    for uid in ids:
        gp2 = fetch_ncbi_protein_genpept(uid)  # UID ile de efetch olur
        try:
            nuc_acc, segs, strand = parse_coded_by(gp2)
            # UID yerine mümkünse ACCESSION'ı yazmak için VERSION satırını çekelim:
            mver = re.search(r"^VERSION\s+(\S+)", gp2, flags=re.MULTILINE)
            chosen = mver.group(1) if mver else uid
            return nuc_acc, segs, strand, chosen
        except Exception:
            continue

    raise ValueError("coded_by not found even after NCBI fallback search.")

def main():
    out_fna = open("cds_from_ncbi.fna", "w", encoding="utf-8")
    out_log = open("cds_fetch_log.tsv", "w", encoding="utf-8")
    out_log.write("uniprot_acc\trefseq_protein_used\tnuc_acc\tsegments\tstrand\tcds_len_nt\ttranslation_match\n")

    for acc in ACCESSIONS:
        try:
            up = fetch_uniprot_json(acc)
            uniprot_prot = (up.get("sequence") or {}).get("value", "")
            uniprot_prot = uniprot_prot.replace("\n", "").strip()

            refseq_prot = pick_refseq_protein_xref(up)
            if not refseq_prot:
                raise ValueError("No RefSeq xref found in UniProt entry.")

            nuc_acc, segs, strand, protein_used = get_coded_by_with_fallback(acc, refseq_prot)

            cds_seq = ""
            for (s, e) in segs:
                fa = fetch_nuccore_segment_fasta(nuc_acc, s, e, strand)
                cds_seq += fasta_to_seq(fa)
                time.sleep(SLEEP_S)

            translated = translate_cds_bacteria(cds_seq)
            match = (translated == uniprot_prot)

            header = f">{acc}|ProteinUsed:{protein_used}|Nuccore:{nuc_acc}|strand:{strand}|segs:{segs}"
            out_fna.write(header + "\n")
            for i in range(0, len(cds_seq), 70):
                out_fna.write(cds_seq[i:i+70] + "\n")

            out_log.write(f"{acc}\t{protein_used}\t{nuc_acc}\t{segs}\t{strand}\t{len(cds_seq)}\t{match}\n")
            print(f"[OK] {acc} -> {protein_used} -> {nuc_acc} len={len(cds_seq)} match={match}")

            time.sleep(SLEEP_S)

        except Exception as e:
            out_log.write(f"{acc}\tNA\tNA\tNA\tNA\tNA\tERROR:{e}\n")
            print(f"[ERR] {acc}: {e}")

    out_fna.close()
    out_log.close()
    print("Wrote: cds_from_ncbi.fna, cds_fetch_log.tsv")

if __name__ == "__main__":
    main()
