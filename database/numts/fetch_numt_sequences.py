#!/usr/bin/env python3
"""
Fetch all 766 Human Nuclear Mitochondrial DNA Segments (NUMTs) with +/-150bp
nuclear genomic flanks from UCSC DAS server and construct:
  1. database/numts/human_nuclear_NUMTs.fasta (pure nuclear NUMT sequences)
  2. reference/rCRS_plus_NUMTs.fasta (composite rCRS + Nuclear NUMT competitive reference)
  3. database/numts/numt_loci_catalog.json (structured lookup mapping nuclear & mtDNA coordinates)
"""
import json
import os
import re
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed

BASE_DIR = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
NUMT_DIR = os.path.join(BASE_DIR, "database", "numts")
REF_DIR = os.path.join(BASE_DIR, "reference")

NUC_TSV = os.path.join(NUMT_DIR, "ucsc_numtS_nuclear.tsv")
MITO_TSV = os.path.join(NUMT_DIR, "ucsc_numtS_mitochondrion.tsv")
ASM_TSV = os.path.join(NUMT_DIR, "ucsc_numtS_assembled.tsv")
GNOMAD_TSV = os.path.join(BASE_DIR, "gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv")
RCRS_FASTA = os.path.join(REF_DIR, "rCRS.fasta")


def load_numt_tables():
    # Load mitochondrial mappings keyed by HSA_NumtS ID
    mito_by_id = {}
    with open(MITO_TSV) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            _, chrom, start, end, numt_id, score, strand = parts[:7]
            mito_by_id.setdefault(numt_id, []).append({
                "mt_start": int(start),
                "mt_end": int(end),
                "identity_pct": round(int(score) / 10.0, 1),
                "mt_strand": strand,
            })

    nuc_entries = []
    with open(NUC_TSV) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 7:
                continue
            _, chrom, start, end, numt_id, score, strand = parts[:7]
            start_i, end_i = int(start), int(end)
            mt_matches = mito_by_id.get(numt_id, [])
            mt_start = mt_matches[0]["mt_start"] if mt_matches else 0
            mt_end = mt_matches[0]["mt_end"] if mt_matches else 0
            nuc_entries.append({
                "numt_id": numt_id,
                "nuc_chrom": chrom,
                "nuc_start": start_i,
                "nuc_end": end_i,
                "nuc_length": end_i - start_i,
                "identity_pct": round(int(score) / 10.0, 1),
                "strand": strand,
                "mt_start": mt_start,
                "mt_end": mt_end,
                "mt_length": (mt_end - mt_start) if mt_end > mt_start else (end_i - start_i),
            })

    return nuc_entries


def build_merged_clusters(nuc_entries, flank=150):
    sorted_entries = sorted(nuc_entries, key=lambda x: (x["nuc_chrom"], x["nuc_start"], x["nuc_end"]))
    clusters = []
    for e in sorted_entries:
        f_start = max(1, e["nuc_start"] - flank)
        f_end = e["nuc_end"] + flank
        if clusters and clusters[-1]["chrom"] == e["nuc_chrom"] and f_start <= clusters[-1]["end"] + 200:
            c = clusters[-1]
            c["end"] = max(c["end"], f_end)
            c["numt_ids"].append(e["numt_id"])
            c["max_identity"] = max(c["max_identity"], e["identity_pct"])
            if e["mt_start"] > 0:
                c["mt_ranges"].append(f"{e['mt_start']}-{e['mt_end']}")
        else:
            clusters.append({
                "chrom": e["nuc_chrom"],
                "start": f_start,
                "end": f_end,
                "core_start": e["nuc_start"],
                "core_end": e["nuc_end"],
                "numt_ids": [e["numt_id"]],
                "max_identity": e["identity_pct"],
                "mt_ranges": [f"{e['mt_start']}-{e['mt_end']}"] if e["mt_start"] > 0 else [],
            })
    return clusters


def fetch_batch(batch_idx, batch):
    segs = "&".join(f"segment={c['chrom']}:{c['start']},{c['end']}" for c in batch)
    url = f"https://genome.ucsc.edu/cgi-bin/das/hg19/dna?{segs}"
    req = urllib.request.Request(url, headers={"User-Agent": "MITOCLIN-NUMT-Builder/1.0"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        xml_data = resp.read()
    root = ET.fromstring(xml_data)
    results = []
    for seq_el, c in zip(root.findall("SEQUENCE"), batch):
        dna_el = seq_el.find("DNA")
        raw_seq = re.sub(r"\s+", "", dna_el.text or "").upper()
        results.append((c, raw_seq))
    return batch_idx, results


def main():
    print("[1/4] Loading UCSC / RHNumtS.2 nuclear and mitochondrial coordinates...")
    nuc_entries = load_numt_tables()
    clusters = build_merged_clusters(nuc_entries, flank=150)
    print(f"      Loaded {len(nuc_entries)} nuclear NUMT segments -> {len(clusters)} genomic decoy contigs.")

    batch_size = 45
    batches = [clusters[i:i + batch_size] for i in range(0, len(clusters), batch_size)]
    print(f"[2/4] Fetching {len(clusters)} nuclear NUMT genomic sequences in {len(batches)} parallel batches...")

    all_results = [None] * len(batches)
    with ThreadPoolExecutor(max_workers=6) as ex:
        futs = [ex.submit(fetch_batch, idx, b) for idx, b in enumerate(batches)]
        for fut in as_completed(futs):
            idx, res = fut.result()
            all_results[idx] = res
            print(f"      Completed batch {idx + 1}/{len(batches)} ({len(res)} contigs)")

    numt_fasta_path = os.path.join(NUMT_DIR, "human_nuclear_NUMTs.fasta")
    comp_fasta_path = os.path.join(REF_DIR, "rCRS_plus_NUMTs.fasta")

    print(f"[3/4] Writing {numt_fasta_path} and {comp_fasta_path}...")
    with open(RCRS_FASTA) as rf:
        rcrs_content = rf.read().strip() + "\n"

    total_bases = 0
    contig_count = 0
    with open(numt_fasta_path, "w") as nf, open(comp_fasta_path, "w") as cf:
        cf.write(rcrs_content)
        for batch_res in all_results:
            for c, seq in batch_res:
                if not seq:
                    continue
                contig_count += 1
                total_bases += len(seq)
                primary_id = c["numt_ids"][0]
                mt_str = ",".join(c["mt_ranges"][:3]) if c["mt_ranges"] else "mtDNA"
                contig_name = f"NUMT_{primary_id}_{c['chrom']}_{c['start']}_{c['end']}"
                header = (
                    f">{contig_name} "
                    f"locus={c['chrom']}:{c['core_start']}-{c['core_end']} "
                    f"mt_region={mt_str} max_identity={c['max_identity']}% "
                    f"segments={','.join(c['numt_ids'][:5])}\n"
                )
                wrapped = "\n".join(seq[i:i + 70] for i in range(0, len(seq), 70)) + "\n"
                nf.write(header + wrapped)
                cf.write(header + wrapped)

    print(f"      Wrote {contig_count} nuclear NUMT contigs ({total_bases:,} bp) + rCRS (16,569 bp).")

    # Save JSON catalog for fast variant & interval lookup
    catalog_path = os.path.join(NUMT_DIR, "numt_loci_catalog.json")
    with open(catalog_path, "w") as jf:
        json.dump({
            "summary": {
                "total_numt_segments": len(nuc_entries),
                "merged_nuclear_contigs": contig_count,
                "total_nuclear_numt_bp": total_bases,
            },
            "segments": nuc_entries,
        }, jf, indent=2)
    print(f"[4/4] Saved structured NUMT catalog to {catalog_path}")


if __name__ == "__main__":
    main()
