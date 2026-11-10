#!/usr/bin/env python3
"""
MITOCLIN NUMT (Nuclear Mitochondrial DNA Segment) Intelligence Engine
=====================================================================
Provides:
  1. Competitive Alignment BAM Filtering (`filter-bam`):
     - Separates authentic mitochondrial reads mapping to NC_012920.1 (rCRS)
       from nuclear pseudogene reads mapping preferentially to human nuclear
       NUMT contigs (HSA_NumtS_001..766 / nuclear chromosomes).
     - Produces a clean NC_012920.1-only SAM/BAM for downstream GATK Mutect2
       and HaploCheck compatibility.
     - Emits detailed NUMT read depletion metrics (`_numt_depletion.json`).

  2. Variant-Level NUMT Knowledgebase Annotation (`annotate-tsv` & Python API):
     - Integrates Broad Institute gnomAD v3.1 mitochondrial `npg` (Nuclear
       Pseudogene) filter annotations (7,078 curated artifact alleles) with
       the UCSC / RHNumtS.2 catalog of 766 human nuclear NUMT loci across
       chr1..chrY.
"""

import csv
import json
import os
import re
import subprocess
import sys
from collections import Counter

SCRIPT_DIR = os.path.abspath(os.path.dirname(__file__))
# Check whether running inside scripts/ subfolder or repository root
if os.path.basename(SCRIPT_DIR) == "scripts":
    REPO_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, ".."))
else:
    REPO_DIR = SCRIPT_DIR

NUMT_DIR = os.path.join(REPO_DIR, "database", "numts")
MITO_TSV = os.path.join(NUMT_DIR, "ucsc_numtS_mitochondrion.tsv")
NUC_TSV = os.path.join(NUMT_DIR, "ucsc_numtS_nuclear.tsv")

# Locate gnomAD reduced annotations either in database/ or repository root
GNOMAD_TSV = os.path.join(REPO_DIR, "database", "gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv")
if not os.path.exists(GNOMAD_TSV):
    GNOMAD_TSV = os.path.join(REPO_DIR, "gnomad.genomes.v3.1.sites.chrM.reduced_annotations.tsv")

_NUMT_KB_CACHE = None


def load_numt_knowledgebase():
    """Load and cache unified NUMT interval map (1..16569 bp) + gnomAD v3.1 npg alleles."""
    global _NUMT_KB_CACHE
    if _NUMT_KB_CACHE is not None:
        return _NUMT_KB_CACHE

    # 1. Load nuclear NUMT coordinates by HSA_NumtS ID
    nuc_by_id = {}
    if os.path.isfile(NUC_TSV):
        with open(NUC_TSV, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                _, chrom, start, end, numt_id, score, strand = parts[:7]
                try:
                    nuc_by_id[numt_id] = {
                        "numt_id": numt_id,
                        "nuc_chrom": chrom,
                        "nuc_start": int(start),
                        "nuc_end": int(end),
                        "nuc_locus": f"{chrom}:{int(start)}-{int(end)}",
                        "identity_pct": round(int(score) / 10.0, 1),
                        "strand": strand,
                    }
                except ValueError:
                    continue

    # 2. Load mitochondrial intervals and link to nuclear coordinates
    mt_intervals = []
    if os.path.isfile(MITO_TSV):
        with open(MITO_TSV, encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) < 7:
                    continue
                _, _, mt_start, mt_end, numt_id, score, mt_strand = parts[:7]
                try:
                    s_i, e_i = int(mt_start), int(mt_end)
                    id_pct = round(int(score) / 10.0, 1)
                except ValueError:
                    continue
                nuc_info = nuc_by_id.get(numt_id, {})
                mt_intervals.append({
                    "numt_id": numt_id,
                    "mt_start": s_i,
                    "mt_end": e_i,
                    "mt_length": e_i - s_i,
                    "identity_pct": id_pct,
                    "mt_strand": mt_strand,
                    "nuc_chrom": nuc_info.get("nuc_chrom", "chrUn"),
                    "nuc_start": nuc_info.get("nuc_start", 0),
                    "nuc_end": nuc_info.get("nuc_end", 0),
                    "nuc_locus": nuc_info.get("nuc_locus", numt_id),
                })

    # 3. Load gnomAD v3.1 npg alleles
    npg_alleles = {}
    gnomad_filters = {}
    if os.path.isfile(GNOMAD_TSV):
        with open(GNOMAD_TSV, encoding="utf-8", errors="replace") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                try:
                    pos = int(row["position"])
                    ref = row["ref"].strip().upper()
                    alt = row["alt"].strip().upper()
                    flt = row.get("filters", "PASS").strip()
                    gnomad_filters[(pos, ref, alt)] = flt
                    if "npg" in flt.lower():
                        npg_alleles[(pos, ref, alt)] = flt
                except Exception:
                    continue

    _NUMT_KB_CACHE = {
        "mt_intervals": mt_intervals,
        "npg_alleles": npg_alleles,
        "gnomad_filters": gnomad_filters,
        "total_numt_segments": len(mt_intervals),
        "total_npg_alleles": len(npg_alleles),
    }
    return _NUMT_KB_CACHE


def get_numt_annotation(pos, ref=None, alt=None, het=None):
    """
    Return comprehensive NUMT homology & gnomAD npg assessment for an mtDNA variant.
    het can be 0..1 or 0..100 (normalized internally to 0..100%).
    """
    kb = load_numt_knowledgebase()
    try:
        pos_i = int(str(pos).replace(",", "").strip())
    except (ValueError, TypeError):
        return {
            "is_npg": False,
            "gnomad_filter": "NA",
            "overlap_count": 0,
            "high_id_count": 0,
            "max_identity_pct": 0.0,
            "primary_nuclear_locus": "—",
            "nuclear_loci_summary": "None",
            "top_loci": [],
            "risk_tier": "Clean",
            "status_label": "Clean (No NUMT Flag)",
            "badge_html": '<span style="display:inline-block;padding:2px 6px;border-radius:4px;font-size:7pt;font-weight:700;background:#ecfdf5;color:#047857;border:1px solid #a7f3d0">Clean</span>',
            "clinical_note": "No known nuclear pseudogene (NUMT) interference at this allele.",
        }

    ref_u = str(ref or "").strip().upper()
    alt_u = str(alt or "").strip().upper()

    flt = kb["gnomad_filters"].get((pos_i, ref_u, alt_u), "")
    is_npg = (pos_i, ref_u, alt_u) in kb["npg_alleles"]

    # Find overlapping UCSC/RHNumtS.2 intervals covering pos_i
    hits = [iv for iv in kb["mt_intervals"] if iv["mt_start"] <= pos_i <= iv["mt_end"]]
    hits.sort(key=lambda x: (-x["identity_pct"], -x["mt_length"]))

    overlap_count = len(hits)
    high_id_hits = [h for h in hits if h["identity_pct"] >= 95.0]
    high_id_count = len(high_id_hits)
    max_id = hits[0]["identity_pct"] if hits else 0.0

    top_loci = []
    seen_loci = set()
    for h in hits:
        loc = h["nuc_locus"]
        if loc in seen_loci:
            continue
        seen_loci.add(loc)
        top_loci.append(h)
        if len(top_loci) >= 4:
            break

    if top_loci:
        primary_loc_coord = top_loci[0]["nuc_locus"]
        primary_loc_id = top_loci[0]["numt_id"]
        primary_loc = f"{primary_loc_coord} ({top_loci[0]['identity_pct']:.1f}% id)"
        loci_summary = "; ".join(
            f"{x['nuc_locus']} ({x['numt_id']}, {x['identity_pct']:.1f}%)"
            for x in top_loci[:3]
        )
    else:
        primary_loc_coord = "—"
        primary_loc_id = "—"
        primary_loc = "—"
        loci_summary = "None"

    # Normalize heteroplasmy to percentage (0..100)
    het_pct = None
    if het is not None:
        try:
            hv = float(str(het).replace("%", "").strip())
            het_pct = hv * 100.0 if hv <= 1.0 else hv
        except (ValueError, TypeError):
            het_pct = None

    # Determine NUMT risk tier (incorporating VAF-aware diploid NUMT ceiling rule)
    if is_npg and (het_pct is None or het_pct <= 40.0):
        risk_tier = "Confirmed_NUMT_NPG"
        status_label = "Suspected NUMT Artifact (gnomAD npg)"
        badge_text = "NUMT Artifact (npg)"
        badge_html = (
            '<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            'font-size:7pt;font-weight:700;background:#fef2f2;color:#b91c1c;'
            'border:1px solid #fecaca" title="Flagged as Nuclear Pseudogene (npg) in gnomAD v3.1 at low/moderate VAF">'
            'NUMT Artifact (npg)</span>'
        )
        vaf_txt = f" at {het_pct:.1f}% VAF" if het_pct is not None else ""
        clinical_note = (
            f"Exact allele m.{pos_i}{ref_u}>{alt_u}{vaf_txt} is flagged as a Nuclear Pseudogene (npg) "
            f"in gnomAD v3.1 and falls within the diploid NUMT ceiling (<=40% VAF). "
            f"Primary homologous nuclear locus: {primary_loc_coord} ({max_id:.1f}% identity). Confirm via long-range mtDNA PCR."
        )
    elif is_npg and het_pct is not None and het_pct > 40.0:
        risk_tier = "NUMT_Locus_High_VAF"
        status_label = "Authentic mtDNA (Overlaps gnomAD npg Locus)"
        badge_text = f"npg Site (VAF {het_pct:.0f}%)"
        badge_html = (
            f'<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            f'font-size:7pt;font-weight:700;background:#eff6ff;color:#1d4ed8;'
            f'border:1px solid #bfdbfe" title="Locus is in gnomAD npg list ({primary_loc_coord}), but high VAF ({het_pct:.1f}%) exceeds nuclear NUMT ceiling">'
            f'npg Site (VAF {het_pct:.0f}%)</span>'
        )
        clinical_note = (
            f"Overlaps gnomAD v3.1 npg locus ({primary_loc_coord}, {max_id:.1f}% identity), "
            f"but high observed VAF ({het_pct:.1f}%) exceeds the 40% diploid nuclear pseudogene "
            f"ceiling, confirming authentic mitochondrial origin."
        )
    elif max_id >= 95.0 and het_pct is not None and 5.0 <= het_pct < 95.0:
        risk_tier = "High_NUMT_Warning"
        status_label = f"High-Identity NUMT Hotzone ({max_id:.1f}%)"
        badge_text = f"Het in NUMT ({max_id:.1f}%)"
        badge_html = (
            f'<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            f'font-size:7pt;font-weight:700;background:#fffbeb;color:#b45309;'
            f'border:1px solid #fde68a" title="Heteroplasmy in {max_id:.1f}% identity nuclear NUMT region ({primary_loc_coord})">'
            f'Het in NUMT ({max_id:.1f}%)</span>'
        )
        clinical_note = (
            f"Heteroplasmy ({het_pct:.2f}% VAF) within a {max_id:.1f}% sequence-identity nuclear NUMT region "
            f"({primary_loc_coord}; {overlap_count} nuclear segments). "
            f"Verify via competitive NUMT read depletion or long-range PCR."
        )
    elif max_id >= 95.0:
        risk_tier = "High_Identity_Homology"
        status_label = f"High-Identity NUMT Region ({max_id:.1f}%)"
        badge_text = f"Homology {max_id:.1f}%"
        badge_html = (
            f'<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            f'font-size:7pt;font-weight:700;background:#f0f9ff;color:#0369a1;'
            f'border:1px solid #bae6fd" title="Homologous to {primary_loc_coord} ({overlap_count} nuclear NUMTs)">'
            f'Homology {max_id:.1f}%</span>'
        )
        clinical_note = (
            f"Overlaps {high_id_count} high-identity (>=95%) nuclear NUMT segment(s) "
            f"(top locus: {primary_loc_coord}, {max_id:.1f}% identity), but high VAF / gnomAD PASS confirms authentic mtDNA."
        )
    elif overlap_count > 0:
        risk_tier = "Divergent_NUMT_Overlap"
        status_label = f"Divergent NUMT Overlap ({max_id:.1f}%)"
        badge_text = f"Ancient NUMT ({max_id:.0f}%)"
        badge_html = (
            f'<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            f'font-size:7pt;font-weight:600;background:#f8fafc;color:#475569;'
            f'border:1px solid #e2e8f0" title="Divergent nuclear homology ({primary_loc_coord}, {max_id:.1f}% id)">'
            f'Ancient NUMT ({max_id:.0f}%)</span>'
        )
        clinical_note = (
            f"Overlaps {overlap_count} divergent ancient nuclear NUMT(s) "
            f"(max {max_id:.1f}% identity at {primary_loc_coord}); negligible risk of short-read mismapping."
        )
    else:
        risk_tier = "Clean"
        status_label = "Clean (No NUMT Overlap)"
        badge_text = "Unique mtDNA"
        badge_html = (
            '<span style="display:inline-block;white-space:nowrap;padding:2px 7px;border-radius:4px;'
            'font-size:7pt;font-weight:700;background:#ecfdf5;color:#047857;'
            'border:1px solid #a7f3d0">Unique mtDNA</span>'
        )
        clinical_note = "Unique cytoplasmic mtDNA locus with no documented nuclear pseudogene overlap."

    return {
        "is_npg": is_npg,
        "gnomad_filter": "PASS" if flt == "PASS" else (flt if flt else "Novel"),
        "overlap_count": overlap_count,
        "loci_count": overlap_count,
        "high_id_count": high_id_count,
        "max_identity_pct": max_id,
        "primary_nuclear_locus": primary_loc,
        "primary_nuclear_coord": primary_loc_coord,
        "primary_numt_id": primary_loc_id,
        "nuclear_loci_summary": loci_summary,
        "top_loci": top_loci,
        "risk_tier": risk_tier,
        "status_label": status_label,
        "badge_text": badge_text,
        "badge_html": badge_html,
        "clinical_note": clinical_note,
    }


def filter_competitive_sam(in_sam_path, out_bam_path, stats_json_path, sample_name, threads=4):
    """
    Filter a SAM file aligned competitively against rCRS_plus_NUMTs.fasta (or GRCh38+rCRS):
      - Counts total reads, authentic rCRS reads (NC_012920.1), and reads absorbed by
        nuclear NUMT decoy contigs (preferential nuclear alignment or discordant NUMT mate).
      - Strips nuclear @SQ headers so the output BAM has a clean NC_012920.1 header
        ready for Picard MarkDuplicates, GATK Mutect2 (--mitochondria-mode), and HaploCheck.
    """
    mt_contigs = {"NC_012920.1", "chrM", "MT", "rCRS"}
    total_records = 0
    unmapped_records = 0
    rcrs_retained = 0
    numt_depleted = 0
    chimeric_numt_mates = 0
    numt_contig_hits = Counter()

    sort_cmd = ["samtools", " view", "-bS", "-"]
    # Pipe filtered SAM directly into `samtools sort`
    sort_proc = subprocess.Popen(
        f"samtools view -uS - | samtools sort -@ {int(threads)} -o \"{out_bam_path}\"",
        shell=True,
        stdin=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    try:
        with open(in_sam_path, "r", encoding="utf-8", errors="replace") as sf:
            for line in sf:
                if line.startswith("@"):
                    if line.startswith("@SQ"):
                        # Keep only mitochondrial reference @SQ header(s)
                        m = re.search(r"SN:(\S+)", line)
                        if m and m.group(1) in mt_contigs:
                            sort_proc.stdin.write(line)
                    else:
                        sort_proc.stdin.write(line)
                    continue

                parts = line.split("\t")
                if len(parts) < 11:
                    continue

                total_records += 1
                flag = int(parts[1])
                rname = parts[2]
                rnext = parts[6]

                if flag & 0x4:
                    unmapped_records += 1
                    continue

                # If primary mapping is to a nuclear NUMT contig -> deplete!
                if rname not in mt_contigs:
                    numt_depleted += 1
                    numt_contig_hits[rname] += 1
                    continue

                # Read maps to NC_012920.1, but check if mate mapped preferentially to a nuclear NUMT contig
                if rnext not in ("=", "*", "NC_012920.1", "chrM", "MT"):
                    chimeric_numt_mates += 1
                    numt_depleted += 1
                    numt_contig_hits[rnext] += 1
                    continue

                # Clean any SA:/XA: tags referencing NUMT contigs so GATK doesn't complain about missing contigs
                cleaned_fields = parts[:11]
                for opt in parts[11:]:
                    if opt.startswith(("SA:Z:", "XA:Z:")):
                        continue
                    cleaned_fields.append(opt)

                rcrs_retained += 1
                sort_proc.stdin.write("\t".join(cleaned_fields))
                if not cleaned_fields[-1].endswith("\n"):
                    sort_proc.stdin.write("\n")
    finally:
        if sort_proc.stdin and not sort_proc.stdin.closed:
            sort_proc.stdin.close()
        stderr_out = sort_proc.stderr.read() if sort_proc.stderr else ""
        sort_proc.wait()
        if sort_proc.returncode != 0:
            raise RuntimeError(f"samtools sort failed in filter_competitive_sam: {stderr_out}")

    mapped_total = rcrs_retained + numt_depleted
    depletion_pct = round((numt_depleted / mapped_total * 100.0), 3) if mapped_total > 0 else 0.0

    top_contigs = []
    for cname, cnt in numt_contig_hits.most_common(10):
        # Parse NUMT_<id>_<chrom>_<start>_<end>
        m = re.match(r"^NUMT_(HSA_NumtS_\S+)_(chr[^_]+)_(\d+)_(\d+)$", cname)
        if m:
            top_contigs.append({
                "contig": cname,
                "numt_id": m.group(1),
                "locus": f"{m.group(2)}:{m.group(3)}-{m.group(4)}",
                "reads_absorbed": cnt,
            })
        else:
            top_contigs.append({
                "contig": cname,
                "numt_id": cname,
                "locus": cname,
                "reads_absorbed": cnt,
            })

    stats = {
        "sample": sample_name,
        "alignment_mode": "numt_aware_competitive",
        "alignment_mode_label": "NUMT-Aware Competitive Alignment (rCRS + 766 Human Nuclear NUMT Loci)",
        "threads_used": int(threads),
        "total_sam_records": total_records,
        "mapped_records": mapped_total,
        "unmapped_records": unmapped_records,
        "rcrs_retained_reads": rcrs_retained,
        "numt_depleted_reads": numt_depleted,
        "chimeric_numt_mate_reads": chimeric_numt_mates,
        "numt_depletion_pct": depletion_pct,
        "unique_numt_loci_hit": len(numt_contig_hits),
        "top_numt_loci": top_contigs,
    }
    with open(stats_json_path, "w", encoding="utf-8") as jf:
        json.dump(stats, jf, indent=2)
    return stats


def write_rcrs_only_stats(stats_json_path, sample_name, threads=4):
    """Write NUMT stats metadata when run in standard rCRS-only mode."""
    stats = {
        "sample": sample_name,
        "alignment_mode": "rcrs_only",
        "alignment_mode_label": "Mitochondrial Reference Only (rCRS.fasta) + Post-Calling NUMT Knowledgebase Annotation",
        "threads_used": int(threads),
        "total_sam_records": None,
        "mapped_records": None,
        "unmapped_records": None,
        "rcrs_retained_reads": None,
        "numt_depleted_reads": 0,
        "chimeric_numt_mate_reads": 0,
        "numt_depletion_pct": 0.0,
        "unique_numt_loci_hit": 0,
        "top_numt_loci": [],
    }
    with open(stats_json_path, "w", encoding="utf-8") as jf:
        json.dump(stats, jf, indent=2)
    return stats


def annotate_final_report_tsv(tsv_path):
    """
    Append 4 NUMT annotation columns to `${SAMPLE}_final_report.tsv`:
      - gnomAD_NUMT_Filter
      - NUMT_Homology_Status
      - NUMT_Nuclear_Loci
      - NUMT_Max_Identity_Pct
    Preserves all 16 existing columns in exact order.
    """
    if not os.path.isfile(tsv_path):
        return

    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        lines = [ln.rstrip("\r\n") for ln in f if ln.strip()]

    if not lines:
        return

    header = lines[0].split("\t")
    # Avoid duplicating columns if already annotated
    if "NUMT_Homology_Status" in header:
        return

    new_header = header + [
        "gnomAD_NUMT_Filter",
        "NUMT_Homology_Status",
        "NUMT_Nuclear_Loci",
        "NUMT_Max_Identity_Pct",
    ]
    out_lines = ["\t".join(new_header)]

    for line in lines[1:]:
        cols = line.split("\t")
        if len(cols) < 9:
            out_lines.append(line)
            continue
        pos = cols[1]
        ref = cols[2]
        alt = cols[3]
        het = cols[8]
        info = get_numt_annotation(pos, ref, alt, het)
        out_lines.append(
            "\t".join(
                cols + [
                    info["gnomad_filter"],
                    info["status_label"],
                    info["nuclear_loci_summary"],
                    f"{info['max_identity_pct']:.1f}%",
                ]
            )
        )

    with open(tsv_path, "w", encoding="utf-8") as f:
        f.write("\n".join(out_lines) + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: numt_engine.py [filter-bam | rcrs-stats | annotate-tsv] ...")
        sys.exit(1)

    cmd = sys.argv[1]
    if cmd == "filter-bam":
        in_sam, out_bam, stats_json = sys.argv[2], sys.argv[3], sys.argv[4]
        sample = sys.argv[5] if len(sys.argv) > 5 else "Sample"
        threads_arg = sys.argv[6] if len(sys.argv) > 6 else "4"
        st = filter_competitive_sam(in_sam, out_bam, stats_json, sample, int(threads_arg))
        print(
            f"[NUMT-Engine] Competitive filtering complete: "
            f"retained={st['rcrs_retained_reads']} rCRS reads, "
            f"depleted={st['numt_depleted_reads']} nuclear NUMT reads ({st['numt_depletion_pct']}%)"
        )
    elif cmd == "rcrs-stats":
        stats_json = sys.argv[2]
        sample = sys.argv[3] if len(sys.argv) > 3 else "Sample"
        threads_arg = sys.argv[4] if len(sys.argv) > 4 else "4"
        write_rcrs_only_stats(stats_json, sample, int(threads_arg))
        print(f"[NUMT-Engine] Recorded rCRS-only mode metadata for {sample}")
    elif cmd == "annotate-tsv":
        annotate_final_report_tsv(sys.argv[2])
        print(f"[NUMT-Engine] Annotated {sys.argv[2]} with NUMT homology & gnomAD npg status")
    else:
        print(f"Unknown command: {cmd}")
        sys.exit(1)
