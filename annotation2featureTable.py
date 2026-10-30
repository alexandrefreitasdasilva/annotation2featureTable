#!/usr/bin/env python3
"""
Convert GenBank (.gb/.gbk) or GFF/GFF3 (.gff/.gff3) annotation files into
NCBI 5-column tab-delimited Feature Table (.featuretable / .tbl) files
ready for submission via NCBI BankIt or table2asn.
"""

import argparse
import os
import sys
import unicodedata
from collections import OrderedDict
from urllib.parse import unquote

from Bio import SeqIO
from Bio.SeqFeature import AfterPosition, BeforePosition, CompoundLocation


GENEIOUS_NOTE_PREFIXES = (
    "Derived using Geneious",
    "Geneious type:",
)

REFSEQ_PROTEIN_PREFIXES = (
    "YP_", "NP_", "XP_", "WP_", "AP_", "ZP_",
)


def detect_format(filepath, fmt_arg="auto"):
    """Automatically detect whether the input file is GenBank or GFF."""
    if fmt_arg and fmt_arg != "auto":
        return fmt_arg.lower()

    ext = os.path.splitext(filepath)[1].lower()
    if ext in (".gb", ".gbk", ".genbank", ".gbff"):
        return "gb"
    if ext in (".gff", ".gff3", ".gtf"):
        return "gff"

    with open(filepath, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("LOCUS"):
                return "gb"
            if stripped.startswith("##gff-version") or "\t" in stripped:
                return "gff"
            break

    raise ValueError(
        f"Could not detect file format for '{filepath}'. "
        "Please specify --format gb or --format gff."
    )


def to_ascii(text):
    """Normalize accented or non-ASCII characters into plain 7-bit ASCII (required by BankIt)."""
    normalized = unicodedata.normalize("NFKD", str(text))
    return normalized.encode("ascii", "ignore").decode("ascii")


def sanitize_seqid(raw_id):
    """Strip accidental '>' or 'Feature ' prefixes and enforce plain ASCII on Sequence IDs."""
    s = str(raw_id).strip()
    if s.startswith(">"):
        s = s[1:].strip()
    if s.lower().startswith("feature "):
        s = s[8:].strip()
    return to_ascii(s.split()[0]) if s else "Unknown_Seq"


def load_fasta_records(fasta_path):
    """Load FASTA records into an ordered dictionary {seqid: SeqRecord}."""
    records = OrderedDict()
    if not fasta_path:
        return records
    for rec in SeqIO.parse(fasta_path, "fasta"):
        if not rec.description.isascii():
            print(
                f"[BankIt WARNING] FASTA header for '{rec.id}' contains non-ASCII "
                "characters (accents/symbols). NCBI BankIt requires plain ASCII!",
                file=sys.stderr,
            )
        seqid = sanitize_seqid(rec.id)
        records[seqid] = rec
    return records


def resolve_seqid(candidate_ids, fasta_records, override_seqid=None):
    """
    Resolve the sequence identifier for '>Feature <SeqID>', ensuring it matches
    the corresponding FASTA header when provided.
    """
    if override_seqid:
        seqid = sanitize_seqid(override_seqid)
    elif fasta_records:
        seqid = None
        for cid in candidate_ids:
            if cid and sanitize_seqid(cid) in fasta_records:
                seqid = sanitize_seqid(cid)
                break
        # Fallback when the FASTA contains a single sequence record
        if not seqid and len(fasta_records) == 1:
            seqid = next(iter(fasta_records.keys()))
            print(
                f"[WARNING] Annotation ID ({candidate_ids[0]}) mapped to match "
                f"FASTA header ID: '{seqid}'",
                file=sys.stderr,
            )
        if not seqid:
            seqid = sanitize_seqid(candidate_ids[0]) if candidate_ids else "Unknown_Seq"
    else:
        seqid = None
        for cid in candidate_ids:
            # Skip internal Geneious dummy accessions such as 'urn.local...'
            if cid and not cid.startswith("urn.local") and cid not in (".", "unknown"):
                seqid = sanitize_seqid(cid)
                break
        if not seqid:
            seqid = sanitize_seqid(candidate_ids[0]) if candidate_ids else "Unknown_Seq"

    if len(seqid) > 50:
        print(
            f"[BankIt WARNING] SeqID '{seqid}' is {len(seqid)} characters long "
            "(GenBank maximum is 50; <= 25 is recommended for BankIt). "
            "Use --seqid <short_ID> if you shortened the FASTA header.",
            file=sys.stderr,
        )
    return seqid


def format_location_intervals(location):
    """
    Convert a BioPython FeatureLocation or CompoundLocation into a list of
    1-based (start_str, stop_str) tuples formatted for NCBI Feature Tables.
    """
    parts = location.parts if isinstance(location, CompoundLocation) else [location]
    intervals = []

    num_parts = len(parts)
    for idx, part in enumerate(parts):
        strand = part.strand if part.strand is not None else 1
        raw_start = int(part.start) + 1
        raw_end = int(part.end)

        is_before = isinstance(part.start, BeforePosition)
        is_after = isinstance(part.end, AfterPosition)

        if strand == -1:
            # On the reverse/complement strand, 5' is the higher coordinate (end)
            # and 3' is the lower coordinate (start)
            col1 = f"<{raw_end}" if (idx == 0 and is_after) else str(raw_end)
            col2 = f">{raw_start}" if (idx == num_parts - 1 and is_before) else str(raw_start)
        else:
            col1 = f"<{raw_start}" if (idx == 0 and is_before) else str(raw_start)
            col2 = f">{raw_end}" if (idx == num_parts - 1 and is_after) else str(raw_end)

        intervals.append((col1, col2))

    return intervals


def clean_text(val):
    """Strip newlines, tabs, non-ASCII characters, and redundant whitespace from qualifier values."""
    ascii_val = to_ascii(val)
    return " ".join(ascii_val.replace("\t", " ").splitlines()).strip()


def validate_cds_feature(feature, seq_record, label):
    """
    Verify that a CDS translates cleanly (reading frame length, start codon,
    terminal stop codon, and absence of internal stop codons) to catch potential
    BankIt validation errors early.
    """
    if seq_record is None or seq_record.seq is None or len(seq_record.seq) == 0:
        return []

    warnings = []
    try:
        dna_seq = feature.extract(seq_record.seq)
    except Exception as exc:
        return [f"{label}: failed to extract nucleotide sequence ({exc})"]

    codon_start = int(feature.qualifiers.get("codon_start", ["1"])[0])
    if codon_start in (2, 3):
        dna_seq = dna_seq[codon_start - 1 :]

    transl_table = int(feature.qualifiers.get("transl_table", ["1"])[0])

    parts = feature.location.parts if isinstance(feature.location, CompoundLocation) else [feature.location]
    strand = parts[0].strand if parts[0].strand is not None else 1
    if strand == -1:
        partial_5 = isinstance(parts[0].end, AfterPosition)
        partial_3 = isinstance(parts[-1].start, BeforePosition)
    else:
        partial_5 = isinstance(parts[0].start, BeforePosition)
        partial_3 = isinstance(parts[-1].end, AfterPosition)

    if len(dna_seq) % 3 != 0 and not (partial_5 or partial_3):
        warnings.append(
            f"{label}: CDS length ({len(dna_seq)} nt) is not a multiple of 3."
        )

    usable_len = (len(dna_seq) // 3) * 3
    prot = str(dna_seq[:usable_len].translate(table=transl_table))

    start_codon = str(dna_seq[:3]).upper()
    if not partial_5 and codon_start == 1 and start_codon != "ATG":
        warnings.append(
            f"{label}: does not start with ATG (found '{start_codon}'). "
            "BankIt may require coordinate adjustment or a 5'-partial ('<') marker."
        )

    if not partial_3 and prot and not prot.endswith("*"):
        last_codon = str(dna_seq[usable_len - 3 : usable_len]).upper()
        warnings.append(
            f"{label}: does not end with a stop codon (last codon: '{last_codon}'). "
            "BankIt may require coordinate adjustment or a 3'-partial ('>') marker."
        )

    internal_stops = prot[:-1].count("*") if prot.endswith("*") else prot.count("*")
    if internal_stops > 0:
        stop_positions = [i + 1 for i, aa in enumerate(prot[:-1]) if aa == "*"]
        warnings.append(
            f"{label}: contains {internal_stops} internal stop codon(s) at amino acid position(s) {stop_positions}!"
        )

    return warnings


def extract_gb_qualifiers(feature, include_gene=False, include_notes=False, keep_locus_tag=False):
    """Extract sanitized, BankIt-compliant qualifiers from a BioPython SeqFeature."""
    quals = []
    qdict = feature.qualifiers

    if feature.type == "gene":
        for gv in qdict.get("gene", []):
            quals.append(("gene", clean_text(gv)))
        if keep_locus_tag:
            for lt in qdict.get("locus_tag", []):
                quals.append(("locus_tag", clean_text(lt)))
        return quals

    # Product (mandatory for CDS, rRNA, tRNA, mat_peptide in BankIt)
    products = qdict.get("product", [])
    if products:
        for prod in products:
            quals.append(("product", clean_text(prod)))
    elif feature.type == "CDS":
        quals.append(("product", "hypothetical protein"))

    # Only include /gene on CDS when paired 'gene' features are enabled to prevent SEQ_FEAT.OnlyGeneXrefs
    if include_gene:
        for gv in qdict.get("gene", []):
            quals.append(("gene", clean_text(gv)))

    if keep_locus_tag:
        for lt in qdict.get("locus_tag", []):
            quals.append(("locus_tag", clean_text(lt)))

    # Reading frame offset for CDS
    if feature.type == "CDS":
        for cs in qdict.get("codon_start", ["1"]):
            quals.append(("codon_start", clean_text(cs)))

    # Genetic code table (when present)
    for tt in qdict.get("transl_table", []):
        quals.append(("transl_table", clean_text(tt)))

    # Translational exceptions (e.g., ribosomal slippage required for frameshifted join CDS)
    for exc in qdict.get("exception", []):
        quals.append(("exception", clean_text(exc)))

    if "ribosomal_slippage" in qdict and not qdict.get("exception"):
        quals.append(("exception", "ribosomal slippage"))

    # Skip external RefSeq protein accessions (YP_, NP_, etc.) which BankIt rejects on new submissions
    for pid in qdict.get("protein_id", []):
        pid_clean = clean_text(pid)
        if not pid_clean.startswith(REFSEQ_PROTEIN_PREFIXES):
            quals.append(("protein_id", pid_clean))

    # Biological notes (filtering out automatic Geneious annotation boilerplate)
    if include_notes:
        raw_notes = qdict.get("note", []) + qdict.get("Note", [])
        for note in raw_notes:
            cleaned = clean_text(note)
            if cleaned and not cleaned.startswith(GENEIOUS_NOTE_PREFIXES):
                quals.append(("note", cleaned))

    return quals


def write_feature_entry(outfile, ftype, intervals, qualifiers, qual_style="ncbi"):
    """Write a single feature (including multi-interval spans) and its qualifiers to the output table."""
    for idx, (start_str, stop_str) in enumerate(intervals):
        if idx == 0:
            outfile.write(f"{start_str}\t{stop_str}\t{ftype}\n")
        else:
            outfile.write(f"{start_str}\t{stop_str}\n")

    for key, val in qualifiers:
        if qual_style == "legacy_space" and key == "product":
            outfile.write(f"\t\t\tProduct {val}\n")
        elif val is None or val == "":
            outfile.write(f"\t\t\t{key}\n")
        else:
            outfile.write(f"\t\t\t{key}\t{val}\n")


def process_genbank(
    gb_path,
    output_path,
    fasta_records,
    allowed_types,
    override_seqid=None,
    include_gene=False,
    include_notes=False,
    keep_locus_tag=False,
    qual_style="ncbi",
):
    """Parse a GenBank (.gb/.gbk) file and write an NCBI 5-column Feature Table."""
    all_warnings = []
    feature_count = 0

    with open(output_path, "w", encoding="utf-8") as outfile:
        for rec in SeqIO.parse(gb_path, "genbank"):
            seqid = resolve_seqid([rec.name, rec.id], fasta_records, override_seqid)
            outfile.write(f">Feature {seqid}\n")

            seq_for_check = (
                fasta_records.get(seqid)
                or fasta_records.get(sanitize_seqid(rec.name))
                or (next(iter(fasta_records.values())) if len(fasta_records) == 1 else rec)
            )

            for feat in rec.features:
                # Skip source and full-contig Geneious region wrappers
                if feat.type == "source":
                    continue
                if feat.type == "misc_feature":
                    notes = " ".join(feat.qualifiers.get("note", []))
                    if "Geneious type: region" in notes:
                        continue

                if feat.type == "gene" and not include_gene:
                    continue
                if feat.type not in allowed_types:
                    continue

                # Skip empty gene features lacking both /gene and /locus_tag
                quals = extract_gb_qualifiers(
                    feat,
                    include_gene=include_gene,
                    include_notes=include_notes,
                    keep_locus_tag=keep_locus_tag,
                )
                if feat.type == "gene" and not quals:
                    continue

                intervals = format_location_intervals(feat.location)
                write_feature_entry(outfile, feat.type, intervals, quals, qual_style=qual_style)
                feature_count += 1

                if feat.type == "CDS":
                    prod_name = feat.qualifiers.get("product", ["CDS"])[0]
                    coord_desc = ",".join(f"{s}..{e}" for s, e in intervals)
                    label = f"CDS '{prod_name}' ({coord_desc})"
                    all_warnings.extend(validate_cds_feature(feat, seq_for_check, label))

    return feature_count, all_warnings


def parse_gff_attributes(attr_str):
    """Parse Column 9 attributes from GFF3 (key=value) or GTF (key "value") lines."""
    attrs = OrderedDict()
    if not attr_str or attr_str == ".":
        return attrs

    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[unquote(k.strip())] = unquote(v.strip()).strip('"')
        elif " " in part:
            k, v = part.split(" ", 1)
            attrs[k.strip()] = v.strip().strip('"')
    return attrs


def process_gff(
    gff_path,
    output_path,
    fasta_records,
    allowed_types,
    override_seqid=None,
    include_gene=False,
    include_notes=False,
    keep_locus_tag=False,
    qual_style="ncbi",
):
    """Parse a GFF/GFF3/GTF file and write an NCBI 5-column Feature Table (supporting multi-segment CDS)."""
    contig_features = OrderedDict()

    with open(gff_path, "r", encoding="utf-8", errors="replace") as infile:
        for raw_line in infile:
            line = raw_line.strip()
            if not line or line.startswith("#"):
                if line == "##FASTA":
                    break
                continue

            cols = line.split("\t")
            if len(cols) < 8:
                continue

            contig, _source, ftype, start_s, end_s, _score, strand, phase = cols[:8]
            attr_str = cols[8] if len(cols) > 8 else ""

            if ftype == "region":
                continue
            if ftype == "gene" and not include_gene:
                continue
            if ftype not in allowed_types:
                continue

            attrs = parse_gff_attributes(attr_str)
            start = int(start_s)
            end = int(end_s)

            if contig not in contig_features:
                contig_features[contig] = OrderedDict()

            feat_id = attrs.get("ID")
            # Group multi-exon/multi-segment CDS features sharing the same ID
            group_key = (
                (ftype, feat_id)
                if (ftype == "CDS" and feat_id)
                else (ftype, f"__line_{len(contig_features[contig])}")
            )

            if group_key not in contig_features[contig]:
                contig_features[contig][group_key] = {
                    "type": ftype,
                    "strand": strand,
                    "segments": [],
                    "attrs": attrs,
                }
            contig_features[contig][group_key]["segments"].append((start, end, phase))

    feature_count = 0
    with open(output_path, "w", encoding="utf-8") as outfile:
        for contig, grouped in contig_features.items():
            seqid = resolve_seqid([contig], fasta_records, override_seqid)
            outfile.write(f">Feature {seqid}\n")

            for item in grouped.values():
                ftype = item["type"]
                strand = item["strand"]
                segments = item["segments"]
                attrs = item["attrs"]

                # Sort segments in biological 5' -> 3' transcription order
                if strand == "-":
                    segments_sorted = sorted(segments, key=lambda x: x[1], reverse=True)
                    intervals = [(str(e), str(s)) for s, e, _ in segments_sorted]
                else:
                    segments_sorted = sorted(segments, key=lambda x: x[0])
                    intervals = [(str(s), str(e)) for s, e, _ in segments_sorted]

                first_phase = segments_sorted[0][2]
                quals = []

                if ftype == "gene":
                    gene_name = attrs.get("gene") or attrs.get("Name") or attrs.get("gene_name")
                    if gene_name:
                        quals.append(("gene", clean_text(gene_name)))
                    if keep_locus_tag and attrs.get("locus_tag"):
                        quals.append(("locus_tag", clean_text(attrs["locus_tag"])))
                    if not quals:
                        continue
                else:
                    product = (
                        attrs.get("product")
                        or attrs.get("Product")
                        or attrs.get("description")
                        or ("hypothetical protein" if ftype == "CDS" else None)
                    )
                    if product:
                        quals.append(("product", clean_text(product)))

                    if include_gene:
                        gene_name = attrs.get("gene") or attrs.get("gene_name")
                        if gene_name:
                            quals.append(("gene", clean_text(gene_name)))

                    if keep_locus_tag and attrs.get("locus_tag"):
                        quals.append(("locus_tag", clean_text(attrs["locus_tag"])))

                    if ftype == "CDS":
                        if "codon_start" in attrs:
                            quals.append(("codon_start", clean_text(attrs["codon_start"])))
                        elif first_phase in ("0", "1", "2"):
                            quals.append(("codon_start", str(int(first_phase) + 1)))
                        else:
                            quals.append(("codon_start", "1"))

                    if "transl_table" in attrs:
                        quals.append(("transl_table", clean_text(attrs["transl_table"])))

                    if "exception" in attrs:
                        quals.append(("exception", clean_text(attrs["exception"])))
                    elif attrs.get("ribosomal_slippage") in ("true", "1", "True"):
                        quals.append(("exception", "ribosomal slippage"))

                    pid = attrs.get("protein_id")
                    if pid and not pid.startswith(REFSEQ_PROTEIN_PREFIXES):
                        quals.append(("protein_id", clean_text(pid)))

                    if include_notes:
                        note = attrs.get("Note") or attrs.get("note")
                        if note and not note.startswith(GENEIOUS_NOTE_PREFIXES):
                            quals.append(("note", clean_text(note)))

                write_feature_entry(outfile, ftype, intervals, quals, qual_style=qual_style)
                feature_count += 1

    return feature_count, []


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Convert GenBank (.gb/.gbk) or GFF/GFF3 (.gff/.gff3) annotations into "
            "an NCBI 5-column tab-delimited Feature Table for BankIt / table2asn."
        )
    )
    parser.add_argument(
        "-i", "--input", required=True, help="Input annotation file (.gb, .gbk, .gff, .gff3)."
    )
    parser.add_argument(
        "-o", "--output", required=True, help="Output Feature Table file (.featuretable / .tbl)."
    )
    parser.add_argument(
        "-f",
        "--fasta",
        default=None,
        help="Optional corresponding FASTA file (matches SeqIDs and validates CDS translations).",
    )
    parser.add_argument(
        "--seqid",
        default=None,
        help="Manually override the sequence identifier in '>Feature <SeqID>'.",
    )
    parser.add_argument(
        "--format",
        choices=["auto", "gb", "gff"],
        default="auto",
        help="Input file format (default: auto).",
    )
    parser.add_argument(
        "--types",
        default="CDS,tRNA,rRNA,ncRNA,mat_peptide",
        help="Comma-separated feature types to include (default: CDS,tRNA,rRNA,ncRNA,mat_peptide).",
    )
    parser.add_argument(
        "--include-gene",
        action="store_true",
        help="Also emit 'gene' features and paired '/gene' qualifiers.",
    )
    parser.add_argument(
        "--include-notes",
        action="store_true",
        help="Include biological '/note' qualifiers (ignoring automatic Geneious notes).",
    )
    parser.add_argument(
        "--keep-locus-tag",
        action="store_true",
        help="Retain 'locus_tag' qualifiers from the source file.",
    )
    parser.add_argument(
        "--style",
        choices=["ncbi", "legacy_space"],
        default="ncbi",
        help=(
            "Qualifier formatting style: 'ncbi' uses official 5-column TAB separation "
            "(\\t\\t\\tproduct\\t<value>); 'legacy_space' uses '\\t\\t\\tProduct <value>'."
        ),
    )

    args = parser.parse_args()

    allowed_types = {t.strip() for t in args.types.split(",") if t.strip()}
    if args.include_gene:
        allowed_types.add("gene")

    fmt = detect_format(args.input, args.format)
    fasta_records = load_fasta_records(args.fasta)

    if fmt == "gb":
        count, warnings = process_genbank(
            gb_path=args.input,
            output_path=args.output,
            fasta_records=fasta_records,
            allowed_types=allowed_types,
            override_seqid=args.seqid,
            include_gene=args.include_gene,
            include_notes=args.include_notes,
            keep_locus_tag=args.keep_locus_tag,
            qual_style=args.style,
        )
    else:
        count, warnings = process_gff(
            gff_path=args.input,
            output_path=args.output,
            fasta_records=fasta_records,
            allowed_types=allowed_types,
            override_seqid=args.seqid,
            include_gene=args.include_gene,
            include_notes=args.include_notes,
            keep_locus_tag=args.keep_locus_tag,
            qual_style=args.style,
        )

    print(f"Successfully generated '{args.output}' ({count} features written).")

    if warnings:
        print("\n[WARNING - CDS Translation Check for BankIt]:", file=sys.stderr)
        for w in warnings:
            print(f"  - {w}", file=sys.stderr)


if __name__ == "__main__":
    main()
