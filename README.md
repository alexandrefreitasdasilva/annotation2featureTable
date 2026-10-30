# annotation2featureTable

Convert **GenBank (`.gb`/`.gbk`)** or **GFF/GFF3 (`.gff`/`.gff3`)** sequence annotations into **NCBI 5-column tab-delimited Feature Tables (`.tbl` / `.featuretable`)** ready for submission via **NCBI BankIt** or **`table2asn`**.

## Features

- **Dual Format Auto-Detection:** Supports GenBank (`.gb`, `.gbk`, `.genbank`) and GFF/GFF3/GTF (`.gff`, `.gff3`, `.gtf`) files out of the box.
- **Full Coordinate Support:**
  - Multi-interval features (`join(...)`), including viral **programmed ribosomal slippage** (`exception	ribosomal slippage`) and spliced exons.
  - Reverse/complement strand intervals (`stop	start` ordered `5' -> 3'`).
  - Partial `5'` (`<`) and `3'` (`>`) coordinates.
  - Multi-contig / segmented genome inputs (one `>Feature <SeqID>` block per sequence).
- **BankIt Sanitization & Validation:**
  - Automatically synchronizes `>Feature <SeqID>` with FASTA headers (ignoring local software URIs such as Geneious `urn.local...`).
  - Strips foreign RefSeq protein accessions (`YP_`, `NP_`, `WP_`, etc.) and software boilerplate notes that trigger BankIt rejection.
  - Enforces strict 7-bit ASCII across headers and qualifiers, warning if the companion FASTA contains accented/non-ASCII characters.
  - Built-in **CDS translation check** (`ATG` start codon, terminal stop codon, multiple-of-3 frame length, and internal stop codons).

## Installation

```bash
pip install -r requirements.txt
```

## Usage

### Basic usage (GenBank + FASTA)
```bash
python3 annotation2featureTable.py \
  -i genome_annotation.gb \
  -f genome_sequence.fasta \
  -o genome.featuretable
```

### Overriding Sequence ID (e.g., shortened ID for BankIt)
```bash
python3 annotation2featureTable.py \
  -i genome_annotation.gb \
  -f genome_sequence.fasta \
  --seqid Short_Seq_ID \
  -o genome.featuretable
```

### Converting GFF/GFF3 annotations
```bash
python3 annotation2featureTable.py \
  -i predictions.gff3 \
  -f contigs.fasta \
  -o predictions.tbl
```

### Including `gene` blocks, `UTR`s, or biological `/note` qualifiers
```bash
python3 annotation2featureTable.py \
  -i annotation.gb \
  -f sequence.fasta \
  -o output.tbl \
  --include-gene \
  --include-notes \
  --types "CDS,mRNA,tRNA,rRNA,5'UTR,3'UTR"
```

## Command-Line Options

| Argument | Description |
| :--- | :--- |
| `-i`, `--input` | Input annotation file (`.gb`, `.gbk`, `.gff`, `.gff3`). |
| `-o`, `--output` | Output 5-column Feature Table path (`.featuretable` or `.tbl`). |
| `-f`, `--fasta` | Optional FASTA file to match `SeqID`s and validate CDS translations. |
| `--seqid` | Manually override the sequence identifier in `>Feature <SeqID>`. |
| `--format` | Input format (`auto`, `gb`, `gff`; default: `auto`). |
| `--types` | Comma-separated feature types to export (default: `CDS,tRNA,rRNA,ncRNA,mat_peptide`). |
| `--include-gene` | Emit `gene` features and paired `/gene` qualifiers. |
| `--include-notes` | Include biological `/note` qualifiers (excluding Geneious boilerplate). |
| `--keep-locus-tag` | Retain `locus_tag` qualifiers (for BioProject/PGAP `table2asn` submissions). |
