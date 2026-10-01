# HipHap 

<img align="left" width="300" height="200" alt="hiphap logo" src="https://github.com/user-attachments/assets/794d2471-cdc8-46cb-b63a-f1c11fb5d51e" /> 

(**previously Diplinator**) Diploid genome assemblies are now routinely available, but most read aligners were designed for haploid references. When reads are aligned to a diploid assembly, the aligner sees two nearly identical alignments to either haplotype, and thus reduces the mapping quality (MapQ) score to reflect this ambiguity. This can cause downstream tools to discard reads from easily mappable regions.
HipHap resolves this issue by aligning reads to each haplotype assembly separately and assigning each read to its best-supported haplotype. We also introduce a haplotype assignment quality score (HapQ) in HipHap to quantify confidence in the haplotype of origin of a read. HipHap is implemented in Rust, and supports SAM, BAM, CRAM, and PAF formats


## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Example Workflow](#example-workflow)
  - [Partitioned mode](#partitioned-mode-separate-output-file-per-haplotype--p)
  - [CRAM input files](#cram-input-files)
  - [Comparing different reference genomes](#comparing-different-reference-genomes)
- [Example PAF Usage](#example-paf-usage)
- [Losing-haplotype alignments](#losing-haplotype-alignments---keep-loser)
  - [SAM/BAM/CRAM](#sambamcram)
  - [PAF](#paf)
- [Citation](#citation)


## Installation

### RECOMMENDED: Download precompiled binary:
Prebuilt binaries can be downloaded from the [Releases](https://github.com/jheinz27/hiphap/releases) page.

```bash
wget https://github.com/jheinz27/hiphap/releases/download/v0.1.0/hiphap-v0.1.0-x86_64-unknown-linux-musl.tar.gz
tar -xvzf hiphap-v0.1.0-x86_64-unknown-linux-musl.tar.gz
hiphap-v0.1.0-x86_64-unknown-linux-musl/hiphap --help
``` 

### Build from source:
```bash
git clone https://github.com/jheinz27/hiphap.git
cd hiphap
cargo build --release
./target/release/hiphap --help
```

If Rust is not already installed, install it with [rustup](https://rustup.rs/). This installs the Rust compiler (`rustc`), Cargo package manager (`cargo`), and the `rustup` toolchain manager:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
``` 
#### Dependencies 
- `Rust >= 1.70`, including Cargo
- The following Rust crate dependencies are installed automatically by Cargo during the build:
  - `clap = "4.0"`
  - `rust-htslib = "0.46.0"`
  - `twox-hash = "1.6"`
  - `flate2 = "1"`
  - `rand = "0.8"`
 

## Usage

```
HipHap: Choose the best alignment of a read to each haploid of a diploid assembly

Usage: hiphap [OPTIONS] <ASM1> <ASM2>

Arguments:
  <ASM1>  asm1 alignment file (sam/bam/cram/paf)
  <ASM2>  asm2 alignment file (sam/bam/cram/paf)

Options:
  -1, --s1 <NAME>           label for asm1 sample (used in output names) [default: asm1]
  -2, --s2 <NAME>           label for asm2 sample (used in output names) [default: asm2]
      --paf                 input files are PAF
      --ms                  use ms:i: tag rather than AS:i: for alignment score
  -p, --partition           write one file per haplotype instead of a single merged output file
  -b, --both                write reads with equal alignment scores to both output files (requires -p)
  -k, --keep-loser          also write each read's best alignment to the losing haplotype as a secondary record, tagged hs:A: (merged output only)
      --loser-frac <FLOAT>  score floor for AS score fraction in --keep-loser [default: 0.8]
  -o, --output <FILE>       output file name [default: hiphap_{s1}_{s2}_merged.*]
  -u, --unmapped <DEST>     where to write reads unmapped in both assemblies: asm1, asm2, or discard [default: asm1] [possible values: asm1, asm2, discard]
      --ref-merged <FILE>   combined reference FASTA for merged CRAM output; required for merged CRAM output
      --ref1 <FILE>         reference FASTA for cram file (asm1)
      --ref2 <FILE>         reference FASTA for cram file (asm2)
  -A, --match-sc <FLOAT>    per-base match score from aligner scoring scheme (auto-estimated if omitted)
      --no-hapq             skip HAPQ score calculation and hq tag output (e.g. for comparing GRCh38 vs CHM13)
      --no-span-chrom       disable writing the chromosome-spanning reads file (*_span_chrom.fastq, or .txt for PAF)
  -t, --threads <INT>       number of threads[default: 6; 8 with -p]
  -h, --help                Print help
  -V, --version             Print version
```

Each output record is annotated with two tags: 
1) An `HP:i:` tag carrying the haplotype it was assigned to (`1` = asm1, `2` = asm2),
2) An `hq:i:` tag carrying the HapQ score, unless `--no-hapq` is set.

A HapQ value of 0 (`hq:i:0`) is reserved to mark a true tie between the two haplotypes. Every non-tied alignment is reported as at least a HapQ of 1 (`hq:i:1`), so a HapQ of 0 never results from rounding.

The `--keep-loser` flag can be set to also write a read's alignment to the *losing* haplotype, flagged secondary and marked with an `hs:A:` tag. See [losing-haplotype alignments](#losing-haplotype-alignments---keep-loser).

The per-base match score used by the HAPQ calculation is auto-estimated from the `ms:i:` tags of the input files. It can be set explicitly with `-A`/`--match-sc <FLOAT>`

## Example Workflow

HipHap reads the two inputs as parallel streams and clusters records by read, so **both files must have the reads in the same order**. Default [minimap2](https://github.com/lh3/minimap2) output satisfies this, as does name-sorting both files. **Coordinate-sorted input files will fail**. 

```bash
# If needed, split diploid genome assembly FASTA into respective haplotypes
separate_haps_fasta -1 MATERNAL -2 PATERNAL hg002v1.1.fa
# writes hg002v1.1.MATERNAL.fa and hg002v1.1.PATERNAL.fa

# Align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.PATERNAL.fa reads.fastq

# Run hiphap merged mode (default)
hiphap -1 mat -2 pat -o sample.sam asm1_alignments.sam asm2_alignments.sam
# Output: sample.sam  sample_span_chrom.fastq

# Save as sorted BAM
samtools sort -@ 12 -o sample_sorted.bam sample.sam
```

By default HipHap writes a **single merged output file** with a merged header, so no separate `samtools merge` step is needed. Each record carries an `HP:i:1`/`HP:i:2` tag naming the haplotype it was assigned to. 


HipHap also writes `hiphap_{s1}_{s2}_span_chrom.fastq`, a FASTQ of reads whose alignments span more than one chromosome. These reads are emitted for easy realignment. Pass `--no-span-chrom` to skip writing this file. If inputs are PAF files, a tsv of the chromosome spanning reads rather than a fastq is output, as the read sequences are not stored in the input PAF files. When `-o` is given, this file mimics the naming (`sample_span_chrom.fastq` above) so it lands beside the alignment output.

#### Notes:
- Merged output requires the two assemblies to have **unique contig names**. The merged header concatenates the two inputs `@SQ` lists, so any contig name shared between them is ambiguous; pass `-p` for such inputs.
- The output format always follows the input format. An `-o` extension that disagrees is ignored, with a warning.

### Partitioned mode: Separate output file per haplotype (`-p`)

`-p`/`--partition` writes one output file per haplotype instead of the merged file:

```bash
hiphap -p -1 mat -2 pat asm1_alignments.sam asm2_alignments.sam
# Output: hiphap_mat.sam  hiphap_pat.sam  hiphap_mat_pat_span_chrom.fastq

# Merge them afterwards if desired
samtools merge -@ 12 merged.sam hiphap_mat.sam hiphap_pat.sam
```
With `-p`, `-o` is the stem the two files share:

```bash
hiphap -p -1 mat -2 pat -o sample asm1_alignments.sam asm2_alignments.sam
# Output: sample_mat.sam  sample_pat.sam  sample_span_chrom.fastq
```

#### Notes:
- `--both` requires `-p`: in a merged file a tied read would get two primary records.
- `--keep-loser` cannot be used with `-p`: a per-haplotype file holds one assembly's contigs, so there is nowhere in it to put the other haplotype's alignments.

### CRAM input files

If input files are CRAM format, the original reference genomes must be provided. Merged output additionally needs a combined reference (`--ref-merged`), typically the original diploid reference, since a single CRAM writer needs one reference covering all contigs of both haplotypes:

```bash
hiphap --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta --ref-merged diploid.fasta asm1_alignments.cram asm2_alignments.cram

# per-haplotype partitioned output needs only the two originals
hiphap -p --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta asm1_alignments.cram asm2_alignments.cram
```

### Comparing different reference genomes

HipHap can also be used to select best alignment of a read between different reference genomes (e.g. GRCh38 and CHM13). For this use case, the HAPQ score is generally not meaningful, so pass `--no-hapq` to skip its calculation:

```bash
hiphap --no-hapq -p -1 grch38 -2 chm13 grch38_alignments.sam chm13_alignments.sam
# Output: hiphap_grch38.sam  hiphap_chm13.sam
```

## Example PAF Usage

 #### Notes: 
- It is important to use the `--paf-no-hit` flags when aligning with minimap2 to output unmapped reads to the file
- If a SAM file is converted to a PAF file with `paftools.js sam2paf`, it will **NOT** have the required AS:i: tag and HipHap will fail to run

```bash
minimap2 -cx map-hifi --paf-no-hit -o asm1_alignments.paf hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -cx map-hifi --paf-no-hit -o asm2_alignments.paf hg002v1.1.PATERNAL.fa reads.fastq

hiphap --paf asm1_alignments.paf asm2_alignments.paf
# Output: hiphap_asm1_asm2_merged.paf  hiphap_asm1_asm2_span_chrom.txt
```

`-p` also works with PAF

```bash
hiphap --paf -p asm1_alignments.paf asm2_alignments.paf
# Output: hiphap_asm1.paf  hiphap_asm2.paf  hiphap_asm1_asm2_span_chrom.txt
```

## Losing-haplotype alignments (`--keep-loser`)

### SAM/BAM/CRAM
By default only the winning haplotype's alignments are written and the losing haplotype's are discarded. However, when the two haplotypes alignments score close to each other, the discarded alignment may also be of interest, so `--keep-loser`  writes it:

```bash
hiphap --keep-loser -1 mat -2 pat asm1_alignments.bam asm2_alignments.bam
```

A losing alignment is only written when its score is at least `--loser-frac` of the winner's — **0.8 by default**. Reads tied between the haplotypes (`hq:i:0`) always qualify. This threshold can be set.  

```bash
hiphap --keep-loser --loser-frac 0.95 asm1.bam asm2.bam   # near-ties only
hiphap --keep-loser --loser-frac 1.0  asm1.bam asm2.bam   # exact ties only
```

For each read, the losing haplotype's **primary and supplementary** alignments are appended to the merged file directly after the winning cluster, rewritten as secondary alignments. Secondary alignments to the losing haplotype are dropped. Each such record carries an — *haplotype status* — `hs:A:` tag naming whether the record was originally a primary or supplemental alignment in that haplotype. 

| Tag | Meaning | Resulting FLAG |
|---|---|---|
| `hs:A:P` | was the **primary** alignment to the losing haplotype | `0x100` (secondary) |
| `hs:A:S` | was a **supplementary** alignment to the losing haplotype | `0x900` (secondary + supplementary) |
| *absent* | was assigned to the **winning** haplotype | unchanged |

`hs` is the only field that differentiates a losing record from a secondary alignment of a winning record. Every read has exactly one `hs:A:P`, which can be used to count, subset, etc.. 

```bash
samtools view -d hs   merged.bam    # every losing-haplotype record
samtools view -d hs:P merged.bam    # primary alignment to the losing haplotype, one per read
samtools view -d hs:S merged.bam    # the losing cluster's split segments
samtools view -e '![hs]' merged.bam # winners only — the file hiphap writes without the --keep-loser flag
```

Details:

- **`SEQ` and `QUAL` are `*`.** The winning record for the same read already holds the read's bases. Writing them a second time would roughly double size of the output.
- **`hq:i:` is written** with the same HapQ as the winning record. 
- **No `HP:i:` tag.** on the records entering the file from the --keep-loser path` 
- **`MAPQ` is unchanged.** It is still the mapping quality of that alignment within the losing haplotype.
- `--keep-loser` records do not appear in the chromosome-spanning reads file, which tracks winning clusters only.
- `--keep-loser` cannot be combined with `-p`/`--partition`.

### PAF

`--keep-loser` works for merged PAF output too. Losing lines get their `tp:A:P` rewritten to `tp:A:S` and the `hq:i:`/`hs:A:` tags appended, and lines already marked `tp:A:S` or with a `*` target are dropped.

In PAF `hs` is always `P`: minimap2 gives every non-secondary chain `tp:A:P`, including the ones that become supplementary in SAM, so a PAF line carries nothing that separates the two. `hs:A:S` only ever appears in SAM/BAM/CRAM output.

## Citation
If HipHap has helped you in your research, please cite our [preprint](https://www.biorxiv.org/content/10.64898/2026.09.24.754153v1). 
