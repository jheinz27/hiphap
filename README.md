# HipHap (previously Diplinator)

Diploid genome assemblies are now routinely available, but most read aligners were designed for haploid references. When reads are aligned to a diploid assembly, the aligner sees two nearly identical alignments to either haplotype, and thus reduces the mapping quality (MapQ) score to reflect this ambiguity. This can cause downstream tools to discard reads from easily mappable regions.
HipHap resolves this issue by aligning reads to each haplotype assembly separately and assigning each read to its best-supported haplotype. We also introduce a haplotype assignment quality score (HapQ) in HipHap to quantify confidence in the haplotype of origin of a read.

HipHap is implemented in Rust, and supports SAM, BAM, CRAM, and PAF formats


## Installation

Recommended: Download precompiled binary:
```bash
# Todo
```

Build from source:
```bash
git clone https://github.com/jheinz27/hiphap.git
cd hiphap
cargo build --release
./target/release/hiphap
```

## Usage

```
HipHap: Choose the best alignment to each haploid of a diploid assembly

Usage: hiphap [OPTIONS] <ASM1> <ASM2>

Arguments:
  <ASM1>  asm1 alignment file (sam/bam/cram/paf)
  <ASM2>  asm2 alignment file (sam/bam/cram/paf)

Options:
  -1, --s1 <NAME>          label for asm1 sample (used in output file names and summary) [default: asm1]
  -2, --s2 <NAME>          label for asm2 sample (used in output file names and summary) [default: asm2]
      --paf                input files are PAF
      --ms                 use ms:i: tag rather than AS:i: for alignment score
  -b, --both               write reads with equal alignment scores to both output files
  -u, --unmapped <DEST>    where to write reads unmapped in both assemblies: asm1, asm2, or discard [default: asm1] [possible values: asm1, asm2, discard]
      --ref1 <FILE>        reference FASTA for cram file (asm1)
      --ref2 <FILE>        reference FASTA for cram file (asm2)
      --match-sc <FLOAT>   per-base match score from aligner scoring scheme (e.g. minimap2 default is 2.0 for long reads) [default: 2.0]
      --no-hapq            skip HAPQ score calculation and hq tag output (e.g. for comparing grch38 vs chm13)
  -t, --threads <INT>      Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance. [default: 8]
  -h, --help               Print help
  -V, --version            Print version
```

Each output record is annotated with an `hq:i:` tag carrying the HAPQ score (see [HAPQ](#hapq-haplotype-assignment-quality)), unless `--no-hapq` is set.

## Example Workflow

HipHap only works on name-sorted files, which is the default [minimap2](https://github.com/lh3/minimap2) output. Therefore, coordinate-sorted files need to b name-sort first, for example:

```bash
samtools sort -n -o name_sort.bam index_sort.bam
```

### Diploid assembly alignment

```bash
# If needed, split diploid genome assembly FASTA into respective haplotypes
separate_haps_fasta -1 MATERNAL -2 PATERNAL hg002v1.1.fa
# writes hg002v1.1.MATERNAL.fa and hg002v1.1.PATERNAL.fa

# Align reads to each haplotype
minimap2 -ax map-hifi -o asm1_alignments.sam hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -ax map-hifi -o asm2_alignments.sam hg002v1.1.PATERNAL.fa reads.fastq

# Run hiphap
hiphap -1 mat -2 pat asm1_alignments.sam asm2_alignments.sam
# Output: hiphap_mat.sam  hiphap_pat.sam

# Merge best alignments into one file (if desired)
samtools merge -@ 12 merged.sam hiphap_mat.sam hiphap_pat.sam

# Save as sorted BAM
samtools sort -@ 12 -o merged.bam merged.sam
```

### Comparing different reference genomes

HipHap can also be used to select best alignments between different reference genomes (e.g. GRCh38 and CHM13). For this use case, the HAPQ score is generally not meaningful, so pass `--no-hapq` to skip its calculation:

```bash
hiphap --no-hapq -1 grch38 -2 chm13 grch38_alignments.sam chm13_alignments.sam
# Output: hiphap_grch38.sam  hiphap_chm13.sam
```

### CRAM input files

If input files are CRAM format, the original reference genomes must be provided, for example:

```bash
hiphap --ref1 asm1_hap.fasta --ref2 asm2_hap.fasta asm1_alignments.cram asm2_alignments.cram
```

## Example PAF Usage

**NOTE:** It is important to use the `--paf-no-hit` flags when aligning with minimap2. If a SAM file is converted to a PAF file with `paftools.js sam2paf`, it will NOT have the required AS:i: tag.

```bash
minimap2 -cx map-hifi -o asm1_alignments.paf hg002v1.1.MATERNAL.fa reads.fastq
minimap2 -cx map-hifi -o asm2_alignments.paf hg002v1.1.PATERNAL.fa reads.fastq

hiphap --paf out_asm1.paf out_asm2.paf
# Output: hiphap_asm1.paf  hiphap_asm2.paf
```

## Weighted Alignment Scoring Mechanism

For each read, HipHap computes a single weighted alignment score per assembly using all primary and supplementary alignments (secondary alignments are passed through to the output but ignored when scoring).

Let $n$ be the number of primary and supplementary alignments for a given read to that reference genome.

Let $L$ be the full read length (sum of query-consuming CIGAR operations on a non-secondary record, or the `qlen` field of the PAF record).

Let $I_i = [r_i^{\mathrm{start}}, r_i^{\mathrm{end}})$ denote the interval on the read covered by alignment $i$.

Let

$$
B = \left| \bigcup_{i=1}^{n} I_i \right|
$$

be the number of **read bases** covered by at least one alignment.

Let $a_i$ be the alignment score for alignment $i$ (`AS:i:` by default, or `ms:i:` if `--ms` is set).

Let $l_i$ be the alignment length in read coordinates for alignment $i$.

$$
S = \frac{\sum_{i=1}^{n} a_i}{\sum_{i=1}^{n} l_i} \cdot B \cdot \frac{B}{L}
$$

The first factor is the average alignment score per aligned base. It is multiplied by the number of unique read bases covered, $B$, and then scaled by the read coverage fraction $B/L$, so that reads which align over a large fraction of their length are weighted more heavily than reads which align only over a small portion.

For each read, the assembly with the higher $S$ wins; its full alignment cluster (including secondary alignments) is written to the corresponding output file. If $S$ is equal in both assemblies, the "better" assignment is determined by a hash of the read name, or the read is written to both output files when `--both` is used.

## HapQ (haplotype assignment quality)

For each read assigned to a winning haplotype, HipHap reports a HapQ score in the `hq:i:` tag of the output record. HapQ is a Phred-like confidence [0-60] that the read was assigned to the correct haplotype. The calculation is modeled on BWA-MEM's `mem_approx_mapq_se`.

Let $S_w$ and $S_l$ be the weighted alignment scores of the winning and losing assemblies, $m$ the per-base match scoreof the alignment software used (`--match-sc`, default `2.0`), and $k$ the number of non-secondary alignments (splits) on the winning side.

HapQ is the product of: 

d is approximately the difference, in matching bases, between the winning and losing alignments.
```math
d = \frac{S_w - S_l}{m}
```

The split penalty $\rho$ down-weights reads with more than three split alignments, which often fall in complex/repetitive regions where haplotype assignment is less reliable.
```math
\rho = \begin{cases} 1 & k \le 3 \\ \frac{3}{k} & k > 3 \end{cases}
```

The final HapQ calculation is then: 
```math
\text{HapQ}\;=\;\lfloor\, 6.02 \cdot \rho \cdot d \,\rfloor
```

Clamped to the range [0,60]. 

Special cases:
- Read mapped in only one assembly: HapQ = 60.
- Read tied between assemblies (winner = `Both`): HapQ = 0.
- Read unmapped in both assemblies: no `hq` tag is written.

If `--no-hapq` is set, HAPQ is not computed and no `hq:i:` tag is added (recommended when the two inputs are not haplotypes of the same sample, e.g. GRCh38 vs CHM13).

## Citation
If HipHap has helped you in your research, please cite our preprint at: TODO
