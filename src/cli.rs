use clap::{Parser, ValueEnum};


#[derive(Parser, Debug)]
#[command( name = "HipHap", about = "HipHap: Choose the best alignment of a read to each haploid of a diploid assembly", version)]
pub struct Cli {
    //
    #[arg(value_name = "ASM1", help="asm1 alignment file (sam/bam/cram/paf)")]
    pub asm1: String,

    #[arg(value_name = "ASM2", help="asm2 alignment file (sam/bam/cram/paf)")]
    pub asm2: String,

    #[arg(short='1', long, value_name = "NAME", default_value = "asm1", help="label for asm1 sample (used in output names)")]
    pub s1: String,

    #[arg(short='2', long, value_name = "NAME", default_value = "asm2", help="label for asm2 sample (used in output names)")]
    pub s2: String,

    // inputs are PAF files
    #[arg(long, default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    //use ms score rather than AS score
    #[arg(long, default_value_t = false, help = "use ms:i: tag rather than AS:i: for alignment score")]
    pub ms: bool,

    // write one file per haplotype rather than a single merged output file (merged is the default)
    #[arg(short = 'p', long, default_value_t = false, help = "write one file per haplotype instead of a single merged output file")]
    pub partition: bool,

    // write tied reads to both output files (requires -p: two primaries cannot share one file)
    #[arg(short, long, default_value_t = false, help = "write reads with equal alignment scores to both output files (requires -p)")]
    pub both: bool,

    // also write each read's best alignment to the losing haplotype, flagged secondary.
    #[arg(short, long, default_value_t = false, conflicts_with = "partition", help = "also write each read's best alignment to the losing haplotype as a secondary record, tagged hs:A: (merged output only)")]
    pub keep_loser: bool,

    // score floor for --keep-loser 
    #[arg(long, value_name = "FLOAT", default_value_t = 0.8, help = "score floor for AS score fraction in --keep-loser")]
    pub loser_frac: f32,

    // output path: the merged file name, or the shared stem for the pair under --partition
    #[arg(short = 'o', long, value_name = "FILE", help = "output file name [default: hiphap_{s1}_{s2}_merged.*]")]
    pub output: Option<String>,

    // where to write reads unmapped in both assemblies
    #[arg(short, long, value_name = "DEST", default_value = "asm1", help="where to write reads unmapped in both assemblies: asm1, asm2, or discard")]
    pub unmapped: UnmappedDest,

    // combined reference FASTA for writing a merged CRAM (must contain all contigs of both haplotypes)
    #[arg(long, value_name = "FILE", required = false, help = "combined reference FASTA for merged CRAM output; required for merged CRAM output")]
    pub ref_merged: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm1)")]
    pub ref1: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm2)")]
    pub ref2: Option<String>,

    // per-base match score from aligner scoring scheme (used in HAPQ calculation)
    #[arg(short = 'A' , long, value_name = "FLOAT", help = "per-base match score from aligner scoring scheme (auto-estimated if omitted)")]
    pub match_sc: Option<f32>,

    // skip HAPQ score calculation and hq tag output (for non-haplotype comparisons)
    #[arg(long, default_value_t = false, help = "skip HAPQ score calculation and hq tag output (e.g. for comparing GRCh38 vs CHM13)")]
    pub no_hapq: bool,

    // disable writing the list of chromosome-spanning reads
    #[arg(long, default_value_t = false, help = "disable writing the chromosome-spanning reads file (*_span_chrom.fastq, or .txt for PAF)")]
    pub no_span_chrom: bool,

    // output format; defaults to the input format
    #[arg(short = 'O', long = "out-fmt", value_name = "FMT", value_enum, help = "output format: sam, bam, or cram [default: same as input]")]
    pub out_fmt: Option<OutFormat>,

    // number of total threads to use;
    #[arg(short, long,value_name = "INT", help = "number of threads[default: 6; 8 with -p]")]
    pub threads: Option<usize>
}

impl Cli {
    pub fn resolved_threads(&self) -> usize {
        self.threads.unwrap_or(if self.partition { 8 } else { 6 })
    }
}

#[derive(Debug, Clone, ValueEnum)]
pub enum UnmappedDest {
    Asm1,
    Asm2,
    Discard,
}

//output alignment format for a run; the alignment output follows the input format unless
//--out-fmt overrides it (e.g. SAM input written as compressed BAM)
#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum OutFormat {
    Sam,
    Bam,
    Cram,
}
