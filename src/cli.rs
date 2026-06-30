use clap::{Parser, ValueEnum};


#[derive(Parser, Debug)]
#[command( name = "HipHap", about = "HipHap: Choose the best alignment to each haploid of a diploid assembly", version)]

pub struct Cli {
    //
    #[arg(value_name = "ASM1", help="asm1 alignment file (sam/bam/cram/paf)")]
    pub asm1: String,

    #[arg(value_name = "ASM2", help="asm2 alignment file (sam/bam/cram/paf)")]
    pub asm2: String,

    #[arg(short='1', long, value_name = "NAME", default_value = "asm1", help="label for asm1 sample (used in output file names and summary)")]
    pub s1: String,

    #[arg(short='2', long, value_name = "NAME", default_value = "asm2", help="label for asm2 sample (used in output file names and summary)")]
    pub s2: String,

    // inputs are PAF files
    #[arg(long, default_value_t = false, help = "input files are PAF")]
    pub paf: bool,

    //use ms score rather than AS score
    #[arg(long, default_value_t = false, help = "use ms:i: tag rather than AS:i: for alignment score")]
    pub ms: bool,

    // write tied reads to both output files
    #[arg(short, long, default_value_t = false, help = "write reads with equal alignment scores to both output files")]
    pub both: bool,

    // write a single merged output file rather than one file per haplotype
    #[arg(short = 'm', long, default_value_t = false, help = "write a single merged output file (hiphap_{s1}_{s2}_merged.*) instead of one file per haplotype")]
    pub merge: bool,

    // combined reference FASTA for writing a merged CRAM (must contain all contigs of both haplotypes)
    #[arg(long, value_name = "FILE", required = false, help = "combined reference FASTA for merged CRAM output (must contain all contigs of both inputs); required with --merge on CRAM input")]
    pub ref_merged: Option<String>,

    // where to write reads unmapped in both assemblies
    #[arg(short, long, value_name = "DEST", default_value = "asm1", help="where to write reads unmapped in both assemblies: asm1, asm2, or discard")]
    pub unmapped: UnmappedDest,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm1)")]
    pub ref1: Option<String>,

    #[arg(long, value_name = "FILE", required = false, help="reference FASTA for cram file (asm2)")]
    pub ref2: Option<String>,

    // per-base match score from aligner scoring scheme (used in HAPQ calculation)
    #[arg(short = 'A' , long, value_name = "FLOAT", help = "per-base match score from aligner scoring scheme (auto-estimated from ms:i tags if omitted)")]
    pub match_sc: Option<f32>,

    // skip HAPQ score calculation and hq tag output (for non-haplotype comparisons)
    #[arg(long, default_value_t = false, help = "skip HAPQ score calculation and hq tag output (e.g. for comparing GRCh38 vs CHM13)")]
    pub no_hapq: bool,

    // disable writing the list of chromosome-spanning reads
    #[arg(long, default_value_t = false, help = "disable writing hiphap_{s1}_{s2}_span_chrom.fastq")]
    pub no_span_chrom: bool,

    // number of total threads to use
    #[arg(short, long,value_name = "INT", default_value_t = 8, help = "Total thread pool size (min 4). Multiples of 8 recommended for optimal read/write balance.")]
    pub threads: usize
}

#[derive(Debug, Clone, ValueEnum)]
pub enum UnmappedDest {
    Asm1,
    Asm2,
    Discard,
}
