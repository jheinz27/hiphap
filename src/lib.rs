pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam;
use std::cmp::max;
use std::hash::{Hash, Hasher};
use twox_hash::XxHash64;

//enum to store best alignment of read
pub enum Winner {
    Asm1,
    Asm2,
    Both,
    Unmapped,
}

//if read has identical alignment to both haps,
//chose which hap to report randomly with equal likelihoods
//use last bit of hash of read ID (as bytes) as random assignment
pub fn choose_random(id: &[u8]) -> Winner {
    //XxHash64 provides reproducible assignment bc is deterministic
    let mut hasher = XxHash64::with_seed(42);
    id.hash(&mut hasher);
    if hasher.finish() & 1 == 0 { Winner::Asm1 } else { Winner::Asm2 }
}

//compute haplotype assignment quality (HapQ) score
//confidence measure that a read is assigned to the correct haplotype
pub fn compute_hapq(score_winner: f32, score_loser: f32, n_splits: u32, match_sc: f32) -> u8 {
    
    //scores are identical between haploypes: hapq = 0 
    if score_winner <= score_loser {
        return 0;
    }
    //approximately the difference in matching bases btwn alignment1 and alignment2
    let diff = (score_winner - score_loser) / match_sc;
    //penalize reads with more that 3 split aligments (likely a complex region)
    let pen_split = if n_splits <= 3 { 1.0 } else { 3.0 / n_splits as f32 };
    let score = 6.02 * diff * pen_split;
    //hapq = 0 means that the AS scores were truly identical, special case 
    //a very low hapq score (<1) arising from the emperical penalty gets hapq=1
    (score.clamp(0.0, 60.0) as u8).max(1)
    
}

// function to determing whether the losing cluster scored close enough to the winner 
// to be worth reporting under --keep-loser 
pub fn loser_is_close(score_winner: f32, score_loser: f32, frac: f32) -> bool {
    score_winner > 0.0 && score_loser >= frac * score_winner
}

//helper function to merge any read alignment segments that overlap in read coordinates
//returns count of unique bps of the read contained in any alignment segment
pub fn merge_intervals(intervals: &mut [(u32, u32)]) -> u32 {
    //sort cluster by read start location of alignment segment
    intervals.sort_unstable_by_key(|k| k.0);
    
    let mut read_bps_aligned = 0; 
    if !intervals.is_empty() {
        //initialize at first interval
        let (mut cur_start, mut cur_end) = intervals[0]; 
        //iterate through intervals and merge adjacent overlapping intervals
        for &(next_start, next_end) in intervals.iter().skip(1) { 
            if next_start < cur_end {
                // intervals overlap, so extend
                if next_end > cur_end {
                    cur_end = next_end;
                }
            } else {
                //no further overlap, add length and start over with next interval grouping
                read_bps_aligned += cur_end - cur_start;
                cur_start = next_start;
                cur_end = next_end;
            }
        
        }
        //add final overlap segment
        read_bps_aligned += cur_end - cur_start
    }
    read_bps_aligned
}

//strip a trailing alignment-file extension so -o can seed derived file names
//(the span-chrom file), and so the extension can be compared against the output format
pub fn strip_aln_ext(path: &str) -> &str {
    for e in [".sam", ".bam", ".cram", ".paf"] {
        if let Some(stem) = path.strip_suffix(e) {
            return stem;
        }
    }
    path
}

//resolve every output path for both modes
pub fn output_paths(args: &Cli, ext: &str, span_ext: &str) -> (String, Option<String>, String) {
    //the span file follows the alignment output, so with -o it lands in the same directory
    let span_stem = match &args.output {
        Some(o) => strip_aln_ext(o).to_string(),
        None => format!("hiphap_{}_{}", args.s1, args.s2),
    };
    let span = format!("{}_span_chrom{}", span_stem, span_ext);

    if args.partition {
        //-o is the stem the two per-haplotype files share; Strip any extension the user supplied
        let stem = match &args.output {
            Some(o) => strip_aln_ext(o).to_string(),
            None => "hiphap".to_string(),
        };
        (
            format!("{}_{}{}", stem, args.s1, ext),
            Some(format!("{}_{}{}", stem, args.s2, ext)),
            span,
        )
    } else {
        //-o names the merged file outright
        let merged = args
            .output
            .clone()
            .unwrap_or_else(|| format!("hiphap_{}_{}_merged{}", args.s1, args.s2, ext));
        (merged, None, span)
    }
}

//the output format follows the input format (or an explicit --out-fmt), so an -o extension that
//disagrees is silently ignored; say so rather than leaving SAM text in a file named .bam
//`overridden` is true when --out-fmt set the format, which changes where it came from
pub fn warn_output_ext_mismatch(args: &Cli, ext: &str, overridden: bool) {
    let Some(o) = &args.output else { return };
    let stem = strip_aln_ext(o);
    //no recognised extension to disagree with
    if stem.len() == o.len() {
        return;
    }
    let given = &o[stem.len()..];
    if !given.eq_ignore_ascii_case(ext) {
        let fmt = ext.trim_start_matches('.').to_uppercase();
        let source = if overridden { "set by --out-fmt" } else { "taken from the input" };
        eprintln!(
            "Warning: output format is {} ({}); the '{}' extension of '{}' is ignored",
            fmt, source, given, o
        );
    }
}

//how a thread budget is divided between the htslib readers and writers
pub struct ThreadPlan {
    //threads for each reader (there are always two, one per assembly)
    pub reader: usize,
    //threads for each writer (the single merged writer, or one per haplotype with -p)
    pub writer: usize,
}

//divide the --threads budget between the readers and the writers
//claude assisted (checked)
pub fn plan_threads(
    requested: usize,
    n_readers: usize,
    n_writers: usize,
    writer_weight: usize,
) -> ThreadPlan {
    debug_assert!(n_readers > 0 && n_writers > 0 && writer_weight > 0);

    //every reader and every writer needs a thread of its own
    let total = max(requested, n_readers + n_writers);
    let unit = n_readers + n_writers * writer_weight;

    //the reader share of the budget, rounded half up so readers grow smoothly with the budget
    //instead of only stepping at exact multiples of unit
    let mut reader = (2 * total + unit) / (2 * unit);
    //but never so much that a writer is left with nothing
    reader = reader.min((total - n_writers) / n_readers).max(1);
    let mut writer = (total - n_readers * reader) / n_writers;

    //for compressed output the writer is the expensive side, so it must never end up with fewer
    //threads than a reader Plain SAM is left alone
    while writer_weight > 1 && writer < reader && reader > 1 {
        reader -= 1;
        writer = (total - n_readers * reader) / n_writers;
    }

    ThreadPlan { reader, writer }
}

//print end of run summary statistics to stderr, shared by the SAM and PAF paths

pub fn print_summary(s1: &str, s2: &str, counts: [u64; 4], bases: [u64; 4]) {
    let total: u64 = counts.iter().sum();
    let total_bases: u64 = bases.iter().sum();
    //avoid NaN% when nothing was parsed (e.g. empty inputs)
    let pct = |n: u64, denom: u64| if denom == 0 { 0.0 } else { n as f64 / denom as f64 * 100.0 };

    //build labels up front so the number columns line up for any sample name length
    let labels = [
        format!("Reads aligned better to {}:", s1),
        format!("Reads aligned better to {}:", s2),
        "Reads with equal scores:".to_string(),
        "Reads unmapped to both:".to_string(),
        "Total reads parsed:".to_string(),
    ];
    let lw = labels.iter().map(|l| l.len()).max().unwrap_or(0);
    //report bases as gigabases
    let gbp = |n: u64| format!("{}.{:09}", n / 1_000_000_000, n % 1_000_000_000);
    let cw = total.to_string().len();
    let bw = gbp(total_bases).len();

    for (i, label) in labels.iter().take(4).enumerate() {
        eprintln!("{:<lw$} {:>cw$} reads ({:>4.1}%) ; {:>bw$} Gbps ({:>4.1}%)",
            label, counts[i], pct(counts[i], total), gbp(bases[i]), pct(bases[i], total_bases));
    }
    //pad where the percentages would be so the totals line up with the rows above
    eprintln!("{:<lw$} {:>cw$} reads          ; {:>bw$} Gbps", labels[4], total, gbp(total_bases));
}

//report how many assigned reads also had their losing haplotype's alignments written under --keep-loser
pub fn print_loser_summary(count_loser: u64, assigned: u64, frac: f32) {
    let pct = if assigned == 0 { 0.0 } else { count_loser as f64 / assigned as f64 * 100.0 };
    eprintln!("Reads with losing-hap alignments kept (>= {:.2} of winner): {} ({:.1}% of assigned reads)",
        frac, count_loser, pct);
}