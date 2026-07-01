pub mod cli;
pub use cli::Cli;
pub mod paf;
pub mod sam;
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

//compute haplotype assignment quality (HAPQ) score
//modeled on BWA-MEM's mem_approx_mapq_se (bwamem.c)
//confidence measure that a read is assigned to the correct haplotype
pub fn compute_hapq(score_winner: f32, score_loser: f32, n_splits: u32, match_sc: f32) -> u8 {
    if score_winner <= 0.0 {
        return 0;
    }
    //approximately the difference in matching bases btwn alignment1 and alignment2
    let diff = (score_winner - score_loser) / match_sc;
    //penalize reads with more that 3 split aligments (likely a complex region)
    let pen_split = if n_splits <= 3 { 1.0 } else { 3.0 / n_splits as f32 };
    let score = 6.02 * diff * pen_split;
    score.clamp(0.0, 60.0) as u8
    
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
 