use clap::Parser;
use hiphap::{Cli, paf, sam};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let start = Instant::now();
    let args = Cli::parse();

    //reject negative or 0 as match score 
    if let Some(v) = args.match_sc {
        if v <= 0.0 || v.is_nan() {
            return Err(format!("--match-sc must be a positive number (got {})", v).into());
        }
    }

    //--loser-frac is a fraction of the winner's score 
    if args.loser_frac < 0.0 || args.loser_frac > 1.0 || args.loser_frac.is_nan() {
        return Err(format!("--loser-frac must be in [0, 1] (got {})", args.loser_frac).into());
    }

    if args.paf {
        if args.threads.is_some() {
            eprintln!("Warning: --threads is ignored in PAF mode");
        }
        if args.ref1.is_some() || args.ref2.is_some() || args.ref_merged.is_some() {
            eprintln!("Warning: --ref1/--ref2/--ref-merged are ignored in PAF mode");
        }
        //PAF output is always plain text; there is no BAM/CRAM equivalent to switch to
        if args.out_fmt.is_some() {
            return Err("--out-fmt is not supported in PAF mode (PAF output is always plain text)".into());
        }
        paf::process_paf(&args)?;
    } else {
        sam::process_sam(&args)?;
    }

    let duration = start.elapsed();
    eprintln!("Time elapsed: {:?}", duration);
    Ok(())
}

