use clap::Parser;
use hiphap::{Cli, paf, sam};
use std::time::Instant;

fn main() -> Result<(), Box<dyn std::error::Error>>  {
    let start = Instant::now();
    let args = Cli::parse(); 
    
    if args.paf {
        if args.threads != 8 {
            eprintln!("Warning: --threads is ignored in PAF mode");
        }
        if args.ref1.is_some() || args.ref2.is_some() {
            eprintln!("Warning: --ref1/--ref2 are ignored in PAF mode");
        }
        paf::process_paf(&args)?;
    } else {
        sam::process_sam(&args)?;
    }

    let duration = start.elapsed();
    eprintln!("Time elapsed: {:?}", duration);
    Ok(())
}

