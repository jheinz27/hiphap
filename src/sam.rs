use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::ffi::CString;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::Peekable;
use std::path::Path;


use rust_htslib::{
    bam::{self, record::Aux, record::Cigar, Read, Record, Writer},
    errors::Error as BamError,
    htslib,
};

use rand::{thread_rng, Rng};

use crate::cli::Cli;

// estimate  the minimap2 `-A` (Match score) parameter from alignment file.
// Samples ~1 in 10,000 primary alignments until 10 reads are sampledd
// returns ceiling of the maximum ms / alignment_length as -A estimate
//claude implemented (checked )
pub fn estimate_minimap2_a(bam_path: &str, reference: Option<&str>, threads: usize) -> Result<i32, Box<dyn std::error::Error>> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| format!("Failed to open '{}' for -A estimation: {}. Set -A/--match_sc explicitly.", bam_path, e))?;

    //runs before full file passes 
    reader.set_threads(max(1, threads))
        .map_err(|e| format!("Failed to set threads for -A estimation on '{}': {}", bam_path, e))?;

    //if reference provided (CRAM), apply
    if let Some(refpath) = reference {
        reader.set_reference(refpath)
            .map_err(|e| format!("Failed to set reference for -A estimation on '{}': {}. Set -A/--match_sc explicitly.", bam_path, e))?;
    }

    let mut record = Record::new();
    let mut rng = thread_rng();
    let mut max_ratio: f64 = 0.0;
    let mut sampled: u32 = 0;

    //parse one record at a time looking for primary reads
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("Error reading '{}' during -A estimation: {}. Set -A/--match_sc explicitly.", bam_path, e))?;

        //skip unmapped/secondary/supplementary alignments
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        //~1 in 10,000 random sampling
        if !rng.gen_bool(0.0001) { continue; }

        //alignment length from CIGAR 
        let aln_len = get_alignment_len(&record);
        if aln_len == 0 { continue; }

        //extract ms:i tag (integer width varies across files)
        let ms_score: i64 = match record.aux(b"ms") {
            Ok(Aux::I8(v))  => v as i64,
            Ok(Aux::I16(v)) => v as i64,
            Ok(Aux::I32(v)) => v as i64,
            Ok(Aux::U8(v))  => v as i64,
            Ok(Aux::U16(v)) => v as i64,
            Ok(Aux::U32(v)) => v as i64,
            _ => continue,
        };

        let ratio = ms_score as f64 / aln_len as f64;
        if ratio > max_ratio { max_ratio = ratio; }

        sampled += 1;
        if sampled >= 10 { break; }
    }

    if sampled == 0 || max_ratio <= 0.0 {
        return Err(format!(
            "Could not estimate minimap2 -A from '{}': no informative sampled reads with valid ms:i tags found. \
             Set -A/--match_sc explicitly.", bam_path
        ).into());
    }

    Ok(max_ratio.ceil() as i32)
}

/// Helper function to peek at the file format using c path
fn get_format_from_path<P: AsRef<Path>>(path: P) -> Result<bam::Format, Box<dyn std::error::Error>> {
    let path_str = path.as_ref().to_str().ok_or("Invalid UTF-8 path")?;
    let c_path = CString::new(path_str)
        .map_err(|_| format!("Invalid path (contains null byte): {}", path_str))?;

    unsafe {
        let hts_file = htslib::hts_open(c_path.as_ptr(), c"r".as_ptr());
        if hts_file.is_null() {
            return Err(format!("Could not open file: {}", path_str).into());
        }
        let format_struct = (*hts_file).format;
        htslib::hts_close(hts_file);

        // Map the C format to the Rust enum
        match format_struct.format {
            htslib::htsExactFormat_bam => Ok(bam::Format::Bam),
            htslib::htsExactFormat_cram => Ok(bam::Format::Cram),
            htslib::htsExactFormat_sam => Ok(bam::Format::Sam),
            _ => {
                //a PAF handed in without --paf flag is the common mistake (at least for the author)
                let name = path_str.to_ascii_lowercase();
                let name = name.strip_suffix(".gz").or_else(|| name.strip_suffix(".bgz")).unwrap_or(&name);
                if name.ends_with(".paf") {
                    Err(format!(
                        "Unsupported or unknown file format for: {} (this looks like a PAF file: rerun with --paf)",
                        path_str
                    ).into())
                } else {
                    Err(format!("Unsupported or unknown file format for: {}", path_str).into())
                }
            }
        }
    }

}

//function to check both input files are of same type
fn formats_equal(a: &bam::Format, b: &bam::Format) -> bool {
    matches!(
        (a, b),
        (bam::Format::Bam, bam::Format::Bam)
            | (bam::Format::Cram, bam::Format::Cram)
            | (bam::Format::Sam, bam::Format::Sam)
    )
}

//main logic 
pub fn process_sam(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {
    //detect format of both input files (i.e sam/cram/bam)
    let asm1_format = get_format_from_path(&args.asm1)
        .map_err(|e| format!("Failed to identify asm1 file format: {}", e))?;
    let asm2_format = get_format_from_path(&args.asm2)
        .map_err(|e| format!("Failed to identify asm2 file format: {}", e))?;

    //ensure both input files are of the same format
    if !formats_equal(&asm1_format, &asm2_format) {
        return Err(format!("Input files must have the same format (found {:?} and {:?})", asm1_format, asm2_format).into());
    }

    //resolve the output format: the input format unless --out-fmt asks for another one
    //(e.g. SAM input written as compressed BAM)
    let out_format = match args.out_fmt {
        None => asm1_format,
        Some(crate::cli::OutFormat::Sam) => bam::Format::Sam,
        Some(crate::cli::OutFormat::Bam) => bam::Format::Bam,
        Some(crate::cli::OutFormat::Cram) => bam::Format::Cram,
    };
    //whether either side of the run is CRAM decides which references are needed
    let input_is_cram = matches!(asm1_format, bam::Format::Cram);
    let out_is_cram = matches!(out_format, bam::Format::Cram);

    // read in both files
    let mut asm1_reader = bam::Reader::from_path(&args.asm1)
        .map_err(|e| format!("Failed to open asm1 file '{}': {}", args.asm1, e))?;
    let mut asm2_reader = bam::Reader::from_path(&args.asm2)
        .map_err(|e| format!("Failed to open asm2 file '{}': {}", args.asm2, e))?;

    //owned header views from both inputs (used for tid to name mapping and merged-header construction)
    let asm1_hdr = asm1_reader.header().to_owned();
    let asm2_hdr = asm2_reader.header().to_owned();

    //pre-compute target name slices once
    let asm1_names = asm1_hdr.target_names();
    let asm2_names = asm2_hdr.target_names();

    //number of @SQ entries in asm1
    let n1 = asm1_hdr.target_count() as i32;
    //offset applied to asm2 reference ids when writing (n1 when merging, 0 with -p)
    let asm2_offset = if args.partition { 0 } else { n1 };

    //get proper file extension for output based on the resolved output format
    let extension = match out_format {
        bam::Format::Bam => ".bam",
        bam::Format::Sam => ".sam",
        bam::Format::Cram => ".cram",
    };
    //the output format follows the input (or --out-fmt), so flag an -o extension that says otherwise
    crate::warn_output_ext_mismatch(args, extension, args.out_fmt.is_some());

    //merged mode validation (merged is the default; -p writes one file per haplotype)
    if !args.partition {
        //--both would put two primary records for one read in a single file, throw error
        if args.both {
            return Err("--both requires -p/--partition (in a merged file a read would get two primary records)".into());
        }
        //the merged header concatenates the @SQ lists, so contig names must be different between assemblies
        let names1: HashSet<&[u8]> = asm1_names.iter().copied().collect();
        if let Some(dup) = asm2_names.iter().find(|n| names1.contains(*n)) {
            return Err(format!(
                "merged output requires unique contig names between the two inputs, but '{}' appears in both. \
                 Pass -p/--partition to write one file per haplotype instead",
                String::from_utf8_lossy(dup)
            ).into());
        }
        //a single CRAM writer needs one combined reference covering all contigs of both haplotypes
        if out_is_cram && args.ref_merged.is_none() {
            return Err("Merged CRAM output requires a combined reference FASTA containing all contigs of both inputs. Use --ref-merged <FILE>".into());
        }
    } else if args.ref_merged.is_some() {
        eprintln!("Warning: --ref-merged is ignored with -p/--partition");
    }

    //make command line for the @PG tag space seperated
    let hiphap_cl: String = std::env::args()
        .map(|a| a.replace(['\t', '\n'], " "))
        .collect::<Vec<_>>()
        .join(" ");

    //assign all output paths
    let (primary_path, secondary_path, span_path) =
        crate::output_paths(args, extension, ".fastq");

    //create output sam writers
    let (mut out_asm1, mut out_asm2): (Writer, Option<Writer>) = match &secondary_path {
        //merged: out_asm1 is the single merged writer and out_asm2 is None
        None => {
            let merged_header = build_merged_header(&asm1_hdr, &asm2_hdr, &hiphap_cl);
            let w = Writer::from_path(&primary_path, &merged_header, out_format)
                .map_err(|e| format!("Failed to create output file '{}': {}", primary_path, e))?;
            (w, None)
        }
        //partitioned: out_asm1/out_asm2 are the per-haplotype writers, copy input headers
        Some(asm2_out_path) => {
            let header_asm1 = header_with_pg(&asm1_hdr, &hiphap_cl);
            let header_asm2 = header_with_pg(&asm2_hdr, &hiphap_cl);
            let w1 = Writer::from_path(&primary_path, &header_asm1, out_format)
                .map_err(|e| format!("Failed to create output file '{}': {}", primary_path, e))?;
            let w2 = Writer::from_path(asm2_out_path, &header_asm2, out_format)
                .map_err(|e| format!("Failed to create output file '{}': {}", asm2_out_path, e))?;
            (w1, Some(w2))
        }
    };

    //resolve CRAM references: any CRAM reader or writer needs its reference FASTA. A CRAM input
    //needs --ref1/--ref2 for the two readers; a partitioned CRAM output needs them for the two
    //writers, while a merged CRAM output needs --ref-merged (checked above). SAM/BAM need none.
    if input_is_cram || out_is_cram {
        //the per-haplotype references are needed by the readers, and by the writers under -p
        let need_haplo_refs = input_is_cram || (out_is_cram && args.partition);
        let r1 = if need_haplo_refs {
            Some(args.ref1.as_deref().ok_or(
                "CRAM input or output, but no reference FASTA for asm1 provided. Use --ref1 <FILE>")?)
        } else {
            None
        };
        let r2 = if need_haplo_refs {
            Some(args.ref2.as_deref().ok_or(
                "CRAM input or output, but no reference FASTA for asm2 provided. Use --ref2 <FILE>")?)
        } else {
            None
        };

        //set reader references when the input is CRAM
        if input_is_cram {
            asm1_reader.set_reference(r1.unwrap())
                .map_err(|e| format!("Failed to set reference for asm1 Reader: {}", e))?;
            asm2_reader.set_reference(r2.unwrap())
                .map_err(|e| format!("Failed to set reference for asm2 Reader: {}", e))?;
        }

        //set writer references when the output is CRAM
        if out_is_cram {
            if !args.partition {
                //a single merged CRAM writer needs one combined reference for both haplotypes
                let rm = args.ref_merged.as_deref().unwrap();
                out_asm1.set_reference(rm)
                    .map_err(|e| format!("Failed to set reference for merged Writer: {}", e))?;
            } else {
                //each per-haplotype file uses its own haplotype's reference
                out_asm1.set_reference(r1.unwrap())
                    .map_err(|e| format!("Failed to set reference for asm1 Writer: {}", e))?;
                out_asm2.as_mut().unwrap().set_reference(r2.unwrap())
                    .map_err(|e| format!("Failed to set reference for asm2 Writer: {}", e))?;
            }
        }
    } else {
        if args.ref1.is_some() { eprintln!("Warning: --ref1 is ignored (neither input nor output is CRAM)"); }
        if args.ref2.is_some() { eprintln!("Warning: --ref2 is ignored (neither input nor output is CRAM)"); }
    }

    //autoestimate match score (-A in minimap2) from MS tag 
    let resolved_match_sc: f32 = if args.no_hapq {
        // skip under --no-hapq 
        args.match_sc.unwrap_or(0.0)
    } else {
        match args.match_sc {
            Some(v) => v,
            None => {
                let a1 = estimate_minimap2_a(&args.asm1, args.ref1.as_deref(), args.resolved_threads())?;
                let a2 = estimate_minimap2_a(&args.asm2, args.ref2.as_deref(), args.resolved_threads())?;
                let est = a1.max(a2);
                eprintln!("Auto-estimated minimap2 -A (--match-score) from files: asm1={}, asm2={}, using={}", a1, a2, est);
                est as f32
            }
        }
    };

    //open side writer for reads whose winning cluster spans multiple chromosomes 
    let mut span_writer: Option<BufWriter<File>> = if args.no_span_chrom {
        //disanble write with --no-span-chrom
        None
    } else {
        Some(BufWriter::new(File::create(&span_path)
            .map_err(|e| format!("Failed to create '{}': {}", span_path, e))?))
    };

    //one writer when merging, one per haplotype with -p; always two readers
    let n_writers = if out_asm2.is_some() { 2 } else { 1 };
    //set threads: for compressed output a writer is weighted 4x a reader when merging and 3x
    //when partitioned; plain SAM text is 1x
    let writer_weight = match out_format {
        bam::Format::Sam => 1,
        bam::Format::Bam | bam::Format::Cram => if n_writers == 1 { 4 } else { 3 },
    };
    let plan = crate::plan_threads(args.resolved_threads(), 2, n_writers, writer_weight);

    //assign threads to each reader and writer
    asm1_reader.set_threads(plan.reader)?;
    asm2_reader.set_threads(plan.reader)?;
    out_asm1.set_threads(plan.writer)?;
    if let Some(out2) = out_asm2.as_mut() {
        out2.set_threads(plan.writer)?;
    }

    //create peakable iterators of each file
    let mut asm1_iter = asm1_reader.records().peekable();
    let mut asm2_iter = asm2_reader.records().peekable();

    //vectors that store all alignments of one read (cluster of alignments)
    //initiallize capacity to 10 to account for supplemental and secondary alignments
    let mut cluster_asm1: Vec<Record> = Vec::with_capacity(10);
    let mut cluster_asm2: Vec<Record> = Vec::with_capacity(10);

    //initialize counts for summary statistics printed to terminal
    let mut count_asm1: u64 = 0;
    let mut count_asm2: u64 = 0;
    let mut count_equal: u64 = 0;
    let mut count_unmapped: u64 = 0;
    //reads that also got their losing haplotype's alignments written under --keep-loser
    let mut count_loser: u64 = 0;

    //initialize summed read lengths (bps) per category
    let mut bases_asm1: u64 = 0;
    let mut bases_asm2: u64 = 0;
    let mut bases_equal: u64 = 0;
    let mut bases_unmapped: u64 = 0;

    //iterate thorugh both files until they are both fully exhausted
    while asm1_iter.peek().is_some() || asm2_iter.peek().is_some() {

        //move forward by one read's alignments for both files
        get_clusters(&mut asm1_iter, &mut cluster_asm1)?;
        get_clusters(&mut asm2_iter, &mut cluster_asm2)?;

        // check for possible errors such as:
        //end of file, empty cluster, clusters don't represent same read in both files
        match (cluster_asm1.first(), cluster_asm2.first()) {

            (None, None) => break,           // end of file reached for both, should occur at same iteration
            (Some(_), None) | (None, Some(_)) => {
                //one file has ended earlier than the other- throw error
                return Err("alignment streams out of sync: one file ended earlier".into());
            }
            (Some(m), Some(p)) => {
                //read ID is not the same in both clusters- throw error
                if m.qname() != p.qname() {
                    return Err(format!(
                        "alignment streams out of sync: asm1={} asm2={}",
                        String::from_utf8_lossy(m.qname()),
                        String::from_utf8_lossy(p.qname()),
                    ).into());
                }
            }
        }

        //full read length for the base level summary statistics
        //must check both incase read is unmapped in one asm
        let read_bases = cluster_read_len(&cluster_asm1).max(cluster_read_len(&cluster_asm2));

        //get cluster with the higher alignment score, returns the Winner enum, HAPQ, and whether
        //the losing cluster scored close enough to the winner to be worth keeping
        let (winner, hapq, loser_close) = compare_clusters(&mut cluster_asm1, &mut cluster_asm2, args, resolved_match_sc)?;
        //goes to out_asm1, the merged writer
        let keep_loser = args.keep_loser && loser_close;

        //logic for which file to write read to given weighted AS comparison output
        match winner {
            //asm1 clear winner, write to the asm1 output
            crate::Winner::Asm1 => {
                count_asm1 += 1;
                bases_asm1 += read_bases;
                write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, Some(1), 0, &mut span_writer, &asm1_names, "asm1")?;
                if keep_loser {
                    count_loser += 1;
                    write_loser_cluster(&mut out_asm1, &mut cluster_asm2, hapq, asm2_offset)?;
                }
            }
            //asm2 clear winner, write to the asm2 output (or merged writer with offset)
            crate::Winner::Asm2 => {
                count_asm2 += 1;
                bases_asm2 += read_bases;
                 //in merge mode out_asm2 is None: asm2 records go to out_asm1 (the merged writer)
                let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                write_winner_cluster(w2, &mut cluster_asm2, hapq, Some(2), asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
                if keep_loser {
                    count_loser += 1;
                    write_loser_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0)?;
                }
            }
            crate::Winner::Both => {
                count_equal += 1;
                bases_equal += read_bases;
                //if --both true, write equal scoring reads to both output files, only works in partition mode (guarded above)

                if args.both {
                    write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, Some(1), 0, &mut span_writer, &asm1_names, "asm1")?;
                    let w2 = out_asm2.as_mut().expect("internal error: --both requires -p/--partition");
                    write_winner_cluster(w2, &mut cluster_asm2, hapq, Some(2), asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
                
                //deterministically randomly assign each tied read to one haplotype
                } else {
                    //hash read name and use last bit value to assign to asm1 or asm2
                    //ensures that assignments will be reproducible
                    match crate::choose_random(cluster_asm1[0].qname()) {
                        crate::Winner::Asm1 => {
                            write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, Some(1), 0, &mut span_writer, &asm1_names, "asm1")?;
                            //the side the hash did not pick is the loser here
                            if keep_loser {
                                count_loser += 1;
                                write_loser_cluster(&mut out_asm1, &mut cluster_asm2, hapq, asm2_offset)?;
                            }
                        }
                        _ => {
                            let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                            write_winner_cluster(w2, &mut cluster_asm2, hapq, Some(2), asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
                            if keep_loser {
                                count_loser += 1;
                                write_loser_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0)?;
                            }
                        }
                    }
                }
            }
             //hapq is None for unmapped reads, so no hq tag is added and no span record is emitted.
            crate::Winner::Unmapped => {
                count_unmapped += 1;
                bases_unmapped += read_bases;
                match args.unmapped {
                    crate::cli::UnmappedDest::Asm1 => {
                        write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, None, 0, &mut span_writer, &asm1_names, "asm1")?;
                    }
                    crate::cli::UnmappedDest::Asm2 => {
                        let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                        write_winner_cluster(w2, &mut cluster_asm2, hapq, None, asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
                    }
                    crate::cli::UnmappedDest::Discard => {}
                }
            }
        }

    }
    // flush span writer 
    if let Some(w) = span_writer.as_mut() {
        w.flush().map_err(|e| format!("Failed to flush '{}': {}", span_path, e))?;
    }

    //print summarry statistics to terminal
    crate::print_summary(
        &args.s1, &args.s2,
        [count_asm1, count_asm2, count_equal, count_unmapped],
        [bases_asm1, bases_asm2, bases_equal, bases_unmapped],
    );
    //print how many reads cleared --loser-frac
    if args.keep_loser {
        crate::print_loser_summary(count_loser, count_asm1 + count_asm2 + count_equal, args.loser_frac);
    }
Ok(())
}

//function to move ahead one read group at a time for SAM/BAM/CRAM
fn get_clusters<I>(records: &mut Peekable<I>, cluster: &mut Vec<Record>)-> Result<(), Box<dyn std::error::Error>>
where
    I: Iterator<Item= Result<Record,BamError>>,
{
    //forget previous cluster
    cluster.clear();

    //access alignment record of next line in iterator if it exist
    let first_record = match records.next() {
        Some(Ok(r)) => r,
        Some(Err(e)) => return Err(Box::new(e)), //throw error if file appears corrupted
        None => return Ok(()), // End of file
    };

    //get read ID of record
    //we cluster any records with the same read ID
    let cur_id = first_record.qname().to_vec();
    //store first record
    cluster.push(first_record);
    //look for further lines with same read ID.
    loop {
        //peek at next line
        let peek_result = records.peek();
        match peek_result {
            //Next record is valid
            Some(Ok(next_rec)) => {
                //check if next record has same read ID
                if next_rec.qname() == cur_id.as_slice() {
                    // next record belongs to this cluster, consume and add to cluster
                    let rec = records.next().unwrap().unwrap();
                    cluster.push(rec);
                } else {
                    // Belongs to the next cluster.
                    break;
                }
            },
            // Next record is corrupt
            Some(Err(_)) => {
                let err = records.next().unwrap().unwrap_err();
                return Err(Box::new(err));
            },
            //end of file
            None => break,
        }
    }
    //mutated cluster vector in place
    Ok(())
}

//helper function to get weighted score of reads using a specified tag (AS or ms)
//supplental alignments read segments may have overlapping alignments in read coords
//want to take average alignment score for every base in the read to determine total score
fn get_weighted_score(cur_clust : &mut Vec<Record>, tag: &[u8]) -> Result<(f32, u32), Box<dyn std::error::Error>> {
    //get read name
    let qname = String::from_utf8_lossy(cur_clust[0].qname()).into_owned();
    let mut sum_alignment_lens = 0;
    let mut sum_alignment_scores = 0;
    let mut n_splits: u32 = 0;
    //store all read intervals mapping anywhere to take union of later (filter out overlapping segments)
    let mut read_intervals: Vec<(u32, u32)> = Vec::with_capacity(cur_clust.len());

    //get full read length from the first non-secondary record's CIGAR
    //sum of all query-consuming ops (M/I/=/X/S/H) gives original read length even for supplementaries
    let read_len: u32 = cluster_read_len(cur_clust) as u32;

    for rec in cur_clust {
        //do not factor secondary alignments into choosing best alignment,
        //but still output them with the cluster 
        if rec.is_secondary() {continue};

        n_splits += 1;

        //get alignment length of this record in read (query) coordinates
        let alen = get_alignment_len(rec);

        sum_alignment_lens += alen;

        //extract alignment score as i32, throw error if tag missing
        //is not the same integer type in every sam file so check every possile type to be robust
        let alignment_score: i32 = match rec.aux(tag) {
            Ok(Aux::I8(v))  => v as i32,
            Ok(Aux::I16(v)) => v as i32,
            Ok(Aux::I32(v)) => v,
            Ok(Aux::U8(v))  => v as i32,
            Ok(Aux::U16(v)) => v as i32,
            Ok(Aux::U32(v)) => v as i32,
            _ => return Err(format!("Read '{}' is missing the '{}' tag",
            String::from_utf8_lossy(rec.qname()),
            String::from_utf8_lossy(tag)).into()),
        };

        sum_alignment_scores += alignment_score;

        //get the read (query) coordinates of the start of the alignment
        let read_start = get_query_start(rec);
        read_intervals.push((read_start, read_start + alen))
    }
    //this should not happen, but handle just in case
    if sum_alignment_lens == 0 {
        return Err(format!("Read '{}' has primary alignment length of 0", qname).into());
    }


    //takes the union of read (query) coordinates over all alignment segments for a read
    //returns total read bases aligned in any record, so we can take average over read, without double counting bases
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals);

    //calc weighted alignment score:
    //average alignment score per base across all aligning segments
    // multiplied by unique aligned bases, scaled by coverage fraction of the read
    let cov_fraction = read_bps_aligned as f32 / read_len as f32;
    Ok(((sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32 * cov_fraction, n_splits))

}

//choose which alignment block to keep
//returns the winner, the HAPQ, and whether the losing cluster is close enough to the winner to
//be worth writing under --keep-loser
fn compare_clusters<'a>(clust1:&'a mut Vec<Record>, clust2:&'a mut Vec<Record>, args:&Cli, match_sc: f32) ->  Result<(crate::Winner, Option<u8>, bool), Box<dyn std::error::Error>> {

    //if either cluster is empty there is a file sync issue as every cluster should have at least one record
    if clust1.is_empty() || clust2.is_empty() {
        return Err("Fatal Error: Attempted to compare empty read clusters. This usually indicates a file sync issue.".into());
    }

    //check if read is unmapped in either or both files
    let unmappeds = (clust1[0].is_unmapped(), clust2[0].is_unmapped());

    //handle unmapped read cases

    //the losing side of these cases holds nothing but an unmapped record, so there is never a
    //losing alignment to keep: the trailing false turns --keep-loser off for all three
    match unmappeds {
        (true, true) => { return Ok((crate::Winner::Unmapped, None, false)); }, //unmapped in both
        //if read only maps to one hap then that hap is the winner
        (true, false) => return Ok((crate::Winner::Asm2, if args.no_hapq { None } else { Some(60u8) }, false)), //  mapped in asm2
        (false, true) => return Ok((crate::Winner::Asm1, if args.no_hapq { None } else { Some(60u8) }, false)), //  mapped in asm1
        _ => {} //mapped in both continue to check below
    }

    //determine what field we are using to compare alignment score
    //default is using alignment score (AS:i:) but using ms:i: can be set by user wiht --ms
    let tag: &[u8] = if args.ms { b"ms" } else { b"AS" };

    //get score and number of non-secondary alignment segments for each cluster
    let (score1, n_splits1) = get_weighted_score(clust1, tag)?;
    let (score2, n_splits2) = get_weighted_score(clust2, tag)?;

    //return respective winner depending on which AS is higher,
    //both is a special case that can be determined by user input
    if score1 > score2 {
        let hapq = if args.no_hapq { None } else { Some(crate::compute_hapq(score1, score2, n_splits1, match_sc)) };
        Ok((crate::Winner::Asm1, hapq, crate::loser_is_close(score1, score2, args.loser_frac)))
    } else if score1 < score2 {
        let hapq = if args.no_hapq { None } else { Some(crate::compute_hapq(score2, score1, n_splits2, match_sc)) };
        Ok((crate::Winner::Asm2, hapq, crate::loser_is_close(score2, score1, args.loser_frac)))
    } else {
        let hapq = if args.no_hapq { None } else { Some(0u8) };
        //for a tie keep_loser bool always true
        Ok((crate::Winner::Both, hapq, true))
    }
}

//function to get the full read length for a whole cluster of alignments
//check the CIGAR of the first non-secondary record, which recovers the full length
//unmapped records have no CIGAR, so fall back to the length of the SEQ field
fn cluster_read_len(cluster: &[Record]) -> u64 {
    for rec in cluster.iter() {
        if !rec.is_secondary() {
            let rlen = get_read_len(rec);
            if rlen > 0 { return rlen as u64; }
            return rec.seq_len() as u64;
        }
    }
    0
}

//function to get full original read length from CIGAR string
//sums all query-consuming operations: M/I/=/X/S/H
fn get_read_len(rec: &Record) -> u32 {
    let mut rlen = 0;
    for c in rec.cigar().iter() {
        match c {
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l)
            | Cigar::SoftClip(l) | Cigar::HardClip(l) => { rlen += *l },
            _ => {}
        }
    }
    rlen
}


//function to get query span of aligned seqment
fn get_alignment_len(rec: &Record) -> u32  {
    let mut qlen = 0;
    //parse cigar string to determine total aligned query length
    for c in rec.cigar().iter() {
        match c {
            //these fields consume query (per https://samtools.github.io/hts-specs/SAMv1.pdf)
            //hard clip or soft clip consume query coordinates, but does not count towards alignment of query
            Cigar::Match(l) | Cigar::Ins(l) | Cigar::Equal(l) | Cigar::Diff(l) => {qlen += *l},
            _ => {}
        }
    }
    qlen
}

//function to get start of query span from cigar string
fn get_query_start(rec: &Record) -> u32 {
    let cigar = rec.cigar();
    //if record is aligned in the reverse direction, take right clip
    if rec.is_reverse() {
        let mut right = 0;
        for c in cigar.iter().rev() {
            match *c {
                //hard clip or soft clip consume read_coordinates, but does not count towards alignment
                //so alignemnt begins after we get through clipped sequence
                Cigar::HardClip(l) | Cigar::SoftClip(l) => right += l,
                _ => break,
            }
        }
        right

    //if record is aligned in the forward direction, take left clip
    } else {
        let mut left = 0;
        for c in cigar.iter() {
            match *c {
                //same logic as for right clip
                Cigar::HardClip(l) | Cigar::SoftClip(l) => left += l,
                _ => break,
            }
        }
        left
    }


}


//complement a single base, preserving case; unknown/ambiguity codes -> N
fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T', b'T' => b'A', b'C' => b'G', b'G' => b'C',
        b'a' => b't', b't' => b'a', b'c' => b'g', b'g' => b'c',
        b'N' | b'n' => b,
        _ => b'N',
    }
}

//function to reverse complement sequece and quality vals to have 
//output fastq of chromosome spannign reads be in the original oritnetation 
//as the input fasta file to the aligner
fn oriented_seq_qual(rec: &Record) -> (Vec<u8>, Vec<u8>) {
    let mut seq = rec.seq().as_bytes();
    let mut qual = rec.qual().to_vec();
    if rec.is_reverse() {
        //rev comp sequence to get original read orientation
        seq.reverse();
        for b in seq.iter_mut() { *b = complement_base(*b); }
        //reverse quality score to match reversed bases
        qual.reverse();
    }
    (seq, qual)
}


//write one FASTQ record for a read whose winning cluster spans multiple chromosomes:
//claude code assisted  (checked )
fn emit_span_fastq( w: &mut Option<BufWriter<File>>, qname: &[u8], seq: &[u8], qual: &[u8], tids: &[i32], names: &[&[u8]],label: &str,
) -> std::io::Result<()> {
    if let Some(file) = w {
        let q = std::str::from_utf8(qname).unwrap_or("?");
        //if primary record has no stored sequence, skip with warning
        if seq.is_empty() {
            eprintln!("warning: chrom-spanning read '{}' has no stored sequence; skipping FASTQ record", q);
            return Ok(());
        }

        //get unique chromosomes that the read spans
        let chroms: Vec<&str> = tids.iter()
            .filter(|&&t| t >= 0)
            .map(|&t| names.get(t as usize)
                .and_then(|n| std::str::from_utf8(n).ok())
                .unwrap_or("?"))
            .collect();

        // append the asm label and chrom list to header
        writeln!(file, "@{}\t{}\t{}", q, label, chroms.join(","))?;
        file.write_all(seq)?;
        writeln!(file)?;
        writeln!(file, "+")?;

      
        let missing = qual.is_empty() || qual[0] == 0xFF;
        if missing {
            //write placeholder phred values if quality scores missing
            eprintln!("warning: chrom-spanning read '{}' has no quality scores; writing placeholder Phred-0 qualities", q);
            let placeholder = vec![b'!'; seq.len()];
            file.write_all(&placeholder)?;
        } else {
            //convert phred to ascii representation
            let ascii: Vec<u8> = qual.iter().map(|&p| p + 33).collect();
            file.write_all(&ascii)?;
        }
        writeln!(file)?;
    }
    Ok(())
}

//collect every @PG ID declared in SAM header bytes `text`.
//claude assisted (checked)
fn collect_pg_ids(text: &[u8]) -> HashSet<Vec<u8>> {
    let mut ids = HashSet::new();
    for line in text.split(|&b| b == b'\n') {
        if !line.starts_with(b"@PG\t") { continue; }
        for field in line.split(|&b| b == b'\t') {
            if let Some(v) = field.strip_prefix(b"ID:".as_slice()) {
                ids.insert(v.to_vec());
            }
        }
    }
    ids
}

//leaf of the @PG chain in `text`: the most recently declared ID that no @PG references via PP,
//so a freshly added program links onto the end of the chain. None when `text` has no @PG lines.
//claude assisted (checked)
fn pg_chain_leaf(text: &[u8]) -> Option<Vec<u8>> {
    let mut id_order: Vec<Vec<u8>> = Vec::new();
    let mut referenced: HashSet<Vec<u8>> = HashSet::new();
    for line in text.split(|&b| b == b'\n') {
        if !line.starts_with(b"@PG\t") { continue; }
        for field in line.split(|&b| b == b'\t') {
            if let Some(v) = field.strip_prefix(b"ID:".as_slice()) {
                if !id_order.iter().any(|x| x == v) { id_order.push(v.to_vec()); }
            } else if let Some(v) = field.strip_prefix(b"PP:".as_slice()) {
                referenced.insert(v.to_vec());
            }
        }
    }
    id_order.into_iter().rev().find(|i| !referenced.contains(i))
}

//append asm2's @PG lines to `text`, renaming any IDs that collide with IDs already present
//(asm1's) and rewriting asm2-internal PP references to the renamed IDs. 
//claude assisted (checked)
fn append_asm2_pg(text: &mut Vec<u8>, asm2_hdr: &bam::HeaderView) {
    if text.last().is_some_and(|&b| b != b'\n') {
        text.push(b'\n');
    }
    let asm2_bytes = asm2_hdr.as_bytes();
    let mut used = collect_pg_ids(text);

    //first pass: map each colliding asm2 ID to a fresh "<id>-N" name (keeps non-colliding IDs)
    let mut rename: HashMap<Vec<u8>, Vec<u8>> = HashMap::new();
    for line in asm2_bytes.split(|&b| b == b'\n') {
        if !line.starts_with(b"@PG\t") { continue; }
        for field in line.split(|&b| b == b'\t') {
            if let Some(v) = field.strip_prefix(b"ID:".as_slice()) {
                let old = v.to_vec();
                if used.contains(&old) {
                    let mut k = 1;
                    let mut newid = format!("{}-{}", String::from_utf8_lossy(&old), k).into_bytes();
                    while used.contains(&newid) {
                        k += 1;
                        newid = format!("{}-{}", String::from_utf8_lossy(&old), k).into_bytes();
                    }
                    used.insert(newid.clone());
                    rename.insert(old, newid);
                } else {
                    used.insert(old);
                }
            }
        }
    }

    //second pass: emit each asm2 @PG line, applying the rename map to its ID and PP fields
    for line in asm2_bytes.split(|&b| b == b'\n') {
        if !line.starts_with(b"@PG\t") { continue; }
        let mut first = true;
        for field in line.split(|&b| b == b'\t') {
            if !first { text.push(b'\t'); }
            first = false;
            if let Some(v) = field.strip_prefix(b"ID:".as_slice()) {
                text.extend_from_slice(b"ID:");
                text.extend_from_slice(rename.get(v).map_or(v, |n| n.as_slice()));
            } else if let Some(v) = field.strip_prefix(b"PP:".as_slice()) {
                text.extend_from_slice(b"PP:");
                text.extend_from_slice(rename.get(v).map_or(v, |n| n.as_slice()));
            } else {
                text.extend_from_slice(field);
            }
        }
        text.push(b'\n');
    }
}

//append a @PG line recording this hiphap run to SAM header bytes `text`.
//`pp` is the previous-program ID to chain onto 
//claude assisted (checked)
fn append_hiphap_pg(text: &mut Vec<u8>, cl: &str, pp: Option<&[u8]>) {
    if text.last().is_some_and(|&b| b != b'\n') {
        text.push(b'\n');
    }
    let ids = collect_pg_ids(text);
    //unique ID: hiphap, else hiphap.1, hiphap.2, ...
    let mut id = b"hiphap".to_vec();
    let mut n = 1;
    while ids.contains(&id) {
        id = format!("hiphap.{}", n).into_bytes();
        n += 1;
    }
    let mut pg = format!("@PG\tID:{}\tPN:HipHap\tVN:{}",
        String::from_utf8_lossy(&id), env!("CARGO_PKG_VERSION"));
    if let Some(pp) = pp {
        pg.push_str(&format!("\tPP:{}", String::from_utf8_lossy(pp)));
    }
    pg.push_str(&format!("\tCL:{}\n", cl));
    text.extend_from_slice(pg.as_bytes());
}

//round-trip a single input header through htslib, appending this run's hiphap @PG line.
//claude assisted (checked)
fn header_with_pg(hdr: &bam::HeaderView, cl: &str) -> bam::Header {
    let mut text: Vec<u8> = hdr.as_bytes().to_vec();
    let pp = pg_chain_leaf(&text);
    append_hiphap_pg(&mut text, cl, pp.as_deref());
    let view = bam::HeaderView::from_bytes(&text);
    bam::Header::from_template(&view)
}

//build a merged output header from both inputs,
//asm1 keeps tids 0..n1 and asm2's contigs are appended, taking tids n1..n1+n2
//claude assisted (checked)
fn build_merged_header(asm1_hdr: &bam::HeaderView, asm2_hdr: &bam::HeaderView, cl: &str) -> bam::Header {
    let asm1_bytes = asm1_hdr.as_bytes();

    //bucket asm1's header lines by type so we can re-emit them grouped instead of interleaved
    let mut hd: Option<&[u8]> = None;
    let mut sq_lines: Vec<&[u8]> = Vec::new();
    let mut pg_lines: Vec<&[u8]> = Vec::new();
    let mut other_lines: Vec<&[u8]> = Vec::new(); 
    for line in asm1_bytes.split(|&b| b == b'\n') {
        if line.is_empty() { continue; }
        if line.starts_with(b"@HD") { hd = Some(line); }
        else if line.starts_with(b"@SQ\t") { sq_lines.push(line); }
        else if line.starts_with(b"@PG\t") { pg_lines.push(line); }
        else { other_lines.push(line); }
    }

    let mut text: Vec<u8> = Vec::with_capacity(asm1_bytes.len() + asm2_hdr.as_bytes().len());
    let push_line = |text: &mut Vec<u8>, line: &[u8]| { text.extend_from_slice(line); text.push(b'\n'); };

    //@HD first if present
    if let Some(h) = hd { push_line(&mut text, h); }
    //all @SQ lines grouped: asm1's first (tids 0..n1), then asm2's (tids n1..) — order defines tids
    for l in &sq_lines { push_line(&mut text, l); }
    for line in asm2_hdr.as_bytes().split(|&b| b == b'\n') {
        if line.starts_with(b"@SQ\t") { push_line(&mut text, line); }
    }
    //asm1's non-@HD/@SQ/@PG lines (@RG, @CO, ...) preserved, before the @PG block
    for l in &other_lines { push_line(&mut text, l); }
    //@PG block, grouped at the end: asm1's chain first
    for l in &pg_lines { push_line(&mut text, l); }
    //capture asm1's @PG leaf before adding asm2's chain, so the hiphap @PG links onto asm1's chain
    let asm1_leaf = pg_chain_leaf(&text);
    //carry over asm2's @PG provenance (renaming colliding IDs) so its minimap2 reference is recorded
    append_asm2_pg(&mut text, asm2_hdr);
    //record this hiphap run as a @PG line, linked onto asm1's existing @PG chain
    append_hiphap_pg(&mut text, cl, asm1_leaf.as_deref());
    //round-trip the assembled header text through htslib so all @SQ sub-fields are preserved
    let view = bam::HeaderView::from_bytes(&text);
    bam::Header::from_template(&view)
}

//write every record of a winning cluster to assigned file (merged or partitioned) 
fn write_winner_cluster(writer: &mut Writer, cluster: &mut [Record],hapq: Option<u8>,hp: Option<u8>,tid_offset: i32,span_writer: &mut Option<BufWriter<File>>,names: &[&[u8]],label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut seen_tids: Vec<i32> = Vec::with_capacity(4);
    let mut primary_idx: Option<usize> = None;
    for (i, rec) in cluster.iter_mut().enumerate() {
        //track distinct tids of the read's primary/supplementary alignments
        if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
            let t = rec.tid();
            if !seen_tids.contains(&t) { seen_tids.push(t); }
            if !rec.is_supplementary() { primary_idx = Some(i); }
        }
        //add hq (HapQ) and HP (haplotype assignment) tags to record; drop any pre-existing copy
        if let Some(hq) = hapq {
            let _ = rec.remove_aux(b"hq");
            rec.push_aux(b"hq", Aux::U8(hq))?;
        }
        if let Some(h) = hp {
            let _ = rec.remove_aux(b"HP");
            rec.push_aux(b"HP", Aux::U8(h))?;
        }
        //a winner never carries hs, so drop any left over from a previous run on this file
        let _ = rec.remove_aux(b"hs");
        //shift reference ids into the merged header tid coordinates
        if tid_offset != 0 {
            let t = rec.tid();
            if t >= 0 { rec.set_tid(t + tid_offset); }
            let mt = rec.mtid();
            if mt >= 0 { rec.set_mtid(mt + tid_offset); }
        }
        writer.write(rec)?;
    }
    //if the read's winning alignments span more than one chromosome, write it to the span FASTQ
    if seen_tids.len() > 1 {
        if let Some(idx) = primary_idx {
            let rec = &cluster[idx];
            let (seq, qual) = oriented_seq_qual(rec);
            emit_span_fastq(span_writer, rec.qname(), &seq, &qual, &seen_tids, names, label)?;
        }
    }
    Ok(())
}

//write the losing haplotype's alignments to the merged file as secondary records (--keep-loser).
//only the records that were primary or supplementary in their own cluster are kept
//the hs tag records what each record was before the secondary flag was set
//no HP tag: the read was not assigned to this haplotype.
//no span-chrom handling either
//claude assisted (checked)
fn write_loser_cluster(writer: &mut Writer, cluster: &mut [Record], hapq: Option<u8>, tid_offset: i32,
) -> Result<(), Box<dyn std::error::Error>> {
    for rec in cluster.iter_mut() {
        if rec.is_secondary() || rec.is_unmapped() { continue; }

        //remember what this record was, then flag it secondary. set_secondary only ORs in 0x100,
        //so a supplementary keeps 0x800 and comes out as 0x900
        let hs = if rec.is_supplementary() { b'S' } else { b'P' };
        rec.set_secondary();

        //clear SEQ and QUAL: the winning record for this read already holds the read's bases, and
        //writing them twice would roughly double the merged file. Record::set leaves the aux data
        //untouched, and an empty seq/qual pair is written out as '*' for both fields
        let cig = rec.cigar().take();
        let qname = rec.qname().to_vec();
        rec.set(&qname, Some(&cig), &[], &[]);

        //same drop-then-push idiom as the winner path, so a re-run replaces rather than duplicates
        if let Some(hq) = hapq {
            let _ = rec.remove_aux(b"hq");
            rec.push_aux(b"hq", Aux::U8(hq))?;
        }
        //this record's haplotype lost, so any HP it carried into hiphap must not survive
        let _ = rec.remove_aux(b"HP");
        let _ = rec.remove_aux(b"hs");
        rec.push_aux(b"hs", Aux::Char(hs))?;

        //shift reference ids into the merged header tid coordinates
        if tid_offset != 0 {
            let t = rec.tid();
            if t >= 0 { rec.set_tid(t + tid_offset); }
            let mt = rec.mtid();
            if mt >= 0 { rec.set_mtid(mt + tid_offset); }
        }
        writer.write(rec)?;
    }
    Ok(())
}
