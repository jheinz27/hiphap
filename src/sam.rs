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
pub fn estimate_minimap2_a(bam_path: &str, reference: Option<&str>) -> Result<i32, Box<dyn std::error::Error>> {
    let mut reader = bam::Reader::from_path(bam_path)
        .map_err(|e| format!("Failed to open '{}' for -A estimation: {}. Set -A/--match_sc explicitly.", bam_path, e))?;

    //use min required number of threads 
    reader.set_threads(4)
        .map_err(|e| format!("Failed to set threads for -A estimation on '{}': {}", bam_path, e))?;

    //if reference provided (CRAM), apply it
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
            _ => Err(format!("Unsupported or unknown file format for: {}", path_str).into()),
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

    //number of @SQ entries in asm1; in a merged header asm2's contigs are appended after these,
    let n1 = asm1_hdr.target_count() as i32;
    //offset applied to asm2 reference ids when writing (n1 in merge mode, 0 in partitioned mode)
    let asm2_offset = if args.merge { n1 } else { 0 };

    //get proper file extension for output based on input format
    //previos checked that fommats were the same between files
    let extension = match asm1_format {
        bam::Format::Bam => ".bam",
        bam::Format::Sam => ".sam",
        bam::Format::Cram => ".cram",
    };

    //--merge mode validation
    if args.merge {
        //--both would put two primary records for one read in a single file, throw error
        if args.both {
            return Err("--both cannot be combined with --merge (a read would get two primary records in one file)".into());
        }
        //the merged header concatenates the @SQ lists, so contig names must be different between assemblies
        let names1: HashSet<&[u8]> = asm1_names.iter().copied().collect();
        if let Some(dup) = asm2_names.iter().find(|n| names1.contains(*n)) {
            return Err(format!(
                "--merge requires unique contig names between the two inputs, but '{}' appears in both",
                String::from_utf8_lossy(dup)
            ).into());
        }
        //a single CRAM writer needs one combined reference covering all contigs of both haplotypes
        if let bam::Format::Cram = asm1_format {
            if args.ref_merged.is_none() {
                return Err("Merged CRAM output requires a combined reference FASTA containing all contigs of both inputs. Use --ref-merged <FILE>".into());
            }
        }
    } else if args.ref_merged.is_some() {
        eprintln!("Warning: --ref-merged is ignored without --merge");
    }

    //make command line for the @PG tag space seperated
    let hiphap_cl: String = std::env::args()
        .map(|a| a.replace(['\t', '\n'], " "))
        .collect::<Vec<_>>()
        .join(" ");

    //create output writers
    let (mut out_asm1, mut out_asm2): (Writer, Option<Writer>) = if args.merge {
        //in merge mode, out_asm1 is the single merged writer and out_asm2 is None;
        let merged_header = build_merged_header(&asm1_hdr, &asm2_hdr, &hiphap_cl);
        let merged_path = format!("hiphap_{}_{}_merged{}", args.s1, args.s2, extension);
        let w = Writer::from_path(&merged_path, &merged_header, asm1_format)
            .map_err(|e| format!("Failed to create output file '{}': {}", merged_path, e))?;
        (w, None)
    } else {
        //in partitioned mode, out_asm1/out_asm2 are the per-haplotype writers, copy input headers
        let header_asm1 = header_with_pg(&asm1_hdr, &hiphap_cl);
        let header_asm2 = header_with_pg(&asm2_hdr, &hiphap_cl);
        let asm1_out_path = format!("hiphap_{}{}", args.s1, extension);
        let asm2_out_path = format!("hiphap_{}{}", args.s2, extension);
        let w1 = Writer::from_path(&asm1_out_path, &header_asm1, asm1_format)
            .map_err(|e| format!("Failed to create output file '{}': {}", asm1_out_path, e))?;
        let w2 = Writer::from_path(&asm2_out_path, &header_asm2, asm2_format)
            .map_err(|e| format!("Failed to create output file '{}': {}", asm2_out_path, e))?;
        (w1, Some(w2))
    };

    // set cram reference for readers and writers if cram input 
    if let bam::Format::Cram = asm1_format {
        //set readers
        let r1 = args.ref1.as_deref()
            .ok_or("Input format is CRAM, but no reference FASTA for asm1 provided. Use --ref1 <FILE>")?;
        asm1_reader.set_reference(r1)
            .map_err(|e| format!("Failed to set reference for asm1 Reader: {}", e))?;
        let r2 = args.ref2.as_deref()
            .ok_or("Input format is CRAM, but no reference FASTA for asm2 provided. Use --ref2 <FILE>")?;
        asm2_reader.set_reference(r2)
            .map_err(|e| format!("Failed to set reference for asm2 Reader: {}", e))?;

        //set to diploid reference writer if --merge setting 
        if args.merge {
            let rm = args.ref_merged.as_deref().unwrap();
            out_asm1.set_reference(rm)
                .map_err(|e| format!("Failed to set reference for merged Writer: {}", e))?;
        //else set each partitioned file to respectinve reference fasta of input
        } else {
            out_asm1.set_reference(r1)
                .map_err(|e| format!("Failed to set reference for asm1 Writer: {}", e))?;
            out_asm2.as_mut().unwrap().set_reference(r2)
                .map_err(|e| format!("Failed to set reference for asm2 Writer: {}", e))?;
        }
    } else {
        if args.ref1.is_some() { eprintln!("Warning: --ref1 is ignored for non-CRAM input"); }
        if args.ref2.is_some() { eprintln!("Warning: --ref2 is ignored for non-CRAM input"); }
    }

    //autoestimate match score (-A in minimap2) from MS tag 
    let resolved_match_sc: f32 = if args.no_hapq {
        // skip under --no-hapq 
        args.match_sc.unwrap_or(0.0)
    } else {
        match args.match_sc {
            Some(v) => v,
            None => {
                let a1 = estimate_minimap2_a(&args.asm1, args.ref1.as_deref())?;
                let a2 = estimate_minimap2_a(&args.asm2, args.ref2.as_deref())?;
                let est = a1.max(a2);
                eprintln!("Auto-estimated minimap2 -A (--match-score) from files: asm1={}, asm2={}, using={}", a1, a2, est);
                est as f32
            }
        }
    };

    //open side writer for reads whose winning cluster spans multiple chromosomes (unless disabled)
    let span_path = format!("hiphap_{}_{}_span_chrom.fastq", args.s1, args.s2);
    let mut span_writer: Option<BufWriter<File>> = if args.no_span_chrom {
        None
    } else {
        Some(BufWriter::new(File::create(&span_path)
            .map_err(|e| format!("Failed to create '{}': {}", span_path, e))?))
    };

    //set threads
    //if user specifies less than 4, set to 4 (1 thread for each reader and each writer is needed)
    let avail_threads = max(4, args.threads);
    //assign write:reader threads (ideally) 3:1
    let r = max(1, avail_threads / 8);
    //if any additional threads available, assign to writers
    //if num threads is odd, leave one idle
    let w = (avail_threads - (2 * r)) / 2;

    //assign threads to each reader/writer pair
    asm1_reader.set_threads(r)?;
    asm2_reader.set_threads(r)?;
    if let Some(out2) = out_asm2.as_mut() {
        out_asm1.set_threads(w)?;
        out2.set_threads(w)?;
    } else {
        //single merged writer gets the full writer thread budget
        out_asm1.set_threads(2 * w)?;
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

    //iterate thorugh both files until they are both fully exhaused
    while asm1_iter.peek().is_some() || asm2_iter.peek().is_some() {

        //move forward by one read for both files
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

        //get cluster with the higher alignment score, returns the Winner enum and HAPQ
        let (winner, hapq) = compare_clusters(&mut cluster_asm1, &mut cluster_asm2, args, resolved_match_sc)?;

        //logic for which file to write read to given weighted AS comparison output
        match winner {
            //asm1 clear winner, write to the asm1 output
            crate::Winner::Asm1 => {
                count_asm1 += 1;
                write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0, &mut span_writer, &asm1_names, "asm1")?;
            }
            //asm2 clear winner, write to the asm2 output (or merged writer with offset)
            crate::Winner::Asm2 => {
                count_asm2 += 1; 
                 //in merge mode out_asm2 is None: asm2 records go to out_asm1 (the merged writer)
                let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                write_winner_cluster(w2, &mut cluster_asm2, hapq, asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
            }
            crate::Winner::Both => {
                count_equal += 1;
                //if user specifies --both, write equal scoring reads to both output files
                //(--both is rejected together with --merge)
                if args.both {
                    write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0, &mut span_writer, &asm1_names, "asm1")?;
                    let w2 = out_asm2.as_mut().expect("internal error: --both requires partitioned mode");
                    write_winner_cluster(w2, &mut cluster_asm2, hapq, 0, &mut span_writer, &asm2_names, "asm2")?;
                //default behavior is to deterministically randomly assign each tied read to one haplotype
                } else {
                    //hash read name and use last bit value to assign to asm1 or asm2
                    //ensures that assignments will be reproducible
                    match crate::choose_random(cluster_asm1[0].qname()) {
                        crate::Winner::Asm1 => {
                            write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0, &mut span_writer, &asm1_names, "asm1")?;
                        }
                        _ => {
                            let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                            write_winner_cluster(w2, &mut cluster_asm2, hapq, asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
                        }
                    }
                }
            }
            crate::Winner::Unmapped => {
                count_unmapped += 1;
                //hapq is None for unmapped reads, so no hq tag is added and no span record is emitted
                match args.unmapped {
                    crate::cli::UnmappedDest::Asm1 => {
                        write_winner_cluster(&mut out_asm1, &mut cluster_asm1, hapq, 0, &mut span_writer, &asm1_names, "asm1")?;
                    }
                    crate::cli::UnmappedDest::Asm2 => {
                        let w2: &mut Writer = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                        write_winner_cluster(w2, &mut cluster_asm2, hapq, asm2_offset, &mut span_writer, &asm2_names, "asm2")?;
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
    let total = count_asm1 + count_asm2 + count_equal + count_unmapped;
    //avoid NaN% when no reads were parsed (e.g. empty inputs)
    let pct = |n: u64| if total == 0 { 0.0 } else { n as f64 / total as f64 * 100.0 };
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s1, count_asm1, pct(count_asm1));
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s2, count_asm2, pct(count_asm2));
    eprintln!("Reads with equal scores:     {} ({:.1}%)", count_equal, pct(count_equal));
    eprintln!("Reads unmapped to both:      {} ({:.1}%)", count_unmapped, pct(count_unmapped));
    eprintln!("Total reads parsed:                 {}", total);
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
    //mutated cluster vector in place, only need to return result Ok
    Ok(())
}

//helper function to get weighted score of reads using a specified tag (AS or ms)
//for supplental alignments read segments may have overlapping alignments in read coords
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
    let mut read_len: u32 = 0;
    for rec in cur_clust.iter() {
        if !rec.is_secondary() {
            read_len = get_read_len(rec);
            break;
        }
    }

    for rec in cur_clust {
        //do not factor secondary alignments into choosing best alignment,
        // but still output them with the cluster, we don't want to lose them
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
fn compare_clusters<'a>(clust1:&'a mut Vec<Record>, clust2:&'a mut Vec<Record>, args:&Cli, match_sc: f32) ->  Result<(crate::Winner, Option<u8>), Box<dyn std::error::Error>> {

    //if either cluster is empty there is a file sync issue as every cluster should have at least one record
    if clust1.is_empty() || clust2.is_empty() {
        return Err("Fatal Error: Attempted to compare empty read clusters. This usually indicates a file sync issue.".into());
    }

    //check if read is unmapped in either or both files
    let unmappeds = (clust1[0].is_unmapped(), clust2[0].is_unmapped());

    //handle unmapped read cases

    match unmappeds {
        (true, true) => { return Ok((crate::Winner::Unmapped, None)); }, //unmapped in both
        //if read only maps to one hap then that hap is the winner
        (true, false) => return Ok((crate::Winner::Asm2, if args.no_hapq { None } else { Some(60u8) })), //  mapped in asm2
        (false, true) => return Ok((crate::Winner::Asm1, if args.no_hapq { None } else { Some(60u8) })), //  mapped in asm1
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
        Ok((crate::Winner::Asm1, hapq))
    } else if score1 < score2 {
        let hapq = if args.no_hapq { None } else { Some(crate::compute_hapq(score2, score1, n_splits2, match_sc)) };
        Ok((crate::Winner::Asm2, hapq))
    } else {
        let hapq = if args.no_hapq { None } else { Some(0u8) };
        Ok((crate::Winner::Both, hapq))
    }
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
fn write_winner_cluster(writer: &mut Writer, cluster: &mut [Record],hapq: Option<u8>,tid_offset: i32,span_writer: &mut Option<BufWriter<File>>,names: &[&[u8]],label: &str,
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
        //add hq tag to record 
        if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
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
