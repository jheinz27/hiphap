use std::cmp::max;
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


// estimate  the minimap2 `-A` (Match score) parameter from an aligned BAM.
// Samples ~1 in 10,000 mapped reads until 10 reads are sampledd
// and returns the ceiling of the maximum ms / alignment_length  value.
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

    //zero-allocation iteration: reuse single Record buffer
    while let Some(result) = reader.read(&mut record) {
        result.map_err(|e| format!("Error reading '{}' during -A estimation: {}. Set -A/--match_sc explicitly.", bam_path, e))?;

        //skip unmapped/secondary/supplementary alignments
        if record.is_unmapped() || record.is_secondary() || record.is_supplementary() {
            continue;
        }

        //~1 in 10,000 unbiased random sampling
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
        let hts_file = htslib::hts_open(c_path.as_ptr(), b"r\0".as_ptr() as *const i8);
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


fn formats_equal(a: &bam::Format, b: &bam::Format) -> bool {
    match (a, b) {
        (bam::Format::Bam, bam::Format::Bam) => true,
        (bam::Format::Cram, bam::Format::Cram) => true,
        (bam::Format::Sam, bam::Format::Sam) => true,
        _ => false,
    }
}

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

    //Store headers from both input files, as these headers are the same as will be needed in output files
    let header_asm1 = bam::Header::from_template(asm1_reader.header());
    let header_asm2 = bam::Header::from_template(asm2_reader.header());

    //get proper file extension for output based on input format
    let extension = match asm1_format {
        bam::Format::Bam => ".bam",
        bam::Format::Sam => ".sam",
        bam::Format::Cram => ".cram",
    };


    //create writers for both outputs that share user specified prefix
    let asm1_out_path = format!("diplinator_{}{}", args.s1, extension);
    let asm2_out_path = format!("diplinator_{}{}", args.s2, extension);
    //headers are same as in original files, so copy them into output
    let mut out_asm1 = Writer::from_path(&asm1_out_path, &header_asm1, asm1_format)
        .map_err(|e| format!("Failed to create output file '{}': {}", asm1_out_path, e))?;
    let mut out_asm2 = Writer::from_path(&asm2_out_path, &header_asm2, asm2_format)
        .map_err(|e| format!("Failed to create output file '{}': {}", asm2_out_path, e))?;

    //if dealing with a cram file, must set reference fastas and ensure the user provided those
    if let bam::Format::Cram = asm1_format {
        //set fasta reference for both asm1 reader and writer
        if let Some(reference) = &args.ref1 {
            asm1_reader.set_reference(reference)
                .map_err(|e| format!("Failed to set reference for asm1 Reader: {}", e))?;
            out_asm1.set_reference(reference)
                .map_err(|e| format!("Failed to set reference for asm1 Writer: {}", e))?;
        } else {
            //throw error reference fasta was not provided on a cram input
            return Err("Input format is CRAM, but no reference FASTA for asm1 provided. Use --ref1 <FILE>".into());
        }
    } else if args.ref1.is_some() {
        //warn user that asm1 reference will be ignored since the input isn't cram
        eprintln!("Warning: --ref1 is ignored for non-CRAM input");
    }

    //repeat above for the asm2 cram file
    if let bam::Format::Cram = asm2_format {
        if let Some(reference) = &args.ref2 {
            asm2_reader.set_reference(reference)
                .map_err(|e| format!("Failed to set reference for asm2 Reader: {}", e))?;
            out_asm2.set_reference(reference)
                .map_err(|e| format!("Failed to set reference for asm2 Writer: {}", e))?;
        } else {
            return Err("Input format is CRAM, but no reference FASTA for asm2 provided. Use --ref2 <FILE>".into());
        }
    } else if args.ref2.is_some() {
        eprintln!("Warning: --ref2 is ignored for non-CRAM input");
    }

    //resolve match score: user override takes precedence, else auto-estimate from both BAMs
    let resolved_match_sc: f32 = match args.match_sc {
        Some(v) => v,
        None => {
            let a1 = estimate_minimap2_a(&args.asm1, args.ref1.as_deref())?;
            let a2 = estimate_minimap2_a(&args.asm2, args.ref2.as_deref())?;
            let est = a1.max(a2);
            eprintln!("Auto-estimated minimap2 -A (--match-score) from files: asm1={}, asm2={}, using={}", a1, a2, est);
            est as f32
        }
    };

    //open side writer for reads whose winning cluster spans multiple chromosomes (unless disabled)
    //per-run filename mirrors the BAM/SAM output naming so concurrent runs in the same cwd
    //do not overwrite each other's span file
    let span_path = format!("diplinator_{}_{}_span_chrom.txt", args.s1, args.s2);
    let mut span_writer: Option<BufWriter<File>> = if args.no_span_chrom {
        None
    } else {
        Some(BufWriter::new(File::create(&span_path)
            .map_err(|e| format!("Failed to create '{}': {}", span_path, e))?))
    };

    //cache owned header views for tid -> contig name resolution (each reader's records()
    //call below mutably borrows the reader, so we can't reach back to the header inside the loop)
    let asm1_hdr = asm1_reader.header().to_owned();
    let asm2_hdr = asm2_reader.header().to_owned();
    //pre-compute target name slices once; otherwise emit_span would re-walk the header text
    //for every chrom-spanning read
    let asm1_names = asm1_hdr.target_names();
    let asm2_names = asm2_hdr.target_names();

    //set threads
    //if user specifies less than 4, set to 4 (1 thread for each reader and each writer is needed)
    let avail_threads = max(4, args.threads);
    //assign write:reader threads (ideally) 3:1
    let r = max(1, avail_threads / 8);
    //if any additional threads available, assign to writer
    //if num threads is odd, leave one idle
    let w = (avail_threads - (2 * r)) / 2;

    //assign threads to each reader/writer pair
    asm1_reader.set_threads(r)?;
    asm2_reader.set_threads(r)?;
    out_asm1.set_threads(w)?;
    out_asm2.set_threads(w)?;

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

    //iterate thorugh both files until they are both exhaused
    while asm1_iter.peek().is_some() || asm2_iter.peek().is_some() {

        //move forward by one read for both files
        get_clusters(&mut asm1_iter, &mut cluster_asm1)?;
        get_clusters(&mut asm2_iter, &mut cluster_asm2)?;

        // check for possible errors such as:
        //end of file / empty cluster / clusters don't represent same read in both files
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
        let (winner, hapq) = compare_clusters(&mut cluster_asm1, &mut cluster_asm2, &args, resolved_match_sc)?;

        //logic for which file to write read to given score comparison output
        match winner {
            //asm1 clear winner, write to out_asm1
            crate::Winner::Asm1 => {
                count_asm1 += 1; // increment read counter
                let mut seen_tids: Vec<i32> = Vec::with_capacity(4);
                for rec in cluster_asm1.iter_mut() {
                    if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                        let t = rec.tid();
                        if !seen_tids.contains(&t) { seen_tids.push(t); }
                    }
                    if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                    out_asm1.write(rec)?;
                }
                if seen_tids.len() > 1 {
                    emit_span(&mut span_writer, cluster_asm1[0].qname(), &seen_tids, &asm1_names, "asm1")?;
                }
            }
            //asm2 clear winner, write to out_asm2
            crate::Winner::Asm2 => {
                count_asm2 += 1; // increment read counter
                let mut seen_tids: Vec<i32> = Vec::with_capacity(4);
                for rec in cluster_asm2.iter_mut() {
                    if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                        let t = rec.tid();
                        if !seen_tids.contains(&t) { seen_tids.push(t); }
                    }
                    if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                    out_asm2.write(rec)?;
                }
                if seen_tids.len() > 1 {
                    emit_span(&mut span_writer, cluster_asm2[0].qname(), &seen_tids, &asm2_names, "asm2")?;
                }
            }
            crate::Winner::Both => {
                count_equal += 1; // increment read counter
                //if user specifices --both, write equal scoring reads to output files
                if args.both {
                    let mut seen_tids1: Vec<i32> = Vec::with_capacity(4);
                    for rec in cluster_asm1.iter_mut() {
                        if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                            let t = rec.tid();
                            if !seen_tids1.contains(&t) { seen_tids1.push(t); }
                        }
                        if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                        out_asm1.write(rec)?;
                    }
                    if seen_tids1.len() > 1 {
                        emit_span(&mut span_writer, cluster_asm1[0].qname(), &seen_tids1, &asm1_names, "asm1")?;
                    }
                    let mut seen_tids2: Vec<i32> = Vec::with_capacity(4);
                    for rec in cluster_asm2.iter_mut() {
                        if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                            let t = rec.tid();
                            if !seen_tids2.contains(&t) { seen_tids2.push(t); }
                        }
                        if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                        out_asm2.write(rec)?;
                    }
                    if seen_tids2.len() > 1 {
                        emit_span(&mut span_writer, cluster_asm2[0].qname(), &seen_tids2, &asm2_names, "asm2")?;
                    }
                //default behavior is randomly assign equal scoring read to one file
                } else {
                    //hash read name and use last bit value to assign to asm1 or asm2
                    //ensures that assignments will be reproducible
                    match crate::choose_random(cluster_asm1[0].qname()) {
                        crate::Winner::Asm1 => {
                            let mut seen_tids: Vec<i32> = Vec::with_capacity(4);
                            for rec in cluster_asm1.iter_mut() {
                                if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                                    let t = rec.tid();
                                    if !seen_tids.contains(&t) { seen_tids.push(t); }
                                }
                                if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                                out_asm1.write(rec)?;
                            }
                            if seen_tids.len() > 1 {
                                emit_span(&mut span_writer, cluster_asm1[0].qname(), &seen_tids, &asm1_names, "asm1")?;
                            }
                        }
                        _ => {
                            let mut seen_tids: Vec<i32> = Vec::with_capacity(4);
                            for rec in cluster_asm2.iter_mut() {
                                if span_writer.is_some() && !rec.is_secondary() && !rec.is_unmapped() {
                                    let t = rec.tid();
                                    if !seen_tids.contains(&t) { seen_tids.push(t); }
                                }
                                if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                                out_asm2.write(rec)?;
                            }
                            if seen_tids.len() > 1 {
                                emit_span(&mut span_writer, cluster_asm2[0].qname(), &seen_tids, &asm2_names, "asm2")?;
                            }
                        }
                    }
                }
            }
            crate::Winner::Unmapped => {
                count_unmapped += 1;
                match args.unmapped {
                    crate::cli::UnmappedDest::Asm1 => {
                        for rec in cluster_asm1.iter_mut() {
                            if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                            out_asm1.write(rec)?;
                        }
                    }
                    crate::cli::UnmappedDest::Asm2 => {
                        for rec in cluster_asm2.iter_mut() {
                            if let Some(hq) = hapq { rec.push_aux(b"hq", Aux::U8(hq))?; }
                            out_asm2.write(rec)?;
                        }
                    }
                    crate::cli::UnmappedDest::Discard => {}
                }
            }
        }

    }
    //explicitly flush the span writer so disk-full / I/O errors surface as Err instead of
    //being silently swallowed by BufWriter::drop (BAM writers are flushed by htslib at close)
    if let Some(w) = span_writer.as_mut() {
        w.flush().map_err(|e| format!("Failed to flush '{}': {}", span_path, e))?;
    }

    //print summarry statistics to terminal
    let total = count_asm1 + count_asm2 + count_equal + count_unmapped;
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s1, count_asm1, count_asm1 as f64 / total as f64 * 100.0);
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s2, count_asm2, count_asm2 as f64 / total as f64 * 100.0);
    eprintln!("Reads with equal scores:     {} ({:.1}%)", count_equal, count_equal as f64 / total as f64 * 100.0);
    eprintln!("Reads unmapped to both:      {} ({:.1}%)", count_unmapped, count_unmapped as f64 / total as f64 * 100.0);
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
    if sum_alignment_lens <= 0 {
        return Err(format!("Read '{}' has primary alignment length of 0", qname).into());
    }


    //takes the union of read (query) coordinates over all alignment segments for a read
    //returns total read bases aligned in any record, so we can take average over read, without double counting bases
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals);

    //calc weighted alignment score:
    //average alignment score per base across all aligning segments
    // multiplied by unique aligned bases, scaled by coverage fraction of the read
    let cov_fraction = read_bps_aligned as f32 / read_len as f32;
    return Ok(((sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32 * cov_fraction, n_splits));

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
    return qlen;
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
        return right;

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
        return left;
    }


}


//write one line to the chrom-spanning side file: qname \t chrom1,chrom2,... \t asm_label
//chrom names are resolved via the caller-provided target_names slice (computed once,
//not once per spanning event). unmapped tids (-1) and out-of-bounds tids are dropped safely.
//chrom order in the output is insertion order = order of appearance in the cluster
//(i.e. primary first, then supplementaries), which is more informative than alphabetical
fn emit_span(
    w: &mut Option<BufWriter<File>>,
    qname: &[u8],
    tids: &[i32],
    names: &[&[u8]],
    label: &str,
) -> std::io::Result<()> {
    if let Some(file) = w {
        let q = std::str::from_utf8(qname).unwrap_or("?");
        let chroms: Vec<&str> = tids.iter()
            .filter(|&&t| t >= 0)
            .map(|&t| names.get(t as usize)
                .and_then(|n| std::str::from_utf8(n).ok())
                .unwrap_or("?"))
            .collect();
        writeln!(file, "{}\t{}\t{}", q, chroms.join(","), label)?;
    }
    Ok(())
}
