use std::{
    fs::File,
    io::{self, BufRead, BufReader, Write, BufWriter},
    iter::Peekable,
};

use rand::{thread_rng, Rng};

use crate::cli::Cli;


// estimate  the minimap2 `-A` (Match score) parameter from PAF .
// Samples ~1 in 10,000 mapped reads until 20 reads are sampledd
// and returns the ceiling of the maximum ms / alignment_length  value.
//claude implemented (checked) 
pub fn estimate_minimap2_a_paf(paf_path: &str) -> Result<i32, Box<dyn std::error::Error>> {
    let file = File::open(paf_path)
        .map_err(|e| format!("Failed to open '{}' for -A estimation: {}. Set -A/--match_sc explicitly.", paf_path, e))?;
    let mut reader = BufReader::new(file);

    let mut line = String::new();
    let mut rng = thread_rng();
    let mut max_ratio: f64 = 0.0;
    let mut sampled: u32 = 0;

    loop {
        line.clear();
        let n = reader.read_line(&mut line)
            .map_err(|e| format!("Error reading '{}' during -A estimation: {}. Set -A/--match_sc explicitly.", paf_path, e))?;
        if n == 0 { break; } //EOF

        let trimmed = line.trim_end_matches(['\n', '\r']);
        if trimmed.is_empty() { continue; }

        //skip unmapped: target name (col 6, index 5) is "*"
        let mut fields = trimmed.split('\t');
        let _qname = match fields.next() { Some(v) => v, None => continue };
        let _qlen  = match fields.next() { Some(v) => v, None => continue };
        let qstart_s = match fields.next() { Some(v) => v, None => continue };
        let qend_s   = match fields.next() { Some(v) => v, None => continue };
        let _strand  = match fields.next() { Some(v) => v, None => continue };
        let tname    = match fields.next() { Some(v) => v, None => continue };
        if tname == "*" { continue; }

        //skip 6 more fields (tlen, tstart, tend, matches, alnlen, mapq) to reach tags,
        //then walk the optional tags once looking for both tp:A:S (secondary) and ms:i:<n>.
        //skip secondaries 
        let mut is_secondary = false;
        let mut ms_score: i64 = -1;
        for f in fields.skip(6) {
            if f == "tp:A:S" {
                is_secondary = true;
            } else if let Some(v) = f.strip_prefix("ms:i:") {
                ms_score = v.parse::<i64>().unwrap_or(-1);
            }
        }
        if is_secondary { continue; }

        //~1 in 10,000 random sampling of primary alignments
        if !rng.gen_bool(0.0001) { continue; }

        //alignment length: qend - qstart (matches get_alignment_len semantics from sam.rs)
        let qstart: u32 = match qstart_s.parse() { Ok(v) => v, Err(_) => continue };
        let qend: u32   = match qend_s.parse()   { Ok(v) => v, Err(_) => continue };
        if qend <= qstart { continue; }
        let aln_len = qend - qstart;

        if ms_score < 0 { continue; }

        let ratio = ms_score as f64 / aln_len as f64;
        if ratio > max_ratio { max_ratio = ratio; }

        sampled += 1;
        //slightly higher sample size than in sam version since
        //supplementary alignemnts also have tp:A:P" tag 
        //and are being checked here
        if sampled >= 20 { break; }
    }

    if sampled == 0 || max_ratio <= 0.0 {
        return Err(format!(
            "Could not estimate minimap2 -A from '{}': no informative sampled lines with valid ms:i tags found. \
             Ensure the PAF was produced with `minimap2 -cx ...`, or set -A/--match_sc explicitly.", paf_path
        ).into());
    }

    Ok(max_ratio.ceil() as i32)
}


//write one line to the chrom-spanning file: qname \t chrom1,chrom2,... \t asm_label
fn emit_span_paf(
    w: &mut Option<BufWriter<File>>,
    qname: &str,
    chroms: &[String],
    label: &str,
) -> std::io::Result<()> {
    if let Some(file) = w {
        writeln!(file, "{}\t{}\t{}", qname, chroms.join(","), label)?;
    }
    Ok(())
}


//write every record of a winning cluster to respective writer 
fn write_paf_cluster(writer: &mut BufWriter<File>, cluster: &[String],hq_suffix: &str, span_writer: &mut Option<BufWriter<File>>, label: &str,) -> Result<(), Box<dyn std::error::Error>> {
    let mut seen_chroms: Vec<String> = Vec::with_capacity(4);
    for rec in cluster.iter() {

        if span_writer.is_some() {
            if let Some(tname) = paf_get_chrom(rec) {
                if !seen_chroms.iter().any(|s| s == tname) {
                    seen_chroms.push(tname.to_string());
                }
            }
        }
        //append hq tag
        writeln!(writer, "{}{}", rec, hq_suffix)?;
    }
    //if winning read is to multiple chrs, write read to txt file 
    if seen_chroms.len() > 1 {
        let qname = cluster[0].split('\t').next().unwrap_or("");
        emit_span_paf(span_writer, qname, &seen_chroms, label)?;
    }
    Ok(())
}

//get chrom of alignment if mapped and non secondary
fn paf_get_chrom(rec: &str) -> Option<&str> {
    let mut tname: Option<&str> = None;
    for (i, f) in rec.split('\t').enumerate() {
        match i {
            5 => {
                if f == "*" { return None; } // unmapped
                tname = Some(f);
            }
            i if i >= 12
                //ignore secondary alignments
                && f == "tp:A:S" => { return None; }
            _ => {} 
        }
    }
    tname
}


pub fn process_paf(args: &Cli) -> Result<(), Box<dyn std::error::Error>> {

    //--both would write a read's alignments to the merged file twice; reject the combination
    if args.merge && args.both {
        return Err("--both cannot be combined with --merge".into());
    }

    //resolve match score: user override takes precedence, else auto-estimate from both PAFs
    let resolved_match_sc: f32 = if args.no_hapq {
        // HAPQ is skipped under --no-hapq, so the match score is unused; don't auto-estimate
        args.match_sc.unwrap_or(0.0)
    } else {
        match args.match_sc {
            Some(v) => v,
            None => {
                let a1 = estimate_minimap2_a_paf(&args.asm1)?;
                let a2 = estimate_minimap2_a_paf(&args.asm2)?;
                let est = a1.max(a2);
                eprintln!("Auto-estimated minimap2 -A (--match-score) from PAFs: asm1={}, asm2={}, using={}", a1, a2, est);
                est as f32
            }
        }
    };

    // read in both files
    let file1 = File::open(&args.asm1)
        .map_err(|e| format!("Failed to open asm1 file '{}': {}", args.asm1, e))?;
    let file2 = File::open(&args.asm2)
        .map_err(|e| format!("Failed to open asm2 file '{}': {}", args.asm2, e))?;

    //create peekable iterators of each file (line-by-line for PAF)
    let mut asm1_iter = BufReader::new(file1).lines().peekable();
    let mut asm2_iter = BufReader::new(file2).lines().peekable();

    //create output writer(s): a single merged file with --merge (out_asm2 = None), else one per haplotype.
    let asm1_out_path = if args.merge {
        format!("hiphap_{}_{}_merged.paf", args.s1, args.s2)
    } else {
        format!("hiphap_{}.paf", args.s1)
    };
    let mut out_asm1 = BufWriter::new(File::create(&asm1_out_path)
        .map_err(|e| format!("Failed to create output file '{}': {}", asm1_out_path, e))?);
    let mut out_asm2: Option<BufWriter<File>> = if args.merge {
        None
    } else {
        let asm2_out_path = format!("hiphap_{}.paf", args.s2);
        Some(BufWriter::new(File::create(&asm2_out_path)
            .map_err(|e| format!("Failed to create output file '{}': {}", asm2_out_path, e))?))
    };

    //open side writer for reads whose winning cluster spans multiple chromosomes (unless disabled)
    let span_path = format!("hiphap_{}_{}_span_chrom.txt", args.s1, args.s2);
    let mut span_writer: Option<BufWriter<File>> = if args.no_span_chrom {
        None
    } else {
        Some(BufWriter::new(File::create(&span_path)
            .map_err(|e| format!("Failed to create '{}': {}", span_path, e))?))
    };

    //vectors that store all alignments of one read (cluster of alignments)
    //initialize capacity to 10 to account for supplemental and secondary alignments
    let mut cluster_asm1: Vec<String> = Vec::with_capacity(10);
    let mut cluster_asm2: Vec<String> = Vec::with_capacity(10);

    //initialize counts for summary statistics printed to terminal
    let mut count_asm1: u64 = 0;
    let mut count_asm2: u64 = 0;
    let mut count_equal: u64 = 0;
    let mut count_unmapped: u64 = 0;

    //iterate through both files until they are both exhausted
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
                return Err("PAF streams out of sync: one file ended earlier".into());
            }
            (Some(m), Some(p)) => {
                //read ID is not the same in both clusters- throw error
                let id1 = m.split('\t').next().unwrap_or("");
                let id2 = p.split('\t').next().unwrap_or("");
                if id1 != id2 {
                    return Err(format!(
                        "PAF streams out of sync: asm1={} asm2={}", id1, id2
                    ).into());
                }
            }
        }

        //get cluster with the higher alignment score, returns the Winner enum and HAPQ
        let (winner, hapq) = compare_clusters(&cluster_asm1, &cluster_asm2, args, resolved_match_sc)?;

       
        //format hq tag suffix if hapq mode is active
        let hq_suffix = match hapq {
            Some(hq) => format!("\thq:i:{}", hq),
            None => String::new(),
        };

        //logic for which file to write read to given score comparison output
        match winner {
            //asm1 clear winner, write to the asm1 output
            crate::Winner::Asm1 => {
                count_asm1 += 1; // increment read counter
                write_paf_cluster(&mut out_asm1, &cluster_asm1, &hq_suffix, &mut span_writer, "asm1")?;
            }
            //asm2 clear winner, write to the asm2 output (or merged writer)
            crate::Winner::Asm2 => {
                count_asm2 += 1; // increment read counter
                let w2 = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                write_paf_cluster(w2, &cluster_asm2, &hq_suffix, &mut span_writer, "asm2")?;
            }
            crate::Winner::Both => {
                count_equal += 1; // increment read counter
                //if user specifies --both, write equal scoring reads to both output files
                //--both is rejected together with --merge
                if args.both {
                    write_paf_cluster(&mut out_asm1, &cluster_asm1, &hq_suffix, &mut span_writer, "asm1")?;
                    let w2 = out_asm2.as_mut().expect("internal error: --both requires partitioned mode");
                    write_paf_cluster(w2, &cluster_asm2, &hq_suffix, &mut span_writer, "asm2")?;
                //default behavior is to deterministically assign each tied read to one haplotype
                } else {
                    //hash read name and use last bit value to assign to asm1 or asm2
                    //ensures that assignments will be reproducible
                    let qname = cluster_asm1[0].split('\t').next().unwrap().to_string();
                    match crate::choose_random(qname.as_bytes()) {
                        crate::Winner::Asm1 => {
                            write_paf_cluster(&mut out_asm1, &cluster_asm1, &hq_suffix, &mut span_writer, "asm1")?;
                        }
                        _ => {
                            let w2 = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                            write_paf_cluster(w2, &cluster_asm2, &hq_suffix, &mut span_writer, "asm2")?;
                        }
                    }
                }
            }
            crate::Winner::Unmapped => {
                count_unmapped += 1;
                //--unmapped selects which input's records to emit (hapq is None, so no hq tag, no span)
                match args.unmapped {
                    crate::cli::UnmappedDest::Asm1 => {
                        write_paf_cluster(&mut out_asm1, &cluster_asm1, &hq_suffix, &mut span_writer, "asm1")?;
                    }
                    crate::cli::UnmappedDest::Asm2 => {
                        let w2 = match out_asm2 { Some(ref mut w) => w, None => &mut out_asm1 };
                        write_paf_cluster(w2, &cluster_asm2, &hq_suffix, &mut span_writer, "asm2")?;
                    }
                    crate::cli::UnmappedDest::Discard => {}
                }
            }
        }

    }
    //flush all writers
    out_asm1.flush().map_err(|e| format!("Failed to flush '{}': {}", asm1_out_path, e))?;
    if let Some(w) = out_asm2.as_mut() {
        w.flush().map_err(|e| format!("Failed to flush output: {}", e))?;
    }
    if let Some(w) = span_writer.as_mut() {
        w.flush().map_err(|e| format!("Failed to flush '{}': {}", span_path, e))?;
    }

    //print summary statistics to terminal
    let total = count_asm1 + count_asm2 + count_equal + count_unmapped;
    //avoid NaN% when no reads were parsed (e.g. empty inputs)
    let pct = |n: u64| if total == 0 { 0.0 } else { n as f64 / total as f64 * 100.0 };
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s1, count_asm1, pct(count_asm1));
    eprintln!("Reads aligned better to {}: {} ({:.1}%)", args.s2, count_asm2, pct(count_asm2));
    eprintln!("Reads with equal scores:     {} ({:.1}%)", count_equal, pct(count_equal));
    eprintln!("Reads unmapped to both:      {} ({:.1}%)", count_unmapped, pct(count_unmapped));
    eprintln!("Total reads parsed:          {}", total);
    Ok(())
}

//function to move ahead one read group at a time for PAF
fn get_clusters<I>(lines: &mut Peekable<I>, cluster: &mut Vec<String>)-> io::Result<()>
where
    I: Iterator<Item= Result<String, std::io::Error>>,
{
    //forget previous cluster
    cluster.clear();

    //access alignment record of next line in iterator if it exists
    let first_line = match lines.next() {
        Some(Ok(line)) => line,
        Some(Err(e)) => return Err(e), //throw error if file appears corrupted
        None => return Ok(()),         // End of file
    };

    //get read ID of record (first tab-delimited field in PAF)
    //we cluster any records with the same read ID
    let cur_id = first_line.split_once('\t')
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Malformed PAF line: no tab delimiter"))?.0.to_string();
    //store first record
    cluster.push(first_line);

    //look for further lines with same read ID
    loop {
        //peek at next line
        match lines.peek() {
            //Next record is valid
            Some(Ok(next)) => {
                //check if next record has same read ID
                let next_id = match next.split_once('\t') {
                    Some((id, _)) => id,
                    None => break, // malformed line, let caller handle it
                };
                if next_id == cur_id {
                    // next record belongs to this cluster, consume and add to cluster
                    let line = lines.next().unwrap()?;
                    cluster.push(line);
                } else {
                    // Belongs to the next cluster.
                    break;
                }
            },
            // Next record is corrupt
            Some(Err(_)) => {
                let err = lines.next().unwrap().unwrap_err();
                return Err(err);
            },
            //end of file
            None => break,
        }
    }

    //mutated cluster vector in place, only need to return result Ok
    Ok(())
}


//helper function to get weighted score of split reads using a specified tag (AS or ms)
//weighted_score = (SUM(score) / SUM(Alignment_len)) * tot read_bps_aligned
pub fn get_weighted_score(cur_clust : &Vec<String>, tag_prefix: &str) -> Result<(f32, u32), Box<dyn std::error::Error>> {
    let mut sum_alignment_lens = 0;
    let mut sum_alignment_scores = 0;
    let mut n_splits: u32 = 0;
    let mut read_intervals: Vec<(u32, u32)> = Vec::with_capacity(cur_clust.len());

    //get read length from PAF field 1 (query length) of the first record
    let read_len: u32 = cur_clust[0].split('\t').nth(1)
        .ok_or("Malformed PAF line: missing qlen field")?
        .parse().map_err(|e| format!("Invalid qlen in PAF: {}", e))?;

    for alignment in cur_clust {
        let mut fields = alignment.split('\t');

        //skip to relevant columns
        let _qname = fields.next().ok_or("Malformed PAF line: missing qname field")?;
        let _qlen = fields.next().ok_or("Malformed PAF line: missing qlen field")?;
        let qstart = fields.next()
            .ok_or("Malformed PAF line: missing qstart field")?
            .parse::<u32>().map_err(|e| format!("Invalid qstart in PAF: {}", e))?;
        let qend = fields.next()
            .ok_or("Malformed PAF line: missing qend field")?
            .parse::<u32>().map_err(|e| format!("Invalid qend in PAF: {}", e))?;

        let mut is_secondary = false;
        let mut as_score = 0;
        //find score and tp tags (start from field 12)
        for field in fields.skip(8) { // skip to tags
            if field.starts_with("tp:A:S") { is_secondary = true; } // check for secondary alignment tag
            if let Some(val) = field.strip_prefix(tag_prefix) {
                 as_score = val.parse()
                    .map_err(|e| format!("Invalid {} value in PAF: {}", tag_prefix, e))?;
            }
        }

        //do not factor secondary alignments into choosing best alignment
        if is_secondary { continue; }

        n_splits += 1;

        //should not happen if paf is formatted correctly
        if qend < qstart {
            return Err(format!("Malformed PAF: qend ({}) < qstart ({}) for read '{}'",
                qend, qstart, cur_clust[0].split('\t').next().unwrap_or("unknown")).into());
        }

        sum_alignment_lens += qend - qstart;

        //store read alignment coordinates to merge all overlaps at end
        read_intervals.push((qstart, qend));

        sum_alignment_scores += as_score;

    }

    if sum_alignment_lens == 0 {
        return Err(format!("Read '{}' has primary alignment length of 0",
            cur_clust[0].split('\t').next().unwrap_or("unknown")).into());
    }

    //get total read bases aligned in any record
    let read_bps_aligned = crate::merge_intervals(&mut read_intervals);

    //weighted_score = (SUM(Alignment_Score) / SUM(Alignment_len)) * tot read_bps_aligned * cov_fraction
    let cov_fraction = read_bps_aligned as f32 / read_len as f32;
    Ok(((sum_alignment_scores as f32 / sum_alignment_lens as f32) * read_bps_aligned as f32 * cov_fraction, n_splits))

}

pub fn compare_clusters<'a>(clust1:&'a Vec<String>, clust2:&'a Vec<String>, args: &Cli, match_sc: f32) ->  Result<(crate::Winner, Option<u8>), Box<dyn std::error::Error>> {

    match (clust1[0].split('\t').nth(5), clust2[0].split('\t').nth(5)) {
        (Some("*"), Some("*")) => {return Ok((crate::Winner::Unmapped, None));}, // both reads unmapped
        (Some("*"), _) => return Ok((crate::Winner::Asm2, if args.no_hapq { None } else { Some(60u8) })), // asm1 hap unmapped
        (_, Some("*")) => return Ok((crate::Winner::Asm1, if args.no_hapq { None } else { Some(60u8) })), // asm2 hap unmapped
        _ => {} // continue if mapped to both haps
    }

    let tag_prefix = if args.ms { "ms:i:" } else { "AS:i:" };
    //get score and number of non-secondary alignment segments for each cluster
    let (score1, n_splits1) = get_weighted_score(clust1, tag_prefix)?;
    let (score2, n_splits2) = get_weighted_score(clust2, tag_prefix)?;

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
