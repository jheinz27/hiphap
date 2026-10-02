## Summary

### 1. Verification of the original command

I built hiphap and ran the exact command (with a generated 500k-read SAM test set, since the repo has no fixtures of its own):

```
hiphap -p -1 mat -2 pat --no-span-chrom asm1_alignments.sam asm2_alignments.sam
→ hiphap_mat.sam  hiphap_pat.sam   (exit 0)
```

**Result: confirmed.** With `-p` and `--no-span-chrom`, hiphap writes exactly two files: `hiphap_mat.sam` and `hiphap_pat.sam`, with no span-chrom file.

### 2. Added compressed BAM/CRAM output

Previously the output format was **hard-wired to the input format** (`let extension = match asm1_format` and `Writer::from_path(..., asm1_format)`), so SAM input could only ever produce SAM text. I added a new `-O`/`--out-fmt <sam|bam|cram>` option so any input can be written in a chosen format.

**Files changed:**

- **`src/cli.rs`** — added `OutFormat` enum and the `-O`/`--out-fmt` arg (defaults to input format when omitted).
- **`src/sam.rs`** — resolve `out_format` from the flag (falling back to the input format); drive the output extension, both `Writer::from_path` calls (merged and partitioned), the CRAM reference handling, and the writer thread-weighting from it. CRAM logic was refactored so a CRAM *output* now correctly requires `--ref1`/`--ref2` (`-p`) or `--ref-merged` (merged), independently of whether the *input* is CRAM.
- **`src/main.rs`** — reject `--out-fmt` in PAF mode with a clear error (PAF has no BAM/CRAM equivalent).
- **`src/lib.rs`** — `warn_output_ext_mismatch` now takes an `overridden` flag so the warning reads "set by --out-fmt" vs "taken from the input".
- **`src/paf.rs`** — updated the call site.
- **`README.md`** — documented `-O`/`--out-fmt` in the usage block and added a "Compressed output (BAM/CRAM)" section.

### Behavior now

```bash
hiphap -p -1 mat -2 pat --out-fmt bam asm1_alignments.sam asm2_alignments.sam
# Output: hiphap_mat.bam  hiphap_pat.bam
```

Default behavior is **unchanged** — omitting `-O` still gives SAM from SAM. Verified with `samtools quickcheck`/`view -c` for partitioned & merged BAM, CRAM output, `-o` stem naming, and clear error messages for missing references and PAF misuse. All 8 final checks passed; the only remaining clippy warning is pre-existing and unrelated.
