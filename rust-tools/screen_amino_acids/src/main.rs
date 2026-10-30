use clap::{Parser, Subcommand};
use flate2::read::MultiGzDecoder;
use rayon::prelude::*;
use screen_amino_acids::Index;
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

#[derive(Parser)]
#[command(
    version,
    about = "Exact six-frame amino-acid screening; preserves original FASTQ pairs"
)]
struct Args {
    #[command(subcommand)]
    command: Action,
}

#[derive(Subcommand)]
enum Action {
    /// Build an index from masked genome FASTA using standard-code six-frame translation.
    Index {
        #[arg(long)]
        reference: PathBuf,
        #[arg(long)]
        output: PathBuf,
        /// Peptide length, in amino acids (1-12). No scientifically validated default.
        #[arg(long)]
        k: u8,
    },
    /// Partition interleaved FASTQ pairs; either matching mate retains the whole pair.
    Screen {
        #[arg(long)]
        index: PathBuf,
        #[arg(long)]
        reads: PathBuf,
        /// Optional separate R2 FASTQ. Without this, reads must be interleaved.
        #[arg(long)]
        reads2: Option<PathBuf>,
        #[arg(long, default_value_t = 1)]
        threads: usize,
        /// Uncompressed FASTQ output; existing files are never overwritten.
        #[arg(long)]
        matched: PathBuf,
        #[arg(long)]
        unmatched: Option<PathBuf>,
    },
}

fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

fn open_input(path: &Path) -> io::Result<Box<dyn BufRead>> {
    let mut input = BufReader::new(File::open(path)?);
    // Detect gzip by magic, so decompression is not dependent on a filename.
    if input.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(input))))
    } else {
        Ok(Box::new(input))
    }
}

fn contents(line: &[u8]) -> &[u8] {
    let line = line.strip_suffix(b"\n").unwrap_or(line);
    line.strip_suffix(b"\r").unwrap_or(line)
}

fn valid_sequence(sequence: &[u8]) -> bool {
    sequence
        .iter()
        .all(|b| b"ACGTURYSWKMBDHVN".contains(&b.to_ascii_uppercase()))
}

/// Create without overwriting, and remove incomplete files on normal errors.
/// Consumers must wait for successful exit before using these outputs.
struct OutputFile {
    destination: PathBuf,
    writer: BufWriter<File>,
    complete: bool,
}

impl OutputFile {
    fn new(destination: &Path) -> io::Result<Self> {
        if destination.extension().is_some_and(|s| s == "gz") {
            return Err(invalid("outputs are uncompressed; do not use a .gz suffix"));
        }
        let writer = BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create_new(true)
                .open(destination)?,
        );
        Ok(Self {
            destination: destination.to_owned(),
            writer,
            complete: false,
        })
    }
}

impl Drop for OutputFile {
    fn drop(&mut self) {
        if !self.complete {
            let _ = fs::remove_file(&self.destination);
        }
    }
}

fn build_index(reference: &Path, output: &Path, k: u8) -> io::Result<()> {
    let mut index = Index::new(k)?;
    let mut input = open_input(reference)?;
    let mut line = Vec::new();
    let mut sequence = Vec::new();
    let mut records = 0u64;
    while input.read_until(b'\n', &mut line)? != 0 {
        let text = contents(&line);
        if text.starts_with(b">") {
            if text[1..].iter().all(u8::is_ascii_whitespace) {
                return Err(invalid("empty FASTA identifier"));
            }
            if records > 0 {
                if sequence.is_empty() {
                    return Err(invalid("empty FASTA record"));
                }
                index.add_reference(&sequence);
                sequence.clear();
            }
            records += 1;
        } else if !text.is_empty() {
            if records == 0 || !valid_sequence(text) {
                return Err(invalid("invalid FASTA sequence or missing header"));
            }
            sequence.extend_from_slice(text);
        }
        line.clear();
    }
    if sequence.is_empty() {
        return Err(invalid("empty FASTA input or final record"));
    }
    index.add_reference(&sequence);
    let mut output = OutputFile::new(output)?;
    index.write(&mut output.writer)?;
    output.writer.flush()?;
    output.complete = true;
    eprintln!(
        "reference_records={records} amino_k={} unique_kmers={}",
        index.k(),
        index.len()
    );
    Ok(())
}

#[derive(Default)]
struct Record {
    lines: [Vec<u8>; 4],
}

impl Record {
    fn read(&mut self, input: &mut dyn BufRead) -> io::Result<bool> {
        for line in &mut self.lines {
            line.clear();
        }
        if input.read_until(b'\n', &mut self.lines[0])? == 0 {
            return Ok(false);
        }
        for i in 1..4 {
            if input.read_until(b'\n', &mut self.lines[i])? == 0 {
                return Err(invalid("truncated FASTQ record"));
            }
        }
        if !self.lines[0].starts_with(b"@") || !self.lines[2].starts_with(b"+") {
            return Err(invalid("invalid FASTQ header or separator"));
        }
        let sequence = self.sequence();
        let quality = contents(&self.lines[3]);
        if sequence.is_empty()
            || !valid_sequence(sequence)
            || sequence.len() != quality.len()
            || !quality.iter().all(|&q| (33..=126).contains(&q))
        {
            return Err(invalid("invalid FASTQ sequence/quality or unequal lengths"));
        }
        Ok(true)
    }

    fn sequence(&self) -> &[u8] {
        contents(&self.lines[1])
    }

    fn pair_id(&self, expected_mate: u8) -> io::Result<&[u8]> {
        let header = &contents(&self.lines[0])[1..];
        let mut fields = header
            .split(|b| b.is_ascii_whitespace())
            .filter(|s| !s.is_empty());
        let token = fields
            .next()
            .ok_or_else(|| invalid("empty FASTQ identifier"))?;
        let mut id = token;
        if token.ends_with(b"/1") || token.ends_with(b"/2") {
            if token[token.len() - 1] != b'0' + expected_mate {
                return Err(invalid("FASTQ mate suffix does not match pair position"));
            }
            id = &token[..token.len() - 2];
        }
        if id.is_empty() {
            return Err(invalid("empty FASTQ pair identifier"));
        }
        if let Some(comment) = fields.next() {
            // Recognize the standard CASAVA mate field and the ` 1` / ` 2` convention.
            if comment == b"1"
                || comment == b"2"
                || comment.starts_with(b"1:")
                || comment.starts_with(b"2:")
            {
                if comment[0] != b'0' + expected_mate {
                    return Err(invalid("FASTQ mate comment does not match pair position"));
                }
            }
        }
        Ok(id)
    }

    fn write(&self, writer: &mut impl Write) -> io::Result<()> {
        for line in &self.lines {
            writer.write_all(line)?;
        }
        Ok(())
    }
}

fn screen(
    index_path: &Path,
    reads: &Path,
    reads2: Option<&Path>,
    matched: &Path,
    unmatched: Option<&Path>,
    threads: usize,
) -> io::Result<()> {
    if threads == 0 {
        return Err(invalid("threads must be positive"));
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| invalid(e.to_string()))?;
    let index = Index::read(BufReader::new(File::open(index_path)?))?;
    let mut input = open_input(reads)?;
    let mut input2 = reads2.map(open_input).transpose()?;
    let mut out_match = OutputFile::new(matched)?;
    let mut out_nomatch = unmatched.map(OutputFile::new).transpose()?;
    let mut total = 0u64;
    let mut matches = 0u64;
    loop {
        let mut batch = Vec::with_capacity(8192);
        for _ in 0..8192 {
            let mut first = Record::default();
            let mut second = Record::default();
            let has_first = first.read(&mut *input)?;
            let has_second = if let Some(ref mut mate) = input2 {
                second.read(&mut **mate)?
            } else if has_first {
                second.read(&mut *input)?
            } else {
                false
            };
            if has_first != has_second {
                return Err(invalid("FASTQ has an unpaired record"));
            }
            if !has_first {
                break;
            }
            if !first.lines[3].ends_with(b"\n") {
                return Err(invalid(
                    "first mate must end in a newline for lossless interleaving",
                ));
            }
            if first.pair_id(1)? != second.pair_id(2)? {
                return Err(invalid("FASTQ pair identifiers differ"));
            }
            batch.push((first, second));
        }
        if batch.is_empty() {
            break;
        }
        let hits: Vec<bool> = pool.install(|| {
            batch
                .par_iter()
                .map(|(first, second)| {
                    index.matches(first.sequence()) || index.matches(second.sequence())
                })
                .collect()
        });
        for ((first, second), hit) in batch.iter().zip(hits) {
            let output = if hit {
                Some(&mut out_match)
            } else {
                out_nomatch.as_mut()
            };
            if let Some(output) = output {
                first.write(&mut output.writer)?;
                second.write(&mut output.writer)?;
            }
            total += 1;
            matches += u64::from(hit);
        }
    }
    out_match.writer.flush()?;
    if let Some(ref mut output) = out_nomatch {
        output.writer.flush()?;
    }
    out_match.complete = true;
    if let Some(ref mut output) = out_nomatch {
        output.complete = true;
    }
    eprintln!(
        "matched_pairs={matches} total_pairs={total} amino_k={}",
        index.k()
    );
    Ok(())
}

fn main() -> io::Result<()> {
    match Args::parse().command {
        Action::Index {
            reference,
            output,
            k,
        } => build_index(&reference, &output, k),
        Action::Screen {
            index,
            reads,
            matched,
            unmatched,
            reads2,
            threads,
        } => screen(
            &index,
            &reads,
            reads2.as_deref(),
            &matched,
            unmatched.as_deref(),
            threads,
        ),
    }
}
