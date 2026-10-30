//! Exact amino-acid k-mers from all six translations using the standard genetic code.
//! Stops and codons containing ambiguous bases break a run; they are never skipped.

use std::collections::HashSet;
use std::io::{self, Read, Write};

const MAGIC: &[u8; 8] = b"NAOAAIDX";
const VERSION: u32 = 1;
const CODE: u8 = 1;
const FRAMES: u8 = 6;
#[cfg(test)]
const AMINO_ACIDS: &[u8] = b"ACDEFGHIKLMNPQRSTVWY";
// NCBI standard code: T,C,A,G codon order; ACDEFGHIKLMNPQRSTVWY residue codes; 255 = stop.
const CODON_CODES: [u8; 64] = [
    4, 4, 9, 9, 15, 15, 15, 15, 19, 19, 255, 255, 1, 1, 255, 18, 9, 9, 9, 9, 12, 12, 12, 12, 6, 6,
    13, 13, 14, 14, 14, 14, 7, 7, 7, 10, 16, 16, 16, 16, 11, 11, 8, 8, 15, 15, 14, 14, 17, 17, 17,
    17, 0, 0, 0, 0, 2, 2, 3, 3, 5, 5, 5, 5,
];

fn invalid(message: &str) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message)
}

fn base_code(base: u8) -> Option<usize> {
    match base.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

fn amino_code(a: u8, b: u8, c: u8) -> Option<u64> {
    let amino = CODON_CODES[base_code(a)? * 16 + base_code(b)? * 4 + base_code(c)?];
    (amino != 255).then_some(u64::from(amino))
}

fn complement(base: u8) -> u8 {
    match base.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        _ => b'N',
    }
}

/// Visit encoded peptide windows, returning early when `visit` returns true.
/// Reads and references use the same six-frame semantics; peptide reversal is
/// not used as a substitute for translating the reverse-complement DNA.
fn visit_kmers(seq: &[u8], k: u8, mut visit: impl FnMut(u64) -> bool) -> bool {
    let mask = (1u64 << (5 * k)) - 1;
    for reverse in [false, true] {
        let base = |i: usize| {
            if reverse {
                complement(seq[seq.len() - 1 - i])
            } else {
                seq[i]
            }
        };
        for offset in 0..3 {
            let mut word = 0;
            let mut valid = 0;
            let mut position = offset;
            while position + 2 < seq.len() {
                match amino_code(base(position), base(position + 1), base(position + 2)) {
                    Some(amino) => {
                        word = ((word << 5) | amino) & mask;
                        valid += 1;
                        if valid >= usize::from(k) && visit(word) {
                            return true;
                        }
                    }
                    None => {
                        word = 0;
                        valid = 0;
                    }
                }
                position += 3;
            }
        }
    }
    false
}

pub struct Index {
    k: u8,
    kmers: HashSet<u64>,
}

impl Index {
    pub fn new(k: u8) -> io::Result<Self> {
        if !(1..=12).contains(&k) {
            return Err(invalid("amino-acid k must be between 1 and 12"));
        }
        Ok(Self {
            k,
            kmers: HashSet::new(),
        })
    }

    pub fn k(&self) -> u8 {
        self.k
    }
    pub fn len(&self) -> usize {
        self.kmers.len()
    }
    pub fn is_empty(&self) -> bool {
        self.kmers.is_empty()
    }

    /// Each call is one reference record. Windows never cross record boundaries.
    pub fn add_reference(&mut self, sequence: &[u8]) {
        visit_kmers(sequence, self.k, |word| {
            self.kmers.insert(word);
            false
        });
    }

    pub fn matches(&self, sequence: &[u8]) -> bool {
        visit_kmers(sequence, self.k, |word| self.kmers.contains(&word))
    }

    /// Deterministic representation, independent of FASTA order and duplicate records.
    pub fn write(&self, mut output: impl Write) -> io::Result<()> {
        if self.is_empty() {
            return Err(invalid("reference has no usable amino-acid k-mers"));
        }
        output.write_all(MAGIC)?;
        output.write_all(&VERSION.to_le_bytes())?;
        output.write_all(&[self.k, CODE, FRAMES, 0])?;
        output.write_all(&(self.len() as u64).to_le_bytes())?;
        let mut words: Vec<_> = self.kmers.iter().copied().collect();
        words.sort_unstable();
        for word in words {
            output.write_all(&word.to_le_bytes())?;
        }
        Ok(())
    }

    /// Reject incompatible settings, unsorted/duplicate codes and truncated/trailing data.
    pub fn read(mut input: impl Read) -> io::Result<Self> {
        let mut header = [0u8; 24];
        input.read_exact(&mut header)?;
        if &header[..8] != MAGIC || u32::from_le_bytes(header[8..12].try_into().unwrap()) != VERSION
        {
            return Err(invalid("unsupported amino-acid index format"));
        }
        if header[13..16] != [CODE, FRAMES, 0] {
            return Err(invalid("unsupported translation settings in index"));
        }
        let mut index = Self::new(header[12])?;
        let count = u64::from_le_bytes(header[16..24].try_into().unwrap());
        if count == 0 {
            return Err(invalid("empty amino-acid index"));
        }
        // Do not allocate from an untrusted header count before reading the payload.
        let mut previous = None;
        for _ in 0..count {
            let mut bytes = [0u8; 8];
            input.read_exact(&mut bytes)?;
            let word = u64::from_le_bytes(bytes);
            if previous.is_some_and(|p| word <= p) {
                return Err(invalid("index codes are not strictly increasing"));
            }
            let mut remaining = word;
            for _ in 0..index.k {
                if remaining & 31 >= 20 {
                    return Err(invalid("invalid amino-acid code"));
                }
                remaining >>= 5;
            }
            if remaining != 0 {
                return Err(invalid("index code exceeds k"));
            }
            index.kmers.insert(word);
            previous = Some(word);
        }
        let mut trailing = [0u8; 1];
        if input.read(&mut trailing)? != 0 {
            return Err(invalid("trailing index data"));
        }
        Ok(index)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // A deliberately separate oracle: explicit codon groups, full translated
    // strings and peptide substring sets, rather than the production rolling encoding.
    fn oracle(seq: &[u8], k: usize) -> HashSet<Vec<u8>> {
        let groups = [
            (b'A', "GCT GCC GCA GCG"),
            (b'C', "TGT TGC"),
            (b'D', "GAT GAC"),
            (b'E', "GAA GAG"),
            (b'F', "TTT TTC"),
            (b'G', "GGT GGC GGA GGG"),
            (b'H', "CAT CAC"),
            (b'I', "ATT ATC ATA"),
            (b'K', "AAA AAG"),
            (b'L', "TTA TTG CTT CTC CTA CTG"),
            (b'M', "ATG"),
            (b'N', "AAT AAC"),
            (b'P', "CCT CCC CCA CCG"),
            (b'Q', "CAA CAG"),
            (b'R', "CGT CGC CGA CGG AGA AGG"),
            (b'S', "TCT TCC TCA TCG AGT AGC"),
            (b'T', "ACT ACC ACA ACG"),
            (b'V', "GTT GTC GTA GTG"),
            (b'W', "TGG"),
            (b'Y', "TAT TAC"),
        ];
        let forward: Vec<_> = seq
            .iter()
            .map(|b| match b.to_ascii_uppercase() {
                b'U' => b'T',
                other => other,
            })
            .collect();
        let reverse: Vec<_> = forward
            .iter()
            .rev()
            .map(|b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                _ => b'N',
            })
            .collect();
        let mut words = HashSet::new();
        for strand in [&forward, &reverse] {
            for frame in 0..3.min(strand.len()) {
                let translated: Vec<_> = strand[frame..]
                    .chunks_exact(3)
                    .map(|codon| {
                        groups
                            .iter()
                            .find(|(_, codons)| codons.split(' ').any(|c| c.as_bytes() == codon))
                            .map(|(aa, _)| *aa)
                            .unwrap_or(b'*')
                    })
                    .collect();
                for word in translated.windows(k) {
                    if !word.contains(&b'*') {
                        words.insert(word.to_vec());
                    }
                }
            }
        }
        words
    }

    fn decoded(seq: &[u8], k: u8) -> HashSet<Vec<u8>> {
        let mut output = HashSet::new();
        visit_kmers(seq, k, |mut word| {
            let mut peptide = vec![0; usize::from(k)];
            for aa in peptide.iter_mut().rev() {
                *aa = AMINO_ACIDS[(word & 31) as usize];
                word >>= 5;
            }
            output.insert(peptide);
            false
        });
        output
    }

    #[test]
    fn translations_match_independent_codon_oracle() {
        // All codons and adjacent codon pairs cover all amino acids, stops,
        // reverse strands and offset translations, without biological sequences.
        for bits in 0u32..4096 {
            let seq: Vec<_> = (0..6)
                .map(|i| b"ACGT"[((bits >> (2 * i)) & 3) as usize])
                .collect();
            for k in 1..=2 {
                assert_eq!(
                    decoded(&seq, k),
                    oracle(&seq, usize::from(k)),
                    "{seq:?}, k={k}"
                );
            }
        }
        for seq in [
            b"ATGGCTNNNGAATTT".as_slice(),
            b"ATGGCTTAAGAATTT",
            b"ATGGCTTAGGAATTT",
            b"ATGGCTTGAGAATTT",
            b"auggcugaauuu",
            b"ATGRYSWKMBDHVNGAATTT",
            b"",
            b"A",
            b"AT",
            b"ATGGCTGAATTTATGGCTGAATTTATGGCTGAATTTATGGCTGAATTT",
        ] {
            for k in [1, 2, 4, 12] {
                assert_eq!(decoded(seq, k), oracle(seq, usize::from(k)));
            }
        }
    }

    #[test]
    fn windows_do_not_cross_reference_records() {
        let mut index = Index::new(4).unwrap();
        index.add_reference(b"ATGGCT");
        index.add_reference(b"GAATTT");
        assert!(index.is_empty());
        assert!(!index.matches(b"ATGGCTGAATTT"));
    }
}
