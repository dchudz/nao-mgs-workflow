// Exploratory AA12 index subset: omit words containing ten or more copies of one residue.
// Preserve sorted encoded words and the existing validated NAOAAIDX file contract.
use std::fs::OpenOptions;
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = std::env::args().collect();
    assert_eq!(args.len(), 3, "usage: filter-aa-repeat-index INPUT OUTPUT");
    let mut input = BufReader::with_capacity(1 << 20, std::fs::File::open(&args[1])?);
    let mut header = [0u8; 24]; input.read_exact(&mut header)?;
    assert_eq!(&header[..8], b"NAOAAIDX");
    assert_eq!(u32::from_le_bytes(header[8..12].try_into()?), 1);
    assert_eq!(&header[12..16], &[12,1,6,0]);
    let original = u64::from_le_bytes(header[16..24].try_into()?);
    let mut output = BufWriter::with_capacity(1 << 20, OpenOptions::new().write(true).create_new(true).open(&args[2])?);
    output.write_all(&header)?;
    let mut retained = 0u64; let mut previous = None;
    for _ in 0..original {
        let mut bytes = [0u8;8]; input.read_exact(&mut bytes)?;
        let word = u64::from_le_bytes(bytes);
        assert!(previous.map_or(true, |p| word > p)); previous = Some(word);
        let mut remaining = word; let mut counts = [0u8;20];
        for _ in 0..12 {
            let residue = (remaining & 31) as usize; assert!(residue < 20);
            counts[residue] += 1; remaining >>= 5;
        }
        assert_eq!(remaining, 0);
        if counts.into_iter().max().unwrap() <= 9 { output.write_all(&bytes)?; retained += 1; }
    }
    let mut trailing = [0u8;1]; assert_eq!(input.read(&mut trailing)?,0); assert!(retained > 0);
    output.flush()?; output.seek(SeekFrom::Start(16))?; output.write_all(&retained.to_le_bytes())?; output.flush()?;
    println!("{{\"original_words\":{},\"retained_words\":{},\"removed_words\":{},\"rule\":\"No residue occupies ten or more of twelve positions\"}}",original,retained,original-retained);
    Ok(())
}
