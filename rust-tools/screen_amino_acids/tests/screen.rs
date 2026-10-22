use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

static NEXT: AtomicU64 = AtomicU64::new(0);

struct Fixture(PathBuf);

impl Fixture {
    fn new() -> Self {
        let path = std::env::temp_dir().join(format!(
            "amino-screen-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).unwrap();
        Self(path)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }

    fn write(&self, name: &str, data: &str) {
        fs::write(self.path(name), data).unwrap();
    }

    fn run(&self, args: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_screen_amino_acids"))
            .current_dir(&self.0)
            .args(args)
            .output()
            .unwrap()
    }

    fn index(&self) {
        success(self.run(&[
            "index",
            "--reference",
            "ref.fa",
            "--output",
            "ref.aa",
            "--k",
            "4",
        ]));
    }

    fn screen(&self) -> Output {
        self.run(&[
            "screen",
            "--index",
            "ref.aa",
            "--reads",
            "reads.fq",
            "--matched",
            "match.fq",
            "--unmatched",
            "nomatch.fq",
        ])
    }
}

impl Drop for Fixture {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}

fn success(output: Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn record(name: &str, sequence: &str) -> String {
    format!(
        "@{name}\n{sequence}\n+original annotation\n{}\n",
        (0..sequence.len())
            .map(|i| b"!I5H"[i % 4] as char)
            .collect::<String>()
    )
}

fn pair(name: &str, first: &str, second: &str) -> String {
    record(&format!("{name}/1"), first) + &record(&format!("{name}/2"), second)
}

fn read(path: impl AsRef<Path>) -> String {
    fs::read_to_string(path).unwrap()
}

fn gzip(data: &str) -> Vec<u8> {
    let mut writer = flate2::write::GzEncoder::new(Vec::new(), flate2::Compression::fast());
    writer.write_all(data.as_bytes()).unwrap();
    writer.finish().unwrap()
}

#[test]
fn pairs_are_partitioned_without_changing_records() {
    // Toy sequences, not organism-derived. Both coding strings translate to MAEF.
    let f = Fixture::new();
    f.write("ref.fa", ">toy\nATGGCTGAATTT\n");
    f.index();
    let rescued = pair("synonymous", "ATGGCCGAGTTC", "NNNNNNNNNNNN");
    let second_mate = pair("second", "NNNNNNNNNNNN", "ATGGCCGAGTTC");
    let both = pair("both", "ATGGCTGAATTT", "ATGGCCGAGTTC");
    let rejected = pair("negative", "GGTGGTGGTGGT", "GGTGGTGGTGGT");
    let casava =
        record("casava 1:N:0:ATCG", "ATGGCCGAGTTC") + &record("casava 2:N:0:ATCG", "NNNNNNNNNNNN");
    let unmarked = (record("unmarked", "ATGGCCGAGTTC") + &record("unmarked", "NNNNNNNNNNNN"))
        .replace('\n', "\r\n");
    f.write(
        "reads.fq",
        &(rescued.clone() + &rejected + &second_mate + &both + &casava + &unmarked),
    );
    let result = f.screen();
    assert!(String::from_utf8_lossy(&result.stderr).contains("matched_pairs=5 total_pairs=6"));
    success(result);
    assert_eq!(
        read(f.path("match.fq")),
        rescued + &second_mate + &both + &casava + &unmarked
    );
    assert_eq!(read(f.path("nomatch.fq")), rejected);
}

#[test]
fn all_frames_and_strands_are_searched() {
    let f = Fixture::new();
    // Reverse complement of ATGGCCGAGTTC, independently written.
    for reference in [
        "ATGGCTGAATTT",
        "AATGGCTGAATTT",
        "CCATGGCTGAATTT",
        "AAATTCAGCCAT",
        "AAAATTCAGCCAT",
        "CCAAATTCAGCCAT",
    ] {
        f.write("ref.fa", &format!(">toy\n{reference}\n"));
        // Use a new output: index never overwrites an existing file.
        let _ = fs::remove_file(f.path("ref.aa"));
        f.index();
        let input: String = [
            "ATGGCCGAGTTC",
            "AATGGCCGAGTTC",
            "CCATGGCCGAGTTC",
            "GAACTCGGCCAT",
            "AGAACTCGGCCAT",
            "CCGAACTCGGCCAT",
        ]
        .iter()
        .enumerate()
        .map(|(i, seq)| pair(&format!("frame{i}"), seq, "NNNNNNNNNNNN"))
        .collect();
        f.write("reads.fq", &input);
        let _ = fs::remove_file(f.path("match.fq"));
        let _ = fs::remove_file(f.path("nomatch.fq"));
        success(f.screen());
        assert_eq!(read(f.path("match.fq")), input);
        assert!(read(f.path("nomatch.fq")).is_empty());
    }
}

#[test]
fn invalid_pairs_fail_without_leaving_outputs() {
    let f = Fixture::new();
    f.write("ref.fa", ">toy\nATGGCTGAATTT\n");
    f.index();
    for invalid in [
        record("odd/1", "ATGGCTGAATTT"),
        record("one/1", "ATGGCTGAATTT") + &record("two/2", "ATGGCTGAATTT"),
        record("one/2", "ATGGCTGAATTT") + &record("one/1", "ATGGCTGAATTT"),
        record("one 2:N:0:ATCG", "ATGGCTGAATTT") + &record("one 1:N:0:ATCG", "ATGGCTGAATTT"),
        "@bad/1\nATGGCTGAATTT\n+\nI\n".to_string(),
        "@bad/1\nATGGCTGAATTT\n".to_string(),
    ] {
        f.write(
            "reads.fq",
            &(pair("valid", "ATGGCTGAATTT", "NNNNNNNNNNNN") + &invalid),
        );
        assert!(!f.screen().status.success());
        assert!(!f.path("match.fq").exists());
        assert!(!f.path("nomatch.fq").exists());
    }
}

#[test]
fn index_is_deterministic_and_rejects_incompatible_or_truncated_data() {
    let f = Fixture::new();
    for reference in [
        "",
        "ATGGCTGAATTT\n",
        ">bad\nATG-GCT\n",
        ">masked\nNNNNNNNNNNNN\n",
        ">a\nATGGCT\n>b\nGAATTT\n",
    ] {
        f.write("ref.fa", reference);
        assert!(!f
            .run(&[
                "index",
                "--reference",
                "ref.fa",
                "--output",
                "ref.aa",
                "--k",
                "4",
            ])
            .status
            .success());
        assert!(!f.path("ref.aa").exists());
    }
    f.write("ref.fa", ">a\nATGGCTGAATTT\n>b\nATGGCCGAGTTC\n");
    f.index();
    f.write("other.fa", ">b\nATGGCCGAGTTC\n>a\nATGGCTGAATTT\n");
    success(f.run(&[
        "index",
        "--reference",
        "other.fa",
        "--output",
        "other.aa",
        "--k",
        "4",
    ]));
    let original = fs::read(f.path("ref.aa")).unwrap();
    assert_eq!(original, fs::read(f.path("other.aa")).unwrap());
    f.write("reads.fq", "");
    for broken in [
        original[..original.len() - 1].to_vec(),
        b"not an amino acid index".to_vec(),
        {
            let mut v = original.clone();
            v[8] = 255;
            v
        },
        {
            let mut v = original.clone();
            v[13] = 2;
            v
        },
        {
            let mut v = original.clone();
            v[12] = 0;
            v
        },
        {
            let mut v = original.clone();
            v.push(0);
            v
        },
    ] {
        fs::write(f.path("ref.aa"), broken).unwrap();
        assert!(!f.screen().status.success());
        assert!(!f.path("match.fq").exists());
    }
    fs::write(f.path("ref.aa"), original).unwrap();
    success(f.screen());
    assert!(read(f.path("match.fq")).is_empty());
    assert!(read(f.path("nomatch.fq")).is_empty());
}

#[test]
fn existing_outputs_are_never_overwritten() {
    let f = Fixture::new();
    f.write("ref.fa", ">toy\nATGGCTGAATTT\n");
    f.index();
    f.write("reads.fq", &pair("one", "ATGGCTGAATTT", "NNNNNNNNNNNN"));
    f.write("match.fq", "keep me");
    assert!(!f.screen().status.success());
    assert_eq!(read(f.path("match.fq")), "keep me");
    assert!(!f.path("nomatch.fq").exists());
    assert!(!f
        .run(&[
            "index",
            "--reference",
            "ref.fa",
            "--output",
            "ref.aa",
            "--k",
            "4"
        ])
        .status
        .success());
}

#[test]
fn gzip_members_preserve_pairs_and_corrupt_streams_fail() {
    let f = Fixture::new();
    fs::write(f.path("ref.fa"), gzip(">toy\nATGGCTGAATTT\n")).unwrap();
    f.index();
    let first = pair("one", "ATGGCCGAGTTC", "NNNNNNNNNNNN");
    let second = pair("two", "NNNNNNNNNNNN", "NNNNNNNNNNNN");
    let mut compressed = gzip(&first);
    compressed.extend(gzip(&second));
    fs::write(f.path("reads.fq"), &compressed).unwrap();
    success(f.screen());
    assert_eq!(read(f.path("match.fq")), first);
    assert_eq!(read(f.path("nomatch.fq")), second);
    fs::remove_file(f.path("match.fq")).unwrap();
    fs::remove_file(f.path("nomatch.fq")).unwrap();
    for broken in [compressed[..compressed.len() - 5].to_vec(), {
        let mut v = compressed.clone();
        let n = v.len();
        v[n - 8] ^= 1;
        v
    }] {
        fs::write(f.path("reads.fq"), broken).unwrap();
        assert!(!f.screen().status.success());
        assert!(!f.path("match.fq").exists());
        assert!(!f.path("nomatch.fq").exists());
    }
}

#[test]
fn paired_and_parallel_inputs_preserve_the_partition() {
    let f = Fixture::new();
    f.write("ref.fa", ">toy\nATGGCTGAATTT\n");
    f.index();
    let mut r1=String::new(); let mut r2=String::new(); let mut interleaved=String::new();
    for i in 0..8200 {
        let a=record(&format!("pair{i}/1"), if i%3==0 {"ATGGCCGAGTTC"} else {"GGTGGTGGTGGT"});
        let b=record(&format!("pair{i}/2"), if i%5==0 {"ATGGCTGAATTT"} else {"GGTGGTGGTGGT"});
        r1.push_str(&a);r2.push_str(&b);interleaved.push_str(&a);interleaved.push_str(&b);
    }
    f.write("reads.fq", &interleaved); f.write("r1.fq", &r1); f.write("r2.fq", &r2);
    success(f.screen());
    success(f.run(&["screen","--index","ref.aa","--reads","r1.fq","--reads2","r2.fq","--matched","parallel.fq","--threads","4"]));
    assert_eq!(read(f.path("match.fq")),read(f.path("parallel.fq")));
    f.write("r2.fq", &(r2+&record("extra/2","GGTGGTGGTGGT")));
    assert!(!f.run(&["screen","--index","ref.aa","--reads","r1.fq","--reads2","r2.fq","--matched","bad.fq"]).status.success());
    assert!(!f.path("bad.fq").exists());
}
