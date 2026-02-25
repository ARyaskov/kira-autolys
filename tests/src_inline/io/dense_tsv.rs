
use super::*;
use std::fs;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn parse_gene_major_dense() {
    let dir = temp_dir("dense_gene_major");
    let file = dir.join("raw_counts.tsv");
    write_file(&file, "gene\tC1\tC2\nG1\t1\t0\nG2\t0\t2\nG3\t3\t0\n");
    let dense = load_dense(&file).unwrap();
    assert_eq!(dense.barcodes, vec!["C1".to_string(), "C2".to_string()]);
    assert_eq!(dense.genes.len(), 3);
    assert_eq!(dense.matrix.n_obs, 2);
    assert_eq!(dense.matrix.n_vars, 3);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn resolve_prefixed_raw_counts() {
    let dir = temp_dir("dense_prefixed");
    write_file(
        &dir.join("GSM123_raw_counts.tsv"),
        "gene\tC1\tC2\nG1\t1\t0\n",
    );
    let resolved = resolve_dense_input_path(&dir).unwrap();
    assert_eq!(
        resolved.file_name().and_then(|v| v.to_str()),
        Some("GSM123_raw_counts.tsv")
    );
    let _ = fs::remove_dir_all(dir);
}

fn temp_dir(label: &str) -> PathBuf {
    let base = std::env::temp_dir();
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = base.join(format!("kira_autolys_{label}_{ts}"));
    fs::create_dir_all(&path).unwrap();
    path
}

fn write_file(path: &Path, content: &str) {
    let mut file = fs::File::create(path).unwrap();
    file.write_all(content.as_bytes()).unwrap();
}
