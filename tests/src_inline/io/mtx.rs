
use super::*;
use std::fs;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

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

fn touch(path: &Path) {
    let mut file = fs::File::create(path).unwrap();
    file.write_all(b"").unwrap();
}

#[test]
fn detect_prefix_none_for_standard_names() {
    let dir = temp_dir("prefix_none");
    touch(&dir.join("matrix.mtx"));
    touch(&dir.join("features.tsv"));
    touch(&dir.join("barcodes.tsv"));
    let prefix = detect_dataset_prefix(&dir).unwrap();
    assert_eq!(prefix, None);
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn detect_prefix_for_prefixed_names() {
    let dir = temp_dir("prefix_yes");
    touch(&dir.join("ABC_matrix.mtx.gz"));
    touch(&dir.join("ABC_features.tsv.gz"));
    touch(&dir.join("ABC_barcodes.tsv.gz"));
    let prefix = detect_dataset_prefix(&dir).unwrap();
    assert_eq!(prefix.as_deref(), Some("ABC"));
    let bin = resolve_cache_bin_path(&dir).unwrap();
    assert_eq!(
        bin.file_name().and_then(|x| x.to_str()),
        Some("ABC.kira-organelle.bin")
    );
    let _ = fs::remove_dir_all(dir);
}

#[test]
fn detect_cache_bin_name_without_prefix() {
    let dir = temp_dir("cache_name");
    touch(&dir.join("matrix.mtx"));
    touch(&dir.join("features.tsv"));
    touch(&dir.join("barcodes.tsv"));
    let bin = resolve_cache_bin_path(&dir).unwrap();
    assert_eq!(
        bin.file_name().and_then(|x| x.to_str()),
        Some("kira-organelle.bin")
    );
    let _ = fs::remove_dir_all(dir);
}
