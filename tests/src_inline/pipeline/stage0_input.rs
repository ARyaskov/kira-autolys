use super::*;
use crate::model::cli::{CliArgs, RunMode};
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn test_alias_resolution() {
    let symbols = vec!["FIP200".to_string()];
    let canonical = canonicalize_symbols(&symbols);
    let (_symbols, index) = build_gene_index(&canonical).unwrap();
    let idx_alias = index.get("RB1CC1").copied();
    let idx_direct = index.get("FIP200").copied();
    assert_eq!(idx_alias, idx_direct);
}

#[test]
fn test_missing_gene_determinism() {
    let mut gene_index = BTreeMap::new();
    gene_index.insert("ATG5".to_string(), 0);
    gene_index.insert("ATG7".to_string(), 1);

    let (resolved, missing, _mask) = resolve_gene_sets(&gene_index, 2).unwrap();
    assert_eq!(resolved.len(), missing.len());
    for panel in missing {
        let mut sorted = panel.missing.clone();
        sorted.sort();
        assert_eq!(panel.missing, sorted);
    }
}

#[test]
fn test_mode_auto_detection() {
    assert_eq!(auto_mode(InputType::Mtx10x, 1000), Mode::Cell);
    assert_eq!(auto_mode(InputType::Mtx10x, 999), Mode::Sample);
    assert_eq!(auto_mode(InputType::DenseTsv, 5000), Mode::Cell);
    assert_eq!(auto_mode(InputType::H5ad, 5000), Mode::Sample);
}

#[test]
fn test_panel_mask() {
    let mut gene_index = BTreeMap::new();
    gene_index.insert("RB1CC1".to_string(), 0);
    gene_index.insert("ATG5".to_string(), 1);
    gene_index.insert("SQSTM1".to_string(), 2);

    let (_resolved, _missing, mask) = resolve_gene_sets(&gene_index, 3).unwrap();
    assert!(mask[0] != 0);
    assert!(mask[1] != 0);
    assert!(mask[2] != 0);
}

#[test]
fn test_stage0_pipeline_uses_cache_without_mtx_files() {
    let dir = temp_dir("stage0_pipeline_cache");
    write_test_cache(&dir.join("kira-organelle.bin"));

    let mut ctx = Ctx::default();
    let cli = CliArgs {
        input: dir.clone(),
        cache_path: None,
        mode: None,
        manifest: None,
        timecourse: false,
        run_mode: RunMode::Pipeline,
    };

    run_stage0(&mut ctx, &cli).unwrap();
    assert_eq!(ctx.n_obs, 2);
    assert_eq!(ctx.n_vars, 3);
    assert_eq!(ctx.obs_ids, vec!["C1".to_string(), "C2".to_string()]);
    assert_eq!(
        ctx.gene_symbols,
        vec!["G1".to_string(), "G2".to_string(), "G3".to_string()]
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

fn pad_to_64(bytes: &mut Vec<u8>) {
    let rem = bytes.len() % 64;
    if rem != 0 {
        bytes.resize(bytes.len() + (64 - rem), 0);
    }
}

fn encode_string_table(values: &[&str]) -> Vec<u8> {
    let mut table = Vec::new();
    table.extend_from_slice(&(values.len() as u32).to_le_bytes());

    let mut offsets = Vec::with_capacity(values.len() + 1);
    let mut total = 0u32;
    offsets.push(total);
    for value in values {
        total += value.len() as u32;
        offsets.push(total);
    }
    for off in offsets {
        table.extend_from_slice(&off.to_le_bytes());
    }
    for value in values {
        table.extend_from_slice(value.as_bytes());
    }
    table
}

fn crc64_ecma(bytes: &[u8]) -> u64 {
    const POLY: u64 = 0x42F0_E1EB_A9EA_3693;
    let mut crc = 0u64;
    for &byte in bytes {
        crc ^= (byte as u64) << 56;
        for _ in 0..8 {
            if (crc & 0x8000_0000_0000_0000) != 0 {
                crc = (crc << 1) ^ POLY;
            } else {
                crc <<= 1;
            }
        }
    }
    crc
}

fn write_test_cache(path: &Path) {
    let n_genes = 3u64;
    let n_cells = 2u64;
    let nnz = 3u64;
    let header_size = 256usize;
    let magic = b"KORG";
    let endian_tag: u32 = 0x1234_5678;

    let genes = encode_string_table(&["G1", "G2", "G3"]);
    let barcodes = encode_string_table(&["C1", "C2"]);
    let col_ptr: [u64; 3] = [0, 2, 3];
    let row_idx: [u32; 3] = [0, 2, 1];
    let values: [u32; 3] = [1, 3, 2];

    let mut bytes = vec![0u8; header_size];
    pad_to_64(&mut bytes);
    let genes_offset = bytes.len() as u64;
    bytes.extend_from_slice(&genes);
    pad_to_64(&mut bytes);
    let barcodes_offset = bytes.len() as u64;
    bytes.extend_from_slice(&barcodes);
    pad_to_64(&mut bytes);
    let col_ptr_offset = bytes.len() as u64;
    for v in col_ptr {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    pad_to_64(&mut bytes);
    let row_idx_offset = bytes.len() as u64;
    for v in row_idx {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    pad_to_64(&mut bytes);
    let values_offset = bytes.len() as u64;
    for v in values {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    let file_bytes = bytes.len() as u64;

    write_file(path, "");
    let mut f = fs::OpenOptions::new()
        .write(true)
        .truncate(true)
        .open(path)
        .unwrap();
    f.write_all(magic).unwrap();
    f.write_all(&1u16.to_le_bytes()).unwrap();
    f.write_all(&0u16.to_le_bytes()).unwrap();
    f.write_all(&endian_tag.to_le_bytes()).unwrap();
    f.write_all(&(header_size as u32).to_le_bytes()).unwrap();
    f.write_all(&n_genes.to_le_bytes()).unwrap();
    f.write_all(&n_cells.to_le_bytes()).unwrap();
    f.write_all(&nnz.to_le_bytes()).unwrap();
    f.write_all(&genes_offset.to_le_bytes()).unwrap();
    f.write_all(&(genes.len() as u64).to_le_bytes()).unwrap();
    f.write_all(&barcodes_offset.to_le_bytes()).unwrap();
    f.write_all(&(barcodes.len() as u64).to_le_bytes()).unwrap();
    f.write_all(&col_ptr_offset.to_le_bytes()).unwrap();
    f.write_all(&row_idx_offset.to_le_bytes()).unwrap();
    f.write_all(&values_offset.to_le_bytes()).unwrap();
    f.write_all(&0u64.to_le_bytes()).unwrap();
    f.write_all(&0u64.to_le_bytes()).unwrap();
    f.write_all(&file_bytes.to_le_bytes()).unwrap();
    f.write_all(&0u64.to_le_bytes()).unwrap();
    f.write_all(&0u64.to_le_bytes()).unwrap();
    f.write_all(&[0u8; 256 - 136]).unwrap();
    f.flush().unwrap();

    let mut all = fs::read(path).unwrap();
    all.resize(header_size, 0);
    let mut header_copy = all[..header_size].to_vec();
    header_copy[120..128].fill(0);
    let crc = crc64_ecma(&header_copy);
    all[120..128].copy_from_slice(&crc.to_le_bytes());
    all.extend_from_slice(&bytes[header_size..]);
    fs::write(path, all).unwrap();
}
