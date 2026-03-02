use super::*;
use crate::expr::normalized::compute_libsize;
use crate::model::cli::RunMode;
use crate::model::ctx::InputType;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

#[test]
fn test_libsize_and_norm() {
    let matrix = CscMatrix {
        n_vars: 2,
        n_obs: 2,
        indptr: vec![0, 2, 3],
        indices: vec![0, 1, 0],
        data: vec![1.0, 3.0, 2.0],
    };
    let libsize = compute_libsize(&matrix);
    assert!((libsize[0] - 4.0).abs() < 1e-6);
    assert!((libsize[1] - 2.0).abs() < 1e-6);

    let norm = CscNormalized::new(matrix, 1e4).unwrap();
    let mut acc = Vec::new();
    norm.for_each_in_obs(0, &mut |idx, val| acc.push((idx, val)));
    assert_eq!(acc.len(), 2);
    let expected = ((1.0_f32 / 4.0_f32) * 1e4_f32).ln_1p();
    assert!((acc[0].1 - expected).abs() < 1e-6);
}

#[test]
fn test_zero_libsize() {
    let matrix = CscMatrix {
        n_vars: 1,
        n_obs: 1,
        indptr: vec![0, 0],
        indices: vec![],
        data: vec![],
    };
    let norm = CscNormalized::new(matrix, 1e4).unwrap();
    let mut called = false;
    norm.for_each_in_obs(0, &mut |_idx, _val| called = true);
    assert!(!called);
}

#[test]
fn pipeline_mode_uses_mmap_cache() {
    let dir = temp_dir("stage1_pipeline_cache");
    write_test_cache(&dir.join("kira-organelle.bin"));

    let mut ctx = Ctx::default();
    ctx.run_mode = RunMode::Pipeline;
    ctx.input_type = Some(InputType::Mtx10x);
    ctx.input_path = Some(dir.clone());
    ctx.obs_ids = vec!["C1".to_string(), "C2".to_string()];
    ctx.gene_symbols = vec!["G1".to_string(), "G2".to_string(), "G3".to_string()];
    ctx.n_obs = 2;
    ctx.n_vars = 3;

    run_stage1(&mut ctx).unwrap();
    assert!(ctx.normalized.is_some());
    let normalized = ctx.normalized.as_ref().unwrap();
    let mut acc = Vec::new();
    normalized.for_each_in_obs(1, &mut |gene, val| acc.push((gene, val)));
    assert_eq!(acc.len(), 1);
    assert_eq!(acc[0].0, 1);

    let _ = fs::remove_dir_all(dir);
}

#[test]
fn pipeline_mode_missing_cache_errors_after_mtx_fallback() {
    let dir = temp_dir("stage1_pipeline_missing_cache");
    write_file(
        &dir.join("matrix.mtx"),
        "%%MatrixMarket matrix coordinate integer general\n2 2 2\n1 1 1\n2 2 3\n",
    );
    write_file(
        &dir.join("features.tsv"),
        "ENSG1\tG1\tGene Expression\nENSG2\tG2\tGene Expression\n",
    );
    write_file(&dir.join("barcodes.tsv"), "C1\nC2\n");

    let mut ctx = Ctx::default();
    ctx.run_mode = RunMode::Pipeline;
    ctx.input_type = Some(InputType::Mtx10x);
    ctx.input_path = Some(dir.clone());
    ctx.obs_ids = vec!["C1".to_string(), "C2".to_string()];
    ctx.gene_symbols = vec!["G1".to_string(), "G2".to_string()];
    ctx.n_obs = 2;
    ctx.n_vars = 2;

    let err = run_stage1(&mut ctx).unwrap_err();
    let msg = err.to_string();
    assert!(msg.contains("pipeline mode expects shared cache file"));
    assert!(msg.contains("kira-organelle.bin"));

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
    let values: [u32; 3] = [5, 1, 7];

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
    pad_to_64(&mut bytes);

    let file_bytes = bytes.len() as u64;
    bytes[0..4].copy_from_slice(magic);
    bytes[4..6].copy_from_slice(&1u16.to_le_bytes());
    bytes[6..8].copy_from_slice(&0u16.to_le_bytes());
    bytes[8..12].copy_from_slice(&endian_tag.to_le_bytes());
    bytes[12..16].copy_from_slice(&(header_size as u32).to_le_bytes());
    bytes[16..24].copy_from_slice(&n_genes.to_le_bytes());
    bytes[24..32].copy_from_slice(&n_cells.to_le_bytes());
    bytes[32..40].copy_from_slice(&nnz.to_le_bytes());
    bytes[40..48].copy_from_slice(&genes_offset.to_le_bytes());
    bytes[48..56].copy_from_slice(&(genes.len() as u64).to_le_bytes());
    bytes[56..64].copy_from_slice(&barcodes_offset.to_le_bytes());
    bytes[64..72].copy_from_slice(&(barcodes.len() as u64).to_le_bytes());
    bytes[72..80].copy_from_slice(&col_ptr_offset.to_le_bytes());
    bytes[80..88].copy_from_slice(&row_idx_offset.to_le_bytes());
    bytes[88..96].copy_from_slice(&values_offset.to_le_bytes());
    bytes[96..104].copy_from_slice(&0u64.to_le_bytes());
    bytes[104..112].copy_from_slice(&0u64.to_le_bytes());
    bytes[112..120].copy_from_slice(&file_bytes.to_le_bytes());
    bytes[128..136].copy_from_slice(&0u64.to_le_bytes());

    let mut header_copy = bytes[..header_size].to_vec();
    header_copy[120..128].fill(0);
    let crc = crc64_ecma(&header_copy);
    bytes[120..128].copy_from_slice(&crc.to_le_bytes());

    let mut file = fs::File::create(path).unwrap();
    file.write_all(&bytes).unwrap();
}
