
use super::*;
use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::time::{SystemTime, UNIX_EPOCH};

fn temp_dir(label: &str) -> PathBuf {
    let base = std::env::temp_dir();
    let ts = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = base.join(format!("kira_autolys_cache_{label}_{ts}"));
    fs::create_dir_all(&path).unwrap();
    path
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

fn write_test_cache(path: &Path) {
    let n_genes = 3u64;
    let n_cells = 2u64;
    let nnz = 3u64;

    let genes = encode_string_table(&["G1", "G2", "G3"]);
    let barcodes = encode_string_table(&["C1", "C2"]);

    let col_ptr: [u64; 3] = [0, 2, 3];
    let row_idx: [u32; 3] = [0, 2, 1];
    let values: [u32; 3] = [5, 1, 7];

    let mut bytes = vec![0u8; HEADER_SIZE];
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

    bytes[0..4].copy_from_slice(MAGIC);
    bytes[4..6].copy_from_slice(&VERSION_MAJOR.to_le_bytes());
    bytes[6..8].copy_from_slice(&0u16.to_le_bytes());
    bytes[8..12].copy_from_slice(&ENDIAN_TAG.to_le_bytes());
    bytes[12..16].copy_from_slice(&(HEADER_SIZE as u32).to_le_bytes());

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

    let mut header_copy = bytes[..HEADER_SIZE].to_vec();
    header_copy[120..128].fill(0);
    let crc = crc64_ecma(&header_copy);
    bytes[120..128].copy_from_slice(&crc.to_le_bytes());

    let mut file = fs::File::create(path).unwrap();
    file.write_all(&bytes).unwrap();
}

#[test]
fn cache_roundtrip_and_header_crc() {
    let dir = temp_dir("roundtrip");
    let path = dir.join("kira-organelle.bin");
    write_test_cache(&path);

    let cache = SharedCache::open(&path).unwrap();
    assert_eq!(cache.n_genes, 3);
    assert_eq!(cache.n_cells, 2);
    assert_eq!(cache.nnz, 3);
    assert_eq!(cache.genes, vec!["G1", "G2", "G3"]);
    assert_eq!(cache.barcodes, vec!["C1", "C2"]);
    assert_eq!(cache.col_ptr().unwrap(), &[0, 2, 3]);
    assert_eq!(cache.row_idx().unwrap(), &[0, 2, 1]);
    assert_eq!(cache.values_u32().unwrap(), &[5, 1, 7]);

    let norm = MmapCscNormalized::new(cache, 1e4).unwrap();
    let mut acc = Vec::new();
    norm.for_each_in_obs(0, &mut |gene, val| acc.push((gene, val)));
    assert_eq!(acc.len(), 2);
    assert_eq!(acc[0].0, 0);
    assert_eq!(acc[1].0, 2);

    let _ = fs::remove_dir_all(dir);
}
