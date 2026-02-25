use std::collections::BTreeSet;
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::{Path, PathBuf};

use flate2::read::GzDecoder;
use tracing::warn;

use crate::expr::normalized::CscMatrix;
use crate::stage_error::StageError;

#[derive(Debug, Clone)]
pub struct DenseInput {
    pub matrix: CscMatrix,
    pub genes: Vec<String>,
    pub barcodes: Vec<String>,
}

#[derive(Debug)]
struct ParsedDenseRows {
    barcodes: Vec<String>,
    genes: Vec<String>,
    columns: Vec<Vec<(usize, f32)>>,
}

pub fn resolve_dense_input_path(path: &Path) -> Result<PathBuf, StageError> {
    if path.is_file() {
        return Ok(path.to_path_buf());
    }
    if !path.is_dir() {
        return Err(StageError::Validation(format!(
            "input path does not exist: {}",
            path.display()
        )));
    }

    for name in ["raw_counts.tsv.gz", "raw_counts.tsv"] {
        let candidate = path.join(name);
        if candidate.is_file() {
            return Ok(candidate);
        }
    }

    let mut candidates = Vec::new();
    for entry in std::fs::read_dir(path)? {
        let entry = entry?;
        let p = entry.path();
        if !p.is_file() {
            continue;
        }
        let Some(name) = p.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        if name.ends_with("_raw_counts.tsv.gz") || name.ends_with("_raw_counts.tsv") {
            candidates.push(p);
        }
    }
    candidates.sort();
    candidates.first().cloned().ok_or_else(|| {
        StageError::Validation(format!(
            "unsupported input: {} (expected raw_counts.tsv(.gz) or *_raw_counts.tsv(.gz))",
            path.display()
        ))
    })
}

pub fn is_dense_tsv_path(path: &Path) -> bool {
    let Some(name) = path.file_name().and_then(|n| n.to_str()) else {
        return false;
    };
    name.ends_with(".tsv") || name.ends_with(".tsv.gz")
}

pub fn load_dense(path: &Path) -> Result<DenseInput, StageError> {
    let source = resolve_dense_input_path(path)?;
    let parsed = parse_dense_rows(&source, true)?;
    let mut indptr = Vec::with_capacity(parsed.barcodes.len() + 1);
    let mut indices = Vec::new();
    let mut data = Vec::new();
    indptr.push(0);
    for col in parsed.columns {
        for (row_idx, value) in col {
            indices.push(row_idx);
            data.push(value);
        }
        indptr.push(indices.len());
    }

    Ok(DenseInput {
        matrix: CscMatrix {
            n_vars: parsed.genes.len(),
            n_obs: parsed.barcodes.len(),
            indptr,
            indices,
            data,
        },
        genes: parsed.genes,
        barcodes: parsed.barcodes,
    })
}

pub fn load_dense_metadata(path: &Path) -> Result<(Vec<String>, Vec<String>), StageError> {
    let source = resolve_dense_input_path(path)?;
    let parsed = parse_dense_rows(&source, false)?;
    Ok((parsed.barcodes, parsed.genes))
}

fn parse_dense_rows(path: &Path, with_values: bool) -> Result<ParsedDenseRows, StageError> {
    let reader = open_maybe_gz(path)?;

    let mut barcodes = Vec::new();
    let mut genes: Vec<String> = Vec::new();
    let mut columns: Vec<Vec<(usize, f32)>> = Vec::new();
    let mut duplicate_symbols = BTreeSet::new();
    let mut seen_symbols = BTreeSet::new();
    let mut duplicate_barcodes = BTreeSet::new();
    let mut seen_barcodes = BTreeSet::new();
    let mut non_finite_count = 0usize;
    let mut cell_major = false;
    let mut header_parsed = false;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let cleaned = line.trim_end_matches(&['\r', '\n'][..]);
        if cleaned.trim().is_empty() || cleaned.trim_start().starts_with('#') {
            continue;
        }

        if !header_parsed {
            let header = cleaned
                .split('\t')
                .map(|s| s.trim().to_string())
                .collect::<Vec<_>>();
            cell_major = first_col_looks_like_cell_header(header.first().map(String::as_str));
            if cell_major {
                genes = header.into_iter().skip(1).collect();
                if genes.len() <= 1 {
                    return Err(StageError::Format(
                        "dense TSV header must contain at least 2 genes".to_string(),
                    ));
                }
                for g in &genes {
                    if !seen_symbols.insert(g.clone()) {
                        duplicate_symbols.insert(g.clone());
                    }
                }
                if with_values {
                    columns = Vec::new();
                }
            } else {
                barcodes = if first_col_looks_like_gene_header(header.first().map(String::as_str)) {
                    header.into_iter().skip(1).collect()
                } else {
                    header
                };
                if barcodes.len() <= 1 {
                    return Err(StageError::Format(
                        "dense TSV header must contain at least 2 cell IDs".to_string(),
                    ));
                }
                if with_values {
                    columns = vec![Vec::new(); barcodes.len()];
                }
            }
            header_parsed = true;
            continue;
        }

        let mut parts = cleaned.split('\t').collect::<Vec<_>>();
        if cell_major {
            let expected = genes.len() + 1;
            if parts.len() + 1 == expected {
                parts.push("");
            }
            if parts.len() != expected {
                return Err(StageError::Format(format!(
                    "{}: line {}: expected {} columns, got {}",
                    path.display(),
                    line_no + 1,
                    expected,
                    parts.len()
                )));
            }
            let barcode = parts[0].trim().to_string();
            if !seen_barcodes.insert(barcode.clone()) {
                duplicate_barcodes.insert(barcode.clone());
            }
            let col_idx = barcodes.len();
            barcodes.push(barcode);
            if with_values {
                let mut col = Vec::new();
                for (gene_idx, token) in parts[1..].iter().enumerate() {
                    let token = token.trim();
                    let mut value = parse_float_token(path, line_no + 1, gene_idx + 2, token)?;
                    if !value.is_finite() {
                        non_finite_count += 1;
                        value = 0.0;
                    }
                    if value != 0.0 {
                        col.push((gene_idx, value));
                    }
                }
                columns.push(col);
                if columns.len() != col_idx + 1 {
                    return Err(StageError::Format(format!(
                        "{}: internal column index mismatch at line {}",
                        path.display(),
                        line_no + 1
                    )));
                }
            }
        } else {
            let expected = barcodes.len() + 1;
            if parts.len() + 1 == expected {
                parts.push("");
            }
            if parts.len() != expected {
                return Err(StageError::Format(format!(
                    "{}: line {}: expected {} columns, got {}",
                    path.display(),
                    line_no + 1,
                    expected,
                    parts.len()
                )));
            }
            let symbol = parts[0].trim().to_string();
            let row_idx = genes.len();
            genes.push(symbol.clone());
            if !seen_symbols.insert(symbol.clone()) {
                duplicate_symbols.insert(symbol);
            }
            if with_values {
                for (col_idx, token) in parts[1..].iter().enumerate() {
                    let token = token.trim();
                    let mut value = parse_float_token(path, line_no + 1, col_idx + 2, token)?;
                    if !value.is_finite() {
                        non_finite_count += 1;
                        value = 0.0;
                    }
                    if value != 0.0 {
                        columns[col_idx].push((row_idx, value));
                    }
                }
            }
        }
    }

    if barcodes.is_empty() {
        return Err(StageError::Format(format!(
            "empty dense TSV input: {}",
            path.display()
        )));
    }
    if genes.is_empty() {
        return Err(StageError::Format(format!(
            "dense TSV contains no genes: {}",
            path.display()
        )));
    }

    for symbol in duplicate_symbols {
        warn!(
            symbol = symbol.as_str(),
            "duplicate gene symbol in dense input"
        );
    }
    for barcode in duplicate_barcodes {
        warn!(
            barcode = barcode.as_str(),
            "duplicate barcode in dense input"
        );
    }
    if non_finite_count > 0 {
        warn!(
            count = non_finite_count,
            "non-finite values in dense input were replaced with 0"
        );
    }

    Ok(ParsedDenseRows {
        barcodes,
        genes,
        columns,
    })
}

fn parse_float_token(
    path: &Path,
    line_no: usize,
    col_no: usize,
    token: &str,
) -> Result<f32, StageError> {
    if token.is_empty() {
        return Ok(0.0);
    }
    token.parse::<f32>().map_err(|_| {
        StageError::Format(format!(
            "{}: line {} column {}: invalid float value `{}`",
            path.display(),
            line_no,
            col_no,
            token
        ))
    })
}

fn open_maybe_gz(path: &Path) -> Result<BufReader<Box<dyn Read>>, StageError> {
    if path
        .extension()
        .and_then(|v| v.to_str())
        .map(|v| v.eq_ignore_ascii_case("gz"))
        .unwrap_or(false)
    {
        let file = File::open(path)?;
        Ok(BufReader::new(Box::new(GzDecoder::new(file))))
    } else {
        let file = File::open(path)?;
        Ok(BufReader::new(Box::new(file)))
    }
}

fn first_col_looks_like_gene_header(v: Option<&str>) -> bool {
    let Some(v) = v else { return false };
    if v.trim().is_empty() {
        return true;
    }
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "gene" | "genes" | "gene_symbol" | "genesymbol" | "symbol" | "feature" | "features"
    )
}

fn first_col_looks_like_cell_header(v: Option<&str>) -> bool {
    let Some(v) = v else { return false };
    matches!(
        v.trim().to_ascii_lowercase().as_str(),
        "barcode" | "barcodes" | "cell" | "cell_id" | "cellid"
    )
}

#[cfg(test)]
#[path = "../../tests/src_inline/io/dense_tsv.rs"]
mod tests;
