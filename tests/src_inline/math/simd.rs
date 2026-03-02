use super::*;

#[test]
fn test_scalar_matches_formula() {
    let values = vec![1.0_f32, 2.0, 3.0];
    let mut out = vec![0.0_f32; 3];
    compute_z_scores_scalar(&values, 2.0, 1.0, &mut out);
    assert_eq!(out, vec![-1.0, 0.0, 1.0]);
}

#[test]
fn test_dispatch_scalar_when_no_simd() {
    let values = vec![1.0_f32, 2.0, 3.0];
    let mut out = vec![0.0_f32; 3];
    compute_z_scores(&values, 2.0, 1.0, &mut out);
    let expected = [-1.0_f32, 0.0_f32, 1.0_f32];
    for (a, b) in out.iter().zip(expected.iter()) {
        assert!((a - b).abs() < 1e-5);
    }
}

#[cfg(feature = "simd")]
#[test]
fn test_simd_matches_scalar_if_available() {
    let values: Vec<f32> = (0..64).map(|v| v as f32 * 0.5).collect();
    let mut out_scalar = vec![0.0_f32; values.len()];
    let mut out_simd = vec![0.0_f32; values.len()];
    let mean = 7.5;
    let std = 2.5;
    let inv = 1.0 / (std + EPS);
    compute_z_scores_scalar(&values, mean, inv, &mut out_scalar);

    compute_z_scores(&values, mean, std, &mut out_simd);

    for (a, b) in out_scalar.iter().zip(out_simd.iter()) {
        assert!((a - b).abs() < 1e-6);
    }
}
