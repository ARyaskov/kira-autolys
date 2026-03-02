use super::*;

#[test]
fn test_running_stats() {
    let mut stats = RunningStats::new();
    let values = [1.0_f32, 2.0, 3.0, 4.0];
    for v in values {
        stats.update(v);
    }
    assert!((stats.mean() - 2.5).abs() < 1e-6);
    let expected_var = 1.25;
    assert!((stats.variance() - expected_var).abs() < 1e-6);
}
