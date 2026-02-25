//! SIMD-optional math kernels for Stage 5.
//!
//! Numerical semantics mirror scalar computation. SIMD may introduce
//! minor lane-level floating-point differences (tolerance ~1e-6).

use crate::config::thresholds::Thresholds;

pub const EPS: f32 = 1e-6;

pub struct SurvivalFlags<'a> {
    pub ssm_high: &'a mut [bool],
    pub afp_high: &'a mut [bool],
    pub lds_high: &'a mut [bool],
    pub prolif_low: &'a mut [bool],
    pub apop_low: &'a mut [bool],
}

pub fn compute_z_scores(values: &[f32], mean: f32, std: f32, out: &mut [f32]) {
    debug_assert_eq!(values.len(), out.len());
    let inv_std = 1.0 / (std + EPS);

    #[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
    {
        // Compile-time backend selection: AVX2 path is only compiled when enabled for target.
        return unsafe { avx2::compute_z_scores(values, mean, inv_std, out) };
    }

    #[cfg(all(feature = "simd", target_arch = "aarch64", target_feature = "neon"))]
    {
        // Compile-time backend selection: NEON path is only compiled when enabled for target.
        return unsafe { neon::compute_z_scores(values, mean, inv_std, out) };
    }

    #[cfg(not(any(
        all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"),
        all(feature = "simd", target_arch = "aarch64", target_feature = "neon")
    )))]
    compute_z_scores_scalar(values, mean, inv_std, out);
}

pub fn compute_ssm_and_flags(
    z_afp: &[f32],
    z_lds: &[f32],
    z_prolif: &[f32],
    z_apop: &[f32],
    thresholds: &Thresholds,
    out_ssm: &mut [f32],
    flags: &mut SurvivalFlags<'_>,
) {
    debug_assert_eq!(z_afp.len(), out_ssm.len());

    #[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
    {
        return unsafe {
            avx2::compute_ssm_and_flags(z_afp, z_lds, z_prolif, z_apop, thresholds, out_ssm, flags)
        };
    }

    #[cfg(all(feature = "simd", target_arch = "aarch64", target_feature = "neon"))]
    {
        return unsafe {
            neon::compute_ssm_and_flags(z_afp, z_lds, z_prolif, z_apop, thresholds, out_ssm, flags)
        };
    }

    #[cfg(not(any(
        all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"),
        all(feature = "simd", target_arch = "aarch64", target_feature = "neon")
    )))]
    compute_ssm_and_flags_scalar(z_afp, z_lds, z_prolif, z_apop, thresholds, out_ssm, flags);
}

pub fn compute_z_scores_scalar(values: &[f32], mean: f32, inv_std: f32, out: &mut [f32]) {
    for (idx, value) in values.iter().enumerate() {
        out[idx] = (*value - mean) * inv_std;
    }
}

pub fn compute_ssm_and_flags_scalar(
    z_afp: &[f32],
    z_lds: &[f32],
    z_prolif: &[f32],
    z_apop: &[f32],
    thresholds: &Thresholds,
    out_ssm: &mut [f32],
    flags: &mut SurvivalFlags<'_>,
) {
    for idx in 0..out_ssm.len() {
        let ssm = 0.35 * z_afp[idx] + 0.35 * z_lds[idx] - 0.15 * z_prolif[idx] - 0.15 * z_apop[idx];
        out_ssm[idx] = ssm;
        flags.ssm_high[idx] = ssm >= thresholds.ssm_high;
        flags.afp_high[idx] = z_afp[idx] >= thresholds.afp_high;
        flags.lds_high[idx] = z_lds[idx] >= thresholds.lds_high;
        flags.prolif_low[idx] = z_prolif[idx] <= thresholds.prolif_low;
        flags.apop_low[idx] = z_apop[idx] <= thresholds.apop_low;
    }
}

#[cfg(all(feature = "simd", target_arch = "x86_64", target_feature = "avx2"))]
mod avx2 {
    #![allow(unsafe_op_in_unsafe_fn)]

    use super::*;
    use std::arch::x86_64::*;

    pub unsafe fn compute_z_scores(values: &[f32], mean: f32, inv_std: f32, out: &mut [f32]) {
        let len = values.len();
        let mut idx = 0;
        let v_mean = _mm256_set1_ps(mean);
        let v_inv = _mm256_set1_ps(inv_std);

        while idx + 8 <= len {
            let v = _mm256_loadu_ps(values.as_ptr().add(idx));
            let v = _mm256_mul_ps(_mm256_sub_ps(v, v_mean), v_inv);
            _mm256_storeu_ps(out.as_mut_ptr().add(idx), v);
            idx += 8;
        }

        for i in idx..len {
            out[i] = (values[i] - mean) * inv_std;
        }
    }

    pub unsafe fn compute_ssm_and_flags(
        z_afp: &[f32],
        z_lds: &[f32],
        z_prolif: &[f32],
        z_apop: &[f32],
        thresholds: &Thresholds,
        out_ssm: &mut [f32],
        flags: &mut SurvivalFlags<'_>,
    ) {
        let len = out_ssm.len();
        let mut idx = 0;

        let c_afp = _mm256_set1_ps(0.35);
        let c_lds = _mm256_set1_ps(0.35);
        let c_pro = _mm256_set1_ps(-0.15);
        let c_apop = _mm256_set1_ps(-0.15);

        let t_ssm = _mm256_set1_ps(thresholds.ssm_high);
        let t_afp = _mm256_set1_ps(thresholds.afp_high);
        let t_lds = _mm256_set1_ps(thresholds.lds_high);
        let t_pro = _mm256_set1_ps(thresholds.prolif_low);
        let t_apop = _mm256_set1_ps(thresholds.apop_low);

        while idx + 8 <= len {
            let v_afp = _mm256_loadu_ps(z_afp.as_ptr().add(idx));
            let v_lds = _mm256_loadu_ps(z_lds.as_ptr().add(idx));
            let v_pro = _mm256_loadu_ps(z_prolif.as_ptr().add(idx));
            let v_apop = _mm256_loadu_ps(z_apop.as_ptr().add(idx));

            let ssm = _mm256_add_ps(
                _mm256_add_ps(_mm256_mul_ps(v_afp, c_afp), _mm256_mul_ps(v_lds, c_lds)),
                _mm256_add_ps(_mm256_mul_ps(v_pro, c_pro), _mm256_mul_ps(v_apop, c_apop)),
            );

            _mm256_storeu_ps(out_ssm.as_mut_ptr().add(idx), ssm);

            let m_ssm = _mm256_cmp_ps(ssm, t_ssm, _CMP_GE_OQ);
            let m_afp = _mm256_cmp_ps(v_afp, t_afp, _CMP_GE_OQ);
            let m_lds = _mm256_cmp_ps(v_lds, t_lds, _CMP_GE_OQ);
            let m_pro = _mm256_cmp_ps(v_pro, t_pro, _CMP_LE_OQ);
            let m_apop = _mm256_cmp_ps(v_apop, t_apop, _CMP_LE_OQ);

            let mask_ssm = _mm256_movemask_ps(m_ssm) as u32;
            let mask_afp = _mm256_movemask_ps(m_afp) as u32;
            let mask_lds = _mm256_movemask_ps(m_lds) as u32;
            let mask_pro = _mm256_movemask_ps(m_pro) as u32;
            let mask_apop = _mm256_movemask_ps(m_apop) as u32;

            for lane in 0..8 {
                let bit = 1u32 << lane;
                let out_idx = idx + lane;
                flags.ssm_high[out_idx] = (mask_ssm & bit) != 0;
                flags.afp_high[out_idx] = (mask_afp & bit) != 0;
                flags.lds_high[out_idx] = (mask_lds & bit) != 0;
                flags.prolif_low[out_idx] = (mask_pro & bit) != 0;
                flags.apop_low[out_idx] = (mask_apop & bit) != 0;
            }

            idx += 8;
        }

        for i in idx..len {
            let ssm = 0.35 * z_afp[i] + 0.35 * z_lds[i] - 0.15 * z_prolif[i] - 0.15 * z_apop[i];
            out_ssm[i] = ssm;
            flags.ssm_high[i] = ssm >= thresholds.ssm_high;
            flags.afp_high[i] = z_afp[i] >= thresholds.afp_high;
            flags.lds_high[i] = z_lds[i] >= thresholds.lds_high;
            flags.prolif_low[i] = z_prolif[i] <= thresholds.prolif_low;
            flags.apop_low[i] = z_apop[i] <= thresholds.apop_low;
        }
    }
}

#[cfg(all(feature = "simd", target_arch = "aarch64", target_feature = "neon"))]
mod neon {
    #![allow(unsafe_op_in_unsafe_fn)]

    use super::*;
    use std::arch::aarch64::*;

    pub unsafe fn compute_z_scores(values: &[f32], mean: f32, inv_std: f32, out: &mut [f32]) {
        let len = values.len();
        let mut idx = 0;
        let v_mean = vdupq_n_f32(mean);
        let v_inv = vdupq_n_f32(inv_std);

        while idx + 4 <= len {
            let v = vld1q_f32(values.as_ptr().add(idx));
            let v = vmulq_f32(vsubq_f32(v, v_mean), v_inv);
            vst1q_f32(out.as_mut_ptr().add(idx), v);
            idx += 4;
        }

        for i in idx..len {
            out[i] = (values[i] - mean) * inv_std;
        }
    }

    pub unsafe fn compute_ssm_and_flags(
        z_afp: &[f32],
        z_lds: &[f32],
        z_prolif: &[f32],
        z_apop: &[f32],
        thresholds: &Thresholds,
        out_ssm: &mut [f32],
        flags: &mut SurvivalFlags<'_>,
    ) {
        let len = out_ssm.len();
        let mut idx = 0;

        let c_afp = vdupq_n_f32(0.35);
        let c_lds = vdupq_n_f32(0.35);
        let c_pro = vdupq_n_f32(-0.15);
        let c_apop = vdupq_n_f32(-0.15);

        let t_ssm = vdupq_n_f32(thresholds.ssm_high);
        let t_afp = vdupq_n_f32(thresholds.afp_high);
        let t_lds = vdupq_n_f32(thresholds.lds_high);
        let t_pro = vdupq_n_f32(thresholds.prolif_low);
        let t_apop = vdupq_n_f32(thresholds.apop_low);

        while idx + 4 <= len {
            let v_afp = vld1q_f32(z_afp.as_ptr().add(idx));
            let v_lds = vld1q_f32(z_lds.as_ptr().add(idx));
            let v_pro = vld1q_f32(z_prolif.as_ptr().add(idx));
            let v_apop = vld1q_f32(z_apop.as_ptr().add(idx));

            let ssm = vaddq_f32(
                vaddq_f32(vmulq_f32(v_afp, c_afp), vmulq_f32(v_lds, c_lds)),
                vaddq_f32(vmulq_f32(v_pro, c_pro), vmulq_f32(v_apop, c_apop)),
            );

            vst1q_f32(out_ssm.as_mut_ptr().add(idx), ssm);

            let m_ssm = vcgeq_f32(ssm, t_ssm);
            let m_afp = vcgeq_f32(v_afp, t_afp);
            let m_lds = vcgeq_f32(v_lds, t_lds);
            let m_pro = vcleq_f32(v_pro, t_pro);
            let m_apop = vcleq_f32(v_apop, t_apop);

            let mut lanes_ssm = [0u32; 4];
            let mut lanes_afp = [0u32; 4];
            let mut lanes_lds = [0u32; 4];
            let mut lanes_pro = [0u32; 4];
            let mut lanes_apop = [0u32; 4];
            vst1q_u32(lanes_ssm.as_mut_ptr(), m_ssm);
            vst1q_u32(lanes_afp.as_mut_ptr(), m_afp);
            vst1q_u32(lanes_lds.as_mut_ptr(), m_lds);
            vst1q_u32(lanes_pro.as_mut_ptr(), m_pro);
            vst1q_u32(lanes_apop.as_mut_ptr(), m_apop);
            for lane in 0..4 {
                let out_idx = idx + lane;
                flags.ssm_high[out_idx] = lanes_ssm[lane] != 0;
                flags.afp_high[out_idx] = lanes_afp[lane] != 0;
                flags.lds_high[out_idx] = lanes_lds[lane] != 0;
                flags.prolif_low[out_idx] = lanes_pro[lane] != 0;
                flags.apop_low[out_idx] = lanes_apop[lane] != 0;
            }

            idx += 4;
        }

        for i in idx..len {
            let ssm = 0.35 * z_afp[i] + 0.35 * z_lds[i] - 0.15 * z_prolif[i] - 0.15 * z_apop[i];
            out_ssm[i] = ssm;
            flags.ssm_high[i] = ssm >= thresholds.ssm_high;
            flags.afp_high[i] = z_afp[i] >= thresholds.afp_high;
            flags.lds_high[i] = z_lds[i] >= thresholds.lds_high;
            flags.prolif_low[i] = z_prolif[i] <= thresholds.prolif_low;
            flags.apop_low[i] = z_apop[i] <= thresholds.apop_low;
        }
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/math/simd.rs"]
mod tests;
