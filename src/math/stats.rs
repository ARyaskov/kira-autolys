#[derive(Debug, Clone, Copy)]
pub struct RunningStats {
    n: u32,
    mean: f32,
    m2: f32,
}

impl RunningStats {
    pub fn new() -> Self {
        Self {
            n: 0,
            mean: 0.0,
            m2: 0.0,
        }
    }

    pub fn update(&mut self, value: f32) {
        self.n = self.n.saturating_add(1);
        let n_f = self.n as f32;
        let delta = value - self.mean;
        self.mean += delta / n_f;
        let delta2 = value - self.mean;
        self.m2 += delta * delta2;
    }

    pub fn mean(&self) -> f32 {
        if self.n == 0 { 0.0 } else { self.mean }
    }

    pub fn variance(&self) -> f32 {
        if self.n < 2 {
            0.0
        } else {
            self.m2 / (self.n as f32)
        }
    }

    pub fn std(&self) -> f32 {
        self.variance().sqrt()
    }
}

#[cfg(test)]
#[path = "../../tests/src_inline/math/stats.rs"]
mod tests;
