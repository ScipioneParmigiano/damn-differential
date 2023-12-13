use super::{ODESYS, ODESysSolver};
pub trait VerletODESysSolver<T: ODESYS> {
    fn ve_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

impl<T: ODESYS> VerletODESysSolver<T> for ODESysSolver {
    fn ve_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        let mut x_val = x;

        let half_dt = h / 2.0;

        while x_val < x_target {
            let mut k1 = vec![0.0; y.len()];
            let mut k2 = vec![0.0; y.len()];

            // First half-step
            let dy = ode.eval(&x_val, &y);
            for (idx, val) in dy.iter().enumerate() {
                k1[idx] = val * half_dt;
                y[idx] += k1[idx];
            }

            // Full step using k1
            x_val += half_dt;

            let dy = ode.eval(&x_val, &y);
            for (idx, val) in dy.iter().enumerate() {
                k2[idx] = val * h;
                y[idx] += k2[idx];
            }

            // Second half-step using k2
            x_val += half_dt;

            let dy = ode.eval(&x_val, &y);
            for (idx, val) in dy.iter().enumerate() {
                k1[idx] += val * half_dt;
                y[idx] += k1[idx];
            }

            x_val += half_dt;
        }

        y
    }
}
