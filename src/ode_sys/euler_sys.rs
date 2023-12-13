use super::{ODESYS, ODESysSolver, vec_scalar_mul};

pub trait EulerODESysSolver<T: ODESYS> {
    fn eu_solve(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
}

impl<T: ODESYS> EulerODESysSolver<T> for ODESysSolver {
    fn eu_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
        let h = (b - a) / n as f64;
        let mut x_val = x;

        for _i in 0..n {
            let dy = vec_scalar_mul(&ode.eval(&x_val, &y), h);

            for (idx, val) in dy.iter().enumerate() {
                y[idx] += *val;
            }

            x_val += h;
        }
        y
    }
}
