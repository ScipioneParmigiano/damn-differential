use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

pub trait LeapfrogODESysSolver<T: ODESYS> {
    fn lf_solve(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
}

impl<T: ODESYS> LeapfrogODESysSolver<T> for ODESysSolver {
    fn lf_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
        let h = (b - a) / n as f64;
        let mut x_val = x;

        for _i in 0..n {
            let mut dy = vec![0.0; y.len()];
            ode.eval(&x_val, &y).iter().enumerate().for_each(|(idx, val)| {
                dy[idx] = val * h / 2.0;
            });

            y = add_vec(&y, &vec_scalar_mul(&dy, 1.0));
            x_val += h;

            let next_dy = ode.eval(&x_val, &y);
            y = add_vec(&y, &vec_scalar_mul(&next_dy, h / 2.0));
        }
        y
    }
}
