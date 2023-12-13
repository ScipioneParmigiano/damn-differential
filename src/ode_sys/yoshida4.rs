use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

pub trait Yoshida4thODESysSolver<T: ODESYS> {
    fn y4_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

impl<T: ODESYS> Yoshida4thODESysSolver<T> for ODESysSolver {
    fn y4_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        let c1 = 1.0 / (2.0 - 2.0_f64.powf(1.0 / 3.0));
        let c2 = -2.0_f64.powf(1.0 / 3.0) * c1;
        let d1 = 1.0 - 2.0 * c1;
        let d2 = 1.0 - 2.0 * c2;

        let mut x_val = x;

        while x_val < x_target {
            let mut dy = vec![0.0; y.len()];

            let k1 = ode.eval(&x_val, &y);
            for (idx, val) in k1.iter().enumerate() {
                dy[idx] = val * c1 * h;
            }
            y = add_vec(&y, &vec_scalar_mul(&dy, d1));

            y = add_vec(&y, &vec_scalar_mul(&dy, c2));

            let k2 = ode.eval(&(x_val + c2 * h), &y);
            for (idx, val) in k2.iter().enumerate() {
                dy[idx] = val * c1 * h;
            }
            y = add_vec(&y, &vec_scalar_mul(&dy, d2));

            x_val += h;
        }
        y
    }
}
