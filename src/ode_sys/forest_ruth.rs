use super::{ODESYS, add_vec, vec_scalar_mul};

pub struct FRSysSolver;

pub trait FRODESysSolver<T: ODESYS> {
    fn solve(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
}

impl<T: ODESYS> FRODESysSolver<T> for FRSysSolver {
    fn solve(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
        let h = (b - a) / n as f64;
        let cbrt2 = 2_f64.sqrt().powf(1.0 / 3.0);

        let c1 = 1.0 / (2.0 - cbrt2);
        let c2 = c1 - cbrt2;
        let d1 = 1.0 / (2.0 * (1.0 - cbrt2));
        let d2 = 1.0 - 2.0 * c1;

        let mut x_val = x;

        for _i in 0..n {
            let mut dy = vec![0.0; y.len()];
            let mut dy_temp = vec![0.0; y.len()];

            let k1 = ode.eval(&x_val, &y);
            for (idx, val) in k1.iter().enumerate() {
                dy[idx] = val * h / 2.0;
                dy_temp[idx] = val * d1 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy, c1));

            let k2 = ode.eval(&(x_val + c1 * h), &y);
            for (idx, val) in k2.iter().enumerate() {
                dy_temp[idx] += val * d2 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy_temp, c2));

            let k3 = ode.eval(&(x_val + c2 * h), &y);
            for (idx, val) in k3.iter().enumerate() {
                dy[idx] = val * h / 2.0;
                dy_temp[idx] = val * d1 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy, c1));

            let k4 = ode.eval(&(x_val + c1 * h), &y);
            for (idx, val) in k4.iter().enumerate() {
                dy_temp[idx] += val * d2 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy_temp, c2));

            let k5 = ode.eval(&(x_val + c2 * h), &y);
            for (idx, val) in k5.iter().enumerate() {
                dy[idx] = val * h / 2.0;
                dy_temp[idx] = val * d1 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy, c1));

            let k6 = ode.eval(&(x_val + c1 * h), &y);
            for (idx, val) in k6.iter().enumerate() {
                dy_temp[idx] += val * d2 * h;
            }

            y = add_vec(&y, &vec_scalar_mul(&dy_temp, c2));

            let k7 = ode.eval(&(x_val + c2 * h), &y);
            for (idx, val) in k7.iter().enumerate() {
                dy[idx] = val * h / 2.0;
            }

            x_val += h;
        }
        y
    }
}
