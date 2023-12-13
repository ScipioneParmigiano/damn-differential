use super::{ODESYS, ODESysSolver, add_vec};

pub trait LeapfrogODESysSolver<T: ODESYS> {
    fn lf_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

impl<T: ODESYS> LeapfrogODESysSolver<T> for ODESysSolver {
    fn lf_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>{
        let mut x_val = x;

        while x_val < x_target {
            let dy = ode.eval(&x_val, &y);

            let mut half_dy = vec![0.0; y.len()];
            for (idx, val) in dy.iter().enumerate() {
                half_dy[idx] = val * h / 2.0;
            }

            y = add_vec(&y, &half_dy);
            x_val += h;

            let next_dy = ode.eval(&x_val, &y);

            let mut next_half_dy = vec![0.0; y.len()];
            for (idx, val) in next_dy.iter().enumerate() {
                next_half_dy[idx] = val * h / 2.0;
            }

            y = add_vec(&y, &next_half_dy);
        }
        y
    }
}
