use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

pub trait RadauODESysSolver<T: ODESYS> {
    fn ia_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

impl<T: ODESYS> RadauODESysSolver<T> for ODESysSolver {
    fn ia_solve(&self, ode: &T, mut x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        while x < x_target {

            let k1 = ode.eval(&x, &add_vec(&y, &vec![h; y.len()]));
             
            y = add_vec(&y, &vec_scalar_mul(&k1, h));
            x += h;
        }
        y
    }     
}
