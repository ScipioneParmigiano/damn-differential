use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

pub trait FRODESysSolver<T: ODESYS> {
    fn fr_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
    }

impl<T: ODESYS> FRODESysSolver<T> for ODESysSolver {
    fn fr_solve(&self, ode: &T, mut x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>{
        while x < x_target {
            // Implement Forest-Ruth method for solving the ODE at each step
            
            // Stage 1
            let k1 = vec_scalar_mul( &ode.eval(&x, &y), h); // Evaluate the ODE at (x, y)
            let y_temp_1 = add_vec(&y, &vec_scalar_mul( &k1, 0.5));

            // Stage 2
            let k2 = vec_scalar_mul( &ode.eval(&(x + h*0.5), &y_temp_1), h);
            let y_temp_2 = add_vec(&y, &vec_scalar_mul(&k2, 0.5));

            // Stage 3
            let k3 = vec_scalar_mul( &ode.eval(&(x + h*0.5), &y_temp_2), h);
            let y_temp_3 = add_vec(&y, &vec_scalar_mul( &k3, 2.0));

            // Stage 4
            let k4 = vec_scalar_mul(&ode.eval(&(x + h), &y_temp_3), h);

            // Update y using weighted averages of these stages
            y = add_vec(&y, &vec_scalar_mul(&add_vec(&k1, &add_vec(&vec_scalar_mul(&add_vec(&k2, &k3), 2.0), &k4)), 1.0/6.0));
            x += h;
        }
        y
    }
}
