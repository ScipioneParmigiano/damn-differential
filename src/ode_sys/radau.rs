use super::{ODESYS, add_vec, vec_scalar_mul};

pub struct RadauSolver;

pub trait RadauODESysSolver<T: ODESYS> {
    fn solve_ia(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
    fn solve_iia(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
}

impl<T: ODESYS> RadauODESysSolver<T> for RadauSolver {
    fn solve_ia(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
            let h = (b - a) / n as f64;
            
            let c1 = (4.0 - 4.0_f64.sqrt()) / 10.0;
            let c2 = (4.0 + 4.0_f64.sqrt()) / 10.0;
            let c3 = (1.0 - 4.0_f64.sqrt()) / 10.0;
            let d1 = (88.0 - 7.0 * 4.0_f64.sqrt()) / 360.0;
            let d2 = (296.0 - 169.0 * 4.0_f64.sqrt()) / 1800.0;
            let d3 = (-2.0 + 3.0 * 4.0_f64.sqrt()) / 225.0;
            
            let mut x_val = x;
            
            for _i in 0..n {
                let mut dy = vec![0.0; y.len()];
                let mut dy_temp = vec![0.0; y.len()];
                
                let k1 = ode.eval(&x_val, &y);
                for (idx, val) in k1.iter().enumerate() {
                    dy[idx] = val * h * c1;
                    dy_temp[idx] = val * h * d1;
                }
    
                let y_temp = add_vec(&y, &vec_scalar_mul(&dy, c2));
                
                let k2 = ode.eval(&(x_val + c1 * h), &y_temp);
                for (idx, val) in k2.iter().enumerate() {
                    dy_temp[idx] += val * h * d2;
                }
    
                y = add_vec(&y, &vec_scalar_mul(&dy_temp, c3));
                
                x_val += h;
            }
            y
        }

        fn solve_iia(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
            let h = (b - a) / n as f64;

            let c1 = 1.0 / (4.0 - 4.0_f64.sqrt());
            let c2 = (6.0 - 4.0_f64.sqrt()) / 10.0;
            let c3 = (6.0 + 4.0_f64.sqrt()) / 10.0;
            let d1 = (88.0 - 7.0 * 4.0_f64.sqrt()) / 360.0;
            let d2 = (296.0 - 169.0 * 4.0_f64.sqrt()) / 1800.0;
            let d3 = (-2.0 + 3.0 * 4.0_f64.sqrt()) / 225.0;

            let mut x_val = x;

            for _i in 0..n {
                let mut dy = vec![0.0; y.len()];
                let mut dy_temp = vec![0.0; y.len()];

                let k1 = ode.eval(&x_val, &y);
                for (idx, val) in k1.iter().enumerate() {
                    dy[idx] = val * h * c1;
                    dy_temp[idx] = val * h * d1;
                }

                let y_temp = add_vec(&y, &vec_scalar_mul(&dy, c2));

                let k2 = ode.eval(&(x_val + c1 * h), &y_temp);
                for (idx, val) in k2.iter().enumerate() {
                    dy_temp[idx] += val * h * d2;
                }

                y = add_vec(&y, &vec_scalar_mul(&dy_temp, c3));

                x_val += h;
            }
            y
    }
}
