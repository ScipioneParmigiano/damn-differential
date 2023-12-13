use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

/// Trait for Runge-Kutta method solving for ordinary differential equation systems.
pub trait RungeKuttaODESysSolver<T: ODESYS> {
    /// Performs the Runge-Kutta 4th order method on a system of ordinary differential equations.
    ///
    /// # Arguments
    ///
    /// * `ode` - A reference to the ODESYS trait.
    /// * `x` - The initial value of the independent variable.
    /// * `y` - The initial values of the dependent variables.
    /// * `a` - The start value of the independent variable for the integration.
    /// * `b` - The end value of the independent variable for the integration.
    /// * `n` - The number of steps for the integration.
    ///
    /// # Returns
    ///
    /// A vector containing the solution of the system of ODEs.
    fn rk4_solve(&self, ode: &T, x: f64, y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64>;
}

impl<T: ODESYS> RungeKuttaODESysSolver<T> for ODESysSolver {
    fn rk4_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, a: f64, b: f64, n: usize) -> Vec<f64> {
        let h = (b - a) / n as f64;

        for _i in 0..n as i32 {
            let result = rk4_step(x, &y, h, &|x, y| ode.eval(&x, &y.to_vec()));
            y.clone_from_slice(&result);
        }
        y
    }
}

/// Performs a single step of the Runge-Kutta 4th order method.
fn rk4_step(x: f64, y_n: &[f64], h: f64, f: &dyn Fn(f64, &[f64]) -> Vec<f64>) -> Vec<f64> {
    let k1 = vec_scalar_mul(&f(x, y_n), h);
    let k2 = vec_scalar_mul(&f(x + h / 2.0, &add_vec(y_n, &vec_scalar_mul(&k1, 0.5))), h);
    let k3 = vec_scalar_mul(&f(x + h / 2.0, &add_vec(y_n, &vec_scalar_mul(&k2, 0.5))), h);
    let k4 = vec_scalar_mul(&f(x + h, &add_vec(y_n, &k3)), h);

    let result = add_vec(
        y_n,
        &vec_scalar_mul(
            &add_vec(
                &k1,
                &add_vec(
                    &vec_scalar_mul(&k2, 2.0),
                    &add_vec(&vec_scalar_mul(&k3, 2.0), &k4),
                ),
            ),
            1.0 / 6.0,
        ),
    );

    result
}
