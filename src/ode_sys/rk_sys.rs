//! Runge-Kutta (RK4) method for solving systems of ordinary differential equations (ODEs).
use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

/// Runge-Kutta (RK4) method for solving systems of Ordinary Differential Equations (ODEs).
///
/// This trait defines the [Runge-Kutta (RK4) method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) for solving systems of ordinary differential equations.
pub trait RungeKuttaODESysSolver<T: ODESYS> {
    /// Solve the system of ODEs using the Runge-Kutta (RK4) method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODESYS` trait.
    /// * `x` - The initial x value.
    /// * `y` - The initial vector of y values (corresponding to the initial x).
    /// * `x_target` - The x value where the solution is desired.
    /// * `h` - The step size or increment for x.
    ///
    /// # Returns
    ///
    /// The vector of estimated y values at `x_target`.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODESYS, ODESysSolver, RungeKuttaODESysSolver};
    ///
    /// struct MyODESystem;
    /// impl ODESYS for MyODESystem {
    ///     fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64> {
    ///         // Define the system of ODEs
    ///         // Example: dy/dx = x * y, dz/dx = x + y
    ///         vec![x * y[0], x + y[1]]
    ///     }
    /// }
    ///
    /// let solver = ODESysSolver;
    /// let my_ode_system = MyODESystem;
    /// let x0 = 0.0;
    /// let y0 = vec![1.0, 2.0];
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.rk_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn rk_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

// Implementing the Runge-Kutta (RK4) method for the system of ODEs Solver
impl<T: ODESYS> RungeKuttaODESysSolver<T> for ODESysSolver {
    /// Implementation of the Runge-Kutta (RK4) method to solve a system of ODEs.
    ///
    /// This method approximates the solution to the system of ODEs at the specified x_target.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODESYS` trait.
    /// * `x` - The initial x value.
    /// * `y` - The initial vector of y values (corresponding to the initial x).
    /// * `x_target` - The x value where the solution is desired.
    /// * `h` - The step size or increment for x.
    ///
    /// # Returns
    ///
    /// The vector of estimated y values at `x_target`.
    /// 
    /// # When to Use: 
    /// 
    /// Suitable for solving non-stiff systems of differential equations.
    /// 
    /// # Pros and Cons:
    /// - Pros: Simple implementation, moderate accuracy for non-stiff problems.
    /// - Cons: May not be as accurate for stiff problems, especially with larger step sizes.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Conditionally stable, suitable for non-stiff problems.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODESYS, ODESysSolver, RungeKuttaODESysSolver};
    ///
    /// struct MyODESystem;
    /// impl ODESYS for MyODESystem {
    ///     fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64> {
    ///         // Define the system of ODEs
    ///         // Example: dy/dx = x * y, dz/dx = x + y
    ///         vec![x * y[0], x + y[1]]
    ///     }
    /// }
    ///
    /// let solver = ODESysSolver;
    /// let my_ode_system = MyODESystem;
    /// let x0 = 0.0;
    /// let y0 = vec![1.0, 2.0];
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.rk_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn rk_solve(&self, ode: &T, mut x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>{
        while x < x_target {
            let result = rk4_step(x, &y, h, &|x, y| ode.eval(&x, &y.to_vec()));
            y.clone_from_slice(&result);
            x += h;
        }
        y
    }
}

/// Function to compute a single RK4 step for a system of ODEs.
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
