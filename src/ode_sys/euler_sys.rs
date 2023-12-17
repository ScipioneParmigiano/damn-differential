//! Euler method for solving systems of ordinary differential equations (ODEs).
use super::{ODESYS, ODESysSolver, vec_scalar_mul};

/// Euler method for solving systems of Ordinary Differential Equations (ODEs).
///
/// This trait defines the [Euler method](https://en.wikipedia.org/wiki/Euler_method) for solving systems of ordinary differential equations.
pub trait EulerODESysSolver<T: ODESYS> {
    /// Solve the system of ODEs using the Euler method.
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
    /// let result = solver.eu_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn eu_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

// Implementing the Euler method for the system of ODEs Solver
impl<T: ODESYS> EulerODESysSolver<T> for ODESysSolver {
    /// Implementation of the Euler method to solve a system of ODEs.
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
    /// Suitable for solving simple non-stiff systems of differential equations.
    /// 
    /// # Pros and Cons:
    /// - Pros: Simplicity and ease of implementation.
    /// - Cons: Accuracy may be lower, especially for stiff systems.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Unconditionally stable but less accurate for stiff systems.
    /// 
    /// # Example
    ///
    /// ```
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
    /// let result = solver.eu_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn eu_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        let mut x_val = x;

        while x_val < x_target {
            let dy = vec_scalar_mul(&ode.eval(&x_val, &y), h);

            for (idx, val) in dy.iter().enumerate() {
                y[idx] += *val;
            }

            x_val += h;
        }
        y
    }
}
