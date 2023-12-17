//! Leapfrog method for solving systems of ordinary differential equations (ODEs).
use super::{ODESYS, ODESysSolver, add_vec};

/// Leapfrog method for solving systems of Ordinary Differential Equations (ODEs).
///
/// This trait defines the [Leapfrog method](https://young.physics.ucsc.edu/115/leapfrog.pdf) for solving systems of ordinary differential equations.
pub trait LeapfrogODESysSolver<T: ODESYS> {
    /// Solve the system of ODEs using the Leapfrog method.
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
    /// let result = solver.lf_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn lf_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

// Implementing the Leapfrog method for the system of ODEs Solver
impl<T: ODESYS> LeapfrogODESysSolver<T> for ODESysSolver {
    /// Implementation of the Leapfrog method to solve a system of ODEs.
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
    /// - Pros: Simple implementation, stability in certain scenarios.
    /// - Cons: May not be as accurate as higher-order methods.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Conditionally stable, suitable for non-stiff problems.
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
    /// let result = solver.lf_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn lf_solve(&self, ode: &T, x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
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
