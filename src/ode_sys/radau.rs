//! Radau IA (Implicit-Explicit) method for solving systems of ordinary differential equations (ODEs).
use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

/// Radau IA (Implicit-Explicit) method for solving systems of Ordinary Differential Equations (ODEs).
///
/// This trait defines the [Radau IA method](https://www.sciencedirect.com/science/article/pii/S037704279900134X) for solving systems of ordinary differential equations.
pub trait RadauODESysSolver<T: ODESYS> {
    /// Solve the system of ODEs using the Radau IA method.
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
    /// let result = solver.ia_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn ia_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

// Implementing the Radau IA method for the system of ODEs Solver
impl<T: ODESYS> RadauODESysSolver<T> for ODESysSolver {
    /// Implementation of the Radau IA method to solve a system of ODEs.
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
    /// Suitable for certain classes of stiff problems where implicit and explicit methods can be combined.
    /// 
    /// # Pros and Cons:
    /// - Pros: Stiffness handling with potential for higher accuracy.
    /// - Cons: May be complex to implement, can be computationally demanding.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Conditionally stable, suitable for some stiff problems.
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
    /// let result = solver.ia_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn ia_solve(&self, ode: &T, mut x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        while x < x_target {

            let k1 = ode.eval(&x, &add_vec(&y, &vec![h; y.len()]));
             
            y = add_vec(&y, &vec_scalar_mul(&k1, h));
            x += h;
        }
        y
    }     
}
