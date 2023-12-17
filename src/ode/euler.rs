//! Euler method
use super::{ODE, ODESolver};

/// Euler Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the [Euler method](https://en.wikipedia.org/wiki/Euler_method) for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait EulerODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Euler method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODE` trait.
    /// * `x0` - The initial x value.
    /// * `y0` - The initial y value (corresponding to the initial x).
    /// * `h` - The step size or increment for x.
    /// * `x_target` - The x value where the solution is desired.
    ///
    /// # Returns
    ///
    /// The estimated y value at `x_target`.
    ///
    /// # Example
    ///
    /// ```
    /// struct MyODE;
    /// impl ODE for MyODE {
    ///     fn eval(&self, x: f64, y: f64) -> f64 {
    ///         // Define the ODE equation, for instance: dy/dx = x + y
    ///         x + y
    ///     }
    /// }
    ///
    /// let solver = ODESolver; 
    /// let my_ode = MyODE;
    /// let x0 = 0.0;
    /// let y0 = 1.0;
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.eu_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn eu_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing the Euler method for the ODE Solver
impl EulerODESolver for ODESolver {
    /// Implementation of the Euler method to solve an IVP for an ODE.
    ///
    /// This method approximates the solution to the ODE at the specified x_target.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODE` trait.
    /// * `x0` - The initial x value.
    /// * `y0` - The initial y value (corresponding to the initial x).
    /// * `h` - The step size or increment for x.
    /// * `x_target` - The x value where the solution is desired.
    ///
    /// # Returns
    ///
    /// The estimated y value at `x_target`.
    /// 
    /// # When to Use: 
    /// 
    /// Suitable for educational purposes or quick approximations.
    /// 
    /// # Pros and Cons:
    /// - Pros: Simple to understand.
    /// - Cons: Less accurate, can diverge significantly for certain differential equations.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Unconditionally unstable for stiff or oscillatory problems.
    /// 
    /// # Example
    ///
    /// ```
    /// struct MyODE;
    /// impl ODE for MyODE {
    ///     fn eval(&self, x: f64, y: f64) -> f64 {
    ///         // Define the ODE equation, for instance: dy/dx = x + y
    ///         x + y
    ///     }
    /// }
    ///
    /// let solver = ODESolver; 
    /// let my_ode = MyODE;
    /// let x0 = 0.0;
    /// let y0 = 1.0;
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.eu_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn eu_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let slope = ode.eval(x, y);
            y += h * slope;
            x += h;
        }

        y
    }
}
