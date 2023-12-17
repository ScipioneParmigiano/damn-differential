//! Heun method
use super::{ODE, ODESolver};

/// Heun Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the Heun's method for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait HeunODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Heun's method.
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
    /// use super::{ODE, ODESolver, HeunODESolver};
    ///
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
    /// let result = solver.he_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn he_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing the Heun's method for the ODE Solver
impl HeunODESolver for ODESolver {
    /// Implementation of the Heun's method to solve an IVP for an ODE.
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
    /// A slightly improved version of Euler's method, used when a bit more accuracy is required without much added complexity.
    /// 
    /// # Pros and Cons:
    /// - Pros: Improved accuracy over Euler's method.
    /// - Cons: Still not as accurate as higher-order methods, may need smaller step sizes for accurate results.
    /// 
    /// # Stability Analysis: 
    /// 
    /// More stable than Euler's method but can still suffer from stability issues for certain problems. 
    /// 
    /// # Example
    ///
    /// ```
    /// use super::{ODE, ODESolver, HeunODESolver};
    ///
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
    /// let result = solver.he_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn he_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let slope = ode.eval(x, y);
            let y_ = y + h * slope;
            
            // The Heun's method formula
            y += (ode.eval(x + h, y_) + ode.eval(x, y)) * h / 2.0;
            x += h;
        }

        y
    }
}
