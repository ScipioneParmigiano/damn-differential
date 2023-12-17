//! Adams-Bashforth method
use super::{ODE, ODESolver};

/// Adams-Bashforth Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the [Adams-Bashforth method](https://web.mit.edu/10.001/Web/Course_Notes/Differential_Equations_Notes/node6.html) for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait ABODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Adams-Bashforth method.
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
    /// use damndiff::{ODE, ODESolver, ABODESolver};
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
    /// let result = solver.ab_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn ab_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing the Adams-Bashforth method for the ODE Solver
impl ABODESolver for ODESolver {
    /// Implementation of the Adams-Bashforth method to solve an IVP for an ODE.
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
    /// Suitable for solving non-stiff differential equations.
    /// 
    /// # Pros and Cons:
    /// - Pros: High accuracy for non-stiff problems.
    /// - Cons: Unstable for stiff problems, accuracy depends on the step size.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Conditionally stable, more suited for non-stiff problems.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODE, ODESolver, ABODESolver};
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
    /// let result = solver.ab_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn ab_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let f0 = ode.eval(x, y);
            let f1 = ode.eval(x + h, y + h * f0);
            y += h * (1.5 * f1 - 0.5 * f0);
            x += h;
        }

        y
    }
}
