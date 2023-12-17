//! Bogacki–Shampine method
use super::{ODE, ODESolver};

/// Bogacki–Shampine Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the [Bogacki–Shampine method](https://en.wikipedia.org/wiki/Bogacki%E2%80%93Shampine_method#:~:text=The%20Bogacki%E2%80%93Shampine%20method%20is%20a%20Runge%E2%80%93Kutta%20method%20of,to%20implement%20adaptive%20step%20size.) for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait BShampineODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Bogacki–Shampine method.
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
    /// let result = solver.bs_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn bs_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing the Bogacki–Shampine method for the ODE Solver
impl BShampineODESolver for ODESolver {
    /// Implementation of the Bogacki–Shampine method to solve an IVP for an ODE.
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
    /// # When to Use: 
    /// 
    /// Well-suited for stiff differential equations and situations requiring adaptive step sizes.
    /// 
    /// # Pros and Cons:
    /// Pros: Adaptive step-size control, maintains good accuracy with fewer evaluations.
    /// Cons: Requires more memory and computational effort compared to fixed-step methods.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Stable and efficient for a wide range of problems, especially those with stiffness.
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
    /// let result = solver.bs_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn bs_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x + h / 2.0, y + k1 / 2.0);
            let k3 = h * ode.eval(x + h * 3.0 / 4.0, y + 3.0 / 4.0 * k2);

            let y_ = y + 2.0 / 3.0 * k1 + 1.0 / 9.0 * k2 + 4.0 / 9.0 * k3;

            let k4 = h * ode.eval(x, y_);

            let q = y + 7.0 / 24.0 * k1 + 1.0 / 4.0 * k2 + 1.0 / 3.0 * k3 + 1.0 / 8.0 * k4;

            y = q;
            x += h;
        }

        y
    }
}
