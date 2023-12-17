//! Runge-Kutta method
use super::{ODE, ODESolver};

/// Runge-Kutta Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the [Runge-Kutta methods](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods) for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait RungeKuttaODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the fourth-order Runge-Kutta method (RK4).
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
    /// let result = solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn rk4_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;

    /// Solve the Initial Value Problem (IVP) for an ODE using the second-order Runge-Kutta method (RK2).
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
    /// let result = solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn rk2_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing the Runge-Kutta methods for the ODE Solver
impl RungeKuttaODESolver for ODESolver {
    /// Implementation of the fourth-order Runge-Kutta method (RK4) to solve an IVP for an ODE.
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
    /// Widely used for solving ordinary differential equations with moderate accuracy requirements.
    /// 
    /// # Pros and Cons:
    /// - Pros: Higher accuracy compared to RK2, suitable for a wide range of problems, reasonably efficient.
    /// - Cons: Still less accurate for stiff differential equations, may require smaller step sizes for some functions.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Stable for most problems within its applicable range.
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
    /// let result = solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn rk4_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x + 0.5 * h, y + 0.5 * k1);
            let k3 = h * ode.eval(x + 0.5 * h, y + 0.5 * k2);
            let k4 = h * ode.eval(x + h, y + k3);

            let slope = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

            y += slope;
            x += h;
        }

        y
    }

    /// Implementation of the second-order Runge-Kutta method (RK2) to solve an IVP for an ODE.
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
    /// Suitable for simple problems when higher accuracy is desired than Euler's method but not as computationally expensive as higher-order methods.
    /// 
    /// # Pros and Cons:
    /// - Pros: Fairly accurate, computationally more efficient than higher-order methods.
    /// - Cons: Less accurate than higher-order methods, may require smaller step sizes for accurate results for some functions.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Generally stable for most problems within its range of applicability.
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
    /// let result = solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn rk2_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x + h, y + k1);

            let slope = 0.5 * (k1 + k2);

            y += slope;
            x += h;
        }

        y
    }
}
