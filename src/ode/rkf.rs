//! Runge-Kutta-Fehlberg method
use super::{ODE,ODESolver};

/// Runge-Kutta-Fehlberg Ordinary Differential Equation (ODE) solver trait.
///
/// This trait defines the Runge-Kutta-Fehlberg (RKF) method for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait RKFODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Runge-Kutta-Fehlberg (RKF) method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODE` trait.
    /// * `x0` - The initial x value.
    /// * `y0` - The initial y value (corresponding to the initial x).
    /// * `h` - The initial step size.
    /// * `x_target` - The x value where the solution is desired.
    ///
    /// # Returns
    ///
    /// The estimated y value at `x_target`.
    fn rkf_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;

    /// Perform a single step of the RKF method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODE` trait.
    /// * `x` - The current x value.
    /// * `y` - The current y value.
    /// * `h` - The step size.
    /// * `tolerance` - The error tolerance for adaptive step control.
    ///
    /// # Returns
    ///
    /// A tuple containing the next estimated `y` value and the new step size `h`.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODE, ODESolver, RKFODESolver};
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
    /// let result = solver.rkf_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```

    fn step(ode: &dyn ODE, x: f64, y: f64, h: f64, tolerance: f64) -> (f64, f64);
}

impl RKFODESolver for ODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the Runge-Kutta-Fehlberg (RKF) method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODE` trait.
    /// * `x0` - The initial x value.
    /// * `y0` - The initial y value (corresponding to the initial x).
    /// * `h` - The initial step size.
    /// * `x_target` - The x value where the solution is desired.
    ///
    /// # Returns
    ///
    /// The estimated y value at `x_target`.
    /// 
    /// # When to Use: 
    /// 
    /// Suitable when adaptive step sizes are needed to balance accuracy and efficiency in solving ordinary differential equations.
    /// 
    /// # Pros and Cons:
    /// - Pros: Offers adaptive step-size control for better accuracy, combines higher and lower-order RK methods for efficiency.
    /// - Cons: Computational overhead due to step size adjustments.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Generally stable and suitable for solving stiff equations with appropriate step size control.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODE, ODESolver, RKFODESolver};
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
    /// let result = solver.rkf_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn rkf_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut h = h;
        let tolerance = 1e-6;
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let (y_next, h_new) = Self::step(ode, x, y, h, tolerance);
            y = y_next;
            x += h;

            h = h_new;
        }

        y
    }

    fn step(ode: &dyn ODE, x: f64, y: f64, h: f64, tolerance: f64) -> (f64, f64) {
        let a2 = 1.0 / 4.0;
        let a3 = 3.0 / 8.0;
        let a4 = 12.0 / 13.0;
        let a5 = 1.0;
        let a6 = 1.0 / 2.0;

        let b21 = 1.0 / 4.0;
        let b31 = 3.0 / 32.0;
        let b32 = 9.0 / 32.0;
        let b41 = 1932.0 / 2197.0;
        let b42 = -7200.0 / 2197.0;
        let b43 = 7296.0 / 2197.0;
        let b51 = 439.0 / 216.0;
        let b52 = -8.0;
        let b53 = 3680.0 / 513.0;
        let b54 = -845.0 / 4104.0;
        let b61 = -8.0 / 27.0;
        let b62 = 2.0;
        let b63 = -3544.0 / 2565.0;
        let b64 = 1859.0 / 4104.0;
        let b65 = -11.0 / 40.0;

        let c1 = 25.0 / 216.0;
        let c3 = 1408.0 / 2565.0;
        let c4 = 2197.0 / 4104.0;
        let c5 = -1.0 / 5.0;

        let d1 = 16.0 / 135.0;
        let d3 = 6656.0 / 12825.0;
        let d4 = 28561.0 / 56430.0;
        let d5 = -9.0 / 50.0;
        let d6 = 2.0 / 55.0;

        let k1 = h * ode.eval(x, y);
        let k2 = h * ode.eval(x + a2 * h, y + b21 * k1);
        let k3 = h * ode.eval(x + a3 * h, y + b31 * k1 + b32 * k2);
        let k4 = h * ode.eval(x + a4 * h, y + b41 * k1 + b42 * k2 + b43 * k3);
        let k5 = h * ode.eval(x + a5 * h, y + b51 * k1 + b52 * k2 + b53 * k3 + b54 * k4);
        let k6 = h * ode.eval(x + a6 * h, y + b61 * k1 + b62 * k2 + b63 * k3 + b64 * k4 + b65 * k5);

        let y_next = y + c1 * k1 + c3 * k3 + c4 * k4 + c5 * k5;
        let y_next_star = y + d1 * k1 + d3 * k3 + d4 * k4 + d5 * k5 + d6 * k6;

        let error = (y_next - y_next_star).abs();
        let h_new = h * (tolerance / error).powf(0.2);

        (y_next, h_new)
    }
}