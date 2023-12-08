use super::ODE;

/// Struct representing the Runge-Kutta-Felbergh solver.
pub struct RKFSolver;

/// Trait defining methods for a Runge-Kutta-Felbergh ODE solver.
pub trait RKFODESolver {
    /// Solves an Initial Value Problem (IVP) for the given ODE using the Runge-Kutta-Felbergh method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation to be solved.
    /// * `x0` - Initial value of the independent variable.
    /// * `y0` - Initial value of the dependent variable.
    /// * `h` - Initial step size for numerical integration.
    /// * `x_target` - Target value of the independent variable for which the solution is computed.
    ///
    /// # Returns
    ///
    /// The approximate solution of the ODE at `x_target`.
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;

    /// Computes a single step using the Runge-Kutta-Felbergh method for the given ODE.
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation for which a single step is computed.
    /// * `x` - Current value of the independent variable.
    /// * `y` - Current value of the dependent variable.
    /// * `h` - Step size for the current step.
    /// * `tolerance` - Tolerance for error control.
    ///
    /// # Returns
    ///
    /// A tuple containing the next approximate value of the dependent variable and the adjusted step size.
    fn step(ode: &dyn ODE, x: f64, y: f64, h: f64, tolerance: f64) -> (f64, f64);
}

impl RKFODESolver for RKFSolver {
    /// Implements the Runge-Kutta-Felbergh method for solving an Initial Value Problem (IVP).
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation to be solved.
    /// * `x0` - Initial value of the independent variable.
    /// * `y0` - Initial value of the dependent variable.
    /// * `h` - Initial step size for numerical integration.
    /// * `x_target` - Target value of the independent variable for which the solution is computed.
    ///
    /// # Returns
    ///
    /// The approximate solution of the ODE at `x_target`.
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
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

    /// Computes a single step using the Runge-Kutta-Felbergh method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation for which a single step is computed.
    /// * `x` - Current value of the independent variable.
    /// * `y` - Current value of the dependent variable.
    /// * `h` - Step size for the current step.
    /// * `tolerance` - Tolerance for error control.
    ///
    /// # Returns
    ///
    /// A tuple containing the next approximate value of the dependent variable and the adjusted step size.
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