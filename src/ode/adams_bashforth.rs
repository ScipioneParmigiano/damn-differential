use super::{ODE, ODESolver};

/// Trait defining methods for an Adam-Bashforth ODE solver.
pub trait ABODESolver {
    /// Solves an Initial Value Problem (IVP) for the given ODE.
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation to be solved.
    /// * `x0` - Initial value of the independent variable.
    /// * `y0` - Initial value of the dependent variable.
    /// * `h` - Step size for numerical integration.
    /// * `x_target` - Target value of the independent variable for which the solution is computed.
    ///
    /// # Returns
    ///
    /// The approximate solution of the ODE at `x_target`.
    fn ab_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl ABODESolver for ODESolver {
    /// Implements the Adam-Bashforth method for solving an Initial Value Problem (IVP).
    ///
    /// # Arguments
    ///
    /// * `ode` - The Ordinary Differential Equation to be solved.
    /// * `x0` - Initial value of the independent variable.
    /// * `y0` - Initial value of the dependent variable.
    /// * `h` - Step size for numerical integration.
    /// * `x_target` - Target value of the independent variable for which the solution is computed.
    ///
    /// # Returns
    ///
    /// The approximate solution of the ODE at `x_target`.
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
