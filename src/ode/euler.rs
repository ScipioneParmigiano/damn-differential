use super::ODE;

/// EulerSolver represents a numerical method for solving Ordinary Differential Equations (ODEs) using Euler's method.
pub struct EulerSolver;

/// EulerODESolver defines methods for solving ODEs using Euler's method.
pub trait EulerODESolver {
    /// Performs an initial value problem (IVP) solving using Euler's method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE (Ordinary Differential Equation) to be solved.
    /// * `x0` - The initial value of the independent variable.
    /// * `y0` - The initial value of the dependent variable.
    /// * `h` - The step size for the approximation.
    /// * `x_target` - The target value of the independent variable to approximate.
    ///
    /// # Returns
    ///
    /// The approximate value of the dependent variable at `x_target`.
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl EulerODESolver for EulerSolver {
    /// Performs an initial value problem (IVP) solving using Euler's method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE (Ordinary Differential Equation) to be solved.
    /// * `x0` - The initial value of the independent variable.
    /// * `y0` - The initial value of the dependent variable.
    /// * `h` - The step size for the approximation.
    /// * `x_target` - The target value of the independent variable to approximate.
    ///
    /// # Returns
    ///
    /// The approximate value of the dependent variable at `x_target`.
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
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
