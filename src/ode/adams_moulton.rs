use super::ODE;

/// A struct representing the Adam-Moulton solver.
pub struct AMSolver;

/// A trait defining the interface for an Adam-Moulton Ordinary Differential Equation (ODE) solver.
pub trait AMODESolver {
    /// Solves an initial value problem (IVP) for an ODE using the Adam-Moulton method.
    ///
    /// # Arguments
    ///
    /// * `ode` - A reference to an object implementing the `ODE` trait.
    /// * `x0` - The initial x value.
    /// * `y0` - The initial y value.
    /// * `h` - The step size.
    /// * `x_target` - The target x value at which the solution is sought.
    ///
    /// # Returns
    ///
    /// The approximate value of y at `x_target` using the Adam-Moulton method.
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl AMODESolver for AMSolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        // Perform the Adam-Moulton method
        while x < x_target {
            let f0 = ode.eval(x, y);
            let f1 = ode.eval(x + h, y + h * f0);
            y += h * (5.0 * f1 + 8.0 * f0) / 12.0;
            x += h;
        }

        y // Return the approximate value of y at x_target
    }
}
