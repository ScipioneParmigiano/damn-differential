use super::ODE;

/// Represents a solver for Ordinary Differential Equations (ODEs) using the Runge-Kutta method.
pub struct RungeKuttaSolver;

/// Defines methods for solving ODEs using the Runge-Kutta method.
pub trait RungeKuttaODESolver {
    /// Solves an initial value problem (IVP) using the Runge-Kutta 4th order method.
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
    fn rk4_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;

    /// Solves an initial value problem (IVP) using the Runge-Kutta 2nd order method.
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
    fn rk2_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl RungeKuttaODESolver for RungeKuttaSolver {
    /// Solves an initial value problem (IVP) using the Runge-Kutta 4th order method.
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

    /// Solves an initial value problem (IVP) using the Runge-Kutta 2nd order method.
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
