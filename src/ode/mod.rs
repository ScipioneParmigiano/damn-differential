pub mod euler;
pub mod rk;
pub mod rkf;
pub mod adams_bashforth;
pub mod adams_moulton;

// A trait representing an Ordinary Differential Equation (ODE).
pub trait ODE {
    /// Evaluates the ODE at a given x and y value.
    ///
    /// # Arguments
    ///
    /// * `x` - The value of the independent variable.
    /// * `y` - The value of the dependent variable.
    ///
    /// # Returns
    ///
    /// The value of the derivative dy/dx at the given x and y.
    fn eval(&self, x: f64, y: f64) -> f64;
}
