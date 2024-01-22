//! Quantized state systems method (QSS1)
use super::{ODE, ODESolver};

/// QSS1 solver trait.
///
/// This trait defines the [QSS1](https://www.fceia.unr.edu.ar/~kofman/files/scsc_08_cellier.pdf) for solving initial value problems (IVPs)
/// of ordinary differential equations (ODEs).
pub trait QSSODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the QSS1.
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
    ///     u(x) - q(y)
    ///}
    /// }
    /// 
    /// fn u(t: f64)->f64{
    ///     if t < 1.0 {
    ///         return 0.0
    ///     } else {
    ///         return 10.0
    ///     }
    /// }
    ///
    /// fn q(t: f64)->f64{
    ///     t.trunc() 
    /// }
    ///
    /// let solver = ODESolver; 
    /// let my_ode = MyODE;
    /// let x0 = 0.0;
    /// let y0 = 1.0;
    /// let h = 0.1;
    /// let x_target = 10.0;
    ///
    /// let result = solver.qss1_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn qss1_ivp(&self, mod_ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

// Implementing QSS1 for the ODE Solver
impl QSSODESolver for ODESolver {
    /// Solve the Initial Value Problem (IVP) for an ODE using the QSS1.
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
    /// When the technique is used to solve a stable linear time-invariant (LTI) system, the global error is bounded by a constant that is proportional to the quantum, but independent of the duration of the simulation.
    /// 
    /// # Pros and Cons:
    /// - Pros: Allow for modeling discontinuities in the system due to their discrete-event nature and asynchronous nature, it also has remarkable global stability and error bounds.
    /// - Cons: Complex.
    /// 
    /// # Stability Analysis: 
    /// 
    /// When the technique is used to solve a stable linear time-invariant (LTI) system, the global error is bounded by a constant that is proportional to the quantum, but independent of the duration of the simulation.
    /// 
    /// # Example
    ///
    /// ```
    /// struct MyODE;
    /// impl ODE for MyODE {
    ///     fn eval(&self, x: f64, y: f64) -> f64 {
    ///     u(x) - q(y)
    ///}
    /// }
    /// 
    /// fn u(t: f64)->f64{
    ///     if t < 1.0 {
    ///         return 0.0
    ///     } else {
    ///         return 10.0
    ///     }
    /// }
    ///
    /// fn q(t: f64)->f64{
    ///     t.trunc() 
    /// }
    ///
    /// let solver = ODESolver; 
    /// let my_ode = MyODE;
    /// let x0 = 0.0;
    /// let y0 = 1.0;
    /// let h = 0.1;
    /// let x_target = 10.0;
    ///
    /// let result = solver.qss1_ivp(&my_ode, x0, y0, h, x_target);
    /// println!("Solution at x = {}: {}", x_target, result);
    /// ```
    fn qss1_ivp(&self, mod_ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let slope = mod_ode.eval(x, y);
            y += h * slope;
            x += h;
        }

        y
    }
}