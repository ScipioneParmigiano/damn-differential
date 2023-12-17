//! Forest-Ruth method for solving systems of ordinary differential equations (ODEs).
use super::{ODESYS, ODESysSolver, add_vec, vec_scalar_mul};

/// Forest-Ruth method for solving systems of Ordinary Differential Equations (ODEs).
///
/// This trait defines the [Forest-Ruth method](https://en.wikipedia.org/wiki/Symplectic_integrator) for solving systems of ordinary differential equations.
pub trait FRODESysSolver<T: ODESYS> {
    /// Solve the system of ODEs using the Forest-Ruth method.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODESYS` trait.
    /// * `x` - The initial x value.
    /// * `y` - The initial vector of y values (corresponding to the initial x).
    /// * `x_target` - The x value where the solution is desired.
    /// * `h` - The step size or increment for x.
    ///
    /// # Returns
    ///
    /// The vector of estimated y values at `x_target`.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODESYS, ODESysSolver, FRODESysSolver};
    ///
    /// struct MyODESystem;
    /// impl ODESYS for MyODESystem {
    ///     fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64> {
    ///         // Define the system of ODEs
    ///         // Example: dy/dx = x * y, dz/dx = x + y
    ///         vec![x * y[0], x + y[1]]
    ///     }
    /// }
    ///
    /// let solver = ODESysSolver;
    /// let my_ode_system = MyODESystem;
    /// let x0 = 0.0;
    /// let y0 = vec![1.0, 2.0];
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.fr_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn fr_solve(&self, ode: &T, x: f64, y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64>;
}

// Implementing the Forest-Ruth method for the system of ODEs Solver
impl<T: ODESYS> FRODESysSolver<T> for ODESysSolver {
    /// Implementation of the Forest-Ruth method to solve a system of ODEs.
    ///
    /// This method approximates the solution to the system of ODEs at the specified x_target.
    ///
    /// # Arguments
    ///
    /// * `ode` - The ODE object implementing the `ODESYS` trait.
    /// * `x` - The initial x value.
    /// * `y` - The initial vector of y values (corresponding to the initial x).
    /// * `x_target` - The x value where the solution is desired.
    /// * `h` - The step size or increment for x.
    ///
    /// # Returns
    ///
    /// The vector of estimated y values at `x_target`.
    /// 
    /// # When to Use: 
    /// 
    /// Suitable for solving non-stiff systems of differential equations.
    /// 
    /// # Pros and Cons:
    /// - Pros: Higher accuracy compared to some other methods for non-stiff problems.
    /// - Cons: May require more computational effort due to multiple stages.
    /// 
    /// # Stability Analysis: 
    /// 
    /// Conditionally stable, suitable for non-stiff problems.
    /// 
    /// # Example
    ///
    /// ```
    /// use damndiff::{ODESYS, ODESysSolver, FRODESysSolver};
    ///
    /// struct MyODESystem;
    /// impl ODESYS for MyODESystem {
    ///     fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64> {
    ///         // Define the system of ODEs
    ///         // Example: dy/dx = x * y, dz/dx = x + y
    ///         vec![x * y[0], x + y[1]]
    ///     }
    /// }
    ///
    /// let solver = ODESysSolver;
    /// let my_ode_system = MyODESystem;
    /// let x0 = 0.0;
    /// let y0 = vec![1.0, 2.0];
    /// let h = 0.1;
    /// let x_target = 1.0;
    ///
    /// let result = solver.fr_solve(&my_ode_system, x0, y0, x_target, h);
    /// println!("Solution at x = {}: {:?}", x_target, result);
    /// ```
    fn fr_solve(&self, ode: &T, mut x: f64, mut y: Vec<f64>, x_target: f64, h: f64) -> Vec<f64> {
        while x < x_target {
            // Implement Forest-Ruth method for solving the ODE at each step
            
            // Stage 1
            let k1 = vec_scalar_mul(&ode.eval(&x, &y), h);
            let y_temp_1 = add_vec(&y, &vec_scalar_mul(&k1, 0.5));

            // Stage 2
            let k2 = vec_scalar_mul(&ode.eval(&(x + h * 0.5), &y_temp_1), h);
            let y_temp_2 = add_vec(&y, &vec_scalar_mul(&k2, 0.5));

            // Stage 3
            let k3 = vec_scalar_mul(&ode.eval(&(x + h * 0.5), &y_temp_2), h);
            let y_temp_3 = add_vec(&y, &vec_scalar_mul(&k3, 2.0));

            // Stage 4
            let k4 = vec_scalar_mul(&ode.eval(&(x + h), &y_temp_3), h);

            // Update y using weighted averages of these stages
            y = add_vec(
                &y,
                &vec_scalar_mul(
                    &add_vec(
                        &k1,
                        &add_vec(
                            &vec_scalar_mul(&add_vec(&k2, &k3), 2.0),
                            &k4,
                        ),
                    ),
                    1.0 / 6.0,
                ),
            );
            x += h;
        }
        y
    }
}
