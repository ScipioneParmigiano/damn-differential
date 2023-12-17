/// A module containing implementations related to solving systems of Ordinary Differential Equations (ODEs).
pub mod rk_sys; 
pub mod leapfrog;
pub mod forest_ruth;
pub mod euler_sys;
pub mod radau;

/// A trait representing a system of Ordinary Differential Equations (ODEs).
pub trait ODESYS {
    /// Evaluates the system of ODEs at a given x and y value.
    ///
    /// # Arguments
    ///
    /// * `x` - The value of the independent variable.
    /// * `y` - A vector containing the values of dependent variables.
    ///
    /// # Returns
    ///
    /// A vector representing the derivatives of the ODE system at the given x and y.
    fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64>;
}

pub struct ODESysSolver;

/// Adds two vectors element-wise.
fn add_vec(a: &[f64], b: &[f64]) -> Vec<f64> {
    a.iter().zip(b.iter()).map(|(&x, &y)| x + y).collect()
}

/// Multiplies a vector by a scalar.
fn vec_scalar_mul(a: &[f64], scalar: f64) -> Vec<f64> {
    a.iter().map(|&x| x * scalar).collect()
}
