pub mod euler;
pub mod rk;

pub trait ODE {
    fn eval(&self, x: f64, y: f64) -> f64;
}
