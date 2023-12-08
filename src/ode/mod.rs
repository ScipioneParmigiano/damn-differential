pub mod euler;
pub mod rk;
pub mod rkf;
pub mod adams_bashforth;
pub mod adams_moulton;
pub mod qss;

pub trait ODE {
    fn eval(&self, x: f64, y: f64) -> f64;
}
