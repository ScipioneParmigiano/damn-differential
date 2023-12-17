pub mod euler;
pub mod rk;
pub mod rkf;
pub mod adams_bashforth;
pub mod adams_moulton;
pub mod heun;
pub mod bogacki_shampine;

/// Trait defining the ODE
pub trait ODE {
    fn eval(&self, x: f64, y: f64) -> f64;
}

/// Struct implementing the solver for an ODE. It has various function associated with it, defining the specifict method to use 
pub struct ODESolver;