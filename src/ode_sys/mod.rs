pub mod rk_sys;

pub trait ODESYS {
    fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64>; // param s.a. x, y, z
}
