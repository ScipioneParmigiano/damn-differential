use super::ODE;

pub struct BShampineSolver;

pub trait BShampineODESolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl BShampineODESolver for BShampineSolver {

    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x+h / 2.0, y + k1 / 2.0);
            let k3 = h * ode.eval(x+h * 3.0 / 4.0, y + 3.0 / 4.0 * k2);
            
            let y_ = y + 2.0 / 3.0 * k1 + 1.0 / 9.0 * k2 + 4.0 / 9.0 * k3;
            
            let k4 = h * ode.eval(x, y_);

            let q = y + 7.0 / 24.0 * k1 + 1.0 / 4.0 * k2  + 1.0 / 3.0 * k3 + 1.0 / 8.0 * k4 ;
            
            y = q;
            x += h;
        }

        y
    }
}
