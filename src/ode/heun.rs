use super::ODE;

pub struct HeunSolver;

pub trait HeunODESolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl HeunODESolver for HeunSolver {

    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let slope = ode.eval(x, y);
            let y_ = y + h * slope;
            
            y += (ode.eval(x+h, y_) + ode.eval(x, y)) * h / 2.0;
            x += h;
        }

        y
    }
}
