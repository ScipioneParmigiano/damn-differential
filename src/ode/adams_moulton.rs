use super::ODE;

pub struct AMSolver;

pub trait AMODESolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl AMODESolver for AMSolver{
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
    let mut y = y0;

    while x < x_target {
        let f0 = ode.eval(x, y);
        let f1 = ode.eval(x + h, y + h * f0);
        y += h * (5.0 * f1 + 8.0 * f0) / 12.0;
        x += h;
    }

    y
    }
}