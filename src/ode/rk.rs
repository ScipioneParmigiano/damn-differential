use super::ODE;

pub struct RungeKuttaSolver;
pub trait RungeKuttaODESolver {
    fn rk2_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
    fn rk4_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl RungeKuttaODESolver for RungeKuttaSolver {
    fn rk4_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x + 0.5 * h, y + 0.5 * k1);
            let k3 = h * ode.eval(x + 0.5 * h, y + 0.5 * k2);
            let k4 = h * ode.eval(x + h, y + k3);

            let slope = (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

            y = y + slope;
            x = x + h;
        }

        y
    }

    fn rk2_ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            let k1 = h * ode.eval(x, y);
            let k2 = h * ode.eval(x + h, y + k1);

            let slope = 0.5 * (k1 + k2);

            y = y + slope;
            x = x + h;
        }

        y
    }
}
