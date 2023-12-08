use super::ODE;

pub struct QSSSolver;

pub trait QSSODESolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl QSSODESolver for QSSSolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;
        let mut step_size = h;

        while x < x_target {
            let y_old = y;
            let mut done = false;

            while !done {
                let k1 = step_size * ode.eval(x, y);
                let k2 = step_size * ode.eval(x + step_size / 2.0, y + k1 / 2.0);
                let k3 = step_size * ode.eval(x + step_size / 2.0, y + k2 / 2.0);
                let k4 = step_size * ode.eval(x + step_size, y + k3);

                let y_n_plus_1 = y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;

                let error = (y_n_plus_1 - y_old).abs();
                if error < 1e-6 { // Set your tolerance here
                    x += step_size;
                    y = y_n_plus_1;
                    done = true;
                } else {
                    step_size *= 0.9 * 1e-6 / error.powf(0.2); // Adjust the tolerance here
                    if step_size < 1e-6 { // Set your minimum step size here
                        done = true;
                    }
                }
            }
        }

        y
    }
}
