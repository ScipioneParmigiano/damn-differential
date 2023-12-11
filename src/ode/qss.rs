use super::ODE;

pub trait QSSSolverODE {
    fn ivp_qss1(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64;
    fn ivp_qss2(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64;
    fn ivp_qss3(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64;
}

pub struct QSSSolver {}

impl QSSSolver {
    fn update(&self, ode: &dyn ODE, x: f64, y: f64, h: f64, delta_q: f64) -> f64 {
        let k1 = ode.eval(x, y);
        let q_prev = y;

        let updated_value = q_prev + h * k1;

        let q_current = if (updated_value - q_prev).abs() >= delta_q {
            updated_value
        } else {
            q_prev
        };

        q_current
    }

    fn update2(&self, ode: &dyn ODE, x: f64, y: f64, h: f64, delta_q: f64) -> f64 {
        let k1 = ode.eval(x, y);
        let k2 = ode.eval(x + h, y + h * k1);
        let q_prev = y;
        
        let updated_value = q_prev+ h * (k1 + k2) / 2.0;

        let q_current = if (updated_value - q_prev).abs() >= delta_q {
            updated_value
        } else {
            q_prev
        };

        q_current
    }

    fn update3(&self, ode: &dyn ODE, x: f64, y: f64, h: f64, delta_q: f64) -> f64 {
        let k1 = ode.eval(x, y);
        let k2 = ode.eval(x + h / 2.0, y + (h / 2.0) * k1);
        let k3 = ode.eval(x + h, y + h * (-k1 + 2.0 * k2));
        let q_prev = y;

        let updated_value = q_prev + h * (k1 + 4.0 * k2 + k3) / 6.0;

        let q_current = if (updated_value - q_prev).abs() >= delta_q {
            updated_value
        } else {
            q_prev
        };

        q_current
    }
}

impl QSSSolverODE for QSSSolver {
    fn ivp_qss1(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            y = self.update(ode, x, y, h, delta_q);
            x += h;
        }

        y
    }

    fn ivp_qss2(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            y = self.update2(ode, x, y, h, delta_q);
            x += h;
        }

        y
    }

    fn ivp_qss3(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            y = self.update3(ode, x, y, h, delta_q);
            x += h;
        }

        y
    }

}
