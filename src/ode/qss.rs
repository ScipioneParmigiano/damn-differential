use super::ODE;

pub trait QSSSolverODE {
    fn ivp_qss(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64;

}

pub struct QSSSolver {}

impl QSSSolver {
    fn update(&self, ode: &dyn ODE, x: f64, y: f64, h: f64, delta_q: f64) -> f64 {
        let k1 = ode.eval(x, y);
        let q_prev = y; 

        let q_current = if (y - q_prev).abs() >= delta_q {
            y 
        } else {
            q_prev 
        };
        

        q_current + h * k1 
    }
}

impl QSSSolverODE for QSSSolver {
    fn ivp_qss(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64, delta_q: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;

        while x < x_target {
            y = self.update(ode, x, y, h, delta_q);
            x += h;
        }

        y
    }
}
