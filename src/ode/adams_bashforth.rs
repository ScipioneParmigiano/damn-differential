use super::ODE;

pub struct ABSolver;

pub trait ABODESolver {
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64;
}

impl ABODESolver for ABSolver{
    fn ivp(&self, ode: &dyn ODE, x0: f64, y0: f64, h: f64, x_target: f64) -> f64 {
        let mut x = x0;
        let mut y = y0;
        
        while x < x_target {
            let f0 = ode.eval(x, y);
            let f1 = ode.eval(x + h, y + h * f0);
            y += h * (1.5 * f1 - 0.5 * f0);
            x += h;
        }
    
        y
    }
}