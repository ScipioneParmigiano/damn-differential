mod ode;
use euler::EulerODESolver;
use euler::EulerSolver;
use ode::{euler, rk, ODE};
use rk::{RungeKuttaODESolver, RungeKuttaSolver};

struct MyOde;

impl ODE for MyOde {
    fn eval(&self, x: f64, y: f64) -> f64 {
        x - y
    }
}

pub struct MyODE;

impl ODE for MyODE {
    fn eval(&self, x: f64, y: f64) -> f64 {
        -y
    }
}

fn main() {
    let my_ode = MyODE;

    // euler
    {
        let x0_ivp = 0.0;
        let y0_ivp = 2.0;
        let h_ivp = 0.01;
        let x_target_ivp = 1.0;

        let euler_solver = EulerSolver {};
        let result_ivp_euler = euler_solver.ivp(&my_ode, x0_ivp, y0_ivp, h_ivp, x_target_ivp);
        println!("Euler's method: {}", result_ivp_euler);
    }

    // Runge-Kutta
    {
        let x0 = 0.0;
        let y0 = 2.0;
        let h = 0.01;
        let x_target = 1.0;

        let rk_solver = RungeKuttaSolver {};
        let result_ivp_rk = rk_solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK2 method: {}", result_ivp_rk);

        let result_ivp_rk = rk_solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK4 method: {}", result_ivp_rk);
    }
}
