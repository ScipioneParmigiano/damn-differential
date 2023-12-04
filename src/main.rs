mod ode;
mod ode_sys;
use euler::EulerODESolver;
use euler::EulerSolver;
use ode::{euler, rk, ODE};
use ode_sys::{rk_sys, ODESYS};
use rk::{RungeKuttaODESolver, RungeKuttaSolver};

use rk_sys::{RungeKuttaODESysSolver, RungeKuttaSysSolver};

struct MyODE;

impl ODE for MyODE {
    fn eval(&self, x: f64, y: f64) -> f64 {
        -y
    }
}

struct MyODESYS;

impl ODESYS for MyODESYS {
    fn eval(&self, _x: &f64, y: &Vec<f64>) -> Vec<f64> {
        let dy_dx = y[1]; // dy/dx = z
        let dz_dx = 6.0 * y[0] - y[1]; // dz/dx = 6y - z
        vec![dy_dx, dz_dx]
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

    // rk4 for sys
    {
        let mut y: Vec<f64> = vec![0.0; 2];
        let a: f64 = 0.0;
        let x: f64 = 0.0;
        let x_target: f64 = 1.0;
        let n: usize = 1000;

        // Initial conditions
        y[0] = 3.0; // y(0)
        y[1] = 1.0; // y'(0)

        let rks_solver = RungeKuttaSysSolver {};
        let result_rks = rks_solver.rk4sys(&MyODESYS, x, y, a, x_target, n);
        println!("y(1): {},\nz(1): {}", result_rks[0], result_rks[1]);
    }
}
