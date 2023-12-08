mod ode;
mod ode_sys;

use ode::{euler, rk, ODE, adams_bashforth, adams_moulton,rkf};
use euler::{EulerODESolver,EulerSolver};
use rk::{RungeKuttaODESolver, RungeKuttaSolver};
use rkf::{RKFODESolver,RKFSolver};
use adams_bashforth::{ABSolver, ABODESolver};
use adams_moulton::{AMSolver, AMODESolver};

use ode_sys::{rk_sys, ODESYS};
use rk_sys::{RungeKuttaODESysSolver, RungeKuttaSysSolver};
struct MyODE;

impl ODE for MyODE {
    fn eval(&self, x: f64, y: f64) -> f64 {
        -y +y*x.powi(2)
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

    // Euler
    {
        let x0 = 0.0;
        let y0 = 3.0;
        let h = 0.01;
        let x_target = 1.0;

        let euler_solver = EulerSolver {};
        let result = euler_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Euler's method: {}", result);
    }

    // Runge-Kutta
    {
        let x0 = 0.0;
        let y0 = 3.0;
        let h = 0.01;
        let x_target = 1.0;

        let rk_solver = RungeKuttaSolver {};
        let result = rk_solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK2 method: {}", result);
        
        let result = rk_solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK4 method: {}", result);
    }
    
    // Runge–Kutta–Fehlberg method
    {
        let x0 = 0.0;
        let y0 = 3.0;
        let h = 0.01;
        let x_target = 1.0;

        let rk_solver = RKFSolver {};
        let result = rk_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("RKF method: {}", result);
    }

    // Adams–Bashforth
    {
        let x0 = 0.0;
        let y0 = 3.0;
        let h = 0.01;
        let x_target = 1.0;

        let ab_solver = ABSolver {};
        let result = ab_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Bashforth's method: {}", result);
    }

    // Adams-Moulton
    {
        let x0 = 0.0;
        let y0 = 3.0;
        let h = 0.01;
        let x_target = 1.0;

        let am_solver = AMSolver {};
        let result = am_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Moulton's method: {}", result);
    }

    // rk4 for sys
    {
        // let mut y: Vec<f64> = vec![0.0; 2];
        // let a: f64 = 0.0;
        // let x: f64 = 0.0;
        // let x_target: f64 = 1.0;
        // let n: usize = 1000;

        // // Initial conditions
        // y[0] = 3.0; // y(0)
        // y[1] = 1.0; // y'(0)

        // let solver = RungeKuttaSysSolver {};
        // let result = solver.rk4sys(&MyODESYS, x, y, a, x_target, n);
        // println!("y(1): {},\nz(1): {}", result[0], result[1]);
    }
}
