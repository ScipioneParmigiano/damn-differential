mod ode;
mod ode_sys;

use ode::{euler, rk, ODE, qss, adams_bashforth, adams_moulton, rkf};
use euler::{EulerODESolver,EulerSolver};
use rk::{RungeKuttaODESolver, RungeKuttaSolver};
use rkf::{RKFODESolver,RKFSolver};
use adams_bashforth::{ABSolver, ABODESolver};
use adams_moulton::{AMSolver, AMODESolver};
use qss::{QSSSolver, QSSSolverODE};

use ode_sys::{euler_sys, rk_sys, ODESYS, leapfrog, forest_ruth, yoshida4, qss_sys};
use euler_sys::{EulerODESysSolver, EulerSysSolver};
use rk_sys::{RungeKuttaODESysSolver, RungeKuttaSysSolver};
use leapfrog::{LeapfrogODESysSolver, LeapfrogSysSolver};
use forest_ruth::{FRODESysSolver, FRSysSolver};
use yoshida4::{Yoshida4thODESysSolver, Yoshida4thSysSolver};
use qss_sys::*;

struct MyODE;

impl ODE for MyODE {
    fn eval(&self, x: f64, y: f64) -> f64 {
        f64::cos(x)
    }
}


fn main() {
    {        
        let my_ode = MyODE;
        let euler_solver = EulerSolver {};
        let rk2_solver = RungeKuttaSolver {};
        let rk4_solver = RungeKuttaSolver {};
        let rkf_solver = RKFSolver {};
        let ab_solver = ABSolver {};
        let am_solver = AMSolver {};
        // let qss_solver = QSSSolver {};

        let x0 = 0.00;
        let y0 = 1.0;
        let h = 0.01;
        let x_target = 5.0;
        // let quantum = 1e-2;

        let result = euler_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Euler's method: {}", result);
        let result = rk2_solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK2 method: {}", result);
        let result = rk4_solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK4 method: {}", result);
        let result = rkf_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("RKF method: {}", result);
        let result = ab_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Bashforth's method: {}", result);
        let result = am_solver.ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Moulton's method: {}", result);
        // let result = qss_solver.ivp_qss1(&my_ode, x0, y0, h, x_target, quantum);
        // println!("QSS1 method: {:?}", result);
        // let result = qss_solver.ivp_qss2(&my_ode, x0, y0, h, x_target, quantum);
        // println!("QSS2 method: {}", result);
        // let result = qss_solver.ivp_qss3(&my_ode, x0, y0, h, x_target, quantum);
        // println!("QSS3 method: {}", result);
    }
    

    struct LotkaVolterra;
    impl ODESYS for LotkaVolterra {
        fn eval(&self, _x: &f64, y: &Vec<f64>) -> Vec<f64> {
            let alpha = 0.1;
            let beta = 0.02;
            let gamma = 0.3;
            let delta = 0.01;
    
            let dx_dt = alpha * y[0] - beta * y[0] * y[1];
            let dy_dt = delta * y[0] * y[1] - gamma * y[1];
    
            vec![dx_dt, dy_dt]
        }
    }

    
    {
        let solver_euler= EulerSysSolver;
        let solver_rk = RungeKuttaSysSolver;
        let solver_fr = FRSysSolver;
        let solver_y4 = Yoshida4thSysSolver;
        let solver_leap = LeapfrogSysSolver;
        // let solver_qss_sys: QSSSysSolver= QSSSysSolver;
        
        let lv_equation = LotkaVolterra;
        
        let x_initial = 40.0; 
        let y_initial = 12.0; 
        let initial_conditions = vec![x_initial, y_initial];
        
        let start_time = 0.0;
        let end_time = 1.0;
        let num_steps = 100;
        
        let result = solver_euler.solve(&lv_equation, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("Euler: {:?}", result);
        let result = solver_rk.solve(&lv_equation, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("RK: {:?}", result);
        let result = solver_fr.solve(&lv_equation, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("FR: {:?}", result);
        let result = solver_y4.solve(&lv_equation, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("Y4: {:?}", result);
        let result = solver_leap.solve(&lv_equation, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("Leapfrog: {:?}", result);
        // let result = solver_qss_sys.solve(&lv_equation, 0.0, initial_conditions, start_time, end_time, num_steps, 1e-3);
        // println!("QSS: {:?}", result);
    }

    struct LorentzSystem;

    impl ODESYS for LorentzSystem {
        fn eval(&self, _x: &f64, y: &Vec<f64>) -> Vec<f64> {
            let sigma = 10.0;
            let rho = 28.0;
            let beta = 8.0 / 3.0;

            let dx = sigma * (y[1] - y[0]);
            let dy = y[0] * (rho - y[2]) - y[1];
            let dz = y[0] * y[1] - beta * y[2];

            vec![dx, dy, dz]
        }
    }

    {
        let solver_euler= EulerSysSolver;
        let solver_rk = RungeKuttaSysSolver;
        let solver_fr = FRSysSolver;
        let solver_leap = LeapfrogSysSolver;
        // let solver_qss_sys: QSSSysSolver= QSSSysSolver;        


        let lorentz_system = LorentzSystem;
        let initial_conditions = vec![1.0, 1.0, 1.0];
        
        let start_time = 0.0;
        let end_time = 1.0;
        let num_steps = 1000;

        let result = solver_euler.solve(&lorentz_system, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("Euler: {:?}", result);
        let result = solver_rk.solve(&lorentz_system, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("RK: {:?}", result);
        let result = solver_fr.solve(&lorentz_system, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("FR: {:?}", result);
        let result = solver_leap.solve(&lorentz_system, 0.0, initial_conditions.clone(), start_time, end_time, num_steps);
        println!("Leapfrog: {:?}", result);
        // let result = solver_qss_sys.solve(&lorentz_system, 0.0, initial_conditions, start_time, end_time, num_steps, 1e-3);
        // println!("QSS: {:?}", result);
    }

}


