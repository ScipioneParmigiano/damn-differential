mod ode;
mod ode_sys;

use ode::*;
use ode_sys::*;

use euler::*;
use rk::*;
use rkf::*;
use adams_bashforth::*;
use adams_moulton::*;
use qss::*;
use heun::*;
use bogacki_shampine::*;

use rk_sys::*; 
use leapfrog::*;
use forest_ruth::*;
use yoshida4::*;
use qss_sys::*;
use euler_sys::*;
use radau::*;
use verlet::*;





fn main() {
    println!("\nMy ODE:\n");

    struct MyODE;
    impl ODE for MyODE {
        fn eval(&self, x: f64, y: f64) -> f64 {
            y.sin()-x.cos()
        }
    }

    {        
        let my_ode = MyODE;
        let solver = ODESolver ;

        let x0 = 0.0;
        let y0 = 1.0;
        let h = 0.01;
        let x_target = 2000.0;
        
        let result = solver.eu_ivp(&my_ode, x0, y0, h, x_target);
        println!("Euler's method: {}", result);
        let result = solver.rk2_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK2 method: {}", result);
        let result = solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
        println!("RK4 method: {}", result);
        let result = solver.rkf_ivp(&my_ode, x0, y0, h, x_target);
        println!("RKF method: {}", result);
        let result = solver.ab_ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Bashforth's: {}", result);
        let result = solver.am_ivp(&my_ode, x0, y0, h, x_target);
        println!("Adams-Moulton's: {}", result);
        let result = solver.he_ivp(&my_ode, x0, y0, h, x_target);
        println!("Heun's method: {}", result);
        let result = solver.bs_ivp(&my_ode, x0, y0, h, x_target);
        println!("Bogacki-Shampine method: {}", result);
    
    }

    println!("\nLotka-Volterra:\n");
    
    struct LotkaVolterra;
    impl ODESYS for LotkaVolterra {
        fn eval(&self, x: &f64, y: &Vec<f64>) -> Vec<f64> {
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
        let solver=ODESysSolver;
        
        let lv_equation = LotkaVolterra;
        
        let x_initial = 40.0; 
        let y_initial = 10.0; 
        let initial_conditions = vec![x_initial, y_initial];
        
        let start_time = 0.0;
        let end_time = 40.0;
        let h = 0.01;


        let result = solver.eu_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Euler: {:?}", result);
        let result = solver.rk_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("RK: {:?}", result);
        let result = solver.fr_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("FR: {:?}", result);
        let result = solver.y4_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Y4: {:?}", result);
        let result = solver.lf_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Leapfrog: {:?}", result);
        let result = solver.ia_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IA: {:?}", result);
        let result = solver.iia_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IIA: {:?}", result);
        let result = solver.ve_solve(&lv_equation, start_time, initial_conditions.clone(), end_time, h);
        println!("Verlet: {:?}", result);
    }

    println!("\nLorentz system:\n");

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
        let solver=ODESysSolver;
        
        let lorentz_system = LorentzSystem;
        let initial_conditions = vec![1.0, 1.0, 1.0];
        
        let start_time = 0.0;
        let end_time = 5.0;
        let h = 0.01;

        let result = solver.eu_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Euler: {:?}", result);
        let result = solver.rk_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("RK: {:?}", result);
        let result = solver.fr_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("FR: {:?}", result);
        let result = solver.y4_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Y4: {:?}", result);
        let result = solver.lf_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Leapfrog: {:?}", result);
        let result = solver.ia_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IA: {:?}", result);
        let result = solver.iia_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IIA: {:?}", result);
        let result = solver.ve_solve(&lorentz_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Verlet: {:?}", result);

    }

    println!("\nMy system:\n");


    struct MySys;
    impl ODESYS for MySys {
        fn eval(&self, _x: &f64, y: &Vec<f64>) -> Vec<f64> {

            let dx = -y[0];
            let dy = -y[1];

            vec![dx, dy]
        }
    }

    {
        let solver = ODESysSolver;  

        let my_system = MySys;
        let initial_conditions = vec![1.0, 1.0];
        
        let start_time = 0.0;
        let end_time = 1.0;
        let h = 0.080;

        let result = solver.eu_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Euler: {:?}", result);
        let result = solver.rk_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("RK: {:?}", result);
        let result = solver.fr_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("FR: {:?}", result);
        let result = solver.y4_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Y4: {:?}", result);
        let result = solver.lf_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Leapfrog: {:?}", result);
        let result = solver.ia_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IA: {:?}", result);
        let result = solver.iia_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Radau IIA: {:?}", result);
        let result = solver.ve_solve(&my_system, start_time, initial_conditions.clone(), end_time, h);
        println!("Verlet: {:?}", result);
        
    }
}