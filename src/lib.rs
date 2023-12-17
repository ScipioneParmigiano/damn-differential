//!  Rust library to efficiently handle inital valued problems with numerical methods. Thanks to [April Rains](https://www.youtube.com/watch?v=YjVT80bShYM) for inspiring the name.
//! # Download
//! To build with `damn-differential` you need to add the library to your project as a dependency.
//! In your terminal run:
//! ```
//! cargo add damn-diff
//! ```
//! 
//! # Example
//! To numerically solve an IVP:
//! - implement the appropriate trait to your ODE;
//! - define the solver;
//! - select the adeguate method based on your IVP.
//! 
//! For example:
//! ```
//! struct MyODE;
//! impl ODE for MyODE {
//!     fn eval(&self, x: f64, y: f64) -> f64 {
//!         // Define the ODE equation, for instance: dy/dx = x + y
//!         x + y
//!     }
//! }
//!
//! let solver = ODESolver;
//! let my_ode = MyODE;
//! let x0 = 0.0;
//! let y0 = 1.0;
//! let h = 0.1;
//! let x_target = 1.0;
//!
//! let result = solver.rk4_ivp(&my_ode, x0, y0, h, x_target);
//! println!("Solution at x = {}: {}", x_target, result);
//! ```
/// Ordinary differential equation
pub mod ode;
/// Systems of ordinary differential equations
pub mod ode_sys;
