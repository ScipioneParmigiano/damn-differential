<div align="center">
  <img src="https://img.shields.io/crates/d/damndiff.svg" alt="Downloads">
   <img src="https://img.shields.io/crates/v/damndiff.svg" alt="Version">
</div>


# Damn-Differential
Rust crate containing some numerical methods for ordinary differential equations (and systems of ordinary differential equations).
Thanks to [April Rains](https://www.youtube.com/watch?v=YjVT80bShYM) for inspiring the name.

## Getting Started
- First, create a project:
    ```
    cargo init your-project-name
    ```
 - Add `damn-diff` as a dependency:
    ```
    cargo add damndiff
    ```
 - Then, in the main file, add
    ```
    use damndiff::*
    ```

## The state of the art
### ODE
 - Adams-Bashforth method;
 - Adams-Moulton method;
 - Bogacki-Shampine method; 
 - Euler method;
 - Heun method; 
 - 2nd order Runge-Kutta method;
 - 4th order Runge-Kutta method;
 - Runge-Kutta-Fehlberg method;
 - Quantize state systems method (QSS1); 

### Systems of ODE
 - Euler method;
 - Forest-Ruth integrator;
 - Leapfrog integration;
 - Radau methods IA;
 - 4th order Runge-Kutta method;

## Future features
We plan to incorporate a wider range of numerical methods to enhance the versatility and robustness of the library and to extend the various equation types to include:
 - Partial Differential Equations (PDE);
 - Stochastic Differential Equations (SDE);
 - Fractional Differential Equations;
 - Variable Order Differential Equations.

## Contributions
Damn-differential welcomes contributions from the community to enhance its features, improve performance, and fix bugs. If you're interested in contributing, feel free to submit pull requests with your improvements.

