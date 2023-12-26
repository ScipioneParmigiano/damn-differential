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
