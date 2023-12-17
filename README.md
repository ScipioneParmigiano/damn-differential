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
    use rustonomicon_optima::*
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

### Systems of ODE
 - Euler method;
 - Forest-Ruth integrator;
 - Leapfrog integration;
 - Radau methods IA;
 - 4th order Runge-Kutta method;

## Disclaimer
All claims, content, designs, algorithms, and specifications described in this project are done with the author's best effort. It is up to the reader to check and validate their accuracy and truthfulness. The author is not responsible for any damage that the use of the library content may cause to the user.
