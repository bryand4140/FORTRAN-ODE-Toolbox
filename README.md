![FORTRAN ODE Toolbox](media/logo.png)

## General Description:
A collection of codes for solving linear and nonlinear ordinary differential equations (ODEs).  
A work in progress.  

## Compiling and Executing
A Makefile is provided to build the code and has been written to be compatible with Windows and Linux OS. The Makefile has been  
configured to use `main` or `test` as the default executable. To build the code, use:
```
make main.exe
```
By default, the code is built with double precision (`real64`). This can be changed within the module `MOD_Select_Real_Kind.f90`.  
After the code is built, module and executable files will be placed in a build directory. Use:
````
./main
````
to run the code with your main program. Additionally, an `fpm.toml` file has been provided for users of the Fortran Package  
Manager (FPM) [Fortran Package Manager](https://github.com/fortran-lang/fpm).  

## Repository Structure
The following module files are used:
1. `MOD_IO_Toolbox.f90`
2. `MOD_ODE_Systems.f90`
3. `MOD_ODE_Toolbox.f90`
4. `MOD_Select_Kind.f90`

The primary file is `MOD_ODE_Toolbox.f90`. All ODE solvers and necessary helper subroutines/functions are contained here. A  
separate module, `MOD_ODE_Systems.f90`, contains examples of various ODE systems that may be useful. Using this module is not required  
since the ODE system may be defined in any convenient location that has access to the abstract interface shown in `MOD_ODE_Toolbox.f90`.  
A module file with data input/export subroutines, `MOD_IO_Toolbox.f90`, is also provided for convenience.  

## Methods
The code currently contains the following solvers:
- `ODE_Numerical_Solve_RK4` - Implements a fixed-step classical Runge-Kutta (RK4) method for solving a system of first-order ODEs.
- `ODE_Numerical_Solve_RK4_Adaptive` - Adaptive RK4 solver with error control for variable step-size integration.
- `ODE_Numerical_Solve_VSS` - Variable Step Size (VSS) ODE solver using the Dormand-Prince RK5 method with adaptive step control for handling moderately stiff ODEs.  

Detailed documentation on how to use the codes is provided in the PDF *FORTRAN_ODE_Toolbox.pdf*.


