# **Documentation for the Fiber Electrodynamics project**

*by* **Ivan Toftul**

toftul.ivan@gmail.com


*2018/12/10*

# **Files structure**



## File types

There are two main types of files in this project:
1. `*.jl` -- the core files. Main (=core) functions are located in those files. Ideally, you do not need to change anything there, reading the input arguments should be enough.  
    All of them are located in `./lib/` folder.  
    Structure of functions approximatly is the following:
```Julia
function functionName(argument1, argument2)
    #=
        This function does ....
    
        Input
        ------
        argument1 - ...
        argument1 - ...
        
        Output
        ------
        output1 - ...
        output2 - ...   
    =#
    
    # function body
    
    return output1, output2
end
```
Each file is independent and contains functions of a quite narrow topic. 



## About `*.jl` files

* `consts.jl` - some physical consts in SI units.
* `fiberDispersion.jl` - functions to solve dispersion equation $\beta = \beta(\omega)$ for the step-index fiber.
* `fiberModeFields.jl` - everthing which is conneted with fiber modes (mode profiles, intensity, etc). Based on [Fam's paper](https://doi.org/10.1103/PhysRevA.96.023835). For more info look at the Manuscript.
* `fiberScatteredFields.jl` - Mie solutions for the scattered filed from the incident plane wave on a cylinder.
* `GFfiber.jl` - everything which is connected with Green's tensor of fiber. Functions to calculate it in different coordinates.
* `GFvacuum.jl` - everything which is connected with vacuum Green's tensor. Functions to calculate it in different coordinates.
* `MieTheory.jl` - standard Mie solutions and extention to the generalized Mie theory for scattering on a sphere. Calculation of Mie coefficients, generalized Mie coefficients, force and torque. For details see the Manuscript.
* `MieTheoryCyl.jl` - Mie solutions for the scattered filed from the incident plane wave on a cylinder. Cylindrical Mie coefficients. 
* `permittivity.jl` - function to interpolate data from [refractiveindex.info](https://refractiveindex.info/).
* folder `refractiveIndexData` - data files for the `permittivity.jl`.
* `polarizability.jl` - all approximations for the dipole polarizability of a spherical particle and atom. For details see the Manuscript.
* `supportFunction.jl` - all basic support functions which julia was lacking when I was developing code. In particular, there are rewritten functions from Fortran for the spherical functions which work extreamly fast (much faster then calling python library using PyCall). 
