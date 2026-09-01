# Enabling Large-Scale Simulation of Axonal Delays in Neurobiological Networks

<p align="center">
  <img src="https://github.com/user-attachments/assets/547b1c13-8d30-44c8-ae0c-032430d00b7b" width="49%" />
  <img src="https://github.com/user-attachments/assets/6257ee86-c3b5-40e2-9f76-c4e976b3dd12" width="49%" />
</p>

# Reproducing the results

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> neural_field

To (locally) reproduce this project, do the following:

0. Download this code base. 
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.
