"""
All equations used in the neural field model are defined in this file 
"""

using LinearAlgebra

function w(x,y)
    WAε(norm(x-y))
end

function f(u)
    μ=5.5
    θ=5.6
    1/(1+exp(-μ*u+θ)) - 1/(1+exp(θ))
end

function φ(r;)
    α=20
    β=1/α
    α / (cosh(β*norm(r)))^2
end
    
function A(x)
    b=0.4
    exp(-b*x)*(b*sin(x)+cos(x))
end

function WAε(x)
    ε=1.0e-3    
    Ax = A(x)
    abs(Ax) >= ε ? Ax : zero(Ax)
end
