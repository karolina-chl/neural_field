"""
All equations used in the neural field model are defined in this file 
"""

using LinearAlgebra
using SpecialFunctions

function w(x,y)
    WAε(norm(x-y))
end

function f(u)
    θ=5.6
    α=10
    V=1
    ϕ(α*(u-θ)/sqrt(1+α^2*V))
end

function ϕ(x)
    1/2*(1+erf(x/sqrt(2)))
end    

function φ(r;)
    α=5
    β=1/α
    α / (cosh(β*norm(r)))^2
end
    
function A(x)
    b = 0.4
    exp(-b*x)*(b*sin(x)+cos(x))
end

function WAε(x)
    ε=0.05
    Ax=A(x)
    abs(Ax) >= ε ? Ax : zero(Ax)
end

function τ(x,y) 
    τ_0=0.1
    τ_1=0.01 
    diff=(x .- y).^2
    distance_xy = sqrt(sum(diff))
    τ_0 + τ_1*distance_xy
end 

function α(x,y, num_layer)
    #this function accepts node coordinates
    num_layer/τ(x,y)
end

# function one_layer(dim_u)
#     z_layer=zeros(dim_u,dim_u) 
#     for i in 1:dim_u
#         for j in 1:dim_u
#             z_layer[i,j] = 0.1* φ(i)*φ(j)
#         end 
#     end
#     return z_layer           
# end  

# function one_layer(dim_u)
#     xs = range(-(dim_u-1)/2, (dim_u-1)/2; length=dim_u)
#     v = [φ(x) for x in xs]
#     return 0.1 .* (v * v')
# end

# function z_initial(num_layer, dim_u)
#     return [one_layer(dim_u) for _ in 1:num_layer] 
# end 

# function one_layer(node_u)
#     n = length(node_u)
#     repeat(reshape(node_u, n, 1), 1, n)
# end

# function z_initial(num_layer, node_u)
#     z0 = one_layer(node_u)
#     return [copy(z0) for _ in 1:num_layer]
# end

function one_layer(dim_u; A=6.2, σ=100)
    c = (dim_u + 1) / 2
    Z = zeros(dim_u, dim_u)
    for i in 1:dim_u
        for j in 1:dim_u
            r2 = (i - c)^2 + (j - c)^2
            Z[j, i] = A * exp(-r2 / (2σ^2))
        end
    end
    return Z
end

function z_initial(num_layer, dim_u)
    z0 = one_layer(dim_u)
    [copy(z0) for _ in 1:num_layer]
end

