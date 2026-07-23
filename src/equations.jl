"""
All equations used in the neural field model are defined in this file 
"""

using LinearAlgebra
using SpecialFunctions

function f(u)
    θ=0.9
    a=10
    V=0.1
    ϕ(a*(u-θ)/sqrt(1+a^2*V))
end

function ϕ(x)
    1/2*(1+erf(x/sqrt(2)))
end    

function φ(r)
    γ=5
    β=1/4
    γ / (cosh(β*norm(r)))
end

function w(x,y)
    W(norm(x-y))
end
    
function W_tilde(x)
    b = 0.4
    exp(-b*x)*(b*sin(x)+cos(x))
end

function W(x)
    ε=0.05
    Wx=W_tilde(x)
    abs(Wx) >= ε ? Wx : zero(Wx)
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

function one_layer(node_x)
    dim_u=length(node_x)
    z_layer=zeros(dim_u,dim_u) 
    for i in 1:dim_u
        for j in 1:dim_u
            z_layer[i,j] = 1* φ(node_x[i])*φ(node_x[j])
        end 
    end
    return z_layer           
end  

function z_initial(num_layer, node_x)
    z0 = one_layer(node_x)
    [copy(z0) for _ in 1:num_layer]
end

function z_initial_one_entry(nodex_i, nodex_j)
    """ This function generates initial condition
    for a given entry of z, base on the coresponding coordinates of u
    Arguments: 
        nodex_i: ith coordinate of node_x 
        nodex_j: jth coordinate of node_x 
    """
    z_ij = 1* φ(nodex_i)*φ(nodex_j)
    return z_ij
end     




