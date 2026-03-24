"""
All equations used in the neural field model are defined in this file 
"""

using LinearAlgebra
using SpecialFunctions
using MatrixEquations

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

# function τ(x,y) 
#     τ_0=0.1
#     τ_1=0.01 
#     diff=(x .- y).^2
#     distance_xy = sqrt(sum(diff))
#     τ_0 + τ_1*distance_xy
# end 

# function α(x,y, num_layer)
#     #this function accepts node coordinates
#     num_layer/τ(x,y)
# end

### Daniele functions

τ_s_fun(x, τs0, τs1) = τs0 * (1 + τs1*abs(x))

function compute_delay_matrix(x, τ_s_fun)

  n = length(x)
  M_τ_s = zeros(n, n)

  y_τ_s = τ_s_fun.(x)

  i_rows = 1:n
  i_shift = -(n ÷ 2):(n ÷ 2)

  for i in 1:n
    M_τ_s[i_rows[i], :] = circshift(y_τ_s, i_shift[i])
  end

  return M_τ_s
end

function compute_equilibrium_variance(τ_u, M_τ_s, σ_u, σ_s1, σ_s2)

  # Parameters
  n = size(M_τ_s, 1)
  ξ_u = 1/τ_u

  # Allocate matrices
  A = zeros(3, 3)
  C_σ = similar(A)
  X = similar(A)
  V_star = similar(M_τ_s)

  # For each spatial point
  for j in 1:n
    for i in 1:n 
      
      # Form matrices
      ξ_s = 1/M_τ_s[i,j]

      B = [-ξ_u     0     0;
            ξ_s  -ξ_s     0;
              0   ξ_s  -ξ_s]

      Λ = [σ_u^2       0      0;
               0  σ_s1^2      0;
               0       0 σ_s2^2]

      # Solve the Sylvester equation and store result
      X = sylvc(-B, -B', Λ)
      V_star[i,j] = X[1,1]

    end
  end
  
  return V_star
end

### Initialize z like Daniele 
function z_initial(x,num_layer)
    τs0 = 0.1
    τs1 = 0.01
    τ = 1
    σu = √(2 * 0.06) # 1.0 
    σs1 = √(2 * 0.06) # 1.0 
    σs2 = √(2 * 0.06) # 1.0 

    M_τ = compute_delay_matrix(x, x -> τ_s_fun(x, τs0, τs1))
    V_s = compute_equilibrium_variance(τ, M_τ, σu, σs1, σs2)
    z0 = V_s
    return [copy(z0) for _ in 1:num_layer]
end     