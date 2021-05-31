using Parameters 
using LinearAlgebra 
using Statistics
using Distributions
using Optim
using SparseArrays 
using QuantEcon
using SpecialFunctions

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize parameters & asset grids

"""
    ModelParams
"""
@with_kw mutable struct ModelParams
    A::Float64              = 1.0                           # output productivity factor
    r::Float64              = 0.06                          # annual interest rate 
    r_mon::Float64          = (1.0 + r)^(1.0/12.0) - 1.0    # monthly interest rate 
    worker_surplus::Float64 = 0.5                           # worker's share of surplus 
    β::Float64              = 0.995133                      # discount factor 
    β_mon::Float64          = 1.0/(1.0 + r_mon)             # monthly depreciation rate 
    γ::Float64              = 1.0                           # curvature of utility function 
    λ::Float64              = 0.02                          # exogenous separation rate 
    b::Float64              = 0.4                           # UI benefits 
    B::Float64              = 0.0                           # utility from leisure 
    κ::Float64              = 0.522                         # vacancy posting cost with log utility 
    θ::Float64              = 1.0                           # labor market tightness, v/u, normalized to 1
    α::Float64              = 0.5                           # power term in the matching function 
    η::Float64              = 0.3133                        # scale parameter in the matching function
    p_θ::Float64            = η*θ^α                         # unemployed worker matching rate  
    q_θ::Float64            = η*θ^(α-1)                     # vacancy matching probability 
    amin::Float64           = -6.0                          # minimum asset 
    amax::Float64           = 120.0                         # maximum asset 
    agrid::Vector{Float64}  = LinRange(amin, amax, 100)     # asset grid 
    ρ_x::Float64            = 0.97                          # persistence of idiosyncratic productivity ln(x)
    σ_x::Float64            = 0.13                          # std. dev. of innovation to ln(x) 
    π_x                     = TauchenApprox(9, ρ_x, σ_x).distribution           # probability weights on ln(x) 
    xgrid                   = exp.(TauchenApprox(9, ρ_x, σ_x).state_values)     # match productivity grid 
end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
Household Bernoulli utility function
"""
function utility(para::ModelParams, c::Real, l::Real)
    
    γ = para.γ # Store risk-aversion parameter
    B = para.B # Store the leisure utility coefficient

    # Define the utility mapping 
    # condition on the input `c`
    if c > 0 
        if γ == 1
          return log(c) + B*l # if σ = 1, CRRA utility is log
        else
          return (c^(1-γ)-1)/(1-γ) + B*l # if σ != 1, regular specification
        end
    else
        return -Inf # infinite disutility if c∈(-∞,0]
    end
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
    TauchenApprox
"""
function TauchenApprox(N::Integer, ρ::T1, σ::T2, μ=zero(promote_type(T1, T2)), n_std::Integer=3) where {T1 <: Real, T2 <: Real}
    
    std_norm_cdf(x::T) where {T <: Real} = 0.5 * SpecialFunctions.erfc(-x/sqrt(2))
    std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* SpecialFunctions.erfc(-x./sqrt(2))

    # Get discretized space
    a_bar = n_std * sqrt(σ^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilities
    Π = zeros(promote_type(T1, T2), N, N)
    for row = 1:N
        # Do end points first
        Π[row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ)
        Π[row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ)

        # fill in the middle columns
        for col = 2:N-1
            Π[row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ))
        end
    end

    # state values 
    ȳ = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable
    
    # renormalize transition matrix 
    Π = Π./sum(Π, dims = 2)

    # Generate stationary distribution from Π
    N = size(Π)[1]
    π0 = ones(1,N)/N
    diff = 1.
    while diff > 1e-10
        π1 = π0*Π
        diff = norm(π1-π0,Inf)
        π0 = π1
    end
    π0 = copy(π0') 

    return (transition = Π, distribution = π0, state_values = ȳ)
end

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
  iterateBellman()

Iterates on the bellman equation using continuation value function.
"""
function iterateBellman(para::ModelParams, W_old, U_old, J_old, V_old)

    @unpack agrid, xgrid, π_x = para #unpack parameters

    u(c,l) = utility(para,c,l) #shorthand for utility

    N_x = length(xgrid)
    N_a = length(agrid)

    #solve for each state
    #preallocate space for V and a_policy

    for x_i  in 1:N_x
        for a_i_ in 1:N_a
            for a′_i in 1:N_a 
                
            end 

        end
    end

    return W_new, U_new, J_new, V_new 
end;