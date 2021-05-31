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
function utility(para::ModelParams, c::Real)
    
    γ = para.γ # Store risk-aversion parameter
    B = para.B # Store the leisure utility coefficient

    # Define the utility mapping 
    # condition on the input `c`
    if c > 0 
        if γ == 1
          return log(c) # if σ = 1, CRRA utility is log
        else
          return (c^(1-γ)-1)/(1-γ) # if σ != 1, regular specification
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
Household Bernoulli utility function applied to a vector
"""
function utility(para::ModelParams, cvec)

  u(c) = utility(para,c)

  return u.(cvec)
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
    TauchenApprox()
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
  HHBellmanMap()

Iterates on the bellman equation using continuation value function.
"""
function HHBellmanMap(para::ModelParams, wage, W_old, U_old)

    @unpack agrid, xgrid, π_x, β, b, B, p_θ, r   = para 
    
    u(c) = utility(para,c) #shorthand for utility

    W_new           = similar(W_old)
    U_new           = similar(U_old)
    emp_policy      = similar(W_new, Int)
    unemp_policy    = similar(U_new, Int) 

    N_x = length(xgrid)
    N_a = length(agrid)

    EU      = U_old 
    EWx̄     = W_old[:,N_x]

    for x_i  in 1:N_x
        for a_i in 1:N_a
            
            cvec_emp = (1+r)*agrid[a_i] + wage[a_i,x_i] .- agrid
            EmaxWU =  max.(W_old*π_x, U_old)
            obj_emp     = u(cvec_emp) .+ B .+ β * EmaxWU
            obj_emp     = vec(obj_emp)
            W_new[a_i,x_i], emp_policy[a_i,x_i] = findmax(obj_emp)

            if x_i == N_x 
                cvec_unemp = (1+r)*agrid[a_i] + b .- agrid
                obj_unemp   = u(cvec_unemp) .+ β*(1-p_θ)*EU .+ β*(p_θ)*EWx̄
                obj_unemp   = vec(obj_unemp) 
                U_new[a_i], unemp_policy[a_i] = findmax(obj_unemp)
            end     
        end
    end

    return W_new, U_new, emp_policy, unemp_policy
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out HHBellmanMap()  

para  = ModelParams() 
W_old = ones(length(para.agrid),length(para.xgrid))
U_old = ones(length(para.agrid))
wage  = ones(length(para.agrid),length(para.xgrid))
W_new, U_new, emp_policy, unemp_policy = HHBellmanMap(para, wage, W_old, U_old)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
    SolveHHBellman()
"""
function SolveHHBellman(para::ModelParams, wage, W0, U0, ϵ=1e-6)

    @unpack xgrid = para 

    W_old = W0
    U_old = U0 
    emp_policy      = similar(W_old, Int)
    unemp_policy    = similar(U_old, Int) 

    diff = 1.

    while diff > ϵ
        W_new, U_new = HHBellmanMap(para, wage, W_old, U_old)
        diff = norm((U_new-U_old)[:],Inf)
        W_old = W_new 
        U_old = U_new 
    end
    
    W, U, emp_policy, unemp_policy = HHBellmanMap(para, wage, W_old, U_old)

    x_star = zeros(length(U))

    for i in 1:length(U)
        acc = zeros(length(W[1,:]))
        for j in 1:length(W[1,:])
            if U[i] > W[i,j]
                acc[j] = 1
            end 
        end 
        x_star[i] = floor(Int, sum(acc))
    end 
    x_star = floor.(Int, x_star)
    xstar = zeros(length(x_star))
    for i in 1:3
        if x_star[i] > 0
            xstar[i] = xgrid[x_star[i]]
        end 
    end 
    x_star = copy(xstar) 

    return W, U, emp_policy, unemp_policy, x_star 
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out SolveHHBellman()  

para  = ModelParams() 
W_old = ones(length(para.agrid),length(para.xgrid))
U_old = ones(length(para.agrid))
wage  = randn(length(para.agrid),length(para.xgrid)).^2
W, U, emp_policy, unemp_policy, x_star = SolveHHBellman(para, wage, W_old, U_old)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
    UpdateWage()
"""
function UpdateWage(para::ModelParams, W, U, emp_policy, wage_old)

    @unpack agrid, xgrid, β, λ, α, r, π_x = para 
    
    N_x = length(xgrid)
    N_a = length(agrid)
    
    J = zeros(N_a, N_x)
    wage_new = similar(wage_old)

    for a_i in N_a
        for x_i in N_x
            c_e = (1+r)*agrid[a_i] + wage_old[a_i,x_i] - emp_policy[a_i,x_i]
            J[a_i,x_i] = ((1-α)/α)*(W[a_i,x_i] - U[a_i])*c_e 
        end 
    end 

    for a_i in N_a
        for x_i in N_x
            wage_new[a_i,x_i] = 1*xgrid[x_i] - J[a_i,x_i] + β*(1-λ)*max(reduce(vcat,J[emp_policy[a_i,x_i],:]'*π_x),0.0)
        end 
    end 

    return wage_new, J 
end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out UpdateWage() 

wage_new, J = UpdateWage(para, W, U, emp_policy, wage)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
"""
    SolveWage()
""" 
function SolveWage(para::ModelParams, ϵ=1e-6)
    
    W_old       = ones(length(para.agrid),length(para.xgrid))
    U_old       = ones(length(para.agrid))
    wage_old    = ones(length(para.agrid),length(para.xgrid)).^2
    difff        = 10
    counter     = 0

    while difff > ϵ 
        W_new, U_new, emp_policy, unemp_policy, x_star = SolveHHBellman(para, wage_old, W_old, U_old)
        wage_new, J             = UpdateWage(para, W_new, U_new, emp_policy, wage_old)
        difff                   = norm(wage_new - wage_old)
        wage_old                = 0.5*wage_new + 0.5*wage_old 
        W_old                   = W_new
        U_old                   = U_new 
        counter = counter + 1
        println(counter) 
    end 
    
    W, U, emp_policy, unemp_policy, x_star      = SolveHHBellman(para, wage_old, W_old, U_old)
    wage, J                                     = UpdateWage(para, W, U, emp_policy, wage_old)

    return W, U, J, emp_policy, unemp_policy, x_star, wage 
end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out SolveWage() 

W, U, J, emp_policy, unemp_policy, x_star, wage = SolveWage(ModelParams())
