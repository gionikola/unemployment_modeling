
#module SSEquilibrium

#export ModelParams, utility, TauchenApprox, HHBellmanMap, SolveHHBellman, UpdateWage, SolveWage 

using Parameters 
using LinearAlgebra 
using Distributions
using Optim
using QuantEcon
using SpecialFunctions
using SparseArrays
using ThreadTools 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize parameters & asset grids

@doc """
    ModelParams

A structure containing all model parameters, including asset and matching productivity shock grids. 
"""
@with_kw mutable struct ModelParams
    A::Float64              = 1.0                           # output productivity factor
    r::Float64              = 0.06                          # annual interest rate 
    r_mon::Float64          = (1.0 + r)^(1.0/12.0) - 1.0    # monthly interest rate 
    worker_surplus::Float64 = 0.5                           # worker's share of surplus 
    β::Float64              = 0.9948                        # discount factor 
    β_mon::Float64          = 1.0/(1.0 + r_mon)             # monthly depreciation rate 
    γ::Float64              = 1.0                           # curvature of utility function 
    λ::Float64              = 0.02                          # exogenous separation rate 
    b::Float64              = 0.4                           # UI benefits 
    B::Float64              = 0.15                          # utility from leisure 
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
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    utility(para::ModelParams, c::Real)

Household Bernoulli utility function.
"""
function utility(para::ModelParams, c::Real)
    
    γ = para.γ # Store risk-aversion parameter

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
@doc """
    utility(para::ModelParams, cvec)

Household Bernoulli utility function applied to a vector.
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
@doc """
    TauchenApprox(N, ρ, σ, μ, n_std) 

Tauchen's (1996) method for approximating AR(1) process with finite markov chain.

**Output:**
⋅ transition    -- Estimated transition matrix
⋅ distribution  -- Estimated stationary distribution (vector)
⋅ state_values  -- Etimated state values (vector) 
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
    Threads.@threads for row = 1:N
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
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    HHBellmanMap(para::ModelParams, wage, W_old, U_old)

Iterates on the bellman equation using continuation value function.

**Output:**
⋅ W_new         -- Updated value of employment (matrix) 
⋅ U_new         -- Updated value of unemployment (vector) 
⋅ emp_policy    -- Asset policy function in case of employment (matrix)
⋅ unemp_policy  -- Asset policy function in case of unemployment (vector) 
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

    Threads.@threads for x_i  in 1:N_x
        
        if x_i == 5 
            for a_i in 1:N_a
                cvec_emp = (1+r)*agrid[a_i] + wage[a_i,x_i] .- agrid
                obj_emp     = u(cvec_emp) .+ β * sum( π_x[i]*max.(W_old[:,i], U_old) for i in 1:N_x)
                obj_emp     = vec(obj_emp)
                W_new[a_i,x_i], emp_policy[a_i,x_i] = findmax(obj_emp)  
            
                cvec_unemp = (1+r)*agrid[a_i] + b .- agrid
                obj_unemp   = u(cvec_unemp) .+ B .+ β*(1-p_θ)*U_old .+ β*(p_θ)*W_old[:,5]
                obj_unemp   = vec(obj_unemp) 
                U_new[a_i], unemp_policy[a_i] = findmax(obj_unemp)
            end 
        else
            for a_i in 1:N_a
                cvec_emp = (1+r)*agrid[a_i] + wage[a_i,x_i] .- agrid
                obj_emp     = u(cvec_emp) .+ β * sum( π_x[i]*max.(W_old[:,i], U_old) for i in 1:N_x)
                obj_emp     = vec(obj_emp)
                W_new[a_i,x_i], emp_policy[a_i,x_i] = findmax(obj_emp)  
            end 
        end 
        """
        for a_i in 1:N_a
            
            cvec_emp = (1+r)*agrid[a_i] + wage[a_i,x_i] .- agrid
            obj_emp     = u(cvec_emp) .+ β * sum( π_x[i]*max.(W_old[:,i], U_old) for i in 1:N_x)
            obj_emp     = vec(obj_emp)
            W_new[a_i,x_i], emp_policy[a_i,x_i] = findmax(obj_emp)

            if x_i == 5 
                cvec_unemp = (1+r)*agrid[a_i] + b .- agrid
                obj_unemp   = u(cvec_unemp) .+ B .+ β*(1-p_θ)*U_old .+ β*(p_θ)*W_old[:,5]
                obj_unemp   = vec(obj_unemp) 
                U_new[a_i], unemp_policy[a_i] = findmax(obj_unemp)
            end     
        end
        """
    end

    return W_new, U_new, emp_policy, unemp_policy
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    SolveHHBellman(para::ModelParams, wage, W0, U0, ϵ=1e-6)

Iterates on the Bellman map until convergence. 

**Output:**
⋅ W             -- Approximated value of employment (matrix)
⋅ U             -- Approximated value of unemployment (vector) 
⋅ emp_policy    -- Approximated asset policy function in case of employment (matrix)
⋅ unemp_policy  -- Approximated asset policy function in case of unemployment (vector)
"""
function SolveHHBellman(para::ModelParams, wage, W0, U0, ϵ=1e-5)

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

    return W, U, emp_policy, unemp_policy
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    UpdateWage(para::ModelParams, W, U, emp_policy, wage_old)

Applies step #3(c) of the Bils, Chang, and Kim (2011) steady state equilibrium algorithm to update the wage, w¹(a,x;θ⁰), object.

**Output:**
⋅ wage_new  -- New iteration of the wage mapping (matrix) 
⋅ J         -- New iteration of the value of a matched job (matrix)
"""
function UpdateWage(para::ModelParams, W, U, emp_policy, wage_old)

    @unpack agrid, xgrid, β, λ, α, r, π_x = para 
    
    N_x = length(xgrid)
    N_a = length(agrid)
    
    J = zeros(N_a, N_x)
    wage_new = similar(wage_old)

    Threads.@threads for x_i in 1:N_x 
        for a_i in 1:N_a
            c_e = (1+r)*agrid[a_i] + wage_old[a_i,x_i] - emp_policy[a_i,x_i]
            J[a_i,x_i] = ((1-α)/α)*(W[a_i,x_i] - U[a_i])*c_e 
        end 
    end 

    Threads.@threads for x_i in 1:N_x 
        for a_i in 1:N_a
            wage_new[a_i,x_i] = 1*xgrid[x_i] - J[a_i,x_i] + β*(1-λ)*sum( π_x[i]*max.( J[emp_policy[a_i,x_i],i] , 0.0 ) for i in 1:N_x)
        end 
    end 

    return wage_new, J 
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize value functions and wage 
@doc """
    InitializeFunctions(para::ModelParams)

Initializes W(a,x), U(a), and wage(a,x).
"""
function InitializeFunctions(para::ModelParams)
    @unpack agrid, xgrid, β, r, b, B = para 
    N_a = length(agrid)
    N_x = length(xgrid) 
    W       = ones(N_a,N_x)
    U       = ones(N_a)
    wage    = ones(N_a,N_x)

    Threads.@threads for x_i  in 1:N_x
        
        if x_i == 5 
            for a_i in 1:N_a
                wage[a_i,x_i] = xgrid[x_i]*2
                W[a_i,x_i] = log((1+r)*agrid[a_i] + wage[a_i,x_i] - agrid[a_i])/(1-β) 
                U[a_i]     = log((1+r)*agrid[a_i] + b - agrid[a_i] + B)/(1-β) 
            end 
        else
            for a_i in 1:N_a
                wage[a_i,x_i] = xgrid[x_i]*2
                W[a_i,x_i] = log((1+r)*agrid[a_i] + wage[a_i,x_i] - agrid[a_i])/(1-β) 
            end 
        end 
    end 

    return W, U, wage 
end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Solve wage 
@doc """
    SolveWage(para::ModelParams, ϵ=1e-6)

Applies steps #1-3 of the Bils, Chang, and Kim (2011) steady state equilibrium algorithm to obtain value functions, policy functions, and the wage mapping.

**Output:**
⋅ W             -- Approximated value of employment (matrix)
⋅ U             -- Approximated value of unemployment (vector) 
⋅ emp_policy    -- Approximated asset policy function in case of employment (matrix)
⋅ unemp_policy  -- Approximated asset policy function in case of unemployment (vector) 
⋅ wage          -- Approximated wage mapping (matrix) 
""" 
function SolveWage(para::ModelParams, ϵ=1e-5)

    W_old, U_old, wage_old = InitializeFunctions(para) 

    """
    W_old       = ones(length(para.agrid),length(para.xgrid))
    U_old       = ones(length(para.agrid))
    c           = 0.9
    wage_old    = ones(length(para.agrid),length(para.xgrid))*(c*θ)
    """

    ζ_x         = 0.01
    difff       = 10
    counter     = 0

    while difff > ϵ 
        W_new, U_new, emp_policy = SolveHHBellman(para, wage_old, W_old, U_old)
        wage_new, J              = UpdateWage(para, W_new, U_new, emp_policy, wage_old)
        difff                    = norm(wage_new - wage_old)
        wage_old                 = ζ_x*wage_new + (1-ζ_x)*wage_old 
        W_old                    = W_new
        U_old                    = U_new 
        counter = counter + 1
        println("Iteration: $(counter). Norm: $(difff).") 
    end 
    
    W, U, emp_policy, unemp_policy              = SolveHHBellman(para, wage_old, W_old, U_old)
    wage, J                                     = UpdateWage(para, W, U, emp_policy, wage_old)

    return W, U, J, emp_policy, unemp_policy, wage 
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Reservation productivity, x(a;w) 
@doc """
    ReservationProductivity(para::ModelParams, W, U)
"""
function ReservationProductivity(para::ModelParams, W, U)

        @unpack xgrid, agrid = para 
        N_a = length(agrid)
        N_x = length(xgrid)
        x_comps = zeros(N_a,N_x) 
        # If for any given (a,x), U(a) > W(a,x)
        # then store x_comps(a,x) = 1.0 (unemployment)
        # other store x_comps(a,x) = 0.0 (employment) 
        for x_i in 1:N_x 
            for a_i in 1:N_a
                x_comps[a_i,x_i] = -1*W[a_i,x_i] > -1*U[a_i] 
            end 
        end 
        # x_grid = position in xgrid for which W surpasses U 
        # indicating entry of xgrid that allows the switch 
        # from unemployment to employment for a given `a`
        x_opt   = zeros(N_a)
        for a_i in 1:N_a
            x_index     = findfirst(x_comps[a_i,:] .== 0.0) 
            if x_index > 1
                x_opt[a_i]  = (xgrid[x_index] - xgrid[x_index-1])/2
            else 
                x_opt[a_i]  = 0.0 
            end 
        end 
        
    return x_opt
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    ConstructTransitionMatrices(para::ModelParams, emp_policy, unemp_policy, x_star)
""" 
function ConstructTransitionMatrices(para::ModelParams, emp_policy, unemp_policy, x_star)

    @unpack agrid, xgrid, π_x, p_θ = para
    N_x = length(xgrid)
    N_a = length(agrid)
    H_emp = zeros(N_a*N_x,N_a*N_x) 
    H_unemp = zeros(N_a*N_x,N_a*N_x) 
    ### Fill out the employment matrix
    Threads.@threads for a_i in 1:N_a
        for x_i in 1:N_x
            j = a_i+N_a*(x_i-1) # index for (a,x)
            a_emp = emp_policy[a_i,x_i] # index for a′ given employment
            a_unemp = unemp_policy[a_i] # index for a′ given unemployment 
            x′_star = x_star[a_emp] # value of x′_star under employment  
            for x_i′ in 1:N_x
                j′ = a_emp+N_a*(x_i′-1) # index for (a′,x′) given employment
                if xgrid[x_i′] > x′_star # if x′ is above reservation x
                    H_emp[j,j′] = H_emp[j,j′] + π_x[x_i′] # prob. of (a,x) → (a′,x′) given employment
                end 
                j′ = a_unemp+N_a*(x_i′-1) # index for (a′,x′) given unemployment 
                if x_i′ == 5 # if x′ = x̄ (x̄ is the middle (5th) entry of xgrid)
                    H_emp[j,j′] = H_emp[j,j′] + p_θ*π_x[x_i′] # prob. of (a,x) → (a′,x′) given employment
                end 
            end
        end
    end
    ### Fill out the unemployment matrix 
    Threads.@threads for a_i in 1:N_a
        for x_i in 1:N_x
            j = a_i+N_a*(x_i-1) # index for (a,x)
            a_emp = emp_policy[a_i,x_i] # index for a′ given employment
            a_unemp = unemp_policy[a_i] # index for a′ given unemployment 
            x′_star = x_star[a_emp] # value of x′_star under employment  
            for x_i′ in 1:N_x
                j′ = a_emp+N_a*(x_i′-1) # index for (a′,x′) given employment
                if xgrid[x_i′] < x′_star # if x′ is above reservation x
                    H_unemp[j,j′] = H_unemp[j,j′] + π_x[x_i′] # prob. of (a,x) → (a′,x′) given employment
                end 
                j′ = a_unemp+N_a*(x_i′-1) # index for (a′,x′) given unemployment 
                H_unemp[j,j′] = H_unemp[j,j′] + (1-p_θ)*π_x[x_i′] # prob. of (a,x) → (a′,x′) given employment
            end
        end
    end
    return H_emp, H_unemp 
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 
@doc """
    FindStationaryMeasures(para::ModelParams, emp_policy, unemp_policy)
""" 
function FindStationaryMeasures(para::ModelParams,  emp_policy, unemp_policy, x_star)

    H_emp, H_unemp = ConstructTransitionMatrices(para, emp_policy, unemp_policy, x_star)
    N_emp = size(H_emp)[1]
    π0_emp = ones(1,N_emp)/N_emp
    N_unemp = size(H_unemp)[1]
    π0_unemp = ones(1,N_unemp)/N_unemp
    diff = 1.
    
    while diff > 1e-10
        π1_emp = π0_emp*H_emp
        π1_unemp = π0_unemp*H_unemp
        diff = norm(π1_emp-π0_emp,Inf)
        diff = norm(π1_unemp-π0_unemp,Inf)
        π0_emp = π1_emp
        π0_unemp = π1_unemp
    end
    π_emp   = π0_emp 
    π_unemp = π0_unemp 

    return π_emp, π_unemp
end;

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# 

#end 