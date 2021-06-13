###########################################################
# Load (old) functions 
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

###########################################################
# Import necessary packages
using BasisMatrices
using Optim 
using Distributions 
using Parameters 
using LinearAlgebra 
using Distributions
using Optim
using QuantEcon
using SpecialFunctions
using SparseArrays
using Plots 
using LaTeXStrings 
using ThreadTools 

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Model parameters 
@with_kw mutable struct ModelParameters
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
    agrid::Vector{Float64}  = LinRange(amin, amax, 20)      # asset grid 
    ρ_x::Float64            = 0.97                          # persistence of idiosyncratic productivity ln(x)
    σ_x::Float64            = 0.13                          # std. dev. of innovation to ln(x) 
    π_x                     = TauchenApprox(9, ρ_x, σ_x).distribution           # probability weights on ln(x) 
    xgrid                   = exp.(TauchenApprox(9, ρ_x, σ_x).state_values)     # match productivity grid 
    spline_A                = SplineParams(agrid,0,1)  
    spline_X                = SplineParams(xgrid,0,1) 
    basis_A                 = Basis(spline_A)
    basis_AX                = Basis(spline_A,spline_X)
    nodes_A                 = nodes(basis_A)[1]
    nodes_AX                = nodes(basis_AX)[1] 
end;

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# utility function 
function utility(para::ModelParameters, c::Real)
    
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

function utility(para::ModelParameters, cvec)

  u(c) = utility(para,c)

  return u.(cvec)
end;

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Initialize interpolated worker value functions
function InitializeValues(para::ModelParameters)
    
    # Retrieve asset and productivity shock grids 
    @unpack basis_A, basis_AX, nodes_A, nodes_AX = para 
    
    # Create functions to map all (a,x) and (a) to 1
    # for both W(a,x) and U(a) 
    g1d = x -> x[1]^0
    gvals = [g1d(nodes_AX[i,:]) for i in 1:size(nodes_A,1)]
    f2d = x -> x[1]^0 * x[2]^0
    fvals = [f2d(nodes_AX[i,:]) for i in 1:size(nodes_AX,1)]

    # Interpolate W(a,x) and U(a) 
    W_init = Interpoland(basis_AX, fvals);
    U_init = Interpoland(basis_A, gvals)

    # Return bases for A and A×X
    # as well as interpolated initial W(a,x) and U(a) 
    return W_init, U_init
end;

###########################################################
# Initialize interpolated worker value functions
W_init, U_init = InitializeValues(ModelParameters());

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Initialize interpolated worker wage schedule
function InitializeWage(para::ModelParameters)
    
    # Retrieve asset and productivity shock grids 
    # as well as the market tightness paramter θ
    @unpack θ, basis_AX, nodes_AX = para
    
    # Map all (a,x) to θ*0.9 = 0.9
    f2d = x -> x[1]^0 * x[2]^0 * θ * 0.9
    fvals = [f2d(nodes_AX[i,:]) for i in 1:size(nodes_AX,1)]

    # Interpolate wage schedule over A×X
    wage_init = Interpoland(basis_AX, fvals);

    # Return basis for A×X 
    # and interpolated wage schedule 
    return wage_init
end;

###########################################################
# Initialize interpolated wage schedule
wage_init = InitializeWage(ModelParameters());

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Interpolated invariant productivity shock process density function
function ShockDensity(para::ModelParameters, xval)

    @unpack ρ_x, σ_x = para  

    distr = LogNormal(0, σ_x/sqrt(1-ρ_x^2))
    density = pdf(distr, xval)

    return density  
end 

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Single iteration of the Bellman map w interpolated objects 
function WorkerBellmanMap(para::ModelParameters, wage, W_old, U_old)

    # Import parameters 
    @unpack xgrid, π_x, β, b, B, p_θ, r, basis_A, basis_AX, nodes_A, nodes_AX = para 
    
    # Shorthand for utility
    u(c) = utility(para,c)

    # Create empty W and U objects 
    W_new           = similar(nodes_AX[:,1])
    U_new           = similar(nodes_A)
    emp_policy      = similar(W_new)
    unemp_policy    = similar(U_new)

    # Initialize index for U(a) 
    k=0

    # Fill out W_new and U_new matrices 
     Threads.@threads for i in 1:length(nodes_AX[:,1])
    #for i in 1:length(nodes_AX[:,1])
        
        # Create (a,x) vector 
        ax = vec(nodes_AX[i,:])

        # Create Bellman mapping under employment 
        c_emp(a′) = (1+r)*ax[1] + wage(ax) - a′ 
        obj_emp     = tmap(a′ -> (1)*(u(c_emp(a′)) + β * sum( π_x[i]*max.(W_old(vec([a′ xgrid[i]])), U_old(xgrid[i])) for i in 1:length(xgrid))), nodes_A)
        # Solve for maximum W(a,x) 
        # and maximizing a′ 
        W_new[i], optimal_a′ = findmax(obj_emp)
        emp_policy[i] = nodes_A[optimal_a′]

        # If x equals E[x], then 
        # start solving for max. U(a) 
        # and maximizing a′
        if ax[2] ≈ 1.0 
            # Update index for U(a) 
            k = k + 1
            # Create Bellman mapping under unemployment 
            c_unemp(a′) = (1+r)*ax[1] + b - a′
            obj_unemp   = tmap( a′ -> u(c_unemp(a′)) + B + β*(1-p_θ)*U_old(ax[1]) .+ β*(p_θ)*W_old(ax), nodes_A)
            # Solve for maximum W(a,x) 
            # and maximizing a′
            U_new[k], optimal_a′ = findmax(obj_unemp)
            unemp_policy[k] = nodes_A[optimal_a′]
        end
    end

    # Interpolate W_new & U_new 
    W_new = Interpoland(basis_AX, W_new)
    U_new = Interpoland(basis_A, U_new)

    return W_new, U_new, emp_policy, unemp_policy
end;

###########################################################
# Iterate Bellman map
para = ModelParameters()
para.agrid = LinRange(para.amin, para.amax, 18)
W_new, U_new, emp_policy, unemp_policy = WorkerBellmanMap(para, wage_init, W_init, U_init)

###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# Value function iteration 
function SolveWorkerBellman(para::ModelParameters, wage, W0, U0, ϵ=1e-5)

    @unpack nodes_A = para 

    W_old = W0
    U_old = U0 

    diff = 1.
    iter = 0
    
    while diff > ϵ
        iter = iter + 1
        W_new, U_new = WorkerBellmanMap(para, wage, W_old, U_old)
        diff = norm(U_new.(nodes_A)-U_old.(nodes_A))
        W_old = W_new 
        U_old = U_new 
        println("Iteration: $(iter); Norm: $(diff)")
    end
    
    W, U, emp_policy, unemp_policy = WorkerBellmanMap(para, wage, W_old, U_old)

    return W, U, emp_policy, unemp_policy
end;

###########################################################
# Iterate Bellman map
W_init, U_init = InitializeValues(ModelParameters());
wage_init = InitializeWage(ModelParameters());
para = ModelParameters()
para.agrid = LinRange(para.amin, para.amax, 100)
para.B = 0.0
W_new, U_new, emp_policy, unemp_policy = SolveWorkerBellman(para, 
                                                            wage_init, 
                                                            W_init, 
                                                            U_init);

plot(a -> U_new(a), LinRange(-6.0,120.0,5000))

U_new(0)
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
###########################################################
# 

using Profile 

@profile InitializeValues(ModelParameters())
Profile.print()
