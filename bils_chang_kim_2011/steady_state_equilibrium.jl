using Parameters 
using LinearAlgebra 
using Statistics
using BasisMatrices
using Distributions
using Optim
using BasisMatrices 
using SparseArrays 
using Plots 
using LaTeXStrings 
using FastGaussQuadrature

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize parameters & asset grids

"""
    InitializeParameters  
"""
@with_kw mutable struct InitializeParameters
    prod::Float64 = 1.0 # output
    irate::Float64 = 0.06 # annual interest rate 
    α::Float64 = 0.5 # worker's share of surplus 
    β::Float64 = 0.995133 # discount factor 
    γ::Float64 = 1.0 # curvature of utility function 
    λ::Float64 = 0.02 # exogenous separation rate 
    UI::Float64 = 0.4 # UI benefits 
    B_leisure::Float64 = 0.0 # utility from leisure 
    κ::Float64 = 0.522 # vacancy posting cost with log utility 
    match_pow::Float64 = 0.5 # power term in the matching function 
    match_scale::Float64 = 0.3133 # scale parameter in the matching function 
    θ::Float64 = 1.0 # labor market tightness, v/u, normalized to 1
    amin::Float64 = -6.0 # minimum asset 
    amax::Float64 = 120.0 # maximum asset 
    agrid::Vector{Float64} = LinRange(amin, amax, 100) # asset grid 
    irate_monthly::Float64 = (1.0 + irate)^(1.0/12.0) - 1.0 # monthly interest rate 
    β_monthly::Float64 = 1.0/(1.0 + irate_monthly) # monthly depreciation rate 
    match_prob::Float64 = match_scale*θ^(match_pow) # matching probability 
end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize value functions and Nash-bargaining wages

"""
    InitializeWageValueDecision()
"""
function InitializeWageValueDecision(para::InitializeParameters)

    @unpack prod, β, λ, UI, B_leisure, match_prob, agrid, irate_monthly, β_monthly = para 

    # Initialize value functions and Nash-bargaining wages 
    # compute Nash-bargaining wage for a assuming that a is fixed 
    βpth = 1.0/(1-β*(1-λ)+β*match_prob) 

    N = length(agrid) 
    AE = zeros(N)
    AU = zeros(N) 
    Wage = zeros(N)
    Util_W = 0.0
    Util_U = 0.0
    W = zeros(N)
    U = zeros(N)
    J = zeros(N)
    
    # Fill out empty vectors 
    for i in 1:N 
        # Initialized asset decision rules 
        AE[i] = agrid[i] 
        AU[i] = agrid[i]
        # Solve Nash-bargaining wages given AE and AU
        Wage[i] = prod 
        # Initialize value functions given AE, AU, and Wage
        Util_W = log(irate_monthly*agrid[i]+Wage[i])
        Util_U = log(irate_monthly*agrid[i]+UI) + B_leisure 
        W[i] = (Util_W - β*λ*βpth*(Util_W-Util_U))/(1-β)
        U[i] = (Util_U - β*match_prob*βpth*(Util_W-Util_U))/(1-β)
        J[i] = (prod-Wage[i])/(1.0-β_monthly*(1.0-λ))
    end 

    return (AE, AU, Wage, Util_W, Util_U, W, U, J) 
end 

#AE, AU, Wage, Util_W, Util_U, W, U, J =  InitializeWageValueDecision(InitializeParameters())

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Find equilibriun Nash-bargaining wages 

"""
    SolveValueFunctions()
"""
function SolveValueFunctions(para::InitializeParameters)

    @unpack prod, β, λ, UI, B_leisure, match_prob, agrid, irate_monthly, β_monthly = para 

    err_value = 1.0
    iter_value = 0

    # Start iteration 
    while err_value > tol_value 

        # Expected values in the next period
        # according to the current employment status 
        EV_W = (1.0 - λ)*W + λ*U
        EV_U = match_prob*W + (1.0 - match_prob)*U 

        # Employed worker's utility maximization 
        

    end 

end 

"""
    EquilibriumWage()
""" 
function EquilibriumWage()

end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Compute time-invariant measures

"""
    InitializeDAgrid()
"""
function InitializeDAgrid()

end 

"""
    InitializeMeasures()
"""
function InitializeMeasures()

end 

"""
    InverseAssetDecisions()
"""
function InverseAssetDecision()

end 

"""
    InvariantMeasures()
"""
function InvariantMeasures()

end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Compute new θ using free entry condition

# Insert code 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Compute aggregate variables 

"""
    AggregateVariables()
"""
function AggregateVariables()

end 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Generate results 

"""
Steps:
(1) Initialize parmameters and asset grids
(2) Initialize value functions and Nash-bargaining wages 
(3) Find equilibrun Nash-bargaining wages
(4) Compute time-invariant measures 
(5) Compute new θ
(6) Compute aggregate variables
"""



# ######################################################################
# ######################################################################