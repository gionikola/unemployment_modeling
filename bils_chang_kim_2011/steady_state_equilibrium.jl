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

# ###################################
# ###################################
# Initialize parameters & asset grids

"""
    InitializeParameters()  
"""
mutable struct InitializeParameters() 
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
end 

# ###################################
# ###################################
# Initialize value functions and Nash-bargaining wages

"""
    InitializeWageValueDecision()
"""
function InitializeWageValueDecision()

end 

# ###################################
# ###################################
# Find equilibriun Nash-bargaining wages 

"""
    EquilibriumWage()
""" 
function EquilibriumWage()

end 

# ###################################
# ###################################
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

# ###################################
# ###################################
# Compute new θ using free entry condition

# Insert code 

# ###################################
# ###################################
# Compute aggregate variables 

"""
    AggregateVariables()
"""
function AggregateVariables()

end 

# ###################################
# ###################################
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



# ###################################
# ###################################