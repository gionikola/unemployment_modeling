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
    Parameters  
"""
@with_kw mutable struct Parameters
    A::Float64              = 1.0 # output productivity factor
    r::Float64              = 0.06 # annual interest rate 
    α::Float64              = 0.5 # worker's share of surplus 
    β::Float64              = 0.995133 # discount factor 
    γ::Float64              = 1.0 # curvature of utility function 
    λ::Float64              = 0.02 # exogenous separation rate 
    UI::Float64             = 0.4 # UI benefits 
    B_leisure::Float64      = 0.0 # utility from leisure 
    κ::Float64              = 0.522 # vacancy posting cost with log utility 
    match_pow::Float64      = 0.5 # power term in the matching function 
    match_scale::Float64    = 0.3133 # scale parameter in the matching function 
    θ::Float64              = 1.0 # labor market tightness, v/u, normalized to 1
    amin::Float64           = -6.0 # minimum asset 
    amax::Float64           = 120.0 # maximum asset 
    agrid::Vector{Float64}  = LinRange(amin, amax, 100) # asset grid 
    r_mon::Float64          = (1.0 + irate)^(1.0/12.0) - 1.0 # monthly interest rate 
    β_mon::Float64          = 1.0/(1.0 + irate_monthly) # monthly depreciation rate 
    match_prob::Float64     = match_scale*θ^(match_pow) # matching probability 
end 