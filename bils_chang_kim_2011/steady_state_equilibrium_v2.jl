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
#using FastGaussQuadrature

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Initialize parameters & asset grids

"""
    ModelParameters
"""
@with_kw mutable struct ModelParameters
    A::Float64              = 1.0                           # output productivity factor
    r::Float64              = 0.06                          # annual interest rate 
    r_mon::Float64          = (1.0 + r)^(1.0/12.0) - 1.0    # monthly interest rate 
    α::Float64              = 0.5                           # worker's share of surplus 
    β::Float64              = 0.995133                      # discount factor 
    β_mon::Float64          = 1.0/(1.0 + r_mon)             # monthly depreciation rate 
    γ::Float64              = 1.0                           # curvature of utility function 
    λ::Float64              = 0.02                          # exogenous separation rate 
    UI::Float64             = 0.4                           # UI benefits 
    B_leisure::Float64      = 0.0                           # utility from leisure 
    κ::Float64              = 0.522                         # vacancy posting cost with log utility 
    θ::Float64              = 1.0                           # labor market tightness, v/u, normalized to 1
    match_pow::Float64      = 0.5                           # power term in the matching function 
    match_scale::Float64    = 0.3133                        # scale parameter in the matching function 
    match_prob::Float64     = match_scale*θ^(match_pow)     # matching probability 
    amin::Float64           = -6.0                          # minimum asset 
    amax::Float64           = 120.0                         # maximum asset 
    agrid::Vector{Float64}  = LinRange(amin, amax, 100)     # asset grid 
    xmin::Float64           = 0.5                           # minimum match productivity 
    xmax::Float64           = 1.5                           # maximum match productivity 
    xgrid::Vector{Float64}  = LinRange(xmin, xmax, 5)       # match productivity grid 
end 

para = ModelParameters() 
