
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

#using Main.SSEquilibrium

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
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
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out SolveHHBellman()  

para  = ModelParams() 
W_old = ones(length(para.agrid),length(para.xgrid))
U_old = ones(length(para.agrid))
#wage  = 10*randn(length(para.agrid),length(para.xgrid)).^2
wage = ones(length(para.agrid),length(para.xgrid))
for i in 1:length(para.agrid)
    for j in 1:length(para.xgrid)
        wage[i,j] = i*j 
    end 
end 

W, U, emp_policy, unemp_policy = SolveHHBellman(para, wage*0.9, W_old, U_old)

using Plots

fig = plot();
for i in 1:9
    plot!(fig, 1:length(W[:,i]), W[:,i]);
end 
plot(fig) 
fig2 = plot();
plot!(fig2, 1:length(U), U);
plot(fig2)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out UpdateWage() 

wage_old = wage 
wage_new, J = UpdateWage(para, W, U, emp_policy, wage_old)
wage_old = wage_new 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out SolveWage() 

# Load functions 
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

# Approximate value, policy, and wage functions
para = ModelParams(θ=1.0)
para.agrid = LinRange(para.amin, para.amax, 50)
W, U, J, emp_policy, unemp_policy, wage = SolveWage(para)

# Plot approximated functions 
using Plots

## Plot value of employment for each x
fig = plot();
for i in 1:9
    plot!(fig, para.agrid, W[:,i]);
end 
plot(fig) 

## Plot value of unemployment 
fig2 = plot();
plot!(fig2, para.agrid, U);
plot(fig2)

## Plot value of vacancy for each x 
fig3 = plot();
for i in 1:9
    plot!(fig3, para.agrid, J[:,i]);
end 
plot(fig3)

## Plot wage for each x
fig4 = plot();
for i in 1:9
    plot!(fig4, para.agrid, wage[:,i]);
end 
plot(fig4) 

# Save objects externally 
using DelimitedFiles
writedlm( "bils_chang_kim_2011/W.csv",  W, ',')
writedlm( "bils_chang_kim_2011/U.csv",  U, ',')
writedlm( "bils_chang_kim_2011/J.csv",  J, ',')
writedlm( "bils_chang_kim_2011/emp_policy.csv",  emp_policy, ',')
writedlm( "bils_chang_kim_2011/unemp_policy.csv",  unemp_policy, ',')
writedlm( "bils_chang_kim_2011/wage.csv",  wage, ',')

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Analyze exported objects from SolveWage()

#####################################
# Import saved objects 
using DataFrames
using CSV

## Value of employment
W = CSV.read("bils_chang_kim_2011/W.csv", DataFrame, header = false) 
W = Matrix(W) 

# Value of unemployment
U = CSV.read("bils_chang_kim_2011/U.csv", DataFrame, header = false)
U = Matrix(U) 

# Value of a matched job
J = CSV.read("bils_chang_kim_2011/J.csv", DataFrame, header = false)
J = Matrix(J) 

# Asset policy function under employment
emp_policy      = CSV.read("bils_chang_kim_2011/emp_policy.csv", DataFrame, header = false)
emp_policy      = Matrix(emp_policy) 

# Asset policy function under unemployment 
unemp_policy    = CSV.read("bils_chang_kim_2011/unemp_policy.csv", DataFrame, header = false)
unemp_policy    = Matrix(unemp_policy)

# Wage function 
wage            = CSV.read("bils_chang_kim_2011/wage.csv", DataFrame, header = false)
wage            = Matrix(wage) 

#####################################
# Plot approximated functions 
using Plots 
using LaTeXStrings 

## Save asset grid for x-axes 
agrid = ModelParams().agrid 
xgrid = ModelParams().xgrid

## Plot value of employment for each x
fig1 = plot();
for i in 1:9
    plot!(fig1, agrid, emp_policy[:,i]);
end 
plot(fig1) 

## Plot value of unemployment 
fig2 = plot();
plot!(fig2, agrid, unemp_policy);
plot(fig2)

## Plot value of vacancy for each x 
fig3 = plot(title = "Matched Firm Value");
for i in 1:9
    plot!(fig3, agrid, J[:,i], label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!(L"J");
plot(fig3) 

## Plot wages for each x 
fig4 = plot(title = "Wage Schedule");
for i in 1:9
    plot!(fig4, agrid, wage[:,i], label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!("Wage");
plot(fig4) 

## Plot W - U for each x 
fig5 = plot(title = "Employment vs. Unemployment");
for i in 1:9
    plot!(fig5, agrid, W[:,i] - U, label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!(L"W-U");
plot(fig5) 

# All plots together 
plot(fig1)
plot(fig2)
plot(fig3)
plot(fig4)
plot(fig5) 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out ReservationProductivity() 

# Load functions 
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

# Test 
x_star = ReservationProductivity(para, W, U)
plot(1:length(x_star),x_star)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out ConstructTransitionMatrices() 

# Load functions 
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

# Test 
para = ModelParams()
para.agrid = LinRange(para.amin, para.amax, 50)
H_ee, H_eu, H_ue, H_uu = InitializeTransitionMatrices(para, emp_policy, unemp_policy)
π_emp, π_unemp = FindStationaryMeasures(para,  emp_policy, unemp_policy, x_star)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out AdjustTheta() 
θ_new = AdjustTheta(para, J, π_unemp)

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test SolveEquilibrium() 
para = ModelParams() 
para, θ, W, U, J, emp_policy, unemp_policy, wage, x_star, π_emp, π_unemp = SolveEquilibrium(para::ModelParams) 

# Save objects externally 
using DelimitedFiles
writedlm( "bils_chang_kim_2011/data/theta.csv",  θ, ',')
writedlm( "bils_chang_kim_2011/data/W.csv",  W, ',')
writedlm( "bils_chang_kim_2011/data/U.csv",  U, ',')
writedlm( "bils_chang_kim_2011/data/J.csv",  J, ',')
writedlm( "bils_chang_kim_2011/data/emp_policy.csv",  emp_policy, ',')
writedlm( "bils_chang_kim_2011/data/unemp_policy.csv",  unemp_policy, ',')
writedlm( "bils_chang_kim_2011/data/wage.csv",  wage, ',')
writedlm( "bils_chang_kim_2011/data/x_star.csv",  x_star, ',')
writedlm( "bils_chang_kim_2011/data/pi_emp.csv",  π_emp, ',')
writedlm( "bils_chang_kim_2011/data/pi_unemp.csv",  π_unemp, ',')