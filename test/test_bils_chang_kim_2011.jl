
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

W, U, emp_policy, unemp_policy, x_star = SolveHHBellman(para, wage, W_old, U_old)

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
W, U, J, emp_policy, unemp_policy, x_star, wage = SolveWage(ModelParams())

# Plot approximated functions 
using Plots

## Plot value of employment for each x
fig = plot();
for i in 1:9
    plot!(fig, 1:length(W[:,i]), W[:,i]);
end 
plot(fig) 

## Plot value of unemployment 
fig2 = plot();
plot!(fig2, 1:length(U), U);
plot(fig2)

## Plot value of vacancy for each x 
fig3 = plot();
for i in 1:9
    plot!(fig3, 1:length(J[:,i]), J[:,i]);
end 
plot(fig3)

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
# 

# Analyze exported objects from SolveWage()

# Import saved objects 
using DataFrames
using CSV
W = CSV.read("bils_chang_kim_2011/W.csv", DataFrame, header = false) 
W = Matrix(W) 
U = CSV.read("bils_chang_kim_2011/U.csv", DataFrame, header = false)
U = Matrix(U) 
J = CSV.read("bils_chang_kim_2011/J.csv", DataFrame, header = false)
J = Matrix(J) 
emp_policy      = CSV.read("bils_chang_kim_2011/emp_policy.csv", DataFrame, header = false)
emp_policy      = Matrix(emp_policy) 
unemp_policy    = CSV.read("bils_chang_kim_2011/unemp_policy.csv", DataFrame, header = false)
unemp_policy    = Matrix(unemp_policy)
wage            = CSV.read("bils_chang_kim_2011/wage.csv", DataFrame, header = false)
wage            = Matrix(wage) 

#

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out FindStationaryMeasures() 