
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

#using Main.SSEquilibrium

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
# Test out SolveHHBellman()  

para  = ModelParams() 
W_old = ones(length(para.agrid),length(para.xgrid))
U_old = ones(length(para.agrid))
wage  = randn(length(para.agrid),length(para.xgrid)).^2
W, U, emp_policy, unemp_policy, x_star = SolveHHBellman(para, wage, W_old, U_old)

using Plots

plot(1:length(U),U)

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
# Test out UpdateWage() 

wage_old = wage 
wage_new, J = UpdateWage(para, W, U, emp_policy, wage_old)
wage_old = wage_new 

# ######################################################################
# ######################################################################
# ######################################################################
# ######################################################################
# Test out SolveWage() 

W, U, J, emp_policy, unemp_policy, x_star, wage = SolveWage(ModelParams())