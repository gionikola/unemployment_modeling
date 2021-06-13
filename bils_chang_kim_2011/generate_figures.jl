
####################################
# Load libraries
using DelimitedFiles
using Plots
using DataFrames
using CSV
using LaTeXStrings 

####################################
# Load functions 
include("bils_chang_kim_2011/steady_state_equilibrium.jl")

####################################
# Approximate value, policy, and wage functions, 
# reservation productivity function, and stationary distributions
# given θ = 1
para            = ModelParams()
W, U, J, emp_policy, unemp_policy, wage = SolveWage(para)
x_star          = ReservationProductivity(para, W, U)
π_emp, π_unemp  = FindStationaryMeasures(para, emp_policy, unemp_policy, x_star)

####################################
# Save all objects externally 
writedlm( "bils_chang_kim_2011/data/W.csv",  W, ',')
writedlm( "bils_chang_kim_2011/data/U.csv",  U, ',')
writedlm( "bils_chang_kim_2011/data/J.csv",  J, ',')
writedlm( "bils_chang_kim_2011/data/emp_policy.csv",  emp_policy, ',')
writedlm( "bils_chang_kim_2011/data/unemp_policy.csv",  unemp_policy, ',')
writedlm( "bils_chang_kim_2011/data/wage.csv",  wage, ',')
writedlm( "bils_chang_kim_2011/data/x_star.csv",  x_star, ',')
writedlm( "bils_chang_kim_2011/data/pi_emp.csv",  π_emp, ',')
writedlm( "bils_chang_kim_2011/data/pi_unemp.csv",  π_unemp, ',') 

####################################
# Import objects (if necessary) 
## Value of employment
W = CSV.read("bils_chang_kim_2011/data/W.csv", DataFrame, header = false) 
W = Matrix(W) 
## Value of unemployment
U = CSV.read("bils_chang_kim_2011/data/U.csv", DataFrame, header = false)
U = Matrix(U) 
## Value of a matched job
J = CSV.read("bils_chang_kim_2011/data/J.csv", DataFrame, header = false)
J = Matrix(J) 
## Asset policy function under employment
emp_policy      = CSV.read("bils_chang_kim_2011/data/emp_policy.csv", DataFrame, header = false)
emp_policy      = Matrix(emp_policy) 
## Asset policy function under unemployment 
unemp_policy    = CSV.read("bils_chang_kim_2011/data/unemp_policy.csv", DataFrame, header = false)
unemp_policy    = Matrix(unemp_policy)
## Wage function 
wage            = CSV.read("bils_chang_kim_2011/data/wage.csv", DataFrame, header = false)
wage            = Matrix(wage) 
## Employed distribution 
pi_emp          = CSV.read("bils_chang_kim_2011/data/pi_emp.csv", DataFrame, header = false)
pi_emp          = Matrix(pi_emp) 
## Unemployed distribution 
pi_unemp        = CSV.read("bils_chang_kim_2011/data/pi_unemp.csv", DataFrame, header = false)
pi_unemp        = Matrix(pi_unemp) 

####################################
# Plot all objects and save externally 
## Save asset grid for x-axes 
agrid = para.agrid 
xgrid = para.xgrid
N_a = length(agrid)
N_x = length(xgrid) 
## Plot value of employment for each x
fig1 = plot();
for i in 1:9
    plot!(fig1, agrid, emp_policy[:,i]);
end 
## Plot value of unemployment for each x 
fig2 = plot();
plot!(fig2, agrid, unemp_policy);
## Plot value of firm-worker match for each x 
fig3 = plot(title = "Matched Firm Value");
for i in 1:9
    plot!(fig3, agrid, J[:,i], label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!(L"J");
## Plot wages for each x 
fig4 = plot(title = "Wage Schedule");
for i in 1:9
    plot!(fig4, agrid, wage[:,i], label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!("Wage");
## Plot W - U for each x 
fig5 = plot(title = "Employment vs. Unemployment");
for i in 1:9
    plot!(fig5, agrid, W[:,i] - U, label="x=$(round(xgrid[i],digits=2))");
end 
xlabel!("Level of Asset Holdings");
ylabel!(L"W-U");
## Plot asset distributions for employed and unemployed  
fig6 = plot(title = "Asset Distributions for Employed and Unemployed");
emp_distr   = zeros(length(agrid))
unemp_distr = zeros(length(agrid))
for i in 1:N_a
    emp_distr[i]    = sum([π_emp[i+N_a*(j-1)] for j in 1:N_x])
    unemp_distr[i]  = sum([π_unemp[i+N_a*(j-1)] for j in 1:N_x])
end 
plot!(fig6, agrid, emp_distr, label = "Employed");
plot!(fig6, agrid, unemp_distr, label = "Unemployed");
xlabel!("Level of Asset Holdings");
## Plot reservation match productivity 
fig7 = plot(title = "Worker Assets and Reservation Match Productivity");
plot!(fig7, agrid, x_star, label = L"x^\ast(a)");
xlabel!("Level of Asset Holdings");
ylabel!("Productivity");
## All plots together 
plot(fig1)
plot(fig2)
plot(fig3)
plot(fig4)
plot(fig5) 
plot(fig6)
plot(fig7) 
## Export plots 
savefig(fig1, "bils_chang_kim_2011/figures/fig1.png")
savefig(fig2, "bils_chang_kim_2011/figures/fig2.png")
savefig(fig3, "bils_chang_kim_2011/figures/fig3.png")
savefig(fig4, "bils_chang_kim_2011/figures/fig4.png")
savefig(fig5, "bils_chang_kim_2011/figures/fig5.png")
savefig(fig6, "bils_chang_kim_2011/figures/fig6.png")
savefig(fig7, "bils_chang_kim_2011/figures/fig7.png")