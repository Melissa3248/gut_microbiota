using DifferentialEquations
using DelimitedFiles
using CSV

function solveODE(f, u0, vars, tspan, t, M, N)
    y = zeros(N,vars) # (simulation, variable)
    for sim in 1:N
        p = M[sim,:]
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob, Rosenbrock23())

        for variable in 1:vars
            y[sim,variable] = sol.u[t][variable]
        end

    end
    return y # matrix of model outputs
end

phi(u, gamma) = u./(gamma+u)

function f(du,u,p0,t)
    B, E, M, a, h, p = u
    betaa, betab, betaE1, betaE2, betah1, betah2, betaM1, betaM2, betap, gammaa, gammaB, gammah, gammap, muaE, muaM, muhM, mupB, mupE, q = p0

    du[1]= ßb.*phi(u[6],gammap).*u[1] - q.*u[1]
    du[2] = (ßE1.*phi(u[4], gammaa) + ßE2.*(1-phi(u[1], gammaB)).*phi(u[6], gammap)).*u[2] - q.*u[2]
    du[3] = (ßM1.*phi(u[4], gammaa) + ßM2.*phi(u[5], gammah)).*u[3] - q.*u[3]

    du[4] = ßa.*phi(u[6], gammap).*u[1] - q.*u[4] - (muaE.*u[2] + muaM.*u[3]).*u[4]
    du[5] = ßh1.*phi(u[4], gammaa).*u[2] + ßh2.*phi(u[6], gammap).*u[1] - q.*u[5] - muhM.*u[5].*u[3]
    du[6] = ßp.*q.*(cos.(t)+1).^3 - q.*u[6] - (mupB.*u[1] + mupE.*u[2]).*u[6]
end

u0 = [0.0004706; 0.0004706; 0.0004706; 9.7079; 7.9551; 32.061]
tspan = (0.0,50.0);
N = 2^18
t = 50
vars = 6

B = CSV.read("B.csv")
B_sol = solveODE(f, u0, vars, tspan, t, B, N)
writedlm("B_sols.csv",  B_sol, ',')

C_1 = CSV.read("C_1.csv")
C_1_sol = solveODE(f, u0, vars, tspan, t, C_1, N)
writedlm("C_1_sols.csv",  C_1_sol, ',')

C_2 = CSV.read("C_2.csv")
C_2_sol = solveODE(f, u0, vars, tspan, t, C_2, N)
writedlm("C_2_sols.csv",  C_2_sol, ',')

C_3 = CSV.read("C_3.csv")
C_3_sol = solveODE(f, u0, vars, tspan, t, C_3, N)
writedlm("C_3_sols.csv",  C_3_sol, ',')

C_4 = CSV.read("C_4.csv")
C_4_sol = solveODE(f, u0, vars, tspan, t, C_4, N)
writedlm("C_4_sols.csv",  C_4_sol, ',')

C_5 = CSV.read("C_5.csv")
C_5_sol = solveODE(f, u0, vars, tspan, t, C_5, N)
writedlm("C_5_sols.csv",  C_5_sol, ',')

C_6 = CSV.read("C_6.csv")
C_6_sol = solveODE(f, u0, vars, tspan, t, C_6, N)
writedlm("C_6_sols.csv",  C_6_sol, ',')
