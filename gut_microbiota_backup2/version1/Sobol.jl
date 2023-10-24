# solveODE - function that returns a matrix of variable outputs (solutions to the differential equation at time t) for
# each simulation N
# M - matrix with rows as parameter values and columns as simulations
using DifferentialEquations
function solveODE1(f, u0, vars, tspan, t, M, N)
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

# f - system of ODEs
# N - number of simulations (can vary from hundreds to thousands)
# u0 - vector of initial conditions for ODE
# tspan - time span
# μ - vector of means from N(μi,σi^2) for all parameters
# σ - vector of standard deviations from N(μi,σi^2) for all parameters
# index - default is the first order index, but can take on the other value "total-order"
# sobol - function that returns the sensitivity indicies for every parameter/variable combination

using Sobol
using Distributions
function sobol(f, N, u0, tspan, α, θ, index = "first-order", t = 5)
    k = length(α)
    vars = length(u0)

    # quasi-randomly selected points
    S = zeros(N,2*k)
    s = SobolSeq(2*k)
    θ = [θ;θ]
    α = [α;α]

    S = zeros(N,2*k)
    for sim in 1:N
        S[sim,:] = next!(s)' # horizontal vector with the simulated values for each parameter
        for parameter in 1:2*k
            S[sim,parameter] = quantile(Gamma(α[parameter], θ[parameter]),S[sim,parameter])
        end
    end

    A = S[:, 1:k] # accounts for half the random sample
    B = S[:, (k+1):2*k] # accounts for second half of the random sample

    sums = zeros(4,vars)
    ya = solveODE1(f, u0, vars, tspan, t, A, N)
    yb = solveODE1(f, u0, vars, tspan, t, B, N)
    sumC = zeros(k,vars) # (parameter, variable)

    for parameter in 1:k # iterates through changing one column of C at a time (parameters)
        C = copy(B)
	C[:,parameter] = A[:,parameter]
        yc = solveODE1(f, u0, vars, tspan, t, C, N)

        for variable in 1:vars
            if index == "first-order"
                sumC[parameter,variable] = yc[:,variable]'*ya[:,variable]/N
            else
                sumC[parameter,variable] = yc[:,variable]'*yb[:,variable]/N
            end
        end
    end

    sensitivities = zeros(k,vars) # (parameter, variable)

    for variable in 1:vars
        sums[1,variable] = ya[:,variable]'*ya[:,variable]/N # sum(ya^2)/N
        sums[2,variable] = (sum(ya[:,variable])/N)^2 # [sum(ya)/N]^2
        sums[3,variable] = sums[1, variable] - sums[2, variable] # denominator for sensitivity calculation
        sums[4,variable] = sum(yb[:,variable])*sum(ya[:,variable])/N.^2 # sum(ya)*sum(yb)/N^2

        for parameter in 1:k
            if index == "first-order"
                sensitivities[parameter,variable] = (sumC[parameter,variable] - sums[4,variable])/sums[3,variable]
            else
                sensitivities[parameter,variable] = 1 - (sumC[parameter,variable] - sums[2,variable])/sums[3,variable]
            end
        end
    end
    return sensitivities
end

ψ(u, γ) = u./(γ+u)
using DifferentialEquations
#q = 0.05
function g(du,u,p0,t)
    B, E, M, a, h, p = u
    βb, γp, q, βE1, γa, γB, βE2, βM1, βM2, γh, βa, μaE, μaM, βh1, βh2, μhM, βp, μpB, μpE = p0

    du[1]= βb.*ψ(u[6],γp).*u[1] - q.*u[1]
    du[2] = (βE1.*ψ(u[4], γa) + βE2.*(1-ψ(u[1], γB)).*ψ(u[6], γp)).*u[2] - q.*u[2]
    du[3] = (βM1.*ψ(u[4], γa) + βM2.*ψ(u[5], γh)).*u[3] - q.*u[3]

    du[4] = βa.*ψ(u[6], γp).*u[1] - q.*u[4] - (μaE.*u[2] + μaM.*u[3]).*u[4]
    du[5] = βh1.*ψ(u[4], γa).*u[2] + βh2.*ψ(u[6], γp).*u[1] - q.*u[5] - μhM.*u[5].*u[3]
    du[6] = βp.*q.*(cos.(t)+1).^3 - q.*u[6] - (μpB.*u[1] + μpE.*u[2]).*u[6]
end
u0 = [0.0004706; 0.0004706; 0.0004706; 9.7079; 7.9551; 32.061]
tspan = (0.0,50.0);

k = [5.0e6;
      0.6;
      0.4;
      0.3;
     75.0;
  16666.67;
      0.40;
      0.25;
   5000.0;
    100.0;
      5.0;
     75.0;
    200.0;
  12500.0;
  25000.0;
     20.83;
 100000.0;
   2500.0;
      0.025]
θ = fill(2.0,19);

firstOrder = sobol(g, 2^18, u0, tspan, k, θ)

F = zeros(6,1)
for i in 1:6
    F[i,1] = sum(firstOrder[:,i])
end
F

using DelimitedFiles
writedlm("First-Order-2-10.csv",  firstOrder, ',')

totalOrder = sobol(g, 2^18, u0, tspan, k, θ, "total-effects")

T = zeros(6,1)
for i in 1:6
    T[i,1] = sum(totalOrder[:,i])
end
T

writedlm("Total-Order-2-10.csv",  totalOrder, ',')
