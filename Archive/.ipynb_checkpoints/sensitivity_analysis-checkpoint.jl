"""
sobol(f, N, u0, tspan, t, k, (*$\theta$*), index = "first-order")

Creates a matrix of Sobol' sensitivities where each row contains 
the sensitivity indices for a parameter across different model 
outputs.

 f - system of ODEs
 N - number of simulations (can vary from hundreds to thousands)
 u0 - vector of initial conditions for ODE
 tspan - time span
 k - vector of scale parameters from Gamma((*$k_i$*),(*$\theta_i$*)) for all 
     model parameters
 (*$\theta$*) - vector of shape parameters from Gamma((*$k_i$*),(*$\theta_i$*)) for all 
     model parameters
 index - default is the first order index, but can take on the
         other value "total-order"
 t - specified time in tspan to evaluate each model variable
"""

using Sobol
using Distributions
function sobol(f, N, u0, tspan, k, (theta), index = "first-order",
               t = 5)   
    K = length(k)
    vars = length(u0)
# quasi-randomly selected points
    S = zeros(N,2*K)
    s = SobolSeq(2*K)
    (theta) = [(theta);(theta)]
    k = [k;k]
    
    S = zeros(N,2*K)
    for sim in 1:N
        S[sim,:] = next!(s)' 
        for parameter in 1:2*K
            S[sim,parameter] = quantile(Gamma(k[parameter], 
                               (theta)[parameter]), S[sim,parameter])
        end
    end
    
    A = S[:, 1:K] # accounts for half of the random sample
    B = S[:, (K+1):2*K] # accounts for second half of the 
                        # random sample
    
    sums = zeros(4,vars)
    ya = solveODE(f, u0, vars, tspan, t, A, N)   
    yb = solveODE(f, u0, vars, tspan, t, B, N)   
    sumC = zeros(K,vars) # (parameter, variable)
    
    for p in 1:K # iterates through changing one column of C at 
                 # a time (parameters)
        C = Ci(p,A,B)
        yc = solveODE(f, u0, vars, tspan, t, C, N)
        
        for v in 1:vars
            if index == "first-order"
                sumC[p,v] = yc[:,v]'*ya[:,v]/N
            else
                sumC[p,v] = yc[:,v]'*yb[:,v]/N

            end
        end
    end  
    
    sens = zeros(K,vars) # (parameter, variable)
    

    for v in 1:vars 
        sums[1,v] = ya[:,v]'*ya[:,v]/N
      
                    # sum(ya^2)/N

        sums[2,v] = (sum(ya[:,v])/N)^2 # [sum(ya)/N]^2 
        sums[3,v] = sums[1, v] - sums[2,v] 
                    # denominator for sensitivity calculation
        sums[4,v] = sum(yb[:,v])*sum(ya[:,v])/N.^2 
                    # sum(ya)*sum(yb)/N^2
        
        for p in 1:K
            if index == "first-order"
                sens[p,v] = (sumC[p,v] - sums[4,v])/sums[3,v] 
            else
                sens[p,v] = 1 - (sumC[p,v] - sums[2,v])/sums[3,v]
            end
        end
    end
    return sens # matrix of sensitivity indices
end

"""
Ci(i, A, B)

Creates a new matrix C by replacing the ith column of matrix B 
with the ith column of matrix A.

 i - column number to be changed in matrix B
 A - matrix containing half of the pseudo-random samples
 B - matrix containing the other half of the pseudo-random 
     samples
"""

function Ci(i, A, B)
    C = zeros(size(B))
    for j in 1:size(A)[2]
        if j == i
            C[:,i] = A[:,i]
        else
            C[:,j] = B[:,j] 
        end
    end
    return C
end

"""
solveODE(f, u0, vars, tspan, t, M, N)

Solves the system of ODEs f and creates a matrix of model 
evaluations for each simulation and output variable. Each row 
corresponds to model solutions at a given time t for one 
simulation.

 f - system of ODEs
 u0 - vector of initial conditions for the ODEs
 vars - number of output variables
 tspan - time span
 t - specified time in tspan to evaluate each model variable
 M - matrix with rows as parameter values and columns as 
     simulations
 N - number of simulations (can vary from hundreds to thousands)
"""

function solveODE(f, u0, vars, tspan, t, M, N)
    y = zeros(N,vars) # (simulation, variable)
    for sim in 1:N
        p = M[sim,:]
        prob = ODEProblem(f,u0,tspan,p)
        sol = solve(prob)
        for variable in 1:vars
            y[sim,variable] = sol.u[t][variable]
        end
    end 
    return y # matrix of model outputs
end