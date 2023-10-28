using DifferentialEquations
using Statistics

Ψ(a, γ_a) = a ./ (a + γ_a)
"""
function f(du,u,p0,t)
    B, E, M, a, h, p = u
    betaa, betab, betaE1, betaE2, betah1, betah2,betah3, betaM1, betaM2, betap, gammaa, gammaB, gammah, gammap, muaE, muaM, muhM, mupB, mupE, q = p0

    du[1]= betab.*phi(p,gammap).*B - q.*B
    du[2] = (betaE1.*phi(a, gammaa) + betaE2.*(1-phi(B, gammaB)).*phi(p, gammap)).*E - q.*E
    du[3] = (betaM1.*phi(a, gammaa) + betaM2.*phi(h, gammah)).*M - q.*M

    du[4] = betaa.*phi(p, gammap).*B - q.*a - (muaE.*E + muaM.*M).*phi(a, gammaa)
    du[5] = betah1.*phi(a, gammaa).*E + betah2.*phi(p, gammap).*B + betah3.*(1-phi(B,gammaB)).*phi(p,gammap)*E- q.*h- muhM.*phi(h, gammah).*M
    du[6] = betap.*q.*(cos.(t)+1).^3 - q.*p - (mupB.*B + mupE.*E).*phi(p, gammap)
end
"""
function f(du,u,p₀,t)
    B, E, M, a, h, p = u
    β_a, β_b, β_E₁, β_E₂, β_h₁, β_h₂, β_h₃, β_M₁, β_M₂, β_p, γ_a, γ_B, γ_h, γ_p, μ_aE, μ_aM, μ_hM, μ_pB, μ_pE, q = p₀

    du[1] = β_b.*Ψ(p,γ_p).*B - q.*B
    du[2] = (β_E₁.*Ψ(a,γ_a) + β_E₂.*(1-Ψ(B,γ_B)).*Ψ(p,γ_p)).*E - q.*E
    du[3] = (β_M₁.*Ψ(a,γ_a) + β_M₂.*Ψ(h,γ_h)).*M - q.*M

    du[4] = β_a.*Ψ(p,γ_p).*B - q.*a - (μ_aE.*E + μ_aM.*M).*Ψ(a,γ_a)
    du[5] = β_h₁.*Ψ(a,γ_a).*E + β_h₂.*Ψ(p,γ_p).*B + β_h₃.*(1-Ψ(B,γ_B)).*Ψ(p,γ_p)*E - q.*h - μ_hM.*Ψ(h,γ_h).*M
    du[6] = β_p.*q.*(cos.(t)+1).^3 - q.*p - (μ_pB.*B + μ_pE.*(1-Ψ(B,γ_B)).*E).*Ψ(p,γ_p)
end

"""
Ci(i, A, B)

Creates a new matrix C by replacing the ith column of matrix B 
with the ith column of matrix A.

 i - column number to be changed in matrix B
 A - matrix containing half of the pseudo-random samples (# params, num_samples)
 B - matrix containing the other half of the pseudo-random (# params, num_samples)
     samples
"""

function Ci(i, A, B) 
    C = copy(B) # super necessary to copy B or else it will get overwritten
    C[i,:] = A[i,:]
    return C
end


function loss(truth, sol, tpts, max_t)
    """ loss computes MAPE (averaged across timepoints length(sol) - tpts  to length(sol)) and ground truth weighted by a weight vector
    truth (Vector): vector containing the ground truth
    pred (Matrix): matrix with shape (variables, timepoints)
    t_start (integer): first timepoint to pull to compute an average across time
    t_end (integer): last timepoint to pull to compute an average across time
    """
    
    # average the losses from the 3 species since we only observed their summed biomass
    s = abs((sum(truth[1:3]) - sum(center(sol, tpts,max_t)[1:3,:])) ./sum(truth[1:3]))
    
    # compute the remaining losses for p, h, and a and return the average
    return (s + sum(broadcast(abs, (truth[4:6]-center(sol, tpts,max_t)[4:6,:]) ./ truth[4:6] )))/4
end

"""
sample: takes in a matrix of parameter values and outputs the correponding output

p: shape (# parameters, # samples)

Output: length # simulations vector with solutions (length vars)
"""
function sample(p)
    # p has shape (params, # samples)
    l = size(p)[2]
    u0 = [0.0004706; 0.0004706; 0.0004706; 9.7079; 7.9551; 32.061]
    max_t = 250.0
    tspan = (0.0,max_t)
    
    # create the ODEProblem for each parameter sample
    prob = [ODEProblem(f,u0,tspan,p[:,n]) for n in 1:l]
    
    # solve the ODEProblem for each parameter sample
    #sol = [solve(prob[n], Tsit5(),reltol=1e-6) for n in 1:l]
    sol = [solve(prob[n],Rodas4()) for n in 1:l]

    # pull only the last 100 timepoints of the ODE solution to compute metrics for
    tpts = 50
    
    #losses = [loss(u0, sol[n], tpts, max_t) for n in 1:l]
    
    return [mean(sol[n], dims = 2) for n in 1:l]#sol #[center(sol[n], tpts, max_t) for n in 1:l] , losses
end


center(sol, tpts, max_t) = (maximum(sol((max_t-tpts):max_t), dims=2) + minimum(sol((max_t-tpts):max_t), dims=2))/2


"""
C: list of (# simulations, # parameters) matrices, this list has length # of parameters
N: number of simulations from A and B
"""
function sobol(A_sol, B_sol, C_sol, K, vars, N, index = "first-order")   
    
    sums = zeros(4,vars)
    ya = A_sol #shape: (N, vars)  #solveODE(f, u0, vars, tspan, t, A, N)   
    yb = B_sol # shape: (N, vars) #solveODE(f, u0, vars, tspan, t, B, N)   
    sumC = zeros(K,vars) # (parameter, variable)
    
    for p in 1:K # iterates through changing one column of C at 
                 # a time (parameters)
        #C = Ci(p,A,B)
        yc = C_sol[p] #index the pth C (replaces B with the pth column of C) #solveODE(f, u0, vars, tspan, t, C, N)
                     # shape: (N, vars)
        
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
