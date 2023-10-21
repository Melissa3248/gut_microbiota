using DifferentialEquations
using Statistics

phi(a, gamma_a) = a ./ (a + gamma_a)

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

Output: 
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