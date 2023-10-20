"""
type in the command line:

julia generate_params.jl --num_samples {} --begin_idx {} --save_directory {}
"""

using GlobalSensitivity, QuasiMonteCarlo
using ArgParse

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
    #C = zeros(size(B))
    #for j in 1:size(A)[2]
    #    if j == i
    #        C[:,i] = A[:,i]
    #    else
    #        C[:,j] = B[:,j] 
    #    end
    #end
    #return C
    
    B[:,i] = A[:,i]
    return B
end


function loss(truth, sol, tpts, max_t)#, truth, t_start, t_end, weight)
    """ loss computes MAPE (averaged across timepoints length(sol) - tpts  to length(sol)) and ground truth weighted by a weight vector
    truth (Vector): vector containing the ground truth
    pred (Matrix): matrix with shape (variables, timepoints)
    t_start (integer): first timepoint to pull to compute an average across time
    t_end (integer): last timepoint to pull to compute an average across time
    """
    
    s = abs((sum(truth[1:3]) - sum(center(sol, tpts,max_t)[1:3,:])) ./sum(truth[1:3]))
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
    max_t = 7500.0
    tspan = (0.0,max_t)
    
    # create the ODEProblem for each parameter sample
    prob = [ODEProblem(f,u0,tspan,p[:,n]) for n in 1:l]
    
    # solve the ODEProblem for each parameter sample
    sol = [solve(prob[n], Tsit5(),reltol=1e-6) for n in 1:l]

    # pull only the last 100 timepoints of the ODE solution to compute metrics for
    tpts = 100
    
    losses = [loss(u0, sol[n], tpts, max_t) for n in 1:l]
    
    return [center(sol[n], tpts, max_t) for n in 1:l], losses
end

center(sol, tpts, max_t) = (maximum(sol((max_t-tpts):max_t), dims=2) + minimum(sol((max_t-tpts):max_t), dims=2))/2

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "num_samples"
            help = "number of parameter samples to generate"
            required = true
        "idx"
            help = "index to save the output with a unique filename"
            required = true
        "save_directory"
            help = "folder where you want to save the output, assumed that the directory already exists"
            required = true
    end

    return parse_args(s)
end

#parse the arguments
args = parse_commandline()
println(args["num_samples"])

#############################
# create the A and B samples#
#############################

# specify the upper and lower limits of the parameter ranges
upper_b =  [1e6, 10,10,10,1000,1e5,1e5,10,10,1e5,1e4,1e4,1e4,1e4,1e6,1e6,1e5,1e8,1e8,10]
lower_b = repeat([0],length(p))

# sample A and B
sampler = SobolSample()
A,B = QuasiMonteCarlo.generate_design_matrices(args["num_samples"],lower_b,upper_b,sampler) # A, B are shapes (# params, num_samples)

# compute the results for A and B for all samples
res_A, loss_A = sample(A)
res_B, loss_B  = sample(B)

################################
# compute C_i and solve the ODE#
################################

# produce the Ci matrices
C = [Ci(i,A,B) for i in 1:args["num_samples"]]

# produce the results for the Ci matrices
res_C = [sample(C[n]) for n in 1:args["num_samples"]]

###################
# save the results#
###################
using JLD2

i = args["idx"]

save_object(args.["save_directory"] + "/A$i.jld2", A)
save_object(args.["save_directory"] + "/B$i.jld2", B)
save_object(args.["save_directory"] + "/loss_A$i.jld2", loss_A)
save_object(args.["save_directory"] + "/loss_B$i.jld2", loss_B)
save_object(args.["save_directory"] + "/res_A$i.jld2", res_A)
save_object(args.["save_directory"] + "/res_B$i.jld2", res_B)

save_object(args.["save_directory"] + "/C$i.jld2", C)
save_object(args.["save_directory"] + "/res_C$i.jld2", res_C) # includes loss