"""
type in the command line:

julia generate_params.jl {num_samples} {begin_idx} {save_directory}
"""
# include modules from utils.jl
const ROOT_DIR = joinpath(homedir(), ".projectflow/profiles")
include("utils.jl")

using QuasiMonteCarlo
using ArgParse
using JLD2
using DifferentialEquations
using Dates


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "num_samples"
            help = "number of parameter samples to generate"
            arg_type = Int
            required = true
        "idx"
            help = "index to save the output with a unique filename"
            required = true
        "save_directory"
            help = "folder where you want to save the output, assumes that the directory already exists"
            required = true
    end

    return parse_args(s)
end

#parse the arguments
println(Dates.format(now(), "HH:MM:SS") )
args = parse_commandline()
#println(args["num_samples"])

#############################
# create the A and B samples#
#############################

# specify the upper and lower limits of the parameter ranges
upper_b =  [1e6, 10,10,10,1000,1e5,1e5,10,10,1e5,1e4,1e4,1e4,1e4,1e6,1e6,1e5,1e8,1e8,10]
lower_b = repeat([0],length(upper_b))

# sample A and B
sampler =LatticeRuleSample() # generates a sample -- seed is not fixed so we will get different draws every time this is called
A = QuasiMonteCarlo.sample(Int(args["num_samples"]/2), lower_b, upper_b, sampler)
B = QuasiMonteCarlo.sample(Int(args["num_samples"]/2), lower_b, upper_b, sampler) # A, B are shapes (# params, num_samples)

i = args["idx"]

# save the parameter matrix
save_object(args["save_directory"] * "/A$i.jld2", A)
save_object(args["save_directory"] * "/B$i.jld2", B)

##############
# compute C_i#
##############

# produce the Ci matrices, i = 1, ..., # of parameters
C = [Ci(i,A,B) for i in 1:length(upper_b)]

# save the C matrices
save_object(args["save_directory"] * "/C$i.jld2", C)

###############################
# save the results for A and B#
###############################
# compute the results for A and B for all samples
res_A = sample(A)

println("starting to solve B")
println(Dates.format(now(), "HH:MM:SS") )
#println("solved for A")
res_B  = sample(B)
#println("solved for B")

# save the value of the loss function for each parameter
#save_object(args["save_directory"] * "/loss_A$i.jld2", loss_A)
#save_object(args["save_directory"] * "/loss_B$i.jld2", loss_B)

# save the results from the ODE solve
save_object(args["save_directory"] * "/res_A$i.jld2", res_A)
save_object(args["save_directory"] * "/res_B$i.jld2", res_B)


#########################
# save the results for C#
#########################

println("starting to solve C")
println(Dates.format(now(), "HH:MM:SS") )
# produce the results for the Ci matrices, i = 1, ..., # of parameters
res_C = [sample(C[n]) for n in 1:length(upper_b)]

# save the ODE solution and the losses
save_object(args["save_directory"] * "/res_C$i.jld2", res_C) 
println(Dates.format(now(), "HH:MM:SS") )