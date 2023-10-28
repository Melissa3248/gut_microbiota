using JLD2 

# function that preprocesses either the res_A or res_B files
function preprocess_AB(r)
    return transpose(reduce(hcat,r))
end

function preprocess_one_C_file(C)
    l = []
    for j in 1:length(C)
        test = C[j][1] # shape: (vars, n_sub_simulations) 
        for i in 2:length(C[j])
            test = hcat(test, C[j][i])
        end
        push!(l,test)
    end
    return l
end

n_files = 512


# preprocessing ODE results from the B parameters
res_B = preprocess_AB(load("/path/to/res_B0.jld2")["single_stored_object"])

for i in 1:(n_files-1)
    file = load("/path/to/res_B"*string(i)*".jld2")["single_stored_object"]
    global res_B = vcat(res_B, preprocess_AB(file))
end



save("/path/to/res_B.jld2", "data", res_B)

print("saved res_B")

# preprocessing ODE results from the A parameters
res_A = preprocess_AB(load("/path/to/res_A0.jld2")["single_stored_object"])

for i in 1:(n_files-1)
    file = load("/path/to/res_A"*string(i)*".jld2")["single_stored_object"]
    global res_A = vcat(res_A, preprocess_AB(file))
end

save("/path/to/res_A.jld2", "data", res_A)
print("saved res_A")

# preprocessing ODE results from the C parameters

#initialize the vector with the first file
res_C = preprocess_one_C_file(load("/path/to/res_C0.jld2")["single_stored_object"])
#print((res_C[1]))
res_C = [transpose(res_C[n]) for n in 1:length(res_C)] # create a 20 long list of (n_simulations, n_params) matrix
for f in 1:(n_files-1) #iterate through all the files
    file = load("/path/to/res_C"*string(f)*".jld2")["single_stored_object"]
    global res_C = [vcat(res_C[n],transpose(preprocess_one_C_file(file)[n])) for n in 1:length(res_C)] # concatenate along the n_simulations dimension
    #if f == 250
    #    print("reached 250")
    #end
end

save("/path/to/res_C.jld2",  "data", res_C)
print("saved res_C")
