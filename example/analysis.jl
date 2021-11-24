# Metropolis Hastings Sampling for BCC Lattice
# with NN interaction

# include("../../src/LatticeMC.jl")
using LatticeMC
using HDF5
using Statistics

# t = parse(Int, ARGS[1])
path = "./result_200/"
files = readdir(path)

covs = Array{Float64,1}[]
pairs = Array{Float64,1}[]
for file in files
    cov, pair, e = analysis_mh(joinpath(path, file), "BCC",
                            [0, 1],
                            [(0, 0), (1, 1), (0, 1)])
    push!(covs, cov)  
    push!(pairs, pair)
end


interactions = read_interaction("./interaction.txt")
dH = analysis_mh_rate("./result_200/mh_output_0.3.h5", "BCC",
                                    [0, 1],
                                    [(1, 1), (1,)], [(0, 0), (0,)],
                                    interactions,
                                    [1.0, 1.0],
                                    [-8.368, 0.0], 500.0)
print(dH)

# lt = load("./result_200/mh_output_0.3.jld")["lattice"]
# print(lt[:, :, 1])



# h5open("cov_pair_$(t).h5", "w") do file
#     write(file, "covs", hcat(covs...))
#     write(file, "pairs", hcat(pairs...))
# end

# print(hcat(pairs...)')
