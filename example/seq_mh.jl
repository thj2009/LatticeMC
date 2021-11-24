# Metropolis Hastings Sampling for BCC Lattice
# with NN interaction

#include("../src/LatticeMC.jl")
using LatticeMC
using LatinHypercubeSampling
using Random
Random.seed!(1)


function mh_one_adsorbate(c1::Float64, tem::Float64)
    fld = "./"
    # # Define BCC Lattice
    bcc = BCCLattice((50, 50), [0, 1])
    init_coverage!(bcc, [1 - c1, c1])

    # println(bcc.site_dict)
    c1 = length(bcc.site_dict[1]) / (bcc.size[1] * bcc.size[2])

    # Read Interactions
    interactions = read_interaction(joinpath(fld, "interaction.txt"))

    int_tem = Int(tem)
    # Metropolis Hastings 

    opts = Dict(
        "tot_step" => 1000000,
        "burn_step" => 20000,
        "write_step" => 10000,
        "save_step" => 1000,
        "report" => joinpath(fld, "result_$(int_tem)/out_$(c1).txt"),
        "savefile" => joinpath(fld, "result_$(int_tem)/mh_output_$(c1).h5")
    )

    mh_sampling(bcc, tem, interactions, opts)
    nothing
end




fld = "./"
t = 500.
mkdir(joinpath(fld, "result_$(Int(t))/"))
# LHS sample
mh_one_adsorbate(0.3, t)



