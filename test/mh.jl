
# Compare with old result

using LatticeMC


bcc = BCCLattice((50, 50), [0, 1, 2])

c1 = 0.6
c2 = 0.1

init_coverage!(bcc, [1 - c1 - c2, c1, c2])


# Read Interactions
interactions = read_interaction("interaction-3.txt")

opts = Dict(
    "tot_step" => 1000000,
    "burn_step" => 20000,
    "write_step" => 10000,
    "save_step" => 1000,
    "report" => "out.txt",
    "savefile" => "mh_output.h5"
)

mh_sampling(bcc, 500.0, interactions, opts)



cov, pair, e = analysis_mh("mh_output.h5", "BCC",
                            [0, 1, 2],
                            [(0, 0), (1, 2)])