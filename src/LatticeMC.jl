__precompile__()
"""
Lattice Monte Carlo Module with Julia
Metropolis Hastings
Wang Landau
"""
module LatticeMC
    export FCCLattice, BCCLattice, OneDLattice
    export mh_sampling, wld_allrange, wld_fixrange
    export init_coverage!, init_load!, update_lattice_energy_range!
    export read_interaction
    export analysis_mh, analysis_mh_deltaH, analysis_mh_rate, analysis_mh_rate_replaceAds
    export site_inter_energy, site_energy, total_energy
    export random_swap, delta_energy_swap, update_lattice_swap!

    include("lattice.jl")
    include("bcc_lattice.jl")
    include("fcc_lattice.jl")
    include("oned_lattice.jl")
    include("model.jl")
    include("post_analysis.jl")
    include("interaction.jl")
end     # module