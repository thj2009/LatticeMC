# Monte Carlo Simulation
# include("lattice.jl")
# include("interaction.jl")
# include("bcc_lattice.jl")
# include("fcc_lattice.jl")


using Random
using Formatting
using HDF5
using Statistics
using JLD
# import Plots

# Declare PHYSICIS Constant
const RG = 8.314        # J.K-1.mol-1


"""
Initialize the Lattice given coverage of certain types
"""
function init_coverage!(l::T, covs::Array{Float64,1}) where T <: Lattice
    ntot = l.ntot
    nx, ny = l.size
    l.lattice = fill(0, (nx, ny))
    shuffled_nums = shuffle(Vector(1:ntot))
    for (cov, name) in zip(covs, l.adsorbate)
        n = trunc(Int, cov * ntot)
        nums = shuffled_nums[1:n]
        shuffled_nums = shuffled_nums[n + 1:end]
        indices = [numsite_to_index(i, l.size) for i in nums]
        for index in indices
            l.lattice[index[1], index[2]] = name
        end
    end
    # Initialize site dictionary
    init_site_dict!(l)
    @assert sum([length(l.site_dict[ads]) for ads in l.adsorbate]) == l.size[1] * l.size[2]
end

"""
Load the lattice given lattice array
"""
function init_load!(l::T, larray::Array{Int,2}, sitedict::Bool=true) where T <: Lattice
    l.lattice = larray
    if sitedict
        # Initialize site dictionary
        init_site_dict!(l)
    end
end

"""
Initialize the site dictionary
"""
function init_site_dict!(l::T) where T <: Lattice
    nx, ny = l.size
    l.ntot = nx * ny
    adsorbate = l.adsorbate
    adsorbate_indices = [Tuple{Int,Int}[] for i in 1:length(adsorbate)]
    for i in 1:nx, j in 1:ny
        site = l.lattice[i, j]
        idx = findall(x->x == site, adsorbate)[1]
        push!(adsorbate_indices[idx], (i, j))
    end
    l.site_dict = Dict(adsorbate .=> adsorbate_indices)
end

"""
Swaping two sites and update the lattice
"""
function update_lattice_swap!(index1::Tuple{Int,Int}, index2::Tuple{Int,Int}, l::T) where T <: Lattice
    site1, site2 = l.lattice[index1[1], index1[2]], l.lattice[index2[1], index2[2]]
    # update the lattice
    l.lattice[index1[1], index1[2]], l.lattice[index2[1], index2[2]] = site2, site1
    # update the site dictionary of site 1
    filter!(x->x != index1, l.site_dict[site1])
    push!(l.site_dict[site1], index2)
    # update the site dictionary of site 2
    filter!(x->x != index2, l.site_dict[site2])
    push!(l.site_dict[site2], index1)
end

"""
Change index to site2 then update the lattice
"""
function update_lattice_change!(index::Tuple{Int,Int}, newsite::Int, l::T) where T <: Lattice
    site = l.lattice[index[1], index[2]]
    # update the lattice
    l.lattice[index[1], index[2]] = newsite
    # update the site dictionary of site 1
    filter!(x->x != index, l.site_dict[site])
    # update the site dictionary of site 2
    push!(l.site_dict[newsite], index)
end

"""



Pick two indices for randomly swap two sites on lattice
"""
function random_swap(l::T) where T <: Lattice
    available_ads = [i for i in l.adsorbate if length(l.site_dict[i]) > 0]
    @assert length(available_ads) >= 2
    shuffled_adsorbate = shuffle(available_ads)
    # Pick First Two Types
    site1 = shuffled_adsorbate[1]
    site2 = shuffled_adsorbate[2]
    # Randomly Pick two indices
    index1 = shuffle(l.site_dict[site1])[1]
    index2 = shuffle(l.site_dict[site2])[1]
    return index1, index2
end

"""
Metropolis-Hastings Sampling
"""
function mh_sampling(l::T1,
                     tem::Float64,
                     inters::Array{T2,1},
                     opts::Dict{String,Any},
                     verbose::Bool=true) where {T1 <: Lattice,T2 <: Interaction}
    # Unpack options
    tot_step = get(opts, "tot_step", 10000)
    burn_step = get(opts, "burn_step", 1000)
    write_step = get(opts, "write_step", 1000)
    save_step = get(opts, "save_step", 100)
    report = get(opts, "report", "")
    savefile = get(opts, "savefile", "mh_output.h5")

    if report != ""
        # open report file
    end

    saved_lattice = zeros(Int64, l.size[1], l.size[2], floor(Int, (tot_step - burn_step) / save_step))
    saved_E = zeros(Float64, floor(Int, (tot_step - burn_step) / save_step))
    if verbose
        printfmt("{1:^15s}   {2:^10s}   {3:^10s}  \n", "step(k)", "Accept", "E(kJ/mol)")
        println("=="^20)
    end
    # Statistical
    jump = 0
    @time begin
        E0 = total_energy(l, inters)
        for i in 1:tot_step
            # Propose swap of two sites
            index1, index2 = random_swap(l)
            # Calculate energy change
            dE = delta_energy_swap(index1, index2, l, inters)
            
            _a = min(1., exp(-dE * 1000 / (RG * tem)))
            _u = rand()
            if _u < _a           # Accept
                update_lattice_swap!(index1, index2, l)
                jump += 1
                E0 += dE
            else                 # Reject
                nothing
            end
            # Write Save and other steps
            if i % write_step == 0 && verbose
                printfmt("{1:^15d}   {2:10.2f}%  {3:10.2e}  \n", floor(Int, i / 1000), 100 * jump / i, E0)
            end
            if i % save_step == 0 && i > burn_step
                j = floor(Int, (i - burn_step) / save_step)
                saved_lattice[:, :, j] = copy(l.lattice)
                saved_E[j] = E0
            end
        end
    end

    if report != ""
        # close report file
    end

    # should also save the computational setting
    h5open(savefile, "w") do file
        write(file, "lattice", saved_lattice)
        write(file, "energy", saved_E)
    end

    # save(savefile, "lattice", saved_lattice, "energy", saved_E)

end







