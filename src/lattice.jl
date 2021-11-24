
include("interaction.jl")

abstract type Lattice end

function numsite_to_index(i::Int, size::Tuple{Int,Int})
    nx, ny = size
    return (mod1(i, nx),
            (i - 1) ÷ nx + 1)
end

function index_to_numsite(index::Tuple{Int,Int}, size::Tuple{Int,Int})
    nx, ny = size
    x, y = index
    return (y - 1) * nx + x
end


function periodic_index(index::Tuple{Int,Int}, size::Tuple{Int,Int})
    nx, ny = size
    ix, iy = index
    return (ix % nx > 0 ? ix % nx : ix % nx + nx,
            iy % ny > 0 ? iy % ny : iy % ny + ny)
end

function total_energy(l::T1, inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    e = 0
    nx, ny = l.size
    for i in 1:nx, j in 1:ny
        index = (i, j)
        # calculate site energy
        es, _ = site_energy(index, l, inters)
        e += es
    end
    return e
end

"""
dE when swap two indices (Enew - Eold)
"""
function delta_energy_swap(index1::Tuple{Int,Int}, index2::Tuple{Int,Int}, l::T1, inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    type1 = l.lattice[index1[1], index1[2]]
    type2 = l.lattice[index2[1], index2[2]]
    if type1 == type2
        return 0
    else
        # dE of change index1 to type2
        dE1 = delta_energy_change(index1, type2, l, inters)
        # update index 1 to type2
        l.lattice[index1[1], index1[2]] = type2
        # dE of change index2 to type1
        dE2 = delta_energy_change(index2, type1, l, inters)
        # restore the index 1
        l.lattice[index1[1], index1[2]] = type1
        return dE2 + dE1
    end
end

"""
dE when change the index to other type (Enew - Eold)
"""
function delta_energy_change(index::Tuple{Int,Int}, newtype::Int, l::T1, inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    oldtype = l.lattice[index[1], index[2]] 
    if oldtype == newtype
        return 0
    else
        oldlattice = copy(l.lattice)
        old_e, old_elist = site_energy(index, l, inters)
        # Update the lattice
        l.lattice[index[1], index[2]] = newtype
        new_e, new_elist = site_energy(index, l, inters)
        # Restore the lattice
        l.lattice[index[1], index[2]] = oldtype
        dE = sum(new_elist) - sum(old_elist)
        return dE
    end
end

"""
dE when change the index 1 to other type 1, and index 2 to other type 2 (Enew - Eold)
"""
function delta_energy_change_pairs(index1::Tuple{Int,Int}, newtype1::Int,
                                   index2::Tuple{Int,Int}, newtype2::Int,
                                   l::T1, inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    oldtype1, oldtype2 = l.lattice[index1[1], index1[2]], l.lattice[index2[1], index2[2]]
    dE1 = delta_energy_change(index1, newtype1, l, inters)
    # Update the first index on the lattice
    l.lattice[index1[1], index1[2]] = newtype1
    dE2 = delta_energy_change(index2, newtype2, l, inters)
    # Restore the lattice
    l.lattice[index1[1], index1[2]] = oldtype1
    return dE1 + dE2
end

"""
Array contain the energy with different degree of interaction
"""
function site_energy(index::Tuple{Int,Int}, l::T1, inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    elist = zeros(5)
    total_e = 0
    for inter in inters
        order, e = site_inter_energy(index, l, inter)
        elist[order] += e
        total_e += e / order
    end
    return total_e, elist
end


"""
EmptySite: Calculate the site energy with specify interaction 
"""
function site_inter_energy(index::Tuple{Int,Int}, l::T, inter::EmptySite) where T <: Lattice
    l.lattice[index[1], index[2]] == 0 ? (1, inter.energy) : (1, 0)
end

"""
OneSite: Calculate the site energy with specify interaction 
"""
function site_inter_energy(index::Tuple{Int,Int}, l::T, inter::OneSite) where T <: Lattice
    l.lattice[index[1], index[2]] == inter.site ? (1, inter.energy) :  (1, 0)
end

"""
TwoSites: Calculate the site energy with specify interaction 
"""
function site_inter_energy(index::Tuple{Int,Int}, l::T, inter::TwoSites) where T <: Lattice
    site = l.lattice[index[1], index[2]]
    count = 0
    site_set_list = Set[]
    if in_pair(site, inter.pair)
        site2 = other_site(site, inter.pair)
        # create neighbor list
        nn_list = nn_enum(index, l, inter.pair.config)
        for nn_index in nn_list
            if l.lattice[nn_index[1], nn_index[2]] == site2
                site_set = Set([index, nn_index])
                if !(site_set in site_set_list)
                    count += 1
                    push!(site_set_list, site_set)
                end
            end
        end
    end
    return (2, count * inter.energy)
end


"""
MultiSites: Calculate the site energy with specify interaction 
Maybe need recursive match 
"""
function site_inter_energy(index1::Tuple{Int,Int}, l::T, inter::MultiSites) where T <: Lattice
    site1 = l.lattice[index1[1], index1[2]]
    count = 0
    site_set_list = Set[]
    for first in inter.sites
        if site1 == first[1]
            # searching the pair with s in it
            for pair in inter.pairs
                if in_pair(first, pair)
                    second = other_site(first, pair)
                    config_12 = pair.config
                    # search site neighbor to check match
                    nn_list_12 = nn_enum(index1, l, config_12)
                    for index2 in nn_list_12
                        site2 = l.lattice[index2[1], index2[2]]
                        if site2 == second[1]
                            # search the common neighbor of those two
                            third = inter.sites[findfirst(x -> x!=first && x!=second, inter.sites)]
                            config_23, config_13 = "", ""
                            for pair in inter.pairs
                                if in_pair(first, pair) && in_pair(third, pair)
                                    config_13 = pair.config
                                end
                                if in_pair(second, pair) && in_pair(third, pair)
                                    config_23 = pair.config
                                end
                            end
                            # search the common neighbor of 1 and 2
                            nn_list_13 = nn_enum(index1, l, config_13)
                            nn_list_23 = nn_enum(index2, l, config_23)
                            for index3 in nn_list_13
                                site3 = l.lattice[index3[1], index3[2]]
                                if index3 in nn_list_23 && site3 == third[1]
                                    site_set = Set([index1, index2, index3])
                                    if !(site_set in site_set_list)
                                        count += 1
                                        push!(site_set_list, site_set)
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    return (3, count * inter.energy)
end

"""

enumerate neighbor list given site and configuration
"""
function nn_enum(index::Tuple{Int,Int}, l::T, config::String) where T <: Lattice
    if config == "NN1"
        nn_list = nn1_list(index, l)
    elseif config == "NN2"
        nn_list = nn2_list(index, l)
    elseif config == "NN3"
        nn_list = nn3_list(index, l)
    elseif config == "NN4"
        nn_list = nn4_list(index, l)
    end
    return nn_list
end