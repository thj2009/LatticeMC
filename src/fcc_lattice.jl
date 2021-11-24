# include("lattice.jl")

mutable struct FCCLattice <: Lattice
    size :: Tuple{Int, Int}
    ntot :: Int
    adsorbate :: Array{Int, 1}
    lattice :: Array{Int, 2}
    site_dict :: Dict{Int, Array{Tuple{Int, Int}, 1}}
    FCCLattice() = new()
end

function FCCLattice(size::Tuple{Int, Int}, ads::Array{Int, 1}=Int[])
    l = FCCLattice()
    l.size = size
    l.ntot = size[1] * size[2]
    l.adsorbate = ads
    return l
end

function nn1_list(index::Tuple{Int, Int}, l::FCCLattice)
    ix, iy = index
    list = [periodic_index((ix+1, iy), l.size),
            periodic_index((ix-1, iy), l.size),
            periodic_index((ix, iy+1), l.size),
            periodic_index((ix, iy-1), l.size),
            periodic_index((ix-1, iy+1), l.size),
            periodic_index((ix+1, iy-1), l.size)]
    return list
end

function nn2_list(index::Tuple{Int, Int}, l::FCCLattice)
    ix, iy = index
    list = [periodic_index((ix+1, iy+1), l.size),
            periodic_index((ix-1, iy+2), l.size),
            periodic_index((ix-2, iy+1), l.size),
            periodic_index((ix-1, iy-1), l.size),
            periodic_index((ix+1, iy-2), l.size),
            periodic_index((ix+2, iy-1), l.size)]
    return list
end

function nn3_list(index::Tuple{Int, Int}, l::FCCLattice)
    ix, iy = index
    list = [periodic_index((ix+2, iy), l.size),
            periodic_index((ix-2, iy), l.size),
            periodic_index((ix, iy+2), l.size),
            periodic_index((ix, iy-2), l.size),
            periodic_index((ix-2, iy+2), l.size),
            periodic_index((ix+2, iy-2), l.size)]
    return list
end

function nn4_list(index::Tuple{Int, Int}, l::FCCLattice)
end

function nn5_list(index::Tuple{Int, Int}, l::FCCLattice)
end
