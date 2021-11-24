"""
One D strap with peridoic condition along y axis
"""
mutable struct OneDLattice <: Lattice
    size :: Tuple{Int, Int}
    ntot :: Int
    adsorbate :: Array{Int, 1}
    lattice :: Array{Int, 2}
    site_dict :: Dict{Int, Array{Tuple{Int, Int}, 1}}
    OneDLattice() = new()
end

function OneDLattice(size::Tuple{Int, Int}, ads::Array{Int, 1}=Int[])
    l = OneDLattice()
    l.size = size
    l.ntot = size[1] * size[2]
    l.adsorbate = ads
    return l
end

# Nearest Neighbor List
function nn1_list(index::Tuple{Int, Int}, l::OneDLattice)
    ix, iy = index
    nx, ny = l.size
    if nx == 1
        list = [periodic_index((ix, iy + 1), l.size),
                periodic_index((ix, iy - 1), l.size)]
    else
        if ix == 1
            list = [periodic_index((ix, iy + 1), l.size),
                    periodic_index((ix, iy - 1), l.size),
                    periodic_index((ix + 1, iy), l.size)]
        elseif ix == nx
            list = [periodic_index((ix, iy + 1), l.size),
                    periodic_index((ix, iy - 1), l.size),
                    periodic_index((ix - 1, iy), l.size)]
        else
            list = [periodic_index((ix + 1, iy), l.size),
                    periodic_index((ix - 1, iy), l.size),
                    periodic_index((ix, iy + 1), l.size),
                    periodic_index((ix, iy - 1), l.size)]
        end
    end
    return list
end
