abstract type Interaction end

# TODO: how to implement multiple site interaction
# The Unit of the energy is kJ/mol
struct Pair
    s1 :: Tuple{Int, String}
    s2 :: Tuple{Int, String}
    config :: String
end

in_pair(i::Int, p::Pair) = i == p.s1[1] || i == p.s2[2]
in_pair(i::Tuple{Int, String}, p::Pair) = i == p.s1 || i == p.s2

other_site(i::Int, p::Pair) = i == p.s1[1] ? p.s2[1] : p.s1[1]
other_site(i::Tuple{Int, String}, p::Pair) = i == p.s1 ? p.s2 : p.s1

mutable struct EmptySite <: Interaction
    energy :: Float64
end

mutable struct OneSite <: Interaction
    name :: String
    energy :: Float64
    site :: Int
end

mutable struct TwoSites <: Interaction
    name :: String
    energy :: Float64
    pair :: Pair

end

mutable struct MultiSites <: Interaction
    name :: String
    energy :: Float64
    sites :: Array{Tuple{Int, String}, 1}
    pairs :: Array{Pair, 1}
end


#TODO Need Reformulate this !!!
"parse interaction from file"
function read_interaction(filename :: String)
    #TODO: need nicer way to read those
    inter_list = Interaction[]
    open(filename) do file
        for ln in eachline(file)
            if length(ln) >= 0 && ln[1] != "#"
                sln = split(ln, ",")
                sln = [strip(s) for s in sln]
                name = sln[2]
                energy = parse.(Float64, sln[3])
                pairs = sln[4]
                if sln[1] == "zero"         # Empty Site
                elseif sln[1] == "one"      # Single Site
                    sites, _ = read_pairs(pairs)[1]
                    s = OneSite(name, energy, sites[1])
                elseif sln[1] == "two"      # Two Sites
                    (site1, site2), config = read_pairs(pairs)[1]
                    # config = read_pairs(pairs)[2]
                    pair = Pair(site1, site2, config)
                    s = TwoSites(name, energy, pair)
                elseif sln[1] == "multi"    # Multi Sites
                    allsite = Tuple{Int, String}[]
                    pair_list = Pair[]
                    for pair in read_pairs(pairs)
                        (site1, site2), config = pair
                        !(site1 in allsite) ? push!(allsite, site1) : nothing
                        !(site2 in allsite) ? push!(allsite, site2) : nothing
                        push!(pair_list, Pair(site1, site2, config))
                    end
                    # reduce all site by half
                    s = MultiSites(name, energy, allsite, pair_list)
                end
                push!(inter_list, s)
            end
        end
    end
    return inter_list
end

function read_pairs(ln::SubString{String})
    # "/" divide pair
    # ";" divide sites and config
    pair_list = Tuple{Any, Any}[]
    pairs = strip.(split(ln, "/"))
    for pair in pairs
        spair = strip.(split(pair, ";"))
        ss = split.(strip.(split(spair[1], "!")))
        if length(ss) == 1
            sites = [parse(Int, ss[1][1])]
        else
            sites = [(parse(Int, ss[1][1]), ss[1][2]),
                     (parse(Int, ss[2][1]), ss[2][2])]
        end
        config = spair[2]
        push!(pair_list, (sites, config))
    end
    return pair_list
end
