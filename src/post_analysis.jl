# Post Analysis of simulation
# include("lattice.jl")

import HDF5
import JLD
import Plots
using Statistics
using Formatting
using StatsBase: harmmean
using Random
using Printf

const RG = 8.314        # J.K-1.mol-1

"Metropolis Hastings Sampling Result"
function analysis_mh(latticeFile::String, type::String,
                     adsorbates::Array{Int, 1}=[],
                     pairs::Array{Tuple{Int, Int}, 1}=Tuple{Int, Int}[],
                     verbose::Bool=true)
    lattice = HDF5.h5read(latticeFile, "lattice")
    energy = HDF5.h5read(latticeFile, "energy")

    normt, skew, kurt = normaltest(energy)

    if verbose
        println("\n")
        println("File: ", latticeFile)
        printfmt("Energy (kJ/mol), mean = {1:6.4e},  std = {2:6.4e} \n", mean(energy), std(energy))
        printfmt("                 min  = {1:6.4e},  max = {2:6.4e} \n", minimum(energy), maximum(energy))
        printfmt("Normal test for energy = {1:4.2f} \n", normt)
        printfmt("Skewness = {1:4.2f} \n", skew)
        printfmt("Kurtosis = {1:4,2f} \n", kurt)
        println("Size of lattice simulation = ", size(lattice))
    end
    nx, ny, nsim = size(lattice)
    # Create fictitous lattice
    if type == "BCC"
        fl = BCCLattice((nx, ny), adsorbates)
    elseif type == "FCC"
        fl = FCCLattice((nx, ny), adsorbates)
    elseif type == "OneD"
        fl = OneDLattice((nx, ny), adsorbates)
    end

    cov_list = Float64[]
    if verbose
        println("Coverage")
    end
    # count surface coverage
    for ad in adsorbates
        n = count(i->(i==ad), lattice[:, :, 1])
        cov = n / (nx * ny)
        push!(cov_list, cov)
        if verbose
            printfmt("Species {1:3d}:  {2:8.2f}\n", ad, cov)
        end
    end
    if verbose
        printfmt("Sum: {1:8.2f}\n\n", sum(cov_list))
    end
    if verbose
        println("Pair Probability Analysis")
        println("=="^20)
        printfmt("Index      {1:10s}   {2:10s}\n", "mean", "std")
    end
    prob_list = Float64[]
    for pair in pairs
        # count pair probability along lattice
        prob = [count_pair(floor.(Int, lattice[:, :, j]), fl, pair) for j in 1:nsim]
        push!(prob_list, mean(prob))
        # Test Convergenc
        if verbose
            printfmt("({1:3d}, {2:3d})   {3:10.2e}   {4:10.2e}\n", pair[1], pair[2], mean(prob), std(prob))
        end
    end
    if verbose
        printfmt("Sum: {1:8.2f}\n\n", sum(prob_list))
    end
    return cov_list, prob_list, energy
end


function analysis_mh_deltaH(latticeFile::String, type::String,
                            adsorbates::Array{Int, 1},
                            pairList1::Any,
                            pairList2::Any,
                            inters::Array{T, 1}) where T <: Interaction
    lattice = HDF5.h5read(latticeFile, "lattice")
    nx, ny, nsim = size(lattice)

    # Create fictitous lattice
    if type == "BCC"
        fl = BCCLattice((nx, ny), adsorbates)
    elseif type == "FCC"
        fl = FCCLattice((nx, ny), adsorbates)
    elseif type == "OneD"
        fl = OneDLattice((nx, ny), adsorbates)
    end

    ddE_list = Float64[]


    
    println("Delta Reaction Enthalpy Analysis")
    println("=="^20)
    for (pair1, pair2) in zip(pairList1, pairList2)
        # count pair probability along lattice
        dE_list = [count_deltaH(lattice[:, :, j], fl, pair1, pair2, inters) for j in 1:nsim]
        println("Pair1: ", pair1, "; Pair2: ", pair2, "; Mean of ddH = ", mean(dE_list))
        push!(ddE_list, mean(dE_list))
    end
    return ddE_list
end


function analysis_mh_rate(latticeFile::String, type::String,
                        adsorbates::Array{Int, 1},
                        pairList1::Any,
                        pairList2::Any,
                        inters::Array{T, 1},
                        omegaList::Array{Float64, 1},
                        dE0List::Array{Float64, 1},
                        tem::Float64) where T <: Interaction
    lattice = HDF5.h5read(latticeFile, "lattice")
    nx, ny, nsim = size(lattice)

    # Create fictitous lattice
    if type == "BCC"
        fl = BCCLattice((nx, ny), adsorbates)
    elseif type == "FCC"
        fl = FCCLattice((nx, ny), adsorbates)
    elseif type == "OneD"
        fl = OneDLattice((nx, ny), adsorbates)
    end

    kk_list = Float64[]
    cov_list = Float64[]
    println("\n")
    println("File: ", latticeFile)
    println("Coverage")
    scov = 0
    # count surface coverage
    for ad in adsorbates
        n = count(i->(i==ad), lattice[:, :, 1])
        cov = n / (nx * ny)
        push!(cov_list, cov)
        printfmt("Species {1:3d}:  {2:8.2f}\n", ad, cov)
        scov += cov
    end
    printfmt("Sum: {1:8.2f}\n\n", scov)
    println("Reaction Rate Analysis")
    println("=="^20)
    for (pair1, pair2, omega, dE) in zip(pairList1, pairList2, omegaList, dE0List)
        mf = 1
        for i in pair1
            mf *= cov_list[i+1]
        end
        # count pair probability along lattice
        k_list = [count_k(lattice[:, :, j], fl, pair1, pair2, inters, omega, dE, tem, false) for j in 1:nsim]
        mmf = mean(k_list)
        # k_list = Float64[]
        # for j in 1:nsim
        #     k, total = count_k(lattice[:, :, j], fl, pair1, pair2, inters, omega, dE, tem)
        #     k_list = vcat(k_list, k)
        # end

        # k, total = count_k(lattice[:, :, 1], fl, pair1, pair2, inters, omega, dE, tem)
        # println(k_list)
        # mmf = mean(k_list) / total

        # println("Pair1: ", pair1, "; Pair2: ", pair2, "; Mean of k = ", mean(k_list), 
        #     "MeanField k = ", mf, "Ratio = ", mean(k_list) / mf)
        print("Pair1: ", pair1, "; Pair2: ", pair2)
        printfmt(";   Mean of k = {1:10.3e}, MF k = {2:10.3e}, Ratio = {3:10.3e}\n",
            mmf, mf, mmf / mf)
        push!(kk_list, mmf)
    end
    return cov_list, kk_list
end

function analysis_mh_rate_replaceAds(latticeFile::String, type::String,
                                     adsorbates::Array{Int, 1},
                                     addads::Array{Int, 1}, addadscov::Array{Float64, 1}, nrep::Int,
                                     pairList1::Any,
                                     pairList2::Any,
                                     inters::Array{T, 1},
                                     omegaList::Array{Float64, 1},
                                     dE0List::Array{Float64, 1},
                                     tem::Float64) where T <: Interaction
    lattice = HDF5.h5read(latticeFile, "lattice")
    nx, ny, nsim = size(lattice)

    # Create fictitous lattice
    if type == "BCC"
        fl = BCCLattice((nx, ny), adsorbates)
    elseif type == "FCC"
        fl = FCCLattice((nx, ny), adsorbates)
    elseif type == "OneD"
        fl = OneDLattice((nx, ny), adsorbates)
    end

    kk_list = Float64[]
    cov_list = Float64[]
    println("\n")
    println("File: ", latticeFile)
    println("Coverage")
    scov = 0
    # count surface coverage
    for ad in adsorbates
        if (ad in addads) == false
            n = count(i->(i==ad), lattice[:, :, 1])
            cov = n / (nx * ny)
        else
            idx = findall(x->x==ad, addads)[1]
            cov = addadscov[idx]
        end
        push!(cov_list, cov)
    end
    cov_list[1] = 1 - sum(cov_list[2:end])
    for (ad, cov) in zip(adsorbates, cov_list)
        printfmt("Species {1:3d}:  {2:8.2e}\n", ad, cov)
        scov += cov
    end
    printfmt("Sum: {1:8.2e}\n\n", scov)
    println("Reaction Rate Analysis")
    println("=="^20)
    for (pair1, pair2, omega, dE) in zip(pairList1, pairList2, omegaList, dE0List)
        mf = 1
        for i in pair1
            mf *= cov_list[i+1]
        end
        # count pair probability along lattice
        k_list = Float64[]
        for i in 1:nrep
            kk = [count_k_replaceAds(lattice[:, :, j], fl, addads, addadscov, pair1, pair2, inters, omega, dE, tem, false) for j in 1:nsim]
            k_list = [k_list; kk]
        end
        mmf = mean(k_list)

        print("Pair1: ", pair1, "; Pair2: ", pair2)
        printfmt(";   Mean of k = {1:10.3e}, MF k = {2:10.3e}, Ratio = {3:10.3e}\n",
                 mmf, mf, mmf / mf)
        push!(kk_list, mmf)
    end
    return cov_list, kk_list
end

function convergence_test(x::Array{Float64, 1})
end

"""Count number of pairs in lattice"""
function count_pair(larray::Array{Int,2}, l::T, pair::Tuple{Int, Int}) where T <: Lattice
    nx, ny = size(larray)

    count::Int = 0
    total::Int = 0
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            # create neighbor list
            nn1 = nn1_list((i, j), l)
            for nn in nn1
                site2 = larray[nn[1], nn[2]]
                if (site1 == site2 && (site1, site2) == pair) || 
                        (site1 != site2 && ((site1, site2) == pair || (site2, site1) == pair))
                    count += 1
                end
                total += 1
            end
        end
    end
    prob = count / total
end

"""
normal test of an array
return skewness^2 + kurtosis^2
skewness tells you the amount and direction of skew (departure from horizontal symmetry)
Kurtosis tells you the height and sharpness of the central peak, relative to that of a standard bell curve.
"""
function normaltest(x::Array{Float64, 1})
    # https://help.gooddata.com/doc/en/reporting-and-dashboards/maql-analytical-query-language/maql-expression-reference/aggregation-functions/statistical-functions/predictive-statistical-use-cases/normality-testing-skewness-and-kurtosis
    n = length(x)
    m = mean(x)
    _var = var(x)
    skew = sum((x .- m) .^ 3 / n) / _var^(1.5)
    kurt = sum((x .- m) .^ 4 / n) / _var^(2) - 3
    return skew*skew + kurt*kurt, skew, kurt
end



function count_deltaH(larray::Array{Int,2}, l::T1,
                      pair::Tuple{Int, Int}, newpair::Tuple{Int, Int},
                      inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    nx, ny = size(larray)
    
    dE_list = Float64[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            # create neighbor list
            nn1 = nn1_list((i, j), l)
            for nn in nn1
                site2 = larray[nn[1], nn[2]]
                if (site1 == site2 && (site1, site2) == pair) || 
                        (site1 != site2 && ((site1, site2) == pair || (site2, site1) == pair))
                    # load array to 
                    init_load!(l, larray, false)
                    dE = delta_energy_change_pairs((i, j), newpair[1],
                                                   (nn[1], nn[2]), newpair[2],
                                                   l, inters)
                    push!(dE_list, dE)
                end
            end
        end
    end
    if size(dE_list, 1) == 0
        return 0
    else
        return mean(dE_list)
    end
end


function count_deltaH(larray::Array{Int,2}, l::T1,
                      ads::Tuple{Int}, newads::Tuple{Int},
                      inters::Array{T2,1}) where {T1 <: Lattice,T2 <: Interaction}
    nx, ny = size(larray)

    dE_list = Float64[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            if site1 == ads[1]
                # load array to 
                init_load!(l, larray, false)
                dE = delta_energy_change((i, j), newads[1], l, inters)
                push!(dE_list, dE)
            end
        end
    end
    if size(dE_list, 1) == 0
        return 0
    else
        return mean(dE_list)
    end
end

function count_k(larray::Array{Int,2}, l::T1,
                 pair::Tuple{Int, Int}, newpair::Tuple{Int, Int},
                 inters::Array{T2,1},
                 omega::Float64, dE0::Float64,
                 tem::Float64,
                 verbose::Bool=false) where {T1 <: Lattice, T2 <: Interaction}
    nx, ny = size(larray)
    # print(nx, ny)
    count::Int = 0
    total::Int = 0

    k_list = Float64[]
    index_list = Tuple{Int,Int}[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            numsite1 = index_to_numsite((i, j), (nx, ny))
            # create neighbor list
            nn1 = nn1_list((i, j), l)
            for nn in nn1
                site2 = larray[nn[1], nn[2]]
                numsite2 = index_to_numsite((nn[1], nn[2]), (nx, ny))
                # if (site1 == site2 && (site1, site2) == pair) || 
                #         (site1 != site2 && ((site1, site2) == pair || (site2, site1) == pair)) && 
                #         !((numsite1, numsite2) in index_list)
                if (site1, site2) == pair && !((numsite1, numsite2) in index_list)
                    # push!(index_list, (numsite1, numsite2))
                    # push!(index_list, (numsite2, numsite1))
                    # println((i, j), (nn[1], nn[2]))
                    # load array to 
                    init_load!(l, larray, false)
                    dE = delta_energy_change_pairs((i, j), newpair[1],
                                                   (nn[1], nn[2]), newpair[2],
                                                   l, inters)
                    # dEa = omega * max(-(dE- dE0), 0)
                    dEa = omega * (-(dE- dE0))
                    k = exp(dEa * 1000 / (RG * tem))
                    # @printf("%d to %d, dE=%.2f, k=%.2f\n", index_to_numsite((i, j), (nx, ny)), index_to_numsite((nn[1], nn[2]), (nx, ny)), dE, k)
                    push!(k_list, k)
                    count += 1
                end
                total += 1
            end
        end
    end
    if verbose
        println(k_list, total)
    end
    # print(sum(k_list) / total)
    # print("\n")
    # return k_list, total
    if size(k_list, 1) == 0
        return 0
    else
        # return harmmean(k_list) * size(k_list, 1) / total
        return sum(k_list) / total
    end
end


function count_k(larray::Array{Int,2}, l::T1,
                 ads::Tuple{Int}, newads::Tuple{Int},
                 inters::Array{T2,1},
                 omega::Float64, dE0::Float64,
                 tem::Float64,
                 verbose::Bool=false) where {T1 <: Lattice, T2 <: Interaction}
    nx, ny = size(larray)
    count::Int = 0
    total::Int = 0

    k_list = Float64[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            if site1 == ads[1]
                # load array to 
                init_load!(l, larray, false)
                dE = delta_energy_change((i, j), newads[1], l, inters)
                # dEa = omega * max(-(dE- dE0), 0)
                dEa = omega * (-(dE- dE0))
                k = exp(dEa * 1000 / (RG * tem))
                push!(k_list, k)
                count += 1
            end
            total += 1
        end
    end
    if verbose
        println(k_list, total)
    end
    # return k_list, total
    if size(k_list, 1) == 0
        return 0
    else
        # return harmmean(k_list) * size(k_list, 1) / total
        #return size(k_list, 1) / total / (sum(1.0 / k_list))
        return sum(k_list) / total
    end
end


function count_k_replaceAds(larray::Array{Int,2}, l::T1, addads::Array{Int, 1}, addadscov::Array{Float64, 1},
                            pair::Tuple{Int, Int}, newpair::Tuple{Int, Int},
                            inters::Array{T2,1},
                            omega::Float64, dE0::Float64,
                            tem::Float64,
                            verbose::Bool=false) where {T1 <: Lattice, T2 <: Interaction}
    nx, ny = size(larray)
    # replace 0 with new adds
    # println(larray)
    # replace 0 with new adds
    index0 = findall(x->x==0, larray)
    for (ads, cov) in zip(addads, addadscov)
        n = floor(Int, cov * nx * ny)
        # println(cov)
        # println(n)
        cov = n / (nx * ny)
        index0 = shuffle(index0)
        # println(index0)
        for idx in index0[1:n]
            larray[idx] = ads
        end
        index0 = index0[n + 1 : end]
    end
    # println(larray)
    count::Int = 0
    total::Int = 0

    k_list = Float64[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            # create neighbor list
            nn1 = nn1_list((i, j), l)
            for nn in nn1
                site2 = larray[nn[1], nn[2]]
                if (site1 == site2 && (site1, site2) == pair) || 
                        (site1 != site2 && ((site1, site2) == pair || (site2, site1) == pair))
                    # load array to 
                    init_load!(l, larray, false)
                    dE = delta_energy_change_pairs((i, j), newpair[1],
                                                    (nn[1], nn[2]), newpair[2],
                                                    l, inters)
                    # dEa = omega * max(-(dE- dE0), 0)
                    dEa = omega * (-(dE- dE0))
                    k = exp(dEa * 1000 / (RG * tem))
                    # print(k)
                    push!(k_list, k)
                    count += 1
                end
                total += 1
            end
        end
    end
    if verbose
        println(k_list, total)
    end
    # print(sum(k_list) / total)
    # print("\n")
    # return k_list, total
    if size(k_list, 1) == 0
        return 0
    else
        # return harmmean(k_list) * size(k_list, 1) / total
        return sum(k_list) / total
    end
end


function count_k_replaceAds(larray::Array{Int,2}, l::T1, addads::Array{Int, 1}, addadscov::Array{Float64, 1},
                            ads::Tuple{Int}, newads::Tuple{Int},
                            inters::Array{T2,1},
                            omega::Float64, dE0::Float64,
                            tem::Float64,
                            verbose::Bool=false) where {T1 <: Lattice, T2 <: Interaction}
    nx, ny = size(larray)
    # println(larray)
    # replace 0 with new adds
    index0 = findall(x->x==0, larray)
    for (ads, cov) in zip(addads, addadscov)
        n = floor(Int, cov * nx * ny)
        # println(cov)
        # println(n)
        cov = n / (nx * ny)
        index0 = shuffle(index0)
        # println(index0)
        larray[index0[1:n]] .= ads
        index0 = index0[n + 1 : end]
    end
    # println(larray)
    count::Int = 0
    total::Int = 0

    k_list = Float64[]
    for i in 1:nx
        for j in 1:ny
            site1 = larray[i, j]
            if site1 == ads[1]
                # load array to 
                init_load!(l, larray, false)
                dE = delta_energy_change((i, j), newads[1], l, inters)
                # dEa = omega * max(-(dE- dE0), 0)
                dEa = omega * (-(dE- dE0))
                k = exp(dEa * 1000 / (RG * tem))
                push!(k_list, k)
                count += 1
            end
            total += 1
        end
    end
    if verbose
        println(k_list, total)
    end
    # return k_list, total
    if size(k_list, 1) == 0
        return 0
    else
        # return harmmean(k_list) * size(k_list, 1) / total
        #return size(k_list, 1) / total / (sum(1.0 / k_list))
        return sum(k_list) / total
    end
end
