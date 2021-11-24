using LatticeMC, Test

bcc = BCCLattice((2, 2), [0, 1])

latt = [0 0;
        1 1]

init_load!(bcc, latt)

@test bcc.site_dict == Dict(0 => [(1, 1), (1, 2)],
                            1 => [(2, 1), (2, 2)])


# interaction
interactions = read_interaction(joinpath("./", "interaction.txt"))

@test total_energy(bcc, interactions) == 2

using LatticeMC: delta_energy_change_pairs, count_k
@test delta_energy_change_pairs((1, 1), 1, (1, 2), 1, bcc, interactions) == 6

@test bcc.site_dict == Dict(0 => [(1, 1), (1, 2)],
                            1 => [(2, 1), (2, 2)])

# new lattice
bcc = BCCLattice((2, 2), [0, 1])

latt = [1 1;
        1 1]

init_load!(bcc, latt)
@test total_energy(bcc, interactions) == 8


bcc = BCCLattice((3, 2), [0, 1])

latt = [1 0;
        1 0
        1 0]

init_load!(bcc, latt)
@test total_energy(bcc, interactions) == 6

using LatticeMC: count_deltaH, count_pair

bcc = BCCLattice((3, 2), [0, 1])

latt = [1 0
        1 0
        1 0]

init_load!(bcc, latt)
@test count_deltaH(latt, bcc, (1, 1), (0, 0), interactions) == -6
@test bcc.lattice == [1 0; 1 0; 1 0]


# another test
bcc = BCCLattice((3, 3), [0, 1])

latt = [1 0 0
        1 0 0
        1 1 0]

init_load!(bcc, latt)
@test count_deltaH(latt, bcc, (1, 1), (0, 0), interactions) == -7
@test bcc.lattice == [1 0 0; 1 0 0; 1 1 0]

# another test
bcc = BCCLattice((3, 3), [0, 1])

latt = [1 0 0
        1 0 0
        1 1 1]


init_load!(bcc, latt)
@test count_deltaH(latt, bcc, (1, 1), (0, 0), interactions) == -52 / 6
@test count_deltaH(latt, bcc, (1, 1), (1, 0), interactions) == -32 / 6
@test count_pair(latt, bcc, (1, 1)) == 6 / 18
@test count_pair(latt, bcc, (1, 0)) == 8 / 18
@test count_pair(latt, bcc, (0, 0)) == 4 / 18



# another test

# interaction
interactions = read_interaction(joinpath("./", "interaction-2.txt"))
bcc = BCCLattice((5, 5), [0, 1, 2])

latt = [1 0 0 0 2
        2 1 2 0 0
        1 1 1 2 0
        0 0 0 1 2
        2 1 0 0 1]

init_load!(bcc, latt)

println(bcc.site_dict)

println(count_k(latt, bcc, (1, 2), (0, 0), interactions, 1.0, 0.0, 500.0, true))



# interaction
interactions = read_interaction(joinpath("./", "interaction-2.txt"))
bcc = BCCLattice((5, 5), [0, 1, 2])

latt = [1 0 0 0 2
        2 1 2 0 0
        1 1 1 2 0
        0 2 1 1 2
        2 1 0 0 1]

init_load!(bcc, latt)

println(bcc.site_dict)

println(count_k(latt, bcc, (1, 2), (0, 0), interactions, 1.0, 0.0, 500.0, true))
