using LatticeMC
using LatticeMC: count_k

# interaction
interactions = read_interaction(joinpath("./", "interaction-2.txt"))

bcc = BCCLattice((3, 3), [0, 1, 2])

latt = [1 0 0
        1 0 0
        2 1 1]

count_k(latt, bcc,
        (1, 2), (0, 0),
        interactions,
        1.0, 0.0,
        500.0, true)



count_k(latt, bcc,
        (1, ), (0, ),
        interactions,
        1.0, 0.0,
        500.0, true)


latt = [1 2 2
        2 2 2
        2 2 1]
println("=========================")
count_k(latt, bcc,
        (1, 2), (0, 0),
        interactions,
        1.0, 0.0,
        500.0, true)



count_k(latt, bcc,
        (1, ), (0, ),
        interactions,
        1.0, 0.0,
        500.0, true)


latt = [1 0 0
        2 2 0
        1 1 0]
println("=========================")
count_k(latt, bcc,
        (1, 2), (0, 0),
        interactions,
        1.0, 0.0,
        500.0, true)

count_k(latt, bcc,
        (0, 0), (0, 0),
        interactions,
        0.0, 0.0,
        500.0, true)


count_k(latt, bcc,
        (1, ), (0, ),
        interactions,
        1.0, 0.0,
        500.0, true)