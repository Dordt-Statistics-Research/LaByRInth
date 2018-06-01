import pdb


def get_progeny_probs(generation):
    if (generation < 1):
        exit("generation must be 1 or greater")

    # F1 probabilities. Must be heterozygous
    probs = [[0,   0, 0, 1/2],
             [0,   0, 0, 0  ],
             [0,   0, 0, 0  ],
             [1/2, 0, 0, 0  ]]

    for gen in range(1, generation):
        probs = _next_probs_helper(probs, gen)

    return probs


def _next_probs_helper(progeny_probs, gen):
    # set up the symbolic variables
    p = var("p_%d" % gen)
    q = var("q_%d" % gen)

    # set up the probabilities of recombination
    ps = [1/2 * (1-p),
          1/2 * p,
          1/2 * p,
          1/2 * (1-p)]

    qs = [1/2 * (1-q),
          1/2 * q,
          1/2 * q,
          1/2 * (1-q)]

    # initiate the next progeny probabilities
    new_probs = [[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0]]

    # for each first possible homolog
    for h1 in range(4):
        # for each second possible homolog
        for h2 in range(4):
            # probability of the progeny defined by the two homologs occurring
            # in the population
            prev_prob = progeny_probs[h1][h2]

            # define possible gametes that can be produced from the progeny
            gametes = [h1,
                      (h1 & 0b10) + (h2 & 0b01),
                      (h2 & 0b10) + (h1 & 0b01),
                      h2]

            # for each possible gamete formed under p probabilities
            for i in range(4):
                g1 = gametes[i]     # new gamete 1
                # for each possible gamete formed under q probabilities
                for j in range(4):
                    g2 = gametes[j] # new gamete 2
                    new_probs[g1][g2] += ps[i]*qs[j]*prev_prob

    return new_probs

