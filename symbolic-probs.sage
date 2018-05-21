import pdb


def xor(a, b):
    return((a).__xor__(b))


def get_symbolic_probs(n, letter="p"):
    return [var(letter + "_%d" % i) for i in xrange(n)]


def get_gamete_probs(n_sites, n_generations):
    n_seps = n_sites - 1  # number of seperations
    n_gametes = 2**n_sites

    #qs = get_symbolic_probs(2**4, "q")
    #qs = qs + [q for q in reversed(qs)]
    general_gamete_probs = get_symbolic_probs(n_gametes, "q")
    recomb_probs = get_symbolic_probs(n_seps, "p")
    recursive_formula = get_next_gamete_probs(general_gamete_probs, recomb_probs)

    gamete_probs = get_F1_gamete_probs(recomb_probs)

    for i in xrange(n_generations - 1):
        gamete_probs = eval_probs(recursive_formula,
                                  general_gamete_probs,
                                  gamete_probs)

    return expand_probs(gamete_probs)


def print_probs(gamete_probs):
    n_sites = log(len(gamete_probs), 2)
    format_str = "{0:0" + str(n_sites) + "b}"
    for i in xrange(len(gamete_probs)):
        print("\n" + format_str.format(i) + ":")
        print(gamete_probs[i])



def get_recomb_config_prob(recomb_config, recomb_probs):
    # Now the crossover locations required to form gamete_k from gamete_i and
    # gamete_j will be determined. By bitshifting k one position to the right,
    # XORing it with k, and masking it to remove any irrelevant leading bits the
    # result will have n_seps meaningful bits where a 0-bit indicates no
    # recombination at the corresponding seperation between sites and a 1-bit
    # indicated a recombination.  Now compute the probability of the gamete
    # being formed in exactly the way specified by k. In sites where a crossover
    # occurred
    prob = 1/4  # because there are 4 possible "base" gametes
    n_seps = len(recomb_probs)
    n_sites = n_seps + 1
    mask = 2**n_seps - 1
    crossovers = xor(recomb_config, recomb_config>>1) & mask
    for i in xrange(n_seps):
        # if crossovers[i] == 1
        if crossovers & 2**i:
            prob *= recomb_probs[i]
        # if crossovers[i] == 0
        else:
            prob *= 1 - recomb_probs[i]

    # If there are no recombinations (recomb is all 0-bits or all 1-bits) then
    # there is an additional 25% chance due to the biology
    if recomb_config == 0 or recomb_config == 2**n_sites - 1:
        prob += 1/4

    return prob
    # TODO(Jason): consider using sum of log(probs) rather than product of probs


def get_F1_gamete_probs(recomb_probs):
    n_sites = len(recomb_probs) + 1
    n_gametes = 2**n_sites
    return [get_recomb_config_prob(config, recomb_probs)
            for config in xrange(n_gametes)]


def get_starting_gamete_probs(n_sites):
    n_sep = n_sites - 1
    n_gametes = 2**n_sites
    recomb_probs = get_symbolic_probs(n_sep, letter="p")
    return [get_recomb_config_prob(config, recomb_probs)
            for config in xrange(n_gametes)]


def get_next_gamete_probs(prev_gamete_probs, recomb_probs):
    n_gametes = len(prev_gamete_probs)
    n_sites = len(recomb_probs) + 1

    # Create an array that will serve as a map to hold the probability of each
    # gamete occuring in the population. Because each gamete's structure will be
    # encoded in one-to-one corresdence with an integer from 0 to (n_gametes-1)
    # an array can be used instead of a python dictionary.
    next_gamete_probs = [0 for i in xrange(n_gametes)]

    # gamete_i and gamete_j are binary encodings of the gamete allele structure
    # where a 0-bit indicates an allele that comes from the first ancestral
    # parent and a 1-bit indicates an allele that comes from the second
    # ancestral parent (ancestral meaning one of the parents of the F1 plant as
    # opposed to the immediate parents of these gametes which may or may not be
    # different).
    for recomb_config in xrange(n_gametes):
        recomb_prob = get_recomb_config_prob(recomb_config, recomb_probs)
        #pdb.set_trace()
        for gamete_i in xrange(n_gametes):
            gamete_i_prob = prev_gamete_probs[gamete_i]
            for gamete_j in xrange(n_gametes):
                progeny_prob = gamete_i_prob * prev_gamete_probs[gamete_j]
                # gamete_i and gamete_j are the two gametes that will become the
                # two homologous chromosomes in the progeny

                # k does not directly represent a gamete in the way that
                # gamete_i and gamete_j do. Rather, every bit of k represents
                # whether the allele at that position comes from gamete_i or
                # gamete_j. Each bit in k is checked to see if it is 1 or 0. If
                # bit m is 0 then bit m in gamete_k will be assigned as bit m in
                # gamete_i and if bit m is 1 then bit m in gamete_k will be
                # assigned as bit m in gamete_j.
                gamete = 0
                for k in xrange(n_sites):
                    if recomb_config & 2**k:
                        gamete += gamete_j & 2**k
                    else:
                        gamete += gamete_i & 2**k

                next_gamete_probs[gamete] += recomb_prob * progeny_prob

    return expand_probs(next_gamete_probs)


def get_recursive_gamete_relation(n_sites):
    n_gametes = 2**n_sites
    n_seps = n_sites - 1

    general_gamete_probs = get_symbolic_probs(n_gametes, "q")
    recomb_probs = get_symbolic_probs(n_seps, "p")
    recursive_formula = get_next_gamete_probs(general_gamete_probs, recomb_probs)

    return recursive_formula


def eval_recursive_relation_with_recomb_probs(relation, recomb_probs):
    if log(len(relation), 2) - 1 != len(recomb_probs):
        exit(1)
    ps  = get_symbolic_probs(len(recomb_probs), "p")
    return eval_probs(relation, ps, recomb_probs)


def get_gamete_probs_rec(recomb_probs, n_gen):
    n_sites = 1 + len(recomb_probs)
    relation = get_recursive_gamete_relation(n_sites)
    f1_probs = get_F1_gamete_probs(recomb_probs)
    relation = eval_recursive_relation_with_recomb_probs(relation, recomb_probs)

    for i in xrange(n_generations - 1):
        gamete_probs = eval_probs(gamete,
                                  general_gamete_probs,
                                  gamete_probs)

    return expand_probs(gamete_probs)


def test():
    n_sites = 5
    n_generations = 5
    n_seps = n_sites - 1
    n_gametes = 2**n_sites

    general_gamete_probs = get_symbolic_probs(n_gametes, "q")
    recomb_probs = get_symbolic_probs(n_seps, "p")
    recursive_formula = get_next_gamete_probs(general_gamete_probs, recomb_probs)

    return recursive_formula


def expand_probs(probs):
    return [prob.expand() for prob in probs]


def eval_probs(probs, symbols, values):
    if (len(symbols) != len(values)):
        exit(1)
    assoc = {}
    for i in xrange(len(symbols)):
        assoc[symbols[i]] = values[i]
    evalled = [prob.subs(assoc) for prob in probs]
    return evalled


def plot_probs(probs):
    points = [(i, probs[i]) for i in range(len(probs))]
    sage.plot.point.point(points)

