#!/sage
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


# The from_to_matrix will be such that if the current state is i, the
# probability of transitioning to state j is from_to_matrix[i][j]
def get_from_to_matrix(generation):
    prog_probs = get_progeny_probs(generation)
    ftmat = [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]]

    prob0 = prog_probs[0][0] + prog_probs[0][1] + prog_probs[1][0] + prog_probs[1][1]
    if (prob0 == 0):
        prob0 = 1  # if 0 all terms are 0 so dividing by 1 prevents error

    ftmat[0][0] = prog_probs[0][0] / prob0
    ftmat[0][1] = prog_probs[0][1] / prob0
    ftmat[0][2] = prog_probs[1][0] / prob0
    ftmat[0][3] = prog_probs[1][1] / prob0


    prob1 = prog_probs[0][2] + prog_probs[0][3] + prog_probs[1][2] + prog_probs[1][3]
    if (prob1 == 0):
        prob1 = 1  # if 0 all terms are 0 so dividing by 1 prevents error

    ftmat[1][0] = prog_probs[0][2] / prob1
    ftmat[1][1] = prog_probs[0][3] / prob1
    ftmat[1][2] = prog_probs[1][2] / prob1
    ftmat[1][3] = prog_probs[1][3] / prob1


    prob2 = prog_probs[2][0] + prog_probs[2][1] + prog_probs[3][0] + prog_probs[3][1]
    if (prob2 == 0):
        prob2 = 1  # if 0 all terms are 0 so dividing by 1 prevents error

    ftmat[2][0] = prog_probs[2][0] / prob2
    ftmat[2][1] = prog_probs[2][1] / prob2
    ftmat[2][2] = prog_probs[3][0] / prob2
    ftmat[2][3] = prog_probs[3][1] / prob2


    prob3 = prog_probs[2][2] + prog_probs[2][3] + prog_probs[3][2] + prog_probs[3][3]
    if (prob3 == 0):
        prob3 = 1  # if 0 all terms are 0 so dividing by 1 prevents error

    ftmat[3][0] = prog_probs[2][2] / prob3
    ftmat[3][1] = prog_probs[2][3] / prob3
    ftmat[3][2] = prog_probs[3][2] / prob3
    ftmat[3][3] = prog_probs[3][3] / prob3

    return ftmat


def fwd_bkw(observations, states, start_prob, trans_prob, emm_prob, end_st):
    # forward part of the algorithm
    fwd = []
    f_prev = {}
    for i, observation_i in enumerate(observations):
        f_curr = {}
        for st in states:
            if i == 0:
                # base case for the forward part
                prev_f_sum = start_prob[st]
            else:
                prev_f_sum = sum(f_prev[k]*trans_prob[k][st] for k in states)

            f_curr[st] = emm_prob[st][observation_i] * prev_f_sum

        fwd.append(f_curr)
        f_prev = f_curr

    p_fwd = sum(f_curr[k] * trans_prob[k][end_st] for k in states)

    # backward part of the algorithm
    bkw = []
    b_prev = {}
    for i, observation_i_plus in enumerate(reversed(observations[1:]+(None,))):
        b_curr = {}
        for st in states:
            if i == 0:
                # base case for backward part
                b_curr[st] = trans_prob[st][end_st]
            else:
                b_curr[st] = sum(trans_prob[st][l] * emm_prob[l][observation_i_plus] * b_prev[l] for l in states)

        bkw.insert(0,b_curr)
        b_prev = b_curr

    p_bkw = sum(start_prob[l] * emm_prob[l][observations[0]] * b_curr[l] for l in states)

    # merging the two parts
    posterior = []
    for i in range(len(observations)):
        posterior.append({st: fwd[i][st] * bkw[i][st] / p_fwd for st in states})

    assert p_fwd - p_bkw < 1e-10
    #return fwd, bkw, posterior
    return posterior


states = ('Healthy', 'Fever')
end_state = 'E'

observations = ('normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy',
                'normal', 'normal', 'cold', 'dizzy', 'normal', 'normal', 'cold', 'dizzy')

start_probability = {'Healthy': 0.6, 'Fever': 0.4}

transition_probability = {
    'Healthy' : {'Healthy': 0.69, 'Fever': 0.3, 'E': 0.01},
    'Fever' : {'Healthy': 0.4, 'Fever': 0.59, 'E': 0.01},
}

emission_probability = {
    'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
    'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
}


def example():
    return fwd_bkw(observations,
                   states,
                   start_probability,
                   transition_probability,
                   emission_probability,
                   end_state)
