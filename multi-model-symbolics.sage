#!/sage
import pdb
import operator as op


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


def eval_mat(mat, ps, qs):
    return [[entry.subs(p_1 = ps[0],
                        p_2 = ps[1],
                        p_3 = ps[2],
                        p_4 = ps[3],
                        q_1 = qs[0],
                        q_2 = qs[1],
                        q_3 = qs[2],
                        q_4 = qs[3]) for entry in row] for row in mat]


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

    # p_fwd = sum(f_curr[k] * trans_prob[k][end_st] for k in states)
    p_fwd = sum(f_curr[k] for k in states)

    # backward part of the algorithm
    bkw = []
    b_prev = {}
    for i, observation_i_plus in enumerate(reversed(observations[1:]+(None,))):
        b_curr = {}
        for st in states:
            if i == 0:
                # base case for backward part
                # b_curr[st] = trans_prob[st][end_st]
                b_curr[st] = 1
            else:
                b_curr[st] = sum(trans_prob[st][l] * emm_prob[l][observation_i_plus] * b_prev[l] for l in states)

        bkw.insert(0,b_curr)
        b_prev = b_curr

    p_bkw = sum(start_prob[l] * emm_prob[l][observations[0]] * b_curr[l] for l in states)

    # merging the two parts
    posterior = []
    for i in range(len(observations)):
        posterior.append({st: fwd[i][st] * bkw[i][st] / p_fwd for st in states})

    print(p_fwd)
    print(p_bkw)
    assert p_fwd - p_bkw < 1e-10
    return fwd, bkw, posterior
    # return posterior


# def fwd_bkw(observations, states, start_prob, trans_prob, emm_prob):
#     # forward part of the algorithm
#     fwd = []
#     f_prev = []
#     for i, observation_i in enumerate(observations):
#         f_curr = []
#         for st in states:
#             if i == 0:
#                 # base case for the forward part
#                 prev_f_sum = start_prob[st]
#             else:
#                 prev_f_sum = sum(f_prev[k]*trans_prob[i-1][k][st] for k in states)

#             f_curr[st] = emm_prob[st][i] * prev_f_sum

#         fwd.append(f_curr)
#         f_prev = f_curr

#     # p_fwd = sum(f_curr[k] * trans_prob[k][end_st] for k in states)
#     p_fwd = sum(f_curr[k] for k in states)

#     # backward part of the algorithm
#     bkw = []
#     b_prev = {}
#     for i, observation_i_plus in enumerate(reversed(observations[1:]+(None,))):
#         b_curr = {}
#         for st in states:
#             if i == 0:
#                 # base case for backward part
#                 b_curr[st] = trans_prob[st][end_st]
#             else:
#                 b_curr[st] = sum(trans_prob[st][l] * emm_prob[l][observation_i_plus] * b_prev[l] for l in states)

#         bkw.insert(0,b_curr)
#         b_prev = b_curr

#     p_bkw = sum(start_prob[l] * emm_prob[l][observations[0]] * b_curr[l] for l in states)

#     # merging the two parts
#     posterior = []
#     for i in range(len(observations)):
#         posterior.append({st: fwd[i][st] * bkw[i][st] / p_fwd for st in states})

#     assert p_fwd - p_bkw < 1e-10
#     #return fwd, bkw, posterior
#     return posterior


def get_trans_probs(positions, generation, p_model, q_model):
    DIST = 1e6
    gen_prob_mat = get_from_to_matrix(generation)
    trans_prob_arr = []

    for i in xrange(1, len(positions)):
        dist = positions[i] - positions[i-1]
        # r is the probability that an odd number of physical recombinations
        # occurrs between sites i-1 and i. This is twice the value that
        # LB-Impute uses because they use it as the probability that a crossover
        # is observed in the progeny which is half of the probability that the
        # crossover actually occurs.
        r = 1 - exp(-1.0 * dist / DIST)
        ps = [r if cross_allowed else 0 for cross_allowed in p_model]
        qs = [r if cross_allowed else 0 for cross_allowed in q_model]

        trans_prob_arr.append(eval_mat(gen_prob_mat, ps, qs))

    return trans_prob_arr


def impute_model(observations, trans_probs, emm_probs):
    states = range(4)
    start_prob = [1.0 / len(states) for i in states]

    def merge_het(state_probs):
        return [state_probs[0],
                state_probs[1] + state_probs[2],
                state_probs[3]]

    def best_state(state_probs):
        best_state = None
        best_prob = 0
        for st in range(len(states)):
            if state_probs[st] >= best_prob:
                best_state = st
                best_prob = state_probs[st]
        return best_state

    posteriors = fwd_bkw(observations,
                         states,
                         start_prob,
                         trans_probs,
                         emm_probs)


    return best_state(merge_het(posteriors))


# return array of differences between sequential entries
def diffs(arr):
    [j-i for i, j in zip(arr[:-1], arr[1:])]


# model_p is a vector of booleans indicating if recombination is allowed in
# each intergeneration
def dist_to_ps_or_qs(dist, model_p_or_q):
    return [dist_to_crossover_prob(dist, model) for model in model_p_or_q]


def dist_to_crossover_prob(dist, model):
    DIST = 1e6
    if model:
        return 1 - exp(-1.0 * dist / DIST)
    else:
        return 0

# def trans_prob(model_p, model_q, dist, k, st):


states = ('Healthy', 'Fever')
end_state = 'E'

observations = ('normal', 'cold', 'normal', 'dizzy')

start_probability = {'Healthy': 0.6, 'Fever': 0.4}

# transition_probability = {
#     'Healthy' : {'Healthy': 0.69, 'Fever': 0.3, 'E': 0.01},
#     'Fever' : {'Healthy': 0.4, 'Fever': 0.59, 'E': 0.01},
# }

transition_probability = {
    'Healthy' : {'Healthy': 0.70, 'Fever': 0.3},
    'Fever' : {'Healthy': 0.4, 'Fever': 0.60},
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


def choose(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, xrange(n, n-r, -1), 1)
    denom = reduce(op.mul, xrange(1, r+1), 1)
    return numer//denom


# def emission(st, obs, rerr):
#     a, b = obs
#     n = a+b
#     if state == "AA":
#         return choose(n, a) * (1-rerr)**a * (rerr)**b
#     elif state == "AB":
#         return choose(n, a) * (rerr)**a * (1-rerr)**b
#     elif state == "BA" || state == "BB":
#         return choose(n, a) * (0.5)**n
#     else:
#         exit("illegal state")

def test():
    ca = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 0, 1, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 2, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 0, 0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 2, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 0, 2, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0]
    cb = [0, 0, 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 1, 3, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 0, 1, 4, 1, 0, 0, 1, 0, 0, 0, 2, 0, 1, 2, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 2, 1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 3, 0, 1, 2, 0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 2, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 3, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    observations = zip(ca, cb)
    positions = [1158042, 1229010, 1865353, 3613321, 3793316, 4090263, 6788188, 7288027, 11181197, 11954091, 12772631, 12778039, 27845033, 28226634, 29994192, 31346621, 31346637, 31670641, 31831318, 31992368, 32859095, 32994613, 33332936, 33571430, 35050210, 35328863, 35491722, 36562581, 38912987, 45592345, 46054117, 46054154, 46143354, 46143369, 46229259, 46346022, 46466865, 46529664, 46529672, 47160734, 47184386, 47404077, 47404128, 47521369, 47664332, 48657588, 48793364, 48904814, 49033878, 49092003, 49111238, 49153584, 49239494, 49245431, 49245501, 49271577, 49281757, 49353979, 49354003, 49832440, 50083670, 50777315, 50844889, 51040853, 51461396, 51478310, 51962196, 52316901, 52440302, 52516438, 52516524, 52763618, 53457986, 54407933, 54410436, 54411554, 54411839, 54500743, 54958838, 55825895, 56556272, 57153353, 57292918, 57640068, 63531481, 78137975, 79028148, 97389558, 110964543, 123091877, 129303850, 133626747, 136657509, 137012910, 140441208, 145703714, 146649557, 152158399, 156564248, 159116444, 167787588, 168414570, 180267240, 180900845, 184140267, 204876492, 207014301, 210301680, 230306199, 240245997, 244510658, 244510662, 244961449, 248318484, 254213953, 266272560, 266459539, 272885718, 274356430, 279021958, 287807698, 301005153, 301024862, 302503926, 302820800, 302820833, 303379904, 303675416, 303728144, 304105702, 304755222, 305770282, 306143327, 306236123, 306379054, 306842493, 307052699, 307052805, 307491855, 307531574, 307531582, 307773877, 307907667, 307944489, 308063506, 308231941, 308334061, 309059445, 309126285, 309126948, 309264487, 309272277, 309572137, 309572797, 309572964, 309890770, 310247417, 310550126, 310794244, 310971341, 311189171, 311284369, 311808383, 311943476, 311979688, 311979705, 312609520, 312761182, 313137890, 313198711, 313200258, 313770990, 313964575, 314336216, 314693135, 314901512, 315115943, 315670094, 316038683, 316043204, 317062384, 317064854, 317225338, 317407330, 317731559, 317731656, 318322481, 319774829, 319774851, 322088264, 322939351, 330025563, 331015800, 333592488, 334508540, 337383836, 337905387, 337938021, 338602551, 338710680, 338857037, 338978551, 339171953, 339335434, 340639504, 340639570, 340654057, 340864086, 341483069, 341983142, 342009830, 342357408, 342437114, 343802182, 344687204, 345709610, 345918789, 345983916, 346871058, 346883790, 346938412, 347242879, 347475568, 347475656, 347648337, 348342378, 348466157, 350328002, 351586819, 352137657, 352993408, 353912420, 354951365, 355278136, 355278162, 355677709, 356507945, 356988685, 357944731, 358404435, 358404450, 358404467, 358509073, 358620990, 358796311, 360438908, 360438921, 361302738, 361466361, 361635441, 361708579, 361723191, 361809047, 362063808, 362336666, 362494504, 362620496, 362775626, 366736468, 368576988, 369738051, 369995560, 371807127, 373799306, 379912791, 379917728, 403654793, 404484678, 405970096, 417443431, 423399512, 425642180, 450925269, 453994150, 462828989, 465733511, 466726565, 468023343, 468234495, 468234507, 469066656, 470375870, 470503339, 470797398, 471114975, 472164667, 472581745, 484602187, 484760825, 484884664, 485290782, 485501809, 485688335, 488233689, 488320055, 488630338, 488818685, 489356628, 503042137, 503188743, 503737526, 503755058, 506042273, 512050179, 512881165, 516520875, 519422843, 520578724, 532487814, 532529494, 532644655, 533670346, 533801756, 535117335, 535282394, 535289851, 535756160, 535878394, 535896242, 536249102, 536249167, 536600194, 536600207, 536973335, 537350548, 538603199, 539617085, 540192855, 540551835, 540592723, 540667747, 540670432, 540712921, 540914885, 541141903, 541263393, 541263419, 541322877, 541729095, 541741174, 545762254, 545888499, 545892585, 546031772, 546044944, 546508579, 546799348, 547315110, 547323163, 547343782, 547419510, 547434794, 547577508, 547598408, 547599643, 547852316, 547958750, 548016972, 548300959, 548305999, 548775451, 548934470, 548956406, 548995077, 549012098, 549295002, 549298905, 549314919, 549371075, 549420673, 549707035, 549851181, 551170221, 551377975, 551436380, 551712446, 551790671, 551802189, 551910214, 551989192, 555565696, 555565700, 555820426, 555878760, 555939835, 557261883, 557295620, 557295935, 557349107, 557726939, 557886121, 557977404, 558031527, 558078298, 558253708, 558307360, 558316677, 558464613, 558512914, 558541135, 558591438, 558598761, 558688020, 559041523, 559114283, 559256562, 559504207, 559613286, 559780257, 559780854, 559780926, 559876073, 559976386, 560369127, 560437999, 561171206, 561704570, 564125293, 564314347, 564570772, 564594337, 564815307, 565121462, 565187582, 565198894, 565389355, 565395781, 565519207, 565665129, 565952956, 566064814, 566160916, 566184384, 566224664, 566625886, 567172345, 567477989, 567478043, 567520175, 567536979, 568003154, 568032856, 568344634, 568357199, 568357210, 568437473, 568459119, 568496490, 568496498, 568650751, 568721487, 569246480, 569961243, 570115371, 570115372, 570226695, 570227644, 571903597, 571933685, 572022462, 572022482, 572053943, 572072836, 572156296, 572680899, 574479469, 576698011, 576898024, 577185921, 577373186, 577377232, 578324269, 578329017, 578785369, 578921601, 578938360, 579406652, 579429644, 579465309, 580421588, 580430996, 580826761, 582175651, 582270895, 582654763, 582655804, 582666606, 582766944, 583451496, 583819106, 583846413, 583868666, 584436785, 585564953, 585645875, 586062476, 586062558, 586095622, 586095643, 586328730, 586597796, 586602130, 586603730, 587142504, 587348951, 587392364, 587703398, 587878836, 587898700, 587994722, 588007962, 588022490, 588333792, 588382900, 588814081, 588882269, 589995287, 590064959, 590237780, 590389811, 591852344, 591970676, 592283151, 592424045]

    res = impute(observations, positions)

