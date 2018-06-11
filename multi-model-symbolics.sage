## Copyright 2017 Jason Vander Woude and Nathan Ryder
##
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
##
##     http://www.apache.org/licenses/LICENSE-2.0
##
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.


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


# def eval_mat(mat, ps, qs):
#     return [[entry.subs(p_1 = ps[0],
#                         p_2 = ps[1],
#                         p_3 = ps[2],
#                         p_4 = ps[3],
#                         q_1 = qs[0],
#                         q_2 = qs[1],
#                         q_3 = qs[2],
#                         q_4 = qs[3]) for entry in row] for row in mat]


# def eval_mat(mat, ps, qs):
#     return [[entry.subs(p_1 = ps[0],
#                         p_2 = ps[1],
#                         q_1 = qs[0],
#                         q_2 = qs[1]) for entry in row] for row in mat]


# var('r')
# [[entry.subs(p_1 = r,
#              p_2 = r,
#              p_3 = r,
#              p_4 = r,
#              q_1 = r,
#              q_2 = r,
#              q_3 = r,
#              q_4 = r).expand() for entry in row] for row in mat]



def get_all_model_probs(generation):
    r = var("r")
    ps = [var("p_%d" % gen) for gen in range(1,generation)]
    qs = [var("q_%d" % gen) for gen in range(1, generation)]

    all_vars = ps + qs  # array concatenation

    n.selfs = generation - 1  # number of selfing events
    n.models = 2^len(all_vars)

    ftmat = get_from_to_matrix(generation)

    def probs(model, all_vars, ftmat):
        dictionary = {all_vars[i]: r if (model & 2^i) else 0 for i in range(len(all_vars))}
        return [[entry.subs(dictionary).expand() for entry in row] for row in ftmat]

    return [probs(model, all_vars, ftmat) for model in range(n.models)]


def parse_all_model_probs_to_R(generation):
    '''This function will print to the console the code for an R file so that the
    transition probabilities for an F{generation} can be calculated under each
    of the recombination models.
    '''

    if (generation < 2):
        exit("can only generate model probs for F2 and beyond")

    model_probs_list = get_all_model_probs(generation)
    last_model = len(model_probs_list) - 1

    print("## Copyright 2017 Jason Vander Woude and Nathan Ryder"                               )
    print("##"                                                                                  )
    print("## Licensed under the Apache License, Version 2.0 (the \"License\");"                )
    print("## you may not use this file except in compliance with the License."                 )
    print("## You may obtain a copy of the License at"                                          )
    print("##"                                                                                  )
    print("##     http://www.apache.org/licenses/LICENSE-2.0"                                   )
    print("##"                                                                                  )
    print("## Unless required by applicable law or agreed to in writing, software"              )
    print("## distributed under the License is distributed on an \"AS IS\" BASIS,"              )
    print("## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied."         )
    print("## See the License for the specific language governing permissions and"              )
    print("## limitations under the License."                                                   )
    print(""                                                                                    )
    print(""                                                                                    )
    print("## This file is auto-generated by multi-model-symbolics.sage for F" + str(generation))
    print("## generation plants. To generate a similar file for any generation of"              )
    print("## recombinant inbred line (RIL) plants, load this SAGE file in a SAGE"              )
    print("## interpreter and run the function 'parse_all_model_probs_to_R' with the single"    )
    print("## parameter being the generation of the plants where a value of 2 corresponds"      )
    print("## to F2, 3 corresponds to F3, etc. Running the function will print the required"    )
    print("## text for the R file to the console where it must be copied and pasted into a"     )
    print("## file called 'F{generation}.R' and saved in the appropriate directory."            )
    print(""                                                                                    )
    print("## To download SAGE, visit https://www.sagemath.org/download.html"                   )
    print(""                                                                                    )
    print(""                                                                                    )
    print("## This function will return a list of transition matrices with one list entry"      )
    print("## for each of the 4^(generation - 1) recomination models."                          )
    print(""                                                                                    )
    print("get.all.model.trans.F" + str(generation) + " <- function() {"                        )
    print(""                                                                                    )
    print("    list(")

    for i, model_probs in enumerate(model_probs_list):
        list_sep = "," if not i==last_model else ""

        print("        function (r) {")
        print("            matrix(c(")

        last_row = len(model_probs) - 1
        for j, row in enumerate(model_probs):
            trailing = "," if not j==last_row else ""

            print("                " + ", ".join([str(x) for x in row]) + trailing)

        print("            ), nrow=4, ncol=4, byrow=T)}" + list_sep + "\n")

    print("    )")
    print("}")
    # print(""                                                                                   )
    # print(""                                                                                   )
    # print("get.trans.mat.F" + str(generation) + " <- function() {"                             )
    # print("    trans.mat.F" + str(generation)                                                  )
    # print("}"                                                                                  )
