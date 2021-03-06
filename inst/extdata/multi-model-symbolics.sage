## Copyright 2018 Jason Vander Woude
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

    # F1 probabilities. Must be heterozygous. Row is binary representation of
    # two sites of homolog 1. Col is binary representation of two sites of
    # homolog 2. This says that the only two options for the F1 generation are
    # in a 2-site genome is that homolog 1 is 00 and homolog 2 is 11 or that
    # homolog 1 is 11 and homolog 2 is 00

    # probs = [[var("s%d" %(4*i + j)) for j in range(4)] for i in range(4)]

    # Initial symbolic probabilities in the F1 generation. Symbolics allows
    # setting the F1 progeny probabilities later. While in theory we know what
    # the F1 progeny should look like, it is possible that the two parents are
    # not truely homozygous within and polymorphic between at all sites
    probs = [[var("s%d%d" %(i, j)) for j in range(4)] for i in range(4)]

    # probs = [[0,   0, 0, 1/2],
    #          [0,   0, 0, 0  ],
    #          [0,   0, 0, 0  ],
    #          [1/2, 0, 0, 0  ]]

    for gen in range(1, generation):
        probs = _next_probs_helper(probs)

    return probs


def breed_parents(p1h1, p1h2, p2h1, p2h2):
    # Set up the symbolic variables. 'ps' is the probability of each of the four
    # gamete types for from parent 1 and 'qs' is the probability of each of the
    # four gamete types from parent 2. When selfing, parent 1 and parent 2 are
    # identical.
    #    ps[0]: probability that gamete 1 is identical to homolog 1 in parent 1
    #    ps[1]: probability that the first allele of gamete 1 is from homolog 1
    #           in parent 1 and the second allele is from homolog 2 of parent 1
    #    ps[2]: probability that the first allele of gamete 1 is from homolog 2
    #           in parent 1 and the second allele is from homolog 1 of parent 1
    #    ps[3]: probability that gamete 1 is identical to homolog 2 in parent 1
    #
    #    qs   : same idea but gamete 2 from parent 2
    r = var("r")

    # set up the probabilities of recombination
    rs = [1/2 * (1-r),
          1/2 * r,
          1/2 * r,
          1/2 * (1-r)]

    # initiate the next progeny probabilities
    progeny_probs = [[0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0]]

    # define possible gametes that can be produced from the progeny
    # using bitwise selection of the alleles
    p1gametes = [p1h1,                        # associated with ps[0] and qs[0]
                (p1h1 & 0b10) + (p1h2 & 0b01),# associated with ps[1] and qs[1]
                (p1h2 & 0b10) + (p1h1 & 0b01),# associated with ps[2] and qs[2]
                 p1h2]                        # associated with ps[3] and qs[3]

    p2gametes = [p2h1,                        # associated with ps[0] and qs[0]
                (p2h1 & 0b10) + (p2h2 & 0b01),# associated with ps[1] and qs[1]
                (p2h2 & 0b10) + (p2h1 & 0b01),# associated with ps[2] and qs[2]
                 p2h2]                        # associated with ps[3] and qs[3]

    # for each possible gamete from parent 1
    for i in range(4):
        g1 = p1gametes[i]  # new gamete 1
        # for each possible gamete from parent 2
        for j in range(4):
            g2 = p2gametes[j]  # new gamete 2
            progeny_probs[g1][g2] += rs[i]*rs[j] / 2
            progeny_probs[g2][g1] += rs[i]*rs[j] / 2

    return progeny_probs


# any of the last four arguments should be between 0 and 3 inclusinve where 0b00 is homozygous reference, 0b01 is heterozygous type I, 0b10 is heterozygous type II, and 0b11 is homozygous alternate
def sub_w_parents(data, p1h1, p1h2, p2h1, p2h2):
    init_progeny_probs = breed_parents(p1h1, p1h2, p2h1, p2h2)
    d = {var("s%d%d" %(i, j)): init_progeny_probs[i][j] for j in range(4) for i in range(4)}

    return rec_apply(lambda x: x.subs(d).expand(), data)


def _next_probs_helper(progeny_probs):

    # Set up the symbolic variables. 'ps' is the probability of each of the four
    # gamete types for from parent 1 and 'qs' is the probability of each of the
    # four gamete types from parent 2. When selfing, parent 1 and parent 2 are
    # identical.
    #    ps[0]: probability that gamete 1 is identical to homolog 1 in parent 1
    #    ps[1]: probability that the first allele of gamete 1 is from homolog 1
    #           in parent 1 and the second allele is from homolog 2 of parent 1
    #    ps[2]: probability that the first allele of gamete 1 is from homolog 2
    #           in parent 1 and the second allele is from homolog 1 of parent 1
    #    ps[3]: probability that gamete 1 is identical to homolog 2 in parent 1
    #
    #    qs   : same idea but gamete 2 from parent 2
    r = var("r")

    # set up the probabilities of recombination
    rs = [1/2 * (1-r),
          1/2 * r,
          1/2 * r,
          1/2 * (1-r)]

    # initiate the next progeny probabilities
    new_probs = [[0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0],
                 [0,0,0,0]]

    # Next, iterate over every possible 2-site genetic makeup, check the
    # probability that it occurs in the current generation, then check what
    # progeny can be produced from such an individual

    # Each genetic individual can be identified by the two homologs where each
    # homolog is encoded as a value in 0..3 inclusive. The binary representation
    # of the homolog indicates if the allele is ancestral reference (0) or
    # ancestral alternate (1).

    # for each first possible homolog
    for h1 in range(4):
        # for each second possible homolog
        for h2 in range(4):
            # probability of the progeny defined by the two homologs occurring
            # in the population
            prev_prob = progeny_probs[h1][h2]

            # define possible gametes that can be produced from the progeny
            # using bitwise selection of the alleles
            gametes = [h1,                      # associated with ps[0] and qs[0]
                      (h1 & 0b10) + (h2 & 0b01),# associated with ps[1] and qs[1]
                      (h2 & 0b10) + (h1 & 0b01),# associated with ps[2] and qs[2]
                       h2]                      # associated with ps[3] and qs[3]

            # for each possible gamete formed under p probabilities
            for i in range(4):
                g1 = gametes[i]     # new gamete 1
                # for each possible gamete formed under q probabilities
                for j in range(4):
                    g2 = gametes[j] # new gamete 2
                    new_probs[g1][g2] += rs[i]*rs[j]*prev_prob

    return new_probs


# The from_to_matrix will be such that if the current state is i, the
# probability of transitioning to state j is from_to_matrix[i][j]

# i and j are not longer binary encodings of homologs, but are each binary
# encodings of the two alleles at a given site where order of alleles at a given
# site matters (because the first site is homolog 1 and the second site is
# homolog 2). i is for the first of the two sites, and j is for the second of
# the two sites.
def get_from_to_matrix(generation, prog_probs):
    # ftmat is from_to_matrix
    ftmat = [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]]

    # all cases where the first state is 00
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


def parse_common_F1_probs_to_R(generation):
    print("## Copyright 2018 Jason Vander Woude"                               )
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
    print("## interpreter and run the function 'parse_common_F1_probs_to_R' with the single"    )
    print("## parameter being the generation of the plants where a value of 2 corresponds"      )
    print("## to F2, 3 corresponds to F3, etc. Running the function will print the required"    )
    print("## text for the R file to the console where it must be copied and pasted into a"     )
    print("## file called 'F{generation}.R' and saved in the appropriate directory."            )
    print(""                                                                                    )
    print("## To download SAGE, visit https://www.sagemath.org/download.html"                   )
    print(""                                                                                    )
    print(""                                                                                    )
    print("## There is a seperate function for each file because at runtime, only one"          )
    print("## function will be needed, so there is no need to source to code for all other"     )
    print("## generations"                                                                      )
    print(""                                                                                    )

    def print_vector(v):
        print("        c(")
        last_elem = len(v) - 1
        for i, elem in enumerate(v):
            trailing = "," if not i==last_elem else ""
            print("            function(r){" + str(elem) + "}" + trailing)
        print("        )")

    def print_matrix(m):
        print("            matrix(c(")

        last_row = len(m) - 1
        for j, row in enumerate(m):
            trailing = "," if not j==last_row else ""

            print("                " + ", ".join([str(x) for x in row]) + trailing)

        print("            ), nrow=" +
              str(len(m)) +
              ", ncol=" +
              str(len(m[0])) +
              ", byrow=T)")

    # Use one fewer generation because we are considering that the F1 is common
    # and not the parents
    site_pair_trans_probs = get_site_pair_trans_probs(generation - 1)

    print("site.pair.transition.probs <- list(")
    for marker1 in range(4):
        print("    list(")
        for marker2 in range(4):
            # p1h1 = ((marker1 & 0b1000) >> 2) + ((marker2 & 0b1000) >> 3)
            # p1h2 = ((marker1 & 0b0100) >> 1) + ((marker2 & 0b0100) >> 2)
            ph1 = ((marker1 & 0b0010) << 0) + ((marker2 & 0b0010) >> 1)
            ph2 = ((marker1 & 0b0001) << 1) + ((marker2 & 0b0001) >> 0)


            print("        function(r) {")
            print_matrix(
                sub_w_parents(
                    site_pair_trans_probs,
                    ph1, ph2, ph1, ph2))

            if marker2 != 3:
                print("        },\n")
            else:
                print("        }")


        if marker1 != 3:
            print("    ),\n")
        else:
            print("    )")
    print(")")


def parse_probs_to_R(generation):
    print("## Copyright 2018 Jason Vander Woude"                                                )
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
    print("## interpreter and run the function 'parse_probs_to_R' with the single"              )
    print("## parameter being the generation of the plants where a value of 2 corresponds"      )
    print("## to F2, 3 corresponds to F3, etc. Running the function will print the required"    )
    print("## text for the R file to the console where it must be copied and pasted into a"     )
    print("## file called 'F{generation}.R' and saved in the appropriate directory."            )
    print(""                                                                                    )
    print("## To download SAGE, visit https://www.sagemath.org/download.html"                   )
    print(""                                                                                    )
    print(""                                                                                    )
    print("## There is a seperate function for each file because at runtime, only one"          )
    print("## function will be needed, so there is no need to source to code for all other"     )
    print("## generations"                                                                      )
    print(""                                                                                    )

    def print_vector(v):
        print("        c(")
        last_elem = len(v) - 1
        for i, elem in enumerate(v):
            trailing = "," if not i==last_elem else ""
            print("            function(r){" + str(elem) + "}" + trailing)
        print("        )")

    def print_matrix(m):
        print("            matrix(c(")

        last_row = len(m) - 1
        for j, row in enumerate(m):
            trailing = "," if not j==last_row else ""

            print("                " + ", ".join([str(x) for x in row]) + trailing)

        print("            ), nrow=" +
              str(len(m)) +
              ", ncol=" +
              str(len(m[0])) +
              ", byrow=T)")

    site_pair_trans_probs = get_site_pair_trans_probs(generation)

    print("site.pair.transition.probs <- list(")
    last_list = 15  # 16 - 1
    for marker1 in range(16):
        print("    list(")
        for marker2 in range(16):
            p1h1 = ((marker1 & 0b1000) >> 2) + ((marker2 & 0b1000) >> 3)
            p1h2 = ((marker1 & 0b0100) >> 1) + ((marker2 & 0b0100) >> 2)
            p2h1 = ((marker1 & 0b0010) << 0) + ((marker2 & 0b0010) >> 1)
            p2h2 = ((marker1 & 0b0001) << 1) + ((marker2 & 0b0001) >> 0)


            print("        function(r) {")
            print_matrix(
                sub_w_parents(
                    site_pair_trans_probs,
                    p1h1, p1h2, p2h1, p2h2))

            if marker2 != 15:
                print("        },\n")
            else:
                print("        }")


        if marker1 != 15:
            print("    ),\n")
        else:
            print("    )")
    print(")")


def rec_apply(fun, data):
    try:
        return [rec_apply(fun, x) for x in data]
    except TypeError:
        return fun(data)


# def verify(generation):
#     cond = condensed_model_probs(generation)

#     # AA = BB
#     assert bool(cond[0][0] == cond[3][3]), "cond[0][0] == cond[3][3]"

#     # AH = HA = BH = HB
#     assert bool(cond[0][1] == cond[0][2]), "cond[0][1] == cond[0][2]"
#     assert bool(cond[0][1] == cond[1][0]), "cond[0][1] == cond[1][0]"
#     assert bool(cond[0][1] == cond[2][0]), "cond[0][1] == cond[2][0]"
#     assert bool(cond[0][1] == cond[1][3]), "cond[0][1] == cond[1][3]"
#     assert bool(cond[0][1] == cond[2][3]), "cond[0][1] == cond[2][3]"
#     assert bool(cond[0][1] == cond[3][1]), "cond[0][1] == cond[3][1]"
#     assert bool(cond[0][1] == cond[3][2]), "cond[0][1] == cond[3][2]"

#     # HH = HH
#     assert bool(cond[1][1] == cond[2][2]), "cond[1][1] == cond[2][2]"
#     assert bool(cond[1][2] == cond[2][1]), "cond[1][1] == cond[1][2]"

#     # AB = BA
#     assert bool(cond[0][3] == cond[3][0]), "cond[0][3] == cond[3][0]"


def get_generalized_trans_probs(generation):

    def prog_probs_to_general_transition_probs(p):
        return [
            # A to A and B to B
            p[0][0] + p[3][3],
            # H to H
            p[0][3] + p[1][2] + p[2][1] + p[3][0],
            # A to H, H to A, B to H, and H to B
            p[0][1] + p[0][2] + p[1][0] + p[1][3] + p[2][0] + p[2][3] + p[3][1] + p[3][2],
            # A to B and B to A
            p[1][1] + p[2][2]
        ]

    def print_vector(v):
        print("        c(")
        last_elem = len(v) - 1
        for i, elem in enumerate(v):
            trailing = "," if not i==last_elem else ""
            print("            function(r){" + str(elem) + "}" + trailing)
        print("        )")

    # result = rec_apply(expand,
    #     prog_probs_to_general_transition_probs(get_progeny_probs(generation)))

    return prog_probs_to_general_transition_probs(get_progeny_probs(generation))



def get_site_pair_trans_probs(generation):

    def prog_probs_to_4_state_transition_probs(p):
        # from row to column
        return [
            #  P1          H1          H2          P2
            [p[0][0],     p[0][1],     p[1][0],     p[1][1]],
            [p[0][2],     p[0][3],     p[1][2],     p[1][3]],
            [p[2][0],     p[2][1],     p[3][0],     p[3][1]],
            [p[2][2],     p[2][3],     p[3][2],     p[3][3]]
        ]

    return prog_probs_to_4_state_transition_probs(get_progeny_probs(generation))


# compute the probability of an odd number of recombination between distant
# sites given the probability of an odd number of recombinations between every
# pair of sequential sites betweeen
def compute_distant_recomb(n):
    rs = [var("r%d" % i) for i in range(1,n+1)]
    expr = 0

    ## for each possible odd-even recombination configurations
    for i in range(2^n):
        ## the binary representation of i will represent which of the n
        ## recombination variables should be odd
        if (bin(i).count("1") % 2 == 1):
            ## There must be an odd number of recombination variables that allow
            ## an odd number of recombinations
            term = 1
            for bit in range(n):
                ## if bit is 1 (i.e. assume odd num recombs for this recomb var)
                if (2**bit & i):
                    term *= rs[bit]
                else:
                    term *= (1 - rs[bit])
            expr += term

    return expr.expand()


# # TODO
# def initial_probs():

#     # Set up the symbolic variables. 'ps' is the probability of each of the four
#     # gamete types for from parent 1 and 'qs' is the probability of each of the
#     # four gamete types from parent 2. When selfing, parent 1 and parent 2 are
#     # identical.
#     #    ps[0]: probability that gamete 1 is identical to homolog 1 in parent 1
#     #    ps[1]: probability that the first allele of gamete 1 is from homolog 1
#     #           in parent 1 and the second allele is from homolog 2 of parent 1
#     #    ps[2]: probability that the first allele of gamete 1 is from homolog 2
#     #           in parent 1 and the second allele is from homolog 1 of parent 1
#     #    ps[3]: probability that gamete 1 is identical to homolog 2 in parent 1
#     #
#     #    qs   : same idea but gamete 2 from parent 2
#     r = var("r")

#     # set up the probabilities of recombination
#     rs = [1/2 * (1-r),
#           1/2 * r,
#           1/2 * r,
#           1/2 * (1-r)]

#     # initiate the next progeny probabilities
#     init_probs = [[0,0,0,0],
#                   [0,0,0,0],
#                   [0,0,0,0],
#                   [0,0,0,0]]

#     # Next, iterate over every possible 2-site genetic makeup, check the
#     # probability that it occurs in the current generation, then check what
#     # progeny can be produced from such an individual

#     # Each genetic individual can be identified by the two homologs where each
#     # homolog is encoded as a value in 0..3 inclusive. The binary representation
#     # of the homolog indicates if the allele is ancestral reference (0) or
#     # ancestral alternate (1).

#     # for each first possible homolog of parent 1
#     for h1 in range(4):
#         # for each second possible homolog of parent 1
#         for h2 in range(4):
#             # for each first possible homolog of parent 2
#             for h3 in range(4):
#                 # for each second possible homolog of parent 2
#                 for h4 in range(4):

#                     # define possible gametes that can be produced from the progeny
#                             # using bitwise selection of the alleles
#                     parent_1_gametes = [h1,                      # associated with ps[0] and qs[0]
#                                         (h1 & 0b10) + (h2 & 0b01),# associated with ps[1] and qs[1]
#                                         (h2 & 0b10) + (h1 & 0b01),# associated with ps[2] and qs[2]
#                                         h2]                      # associated with ps[3] and qs[3]

#                     # define possible gametes that can be produced from the progeny
#                     # using bitwise selection of the alleles
#                     parent_2_gametes = [h3,                      # associated with ps[0] and qs[0]
#                                         (h3 & 0b10) + (h4 & 0b01),# associated with ps[1] and qs[1]
#                                         (h4 & 0b10) + (h3 & 0b01),# associated with ps[2] and qs[2]
#                                         h4]                      # associated with ps[3] and qs[3]

#                     # for each possible parent 1 gamete
#                     for i in range(4):
#                         g1 = parent_1_gametes[i]     # new gamete 1
#                         # for each possible parent 2 gamete
#                         for j in range(4):
#                             g2 = parent_2_gametes[j] # new gamete 2
#                             init_probs[g1][g2] += rs[i]*rs[j]

#                             return init_probs



def print_vector(v):
    print("        c(")
    last_elem = len(v) - 1
    for i, elem in enumerate(v):
        trailing = "," if not i==last_elem else ""
        print("            function(r){" + str(elem) + "}" + trailing)
    print("        )")
