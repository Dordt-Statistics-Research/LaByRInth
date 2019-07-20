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


##               __          ____        ____  _____
##              / /         / __ \      / __ \/_  _/
##             / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __
##            / /   / _  \/ _  _/\ \/ / _  _/ / / /  | / /_  __/ /_/ /
##           / /___/ /_/ / /_\ \  \  / / \ \_/ /_/ /||/ / / / / __  /
##          /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/
##
##                  L O W - C O V E R A G E   B I A L L E L I C
##                    R - P A C K A G E   I M P U T A T I O N
##


import pdb
import operator as op


# To add an additional breeding scheme, uncomment the lines below and follow the
# example. Further breeding schemes can be added by adding additional `elif` blocks

def get_progeny_probs(model, name_1, name_2):
    if (model in ["F1","F2","F3","F4","F5","F6","F7","F9","F10"]):
        # call a special helper function for RIL populations
        generation = int(model[1:len(model)])  # integer value after 'F'
        return get_symbolic_RIL_probs(generation, name_1, name_2)

    elif (model == "F1BC1"):
        p1 = init_symbolic_probabilistic_parent(name_1)
        p2 = init_symbolic_probabilistic_parent(name_2)

        F1 = breed_probabilistic_taxa(p1, p2)
        F1BC1 = breed_probabilistic_taxa(p1, F1)

        return F1BC1

    # UNCOMMENT THE LINES BELOW

    elif (model == "NEW_USER_DEFINED_MODEL"):
        # Initialize symbolic parents for the population
        p1 = init_symbolic_probabilistic_parent(name_1)
        p2 = init_symbolic_probabilistic_parent(name_2)

        # use `breed_probabilistic_taxa(taxa1, taxa2)` to breed two members

        # use `breed_probabilistic_self(taxa)` to self a member. This can also be
        # accomplished by calling `breed_probabilistic_taxa(taxa, taxa)`, but
        # using `breed_probabilistic_self` will be faster

        # for example, to represent an F3 population, do the following:
        F1 = breed_probabilistic_taxa(p1, p2)
        F2 = breed_probabilistic_self(F1)
        F3 = breed_probabilistic_self(F2)
        # etc.
        BC1 = breed_probabilistic_taxa(p1, F1)
        BC2 = breed_probabilistic_taxa(p1, BC1 )

    else:
        raise ValueError("Specified model not found")


def parse_probs_to_R(model):
    print('## Copyright 2018 Jason Vander Woude'                                            )
    print('##'                                                                              )
    print('## Licensed under the Apache License, Version 2.0 (the "License");'              )
    print('## you may not use this file except in compliance with the License.'             )
    print('## You may obtain a copy of the License at'                                      )
    print('##'                                                                              )
    print('##     http://www.apache.org/licenses/LICENSE-2.0'                               )
    print('##'                                                                              )
    print('## Unless required by applicable law or agreed to in writing, software'          )
    print('## distributed under the License is distributed on an "AS IS" BASIS,'            )
    print('## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.'     )
    print('## See the License for the specific language governing permissions and'          )
    print('## limitations under the License.'                                               )
    print(''                                                                                )
    print(''                                                                                )
    print('##               __          ____        ____  _____'                            )
    print('##              / /         / __ \      / __ \/_  _/'                            )
    print('##             / /   ____  / /_/ /_  __/ /_/ / / / __   __________  __'          )
    print('##            / /   / _  \/ _  _/\ \/ / _  _/ / / /  | / /_  __/ /_/ /'          )
    print('##           / /___/ /_/ / /_\ \  \  / / \ \_/ /_/ /||/ / / / / __  /'           )
    print('##          /_____/_/ /_/______/  /_/_/  /_/____/_/ |__/ /_/ /_/ /_/'            )
    print('##'                                                                              )
    print('##                  L O W - C O V E R A G E   B I A L L E L I C'                 )
    print('##                    R - P A C K A G E   I M P U T A T I O N'                   )
    print('##'                                                                              )
    print(''                                                                                )
    print(''                                                                                )
    print('################################################################################')
    print('##                                                                            ##')
    print('##            ATTENTION! THIS FILE IS AUTO-GENERATED. DO NOT EDIT             ##')
    print('##                                                                            ##')
    print('################################################################################')
    print('##                                                                            ##')
    print('## This file is auto-generated by multi-model-symbolics.sage for populations  ##')
    print('## bred under the following breeding scheme:                                  ##')
    print('##                                                                            ##')
    print('## ' +model + ' '*(80 - len(model) - 5) +                                    '##')
    print('##                                                                            ##')
    print('## To generate a similar file for other types of breeding schemes, edit the   ##')
    print('## multi-model-symbolics.sage file according to the comments in that          ##')
    print('## file. Then load multi-model-symbolics.sage in a SAGE interpreter and run   ##')
    print('## the function `parse_probs_to_R` with the model parameter being the string  ##')
    print('## that you use to define your breeding scheme. Running the function will     ##')
    print('## print the required text for the R file to the console where it must be     ##')
    print('## copied and pasted into a file called `{model}.R` where {model} is replaced ##')
    print('## with your breeding scheme model name. Save `{model}.R` in the appropriate  ##')
    print('## directory in the labyrinth package (/inst/extdata/transitions-probs/).     ##')
    print('##                                                                            ##')
    print('## To download SAGE, visit https://www.sagemath.org/download.html             ##')
    print('##                                                                            ##')
    print('## There is a seperate file for each model because at when LaByRInth runs,    ##')
    print('## only one model will be needed, so there is no need to source the code for  ##')
    print('## all other models.                                                          ##')
    print('##                                                                            ##')
    print('################################################################################')
    print(''                                                                                )

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


    ################################ MAIN LOGIC ################################


    # Get the progeny probabilities for this model once, with symbolic parent
    # information. The actual parental types will be filled in in the loops below
    name_1 = "parent_1"
    name_2 = "parent_2"
    progeny_probs = get_progeny_probs(model, name_1, name_2)

    print("site.pair.transition.probs <- list(")
    last_list = 15  # 16 list entries in 0-based indexing
    for marker1 in range(16): # =2^(2*2) = (#alleles)^(#parents * #homologs)
        print("    list(")
        for marker2 in range(16): # =2^(2*2) = (#alleles)^(#parents * #homologs)
            p1h1 = ((marker1 & 0b1000) >> 2) + ((marker2 & 0b1000) >> 3)
            p1h2 = ((marker1 & 0b0100) >> 1) + ((marker2 & 0b0100) >> 2)
            p2h1 = ((marker1 & 0b0010) << 0) + ((marker2 & 0b0010) >> 1)
            p2h2 = ((marker1 & 0b0001) << 1) + ((marker2 & 0b0001) >> 0)

            # convert the parent homologs into the required probabilistic form
            # which will give a probability of 1 in a single entry of the matrix
            p1_probs = convert_parent_to_probabilistic(p1h1, p1h2)
            p2_probs = convert_parent_to_probabilistic(p2h1, p2h2)

            with_parent_1 = numeric_instantiate_parent_in_progeny(
                progeny_probs, name_1, p1_probs)
            with_both_parents = numeric_instantiate_parent_in_progeny(
                with_parent_1, name_2, p2_probs)
            trans_probs = progeny_probs_to_4_state_transition_probs(with_both_parents)


            print("        function(r) {")
            print_matrix(trans_probs)


            ########################## END MAIN LOGIC ##########################


            if marker2 != 15:
                print("        },\n")
            else:
                print("        }")


        if marker1 != 15:
            print("    ),\n")
        else:
            print("    )")
    print(")")


def breed_deterministic_taxa(p1h1, p1h2, p2h1, p2h2):

    r = var("r")
    # The previous line instantiates a symbolic variable which applies
    # symmetrically to both parents because this is a parameter of the
    # species. An organism has 2 homologous chromosomes and we examine two
    # markers on this homologous pair.  Each of those two homologous chromosomes
    # has two sister chromatids. When a cell undergoes meiosis, four gametes
    # will be produced, and for a fixed gamete, the allele at the first marker
    # came from one of the four sister chromatids. Call this fixed gamete G and
    # the corresponding sister chromatid S. The paremeter r represents the
    # probability that there is some kind of crossover event between the two
    # markers such that allele at the second marker of gamete G belongs to a
    # chromatid which is not S or the sister of S. The parameter is specified
    # this way because if there is a crossover event such that the allele of G
    # at the second marker belongs to S or the sister of S, the type of allele
    # in G is the same as if no crossover event had occurred.


    # set up the probabilities of recombination. For a fixed organism, a
    # randomly produced gamete will have a 1/2 probability of having an allele
    # at the first marker from homolog 1 (either chromatid) and a 1/2
    # probability of having an allele at the first marker from homolog 2 (either
    # chromatid). In the array rs below the first two entries together represent
    # the prior, and the last two elements represent the latter. Then, fixing
    # one of these events, there is a probability of r that the allele at the
    # second marker comes from a chromatid of the opposite homolog.
    #
    #     Array element:    binary encoding of index (0 means allele is from
    #                                                 a chromatid of homolog 1
    #                                                 and 1 means allele is from
    #                                                 a chromatid of homolog 2)
    rs = [1/2 * (1-r),      # 00
          1/2 * r,          # 01
          1/2 * r,          # 10
          1/2 * (1-r)]      # 11

    # initiate the next progeny probabilities
    progeny_probs = [[0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0],
                     [0,0,0,0]]

    # define possible gametes that can be produced from the progeny using
    # bitwise selection of the alleles. The code "& 0b01" zeros the first/most
    # significant bit of the value on the left and similarly the code "& 0b10"
    # zeros the second/least significant bit of the value on the left.
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

            # The progeny that would be formed by these gametes would have one
            # homolog which is identically gamete 1 and the other homolog would
            # be identicaly gamete 2. However, because we need to consider fixed
            # homologs in an organism in order to perform the calculations of
            # this script (in the Hidden Markov Model of LaByRInth where these
            # computations are used, it is essential to keep track of the two
            # separate homologs), and it does not make sense to fix the ordering
            # in the progeny yet, (i.e. in the progeny we could consider gamete
            # 1 becoming homolog 1 and gamete 2 becoming homolog 2, but we could
            # also consider gamete 1 becoming homolog 2 and gamete 2 becoming
            # homolog 1) we apply half of the appropriate probability
            # symmetrically to both of these possibilities.
            progeny_probs[g1][g2] += rs[i]*rs[j] / 2
            progeny_probs[g2][g1] += rs[i]*rs[j] / 2

    # The return value is a 4x4 matrix of probabilities. The matrix is
    # indexed by genetic makeup of a potential progeny. More specifically, the
    # rows and columns are indexed by the genetic composition of a homolog. So
    # to find the probability that a progeny would have an reference allele at
    # both markers in the first homolog (represented by 00 which is the binary
    # representation of 0) and a second homolog with an alternate allele at the
    # first marker and a reference allele at the second (represented by 10 which
    # is the binary representation of 2), then use progeny_probs[0][2].
    return progeny_probs


def breed_probabilistic_self(progeny_probs):
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


def breed_probabilistic_taxa(taxon_1_probs, taxon_2_probs):

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

    for taxon_1_h1 in range(4):
        for taxon_1_h2 in range(4):
            taxon_1_prob = taxon_1_probs[taxon_1_h1][taxon_1_h2]

            for taxon_2_h1 in range(4):
                for taxon_2_h2 in range(4):

                    taxon_2_prob = taxon_2_probs[taxon_2_h1][taxon_2_h2]

                    taxon_1_and_2_progeny_probs = breed_deterministic_taxa(
                        taxon_1_h1, taxon_1_h2, taxon_2_h1, taxon_2_h2)

                    new_probs = add_matrices(
                        new_probs,
                        mult_matrix(taxon_1_and_2_progeny_probs,
                                    taxon_1_prob*taxon_2_prob)
                    )

    return new_probs


def rec_apply(fun, data):
    # Apply function fun to all non-iterable portions of data. For example, if
    # data is a matrix (an array of arrays) then fun will be applied to all the
    # elements of the matrix by recursively moving through the arrays.
    try:
        return [rec_apply(fun, x) for x in data]
    except TypeError:
        return fun(data)


def add_arrays(arr_1, arr_2):
    return map(sum, zip(arr_1, arr_2))

def add_matrices(m_1, m_2):
    return map(lambda matrices: add_arrays(
        matrices[0], matrices[1]), zip(m_1, m_2))

def mult_array(arr, const):
    return map(lambda x: const*x, arr)

def mult_matrix(m, const):
    return map(lambda arr: mult_array(arr, const), m)


def progeny_probs_to_4_state_transition_probs(p):
    # This is the matrix specifying the transition probabilities used in
    # LaByRInth. Index i,j of the returned matrix is the probability of
    # transitioning from state i at a marker to state j at the next marker
    # where states are specified by the allele type (reference or alternate)
    # in the first homolog and second homolog. Note that states of the
    # matrix are not dependent on the genetic makeup of the parents.

    return [
        # ref,ref      ref,alt      alt,ref      alt,alt
        [p[0][0],     p[0][1],     p[1][0],     p[1][1]],  # ref,ref
        [p[0][2],     p[0][3],     p[1][2],     p[1][3]],  # ref,alt
        [p[2][0],     p[2][1],     p[3][0],     p[3][1]],  # alt,ref
        [p[2][2],     p[2][3],     p[3][2],     p[3][3]]   # alt,alt
    ]


def print_vector(v):
    print("        c(")
    last_elem = len(v) - 1
    for i, elem in enumerate(v):
        trailing = "," if not i==last_elem else ""
        print("            function(r){" + str(elem) + "}" + trailing)
    print("        )")


def init_symbolic_probabilistic_parent(parent_name):
    # Returns a 4x4 matrix of unique symbolic variables based on
    # parent_name. This is used to create and breed parents prior to knowing the
    # deterministec genetic content of the parents. Once the genetics are known,
    # then the symbolic variables are replaced.

    return [[var("_"+parent_name+"_%d%d" %(i, j)) for j in range(4)] for i in range(4)]


def numeric_instantiate_parent_in_progeny(progeny_probs, parent_name, numeric_parent_probs):
    # Get all the symbolic variables that will be replaced
    symbolics_to_sub = init_symbolic_probabilistic_parent(parent_name)
    # create a dictionary with the symbolic variables and the substitution values
    d = {symbolics_to_sub[i][j]:numeric_parent_probs[i][j]
         for j in range(4) for i in range(4)}
    # make the replacements in progeny_probs
    return rec_apply(lambda x: x.subs(d).expand(), progeny_probs)


def get_symbolic_RIL_probs(generation, par_name_1, par_name_2):
    if (generation < 1):
        raise ValueError("Generation of RIL must be 1 or greater")
    p1 = init_symbolic_probabilistic_parent(par_name_1)
    p2 = init_symbolic_probabilistic_parent(par_name_2)
    F_i = breed_probabilistic_taxa(p1, p2)  # initially F1

    # Execute this loop (generation-1) times, so e.g. if generation=1 then the
    # loop will not run at all.
    for gen in range(1, generation):
        # given that we know what a member of gen looks like probabilistically,
        # we can self it probabilistically to obtain a probabilistic description
        # of the next generation.
        F_i = breed_probabilistic_self(F_i)

    return F_i


def convert_parent_to_probabilistic(homolog_1, homolog_2):
    probs = [[0,0,0,0],
             [0,0,0,0],
             [0,0,0,0],
             [0,0,0,0]]
    probs[homolog_1][homolog_2] = 1
    return probs
