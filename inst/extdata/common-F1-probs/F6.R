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


## This file is auto-generated by multi-model-symbolics.sage for F6
## generation plants. To generate a similar file for any generation of
## recombinant inbred line (RIL) plants, load this SAGE file in a SAGE
## interpreter and run the function 'parse_common_F1_probs_to_R' with the single
## parameter being the generation of the plants where a value of 2 corresponds
## to F2, 3 corresponds to F3, etc. Running the function will print the required
## text for the R file to the console where it must be copied and pasted into a
## file called 'F{generation}.R' and saved in the appropriate directory.

## To download SAGE, visit https://www.sagemath.org/download.html


## There is a seperate function for each file because at runtime, only one
## function will be needed, so there is no need to source to code for all other
## generations

site.pair.transition.probs <- list(
    list(
        function(r) {
            matrix(c(
                1, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                31/64, 1/64, 1/64, 31/64,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                31/64, 1/64, 1/64, 31/64,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 1,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        }
    ),

    list(
        function(r) {
            matrix(c(
                31/64, 0, 0, 0,
                1/64, 0, 0, 0,
                1/64, 0, 0, 0,
                31/64, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 31/64,
                0, 0, 0, 1/64,
                0, 0, 0, 1/64,
                0, 0, 0, 31/64
            ), nrow=4, ncol=4, byrow=T)
        }
    ),

    list(
        function(r) {
            matrix(c(
                31/64, 0, 0, 0,
                1/64, 0, 0, 0,
                1/64, 0, 0, 0,
                31/64, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 35/16*r^4 - 5/8*r^3 + 5/64*r^2, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 55/16*r^4 - 15/8*r^3 + 45/64*r^2 - 5/32*r + 1/64, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r,
                1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 35/8*r^5 + 31/16*r^4 + 1/8*r^3 - 59/64*r^2 + 13/16*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, -1/4*r^10 + 5/4*r^9 - 25/8*r^8 + 5*r^7 - 45/8*r^6 + 37/8*r^5 - 45/16*r^4 + 5/4*r^3 - 25/64*r^2 + 5/64*r, 1/4*r^10 - 5/4*r^9 + 25/8*r^8 - 5*r^7 + 45/8*r^6 - 39/8*r^5 + 59/16*r^4 - 21/8*r^3 + 109/64*r^2 - 31/32*r + 31/64
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 31/64,
                0, 0, 0, 1/64,
                0, 0, 0, 1/64,
                0, 0, 0, 31/64
            ), nrow=4, ncol=4, byrow=T)
        }
    ),

    list(
        function(r) {
            matrix(c(
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                1, 0, 0, 0
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                31/64, 1/64, 1/64, 31/64
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                31/64, 1/64, 1/64, 31/64
            ), nrow=4, ncol=4, byrow=T)
        },

        function(r) {
            matrix(c(
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 0,
                0, 0, 0, 1
            ), nrow=4, ncol=4, byrow=T)
        }
    )
)
