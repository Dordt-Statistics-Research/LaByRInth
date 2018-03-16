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



## This file contains a number of functions that are used throughout the rest of
## the code



ensure_writability <- function(filename) {
  # The '/' is literal.
  # The '[^/]*' will match any character other than '/' any number of time.
  # The '$' only matches the end of the line
  # Basically, this strips any trailing characters and the last '/' to get
  # the directory of the file
  dir <- gsub("/[^/]*$", "", filename)
  if (!dir.exists(dir))
    dir.create(dir, recursive=T)
}


# Function that will call the standard writeLines
# function after first ensuring that the directory
# being saved to exists. If the directory does not
# exist it is created.
force_write_lines <- function(data, filename) {
  ensure_writability(filename)
  writeLines(data, filename)
}


str.to.num <- function(str, sep) {
    as.numeric(str.split(str, sep))
}


## More convenient strsplit if length of vector is 1
str.split <- function(str, sep) {
    if (length(str) != 1) {
        warning("Only splitting the first element of the vector")
    }
    strsplit(str, sep)[[1]]
}


vec.or <- function(vec) {
    Reduce(`|`, vec)
}


## Used for creating a 3-D structure
reorder <- function(x, n) {
    as.vector(sapply(1:n, function(i) {
        v <- rep(FALSE, n)
        v[i] <- TRUE
        x[v]
    }))
}


## A way to check if a string contains only numeric characters. Code comes from
## stackoverflow.com/questions/13301437/how-to-check-if-the-value-is-numeric
check.int <- function(vec){
    sapply(vec, function(N) {
        !length(grep("[^[:digit:]]", as.character(N)))})
}
