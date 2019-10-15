# Settlement violation estimates for proof-of-stake blockchains with the longest-chain rule
Compute the relative margin in time and space O(T^3) where T is the length of the characteristic string. 
See Appendix A of [the paper](https://eprint.iacr.org/2017/241.pdf) for reference. 

Naively computing these probabilities on a set of size 2^T would take time/space exponential in T.

## Author
Saad Quader and Alexander Russell

## Description
The code is in C++. The following files should be built as executables:
* `forkability_table.cpp` 
* `prob_forkable.cpp`

Look at the comments around the `main()` function in each file. They are fairly self-explanatory. 

The file `prob_forkable.cpp`, when built into an executable, is an interactive program which asks for
* `N`, a positive integer, the length of the characteristic string `w = w_1 ... w_N` 
* `R`, a positive integer, the maximum reach allowed in computation, and 
* `eps`, a real between 0 and 1, so that `Pr[w_i = 1] = (1 - eps)/2` independently for each `i = 1, ... , N`.
It outputs an upper bound on `Pr[w is forkable]`.

The file `forkability_table.cpp`, when built into an executable, is an interactive program which requests an output file name and then reproduces the table in Appendix A of [the paper](https://eprint.iacr.org/2017/241.pdf). Look in the `main()` function to see which values of `N`, `eps`, and `R` are being used.


If building on Windows, the header `win_mem_leak.h` helps detect memory leak, if any.
