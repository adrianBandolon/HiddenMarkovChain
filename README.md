# Hidden Markov Chain Implementation in R

### Baum-Welch Algorithm

+ This is my implementation of the Baum-Welch Algorithm and Viterbi Algorithm.
+ The `forward.prob`, `backward.prob`, `myBW` and `BWonestep` functions were provided. I wrote the code that updates `myGamma` in the `BWonestep` function.
+ The structure of `mygamma` could be improved using matrix algebra instead of the 3 `for` loops used in this implementation

### Viterbi Algorithm

+ The Viterbi Algorithm was implemented in the `myViterbi` function. The output of `myViterbi` (`Coding3_HMM_Viterbi_Output.txt`) will be compared with the output from the `viterbi` function (`Coding3_HMM_True_Viterbi_Output.txt`) from the `HMM` library. `Coding3_HMM_Viterbi_Output.txt` will be written in the directory where this `*.rmd` file is executed.

+ For my implementation fo the Viterbi Algorithm, I modified the `viterbi` function from the `HMM` library to accomodate the data provided and the output from the Baum-Welch Algorithm above as inputs . Source code for the `HMM` library is **[here](https://github.com/cran/HMM/blob/bfbe14735b48325b1636a36081f798fdf00199e2/R/HMM.r#L39)**.

+ The inputs for the `myViterbi` function are `(x, A, B, w)`, where `x` is the data, `A` is a 2-by-2 matrix from the Baum-Welch Algorithm as implemented in the previous section, `B` is a 2-by-3 matrix also from the previous section, and `w` is the initial distribution matrix.

+ The code below compares `Coding3_HMM_Viterbi_Output.txt` against `Coding3_HMM_True_Viterbi_Output.txt`.
