# optimization_algorithms

Useful Optimization Algorithms

## L-BFGS

L-BFGS stands for “Limited memory BFGS”. It belongs to the family of quasi-Newton methods that approximates the [BFGS](https://en.wikipedia.org/wiki/Broyden%E2%80%93Fletcher%E2%80%93Goldfarb%E2%80%93Shanno_algorithm) algorithm using limited amount memory. It is particularly useful for parameter estimation in Machine Learning.  In case of big dimensions, the amount of memory required to store a Hessian (N2) is too big, along with the machine time required to process it, hence the need to use L-BFGS . L-BFGS stores only a few vectors that represent the approximation of the dense (N2) matrix implicitly. Due to its moderate memory requirement, L-BFGS method is particularly well suited for optimization problems with a large number of variables.

## Wolfeline Search

[Wolfe line search algorithm](https://en.wikipedia.org/wiki/Wolfe_conditions)
