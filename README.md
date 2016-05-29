PMEM v1.04
Deconvolves dirty maps using the Maximum Entropy Method. Supports polarisation deconvolution.
Requires quickfits v1.1

1.04

Added in two quasi-newton solvers : DFP and BFGS. Also added a conjugate gradient solver. The newton solver (option 0) remains the fastest, though the quasi newton solvers are also quite good. Reordered the header files and some of the code. Made sure multithreading was enabled everywhere. Added comments.

1.03
Added in option to use a steepest descent-based solver. This is usually slower than the Newton-Raphson solver, but can be faster for the first few iterations. Also includes a bugfix for when the user specified a non-standard default map.

1.02
Fixed bugs related to quickfits v1.1
Turned off residual scaling by beam in code, may expose option to user later
Now uses proper masking on output
