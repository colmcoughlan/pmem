PMEM v1.03
Deconvolves dirty maps using the Maximum Entropy Method. Supports polarisation deconvolution.
Requires quickfits v1.1

1.03
Added in option to use a steepest descent-based solver. This is usually slower than the Newton-Raphson solver, but can be faster for the first few iterations. Also includes a bugfix for when the user specified a non-standard default map.

1.02
Fixed bugs related to quickfits v1.1
Turned off residual scaling by beam in code, may expose option to user later
Now uses proper masking on output
