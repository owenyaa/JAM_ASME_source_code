The program in this folder is to obtain the periodic solution of airfoil-store system with cubic nonlinearity through the IHB medhod and TIHB method (example 2 in the paper).

This folder mainly contains the following files,
Main file: main_program.m
The parameters of the algorithm can be modified and set in this document, including the order of series solution retained, the initial value of iteration, etc.

The subroutines related to the ERSA algorithm are as follows:csvd.m; l_corner.m; l_curve.m; lcfun.m; plot_lc.m; tikhonov.m.
These files are subroutines that iteratively solve the optimal solution of the objective function through the ERSA algorithm.