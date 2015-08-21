# CENSO
CENSO is a framework for global optimization of nonconvex, possibly spline-constrained, MINLP problems.

### Dependencies:
* [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) for linear algebra.
* [SPLINTER](https://github.com/bgrimstad/splinter) to compute with splines and other approximations.
* [GUROBI](http://www.gurobi.com/) for solving LP and MIP problems.
* [IPOPT](https://projects.coin-or.org/Ipopt) for solving NLP problems.
* [BONMIN](https://projects.coin-or.org/Bonmin/) for solving convex MINLP problems and for heuristically solving non-convex MINLP problems.

### Notes:
* The code examples in the user guide (Docs/UserGuide/manual.pdf) are currently outdated. Refer to the test problems (TestProblems/*) for examples on how to create and solve optimization problems.

### Reference:
If you use CENSO in a scientific work we kindly ask you to cite it. You can cite it as shown in the bibtex entry below (remember to update the date accessed).
```
@article{Grimstad2015,
author = {Grimstad, B. and Sandnes, A.},
journal = {To appear in Journal of Global Optimization},
title = {{Global optimization with spline constraints: a new branch-and-bound method based on B-splines}},
year = {2015}
}
```
