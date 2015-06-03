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
