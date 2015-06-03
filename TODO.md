##TODO list

###Architectual improvements
- Add exception handling
- Consider adding more flexibility to the B-spline constraint. For example: allowing inequality constraints and selection of relaxation method.
- Extend ConstraintSet to support constraint handlers (use pointers?) and removal of constraints
- With constraint handler support in ConstraintComposite: add support for disjunctive programming in BB (need some way/structure to specify which constraints to remove)

###Algorithmic improvements
- Consider implementing strong branching and pseudo-cost branching: started.
- Look into other knot refinement methods: currently have global and local refinement.
- Add branching variable priority. Important for B-splines in high dimensions since they need more branching. Alternative: let f(x) be a constraint and g(x) be its convex relaxation. Use |f(x) - g(x)| to measure error and prioritize branching using this error.
- Add detection of near-linearity in a subset of the domain of constraint functions: that is, g(x) + A(y-bar(y)) ~= f(x,y), where g(x) = f(x,bar(y)). Use norm(g(x) + A(y-bar(y)) - f(x,y)) to see if A and b are reasonable approximations.

###Speed improvements
- Refactor BranchAndBound and utilize parallell computing
- InterfaceIpopt: point to parent optimizer instead of constraints and objective - make the optimizer responsible for constraints and objective. The constraints and objective may then be updated in BB instead of creating a new optimizer for each class (removes initialization time of Ipopt).
- InterfaceIpopt and BB: implement warm-start using previous lagrange variables.

###Testing and bug-fixing
- Add unit testing framework
- Test bounds tightening procedures (especially for integer variables)
