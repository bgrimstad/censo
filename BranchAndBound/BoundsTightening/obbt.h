/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef OBBT_H
#define OBBT_H

#include "boundstightener.h"
#include "OptimizationProblem/constraint.h"

namespace CENSO
{

namespace BB
{

/*
 * Optimality-Based Bounds Tightening (FBBT)
 * Algorithm for bound tightening that may improve solution speed of non-convex problems with convex relaxations by:
 * 1. Removing branching variables in problem formulation
 * 2. Reducing number of branchings needed in the branch and bound tree
 *
 * Let S be the feasible set defined by a constraint set (including its variable bounds)
 * Let C be the feasible set defined by the variable bounds and linear, convex, and possible convex relaxations, of the constraints
 * Let H be the set described by the variable bounds (hypercube)
 * Then S <= C <= H, and C and H are convex
 *
 * The bound tightening procedure:
 * 1. Retrieve the constraints describing C from S (possible with convex relaxations)
 * 2. Copy the variable bounds to lb and ub (lower and upper variable bounds)
 * 3. for each variable i
 *      lb(i) = arg min i subject to C
 *      ub(i) = arg max i subject to C
 *    end for
 * 4. Set lb and ub as new varaible bounds of S
 *
 * NOTE: The variable bounds (box constraints) contains the convex set that again contains the feasible region.
 * In bound tightening the bounds are pushed inwards to the boundary of the convex set, but it does not
 * cut the convex set. Since the convex set is dependent on the variable bounds it may shrink together and away
 * from the tightened bounds. When bound tightening is run successively the variable bounds will approach
 * the boundary of the convex set.
 *
 */

class OBBT : public BoundsTightener
{

public:
    OBBT(double threshold, unsigned int maxIterations);

    virtual ~OBBT() {}

private:

    // Do parallel computing
    bool doParallelComputing;

    // Do aggressive FBBT
    bool doAggressiveBoundsTightening;

    // Tighten domain bounds of constraint
    bool doTightening(ConstraintPtr constraints) override;

    // Tighten domain bounds of constraint (sequential)
    bool tightenBoundsSequential(ConstraintPtr cs, std::vector<VariablePtr> variables);

    // Tighten domain bounds of constraint (parallel)
    void tightenBoundsParallel(ConstraintPtr cs, std::vector<VariablePtr> variables);

    // Tighten bounds of a variable subset (used by tightenBoundsParallel)
    void tightenVariableBounds(ConstraintPtr cs, std::vector<VariablePtr> variables);

    // Tighten bounds of one variable
    bool tightenVariableBound(ConstraintPtr cs, VariablePtr variable);

//    // Tighten domain bounds of constraint (sequential)
//    void tightenBoundsAggressive(ConstraintPtr cs,
//                                 std::vector<double> z0,
//                                 std::vector<int> variables,
//                                 std::vector<double> &lb,
//                                 std::vector<double> &ub);

};

} // namespace BB

} // namespace CENSO

#endif // OBBT_H
