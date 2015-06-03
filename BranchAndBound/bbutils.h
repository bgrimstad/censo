/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BBUTILS_H
#define BBUTILS_H

#include "OptimizationProblem/constraint.h"

namespace CENSO
{

namespace BB
{

// Check if problem is of type mixed-integer
bool isIntegerProblem(ConstraintPtr constraints);

// Check if a solution (vector) is integer feasible
bool isSolutionIntegerFeasible(std::vector<VariablePtr> variables, std::vector<double> z);

// Get branching variables (union of complicating and integer variables)
std::vector<VariablePtr> findBranchingVariables(ConstraintPtr constraints);

// Construct objective function cut
ConstraintPtr objectiveCut(ConstraintPtr constraints, double objectiveBound, bool upperCut = true);

} // namespace BB

} // namespace CENSO

#endif // BBUTILS_H
