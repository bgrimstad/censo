/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "brancherintegermostfractional.h"

namespace CENSO
{

namespace BB
{

// Implements MOST_FRACTIONAL integer branching rule
bool BrancherIntegerMostFractional::selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const
{
    // Get branching variables and variable types
    auto branchingVariables = node->getBranchingVariables();

    // Search for branching variable
    bool foundBranchingVariable = false;
    double mostFractional = 0.0;

    for (auto var : branchingVariables)
    {
        if (var->getType() != VariableType::INTEGER
              && var->getType() != VariableType::BINARY)
            continue;

        if (var->getUpperBound() == var->getLowerBound())
            continue;

        // Get variable value at optimal point of lower bounding problem

        double zlb = node->getLowerBoundSolution().at(node->getVariableIndex(var));
        // NOTE: All nodes that are branched upon have a feasible lower bound
        //double zlb = node->getLowerBoundSolution().at(var);

        // Calculate fraction
        double fraction = 0.5 - std::abs(zlb - std::floor(zlb) - 0.5); // frac is [0, 0.5], where 0.5 is most fractional

        // This inequality must be non-strict, i.e. >=
        if (fraction >= mostFractional)
        {
            mostFractional = fraction;
            branchingVariable = var;
            foundBranchingVariable = true;
        }
    }

    if (foundBranchingVariable)
        return true;

    // Found none
    return false;
}

} // namespace BB

} // namespace CENSO
