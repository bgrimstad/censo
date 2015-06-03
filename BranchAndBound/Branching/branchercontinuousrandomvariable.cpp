/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "branchercontinuousrandomvariable.h"

namespace CENSO
{

namespace BB
{

// Implements the random variable branching rule for continuous variables
bool BrancherContinuousRandomVariable::selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const
{
    // Get branching variables and variable types
    auto branchingVariables = node->getBranchingVariables();

    // Search for branching variable
    std::vector<VariablePtr> candidates;

    for (auto var : branchingVariables)
    {
        if (var->getType() != VariableType::CONTINUOUS)
            continue;

        if (var->getUpperBound() - var->getLowerBound() <= branchingThreshold)
            continue;

        // Add candidate for branching
        candidates.push_back(var);
    }

    if (candidates.size() > 0)
    {
        int randomIndex = randomInteger(0, candidates.size()-1);
        branchingVariable = candidates.at(randomIndex);
        return true;
    }

    // Found none
    return false;
}

} // namespace BB

} // namespace CENSO
