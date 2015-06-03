/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "branchercontinuousmostpromising.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

// Implements MOST_PROMISING continuous branching rule
// Selects the continuous variable with the best record in increasing the lower bound
bool BrancherContinuousMostPromising::selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const
{
    // Select branching variable based on branching costs
    std::unordered_map<int,double> pCosts = node->getPseudoCostsContinuous();

    bool foundBranchingVariable = false;
    double mostPromising = 0.0;

    for (const auto &it: pCosts)
    {
        auto var = node->getConstraints()->getVariableAt(it.first);

        if (var->getType() != VariableType::CONTINUOUS)
            continue;

        if (var->getUpperBound() - var->getLowerBound() <= branchingThreshold)
            continue;

        if (it.second > mostPromising)
        {
            mostPromising = it.second;
            branchingVariable = var;
            foundBranchingVariable = true;
        }
    }

    if (foundBranchingVariable)
    {
        cout << "The most promising's cost is: " << mostPromising << endl;
        return true;
    }

    // Found none
    return false;
}

} // namespace BB

} // namespace CENSO
