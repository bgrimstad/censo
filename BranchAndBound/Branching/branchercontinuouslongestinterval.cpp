/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "branchercontinuouslongestinterval.h"

namespace CENSO
{

namespace BB
{

// Implements the longest (bound) interval branching rule for continuous variables
bool BrancherContinuousLongestInterval::selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const
{
    // Get branching variables and variable types
    auto branchingVariables = node->getBranchingVariables();

    // Search for branching variable
    bool foundBranchingVariable = false;
    double longestInterval = 0;

    for (auto var : branchingVariables)
    {
        double interval = var->getUpperBound() - var->getLowerBound();

        if (var->getType() != VariableType::CONTINUOUS)
            continue;

        if (interval <= branchingThreshold)
            continue;

        if (interval > longestInterval)
        {
            longestInterval = interval;
            branchingVariable = var;
            foundBranchingVariable = true;
        }
    }

    return foundBranchingVariable;

    /*
     * Testing branching on infeasibility below
     */

//    // Get branching variables and variable types
//    std::vector<int> branchingVariables = node->getBranchingVariables();
//    std::vector<int> variableTypes = node->getVariableTypes();

//    // Get constraints
//    ConstraintPtr cs(node->getConstraints());

//    // Get domain bounds
//    std::vector<double> lbParent, ubParent;
//    node->getDomainBounds(lbParent, ubParent);

//    // Calculate longest interval
//    double longestInterval = 0;
//    for (unsigned int j = 0; j < branchingVariables.size(); j++)
//    {
//        int var = branchingVariables.at(j);
//        double interval = ubParent.at(var) - lbParent.at(var);
//        if (interval > longestInterval)
//            longestInterval = interval;
//    }
//    if (longestInterval == 0)
//        return false;

//    // Get infeasibility of lower bound solution
//    std::vector<double> xLower = node->getLowerBoundSolution();
//    DenseVector xLower2(xLower.size());
//    for (unsigned int i = 0; i < xLower.size(); i++)
//        xLower2(i) = xLower.at(i);

//    std::vector<double> infeas = cs->variableInfeasibility(xLower2);
//    //printVector(infeas);

//    // Search for branching variable
//    bool foundBranchingVariable = false;
//    double mostPromising = 0;

//    // Weight 0 <= alpha < 1
//    // alpha = 0 gives longest interval branching-variable selection rule
//    // alpha*infeasibilityMeasure + (1-alpha)*interval/longestInterval
//    double alpha = 0.0;

//    for (unsigned int j = 0; j < branchingVariables.size(); j++)
//    {
//        int var = branchingVariables.at(j);
//        double interval = ubParent.at(var) - lbParent.at(var);

//        if (variableTypes.at(var) != CONTINUOUS) continue;
//        if (interval <= branchingThreshold) continue;

//        double cost = alpha*infeas.at(var) + (1.0-alpha)*interval/longestInterval;

//        if (cost > mostPromising)
//        {
//            mostPromising = cost;
//            branchingVariable = var;
//            foundBranchingVariable = true;
//        }
//    }

//    return foundBranchingVariable;

}

} // namespace BB

} // namespace CENSO
