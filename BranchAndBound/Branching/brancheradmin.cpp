/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "brancheradmin.h"

// Branchers
#include "brancherintegermostfractional.h"
#include "branchercontinuousmostpromising.h"
#include "branchercontinuouslongestinterval.h"
#include "branchercontinuousrandomvariable.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

BrancherAdmin::BrancherAdmin()
    : BrancherAdmin(MOST_FRACTIONAL, LONGEST_INTERVAL)//MOST_PROMISING, LONGEST_INTERVAL)
{
    // Default branchers
}

BrancherAdmin::BrancherAdmin(BranchingRuleInteger branchingRuleInteger,
                             BranchingRuleContinuous branchingRuleContinuous)
    : branchingRuleInteger(branchingRuleInteger),
      branchingRuleContinuous(branchingRuleContinuous)
{
    brancherInteger = BrancherPtr(createIntegerBrancher(branchingRuleInteger));
    brancherContinuous = BrancherPtr(createContinuousBrancher(branchingRuleContinuous));
}

Brancher* BrancherAdmin::createIntegerBrancher(BranchingRuleInteger branchingRule)
{
    switch (branchingRule)
    {
    case MOST_FRACTIONAL:
        return new BrancherIntegerMostFractional();
        break;

    case PSEUDOCOST:
        cout << "Warning: Pseudo cost branching is not implemented. Using most fractional branching." << endl;
        return new BrancherIntegerMostFractional();
        break;

    case STRONG_BRANCHING:
        cout << "Warning: Strong branching is not implemented. Using most fractional branching." << endl;
        return new BrancherIntegerMostFractional();
        break;

    default:
        // Default is MOST_FRACTIONAL
        return new BrancherIntegerMostFractional();
        break;
    }
}

Brancher* BrancherAdmin::createContinuousBrancher(BranchingRuleContinuous branchingRule)
{
    switch (branchingRule)
    {
    case MOST_PROMISING:
        return new BrancherContinuousMostPromising();
        break;

    case LONGEST_INTERVAL:
        return new BrancherContinuousLongestInterval();
        break;

    case RANDOM_VARIABLE:
        return new BrancherContinuousRandomVariable();
        break;

    default:
        // Default is MOST_PROMISING
        return new BrancherContinuousMostPromising();
        break;
    }
}

NodeList BrancherAdmin::branch(const NodePtr node) const
{
    NodeList children;
    VariablePtr branchingVariable = nullptr;
    bool foundBranchingVariable = false;

    // Prioritization of integer/binary variables over continuous variables
    foundBranchingVariable = brancherInteger->selectBranchingVariable(node, branchingVariable);

    if (!foundBranchingVariable)
    {
        // cout << "No integer branching variable found. Trying continuous variables." << endl;
        foundBranchingVariable = brancherContinuous->selectBranchingVariable(node, branchingVariable);
    }

    if (foundBranchingVariable)
    {
        if (branchingVariable->getType() == VariableType::BINARY)
        {
            children = brancherInteger->branchOnVariable(node, branchingVariable);
        }
        else if (branchingVariable->getType() == VariableType::INTEGER)
        {
            children = brancherInteger->branchOnVariable(node, branchingVariable);
        }
        else if (branchingVariable->getType() == VariableType::CONTINUOUS)
        {
            children = brancherContinuous->branchOnVariable(node, branchingVariable);
        }
        else
        {
            cout << "Variable type not recognized during branching! Exiting..." << endl;
            exit(1); // Missing elses are a scary thing
        }
    }
    else
    {
        cout << "No branching variable found. Node should have converged. Exiting..." << endl;
        exit(1);
    }

    return children;
}

} // namespace BB

} // namespace CENSO
