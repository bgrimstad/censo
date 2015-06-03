/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "branchercontinuous.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

BrancherContinuous::BrancherContinuous()
    : Brancher(),
      branchingThreshold(1e-10)
{
}

NodeList BrancherContinuous::branchOnVariable(const NodePtr node, VariablePtr branchingVariable) const
{
    // Check if branching variable is valid
    assert(validateBranchingVariable(node, branchingVariable));

    // Split at convex combination of variable lower bound solution and midpoint.
    double alpha = 0.25;    // 0 < alpha < 1 (default: 0.25)
    double beta = 0.2;      // 0 < beta < 0.5 (default: 0.2)

    std::vector<double> zL = node->getLowerBoundSolution();
    double lb = branchingVariable->getLowerBound();
    double ub = branchingVariable->getUpperBound();
    double b = beta*(ub - lb);

    double zLi = zL.at(node->getVariableIndex(branchingVariable));
    double zsplit = alpha*(ub + lb)/2.0 + (1-alpha)*zLi;
    if (zsplit < lb + b)
    {
        zsplit = lb + b;
    }
    else if (zsplit > ub - b)
    {
        zsplit = ub - b;
    }

    // Split in middle (old method)
    //double zsplit = (ub + lb)/2.0;

    // Create left child
    branchingVariable->setUpperBound(zsplit);
    NodePtr child1(node->inherit()); // Deep copy (new variables)
    //child1->setParentBranchingVariable(branchingVariable);

    // Create right child
    branchingVariable->setUpperBound(ub);
    branchingVariable->setLowerBound(zsplit);
    NodePtr child2(node->inherit()); // Deep copy (new variables)
    //child2->setParentBranchingVariable(branchingVariable);

    // Add children to node list
    NodeList children;
    children.addNode(child1);
    children.addNode(child2);

    // Reset bounds on variable
    branchingVariable->setLowerBound(lb);
    branchingVariable->setUpperBound(ub);

    return children;
}

bool BrancherContinuous::validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const
{
    if (!Brancher::validateBranchingVariable(node, branchingVariable)) return false;

    // Check that branching variable is continuous
    if (branchingVariable->getType() != VariableType::CONTINUOUS)
    {
        //cout << "BrancherContinuous: Selected branching variable is not continuous!" << endl;
        return false;
    }

    // Check that bounds on variable are not equal
    auto diff = std::abs(branchingVariable->getUpperBound() - branchingVariable->getLowerBound());
    if (diff <= branchingThreshold)
    {
        //cout << "BrancherContinuous: Selected branching variable's bounds within branching threshold!" << endl;
        return false;
    }

    return true;
}

} // namespace BB

} // namespace CENSO
