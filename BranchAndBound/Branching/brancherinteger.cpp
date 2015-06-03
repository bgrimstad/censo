/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "brancherinteger.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

NodeList BrancherInteger::branchOnVariable(const NodePtr node, VariablePtr branchingVariable) const
{
    // Check if branching variable is valid
    assert(validateBranchingVariable(node, branchingVariable));

    double lb = branchingVariable->getLowerBound();
    double ub = branchingVariable->getUpperBound();

    // Calculate split value
    double zsplit = (ub + lb)/2.0;

    // Create left child
    branchingVariable->setUpperBound(std::floor(zsplit));
    NodePtr child1(node->inherit());
    //child1->setParentBranchingVariable(branchingVariable);

    // Create right child
    branchingVariable->setUpperBound(ub);
    branchingVariable->setLowerBound(std::ceil(zsplit));
    NodePtr child2(node->inherit());
    //child2->setParentBranchingVariable(branchingVariable);

    // Add children to node list
    NodeList children;
    children.addNode(child1);
    children.addNode(child2);

    // Reset bounds
    branchingVariable->setLowerBound(lb);
    branchingVariable->setUpperBound(ub);

    return children;
}

bool BrancherInteger::validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const
{
    if (!Brancher::validateBranchingVariable(node, branchingVariable)) return false;

    // Check that branching variable is integer
    if (branchingVariable->getType() != VariableType::INTEGER
        && branchingVariable->getType() != VariableType::BINARY)
    {
        //cout << "BrancherInteger: Selected branching variable is not integer!" << endl;
        return false;
    }

    // Check that bounds on variable are not equal
    if (branchingVariable->getUpperBound() == branchingVariable->getLowerBound())
    {
        //cout << "BrancherInteger: Selected branching variable is already fixed!" << endl;
        return false;
    }

    return true;
}

} // namespace BB

} // namespace CENSO
