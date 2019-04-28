/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "node.h"
#include "Utils/eigen_utils.h"
#include <iostream>

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

int Node::nextId = 0;

// Constructor, intended for root node
Node::Node(ConstraintPtr constraints)
    : constraints(constraints->clone(true)), // Deep copy
      nodeId(getNextNodeId())
{
    objLowerBound = -INF;
    objUpperBound = INF;
    solLowerBound = std::vector<double>(constraints->getNumVariables(), 0.0);
    solUpperBound = std::vector<double>(constraints->getNumVariables(), 0.0);
    depth = 0;
    parentBranchingVariable = -1; // -1 for root node

    auto bvars = findBranchingVariables(this->constraints);

    // Initialize branching/pseudo costs
    double cost = 1.0;
    for (auto var : bvars)
    {
        if (var->getType() == VariableType::CONTINUOUS)
        {
            // Continuous
            pCostsContinuous.emplace(getVariableIndex(var), cost);
        }
        else
        {
            // Integer
            pCostsInteger.emplace(getVariableIndex(var), 1.0);
        }
    }

    assert(bvars.size() == pCostsInteger.size() + pCostsContinuous.size());
}

// Copy constructor - constraints are cloned
Node::Node(const Node &copy)
    : constraints(copy.constraints->clone(true)), // Deep copy
      nodeId(getNextNodeId()),
      depth(copy.depth),
      parentBranchingVariable(copy.parentBranchingVariable),
      objLowerBound(copy.objLowerBound),
      objUpperBound(copy.objUpperBound),
      solLowerBound(copy.solLowerBound),
      solUpperBound(copy.solUpperBound),
      pCostsInteger(copy.pCostsInteger),
      pCostsContinuous(copy.pCostsContinuous)
{
    // Nothing to do here, all values set
}

// Implements node inheritance
Node* Node::inherit() const
{
    // Create a child of this node
    Node* child = new Node(*this);

    // Update child depth
    child->setDepth(depth+1);

    return child;
}

bool Node::hasConverged(double epsilon) const
{
    if (std::abs(objUpperBound - objLowerBound) < epsilon)
        return true;

    if (assertNear(objUpperBound, objLowerBound, epsilon))
        return true;

    return false;
}

void Node::localRefinement()
{
    DenseVector sol = stdToEigenVec(solLowerBound);
    constraints->localRefinement(sol);
}

int Node::getNextNodeId()
{
    return nextId++;
}

void Node::updateLowerBound(double fLowerBound, std::vector<double> zLowerBound)
{
    double oldLowerBound = this->objLowerBound;
    this->objLowerBound = fLowerBound;
    this->solLowerBound = zLowerBound;

    // Update starting point
    //this->startingPoint = zLowerBound;

//        unordered_map<int,double>::iterator it_int = pCostsInteger.find(parentBranchingVariable);
//        if (it_int != pCostsInteger.end())
//        {
//        }

    // This assumes that updateLowerBound() is run once per node
    // OBS: bound tightening may give favorable cost to a variable if run between branching and lower bound solving
    std::unordered_map<int,double>::iterator it_cont = pCostsContinuous.find(parentBranchingVariable);
    if (it_cont != pCostsContinuous.end())
    {
        // This will not be run for root node since parentBranchingVariable is -1
        // Thus, oldLowerBound will be bounded > -INF
        if (true)
        {
            //it_cont->second *= 0.9 + std::abs(fLowerBound - oldLowerBound)/std::abs(fLowerBound);
            if (std::abs(objUpperBound - fLowerBound) < 0.5*std::abs(objUpperBound - oldLowerBound))
            {
                it_cont->second = it_cont->second*1.1; // Praise
            }
            else
            {
                it_cont->second = it_cont->second*0.5; // Bash
            }
        }
    }

}

void Node::updateUpperBound(double fUpperBound, std::vector<double> zUpperBound)
{
    // NOTE: consider testing that bounds are improved
    this->objUpperBound = fUpperBound;
    this->solUpperBound = zUpperBound;
}

void Node::printNode() const
{
    cout << "\n----------------------------------------" << endl;
    cout << "Information about node (" << nodeId << ")" << endl;
    cout << "Depth: " << depth << endl;
    cout << "Lower bound: " << objLowerBound << endl;
    cout << "Upper bound: " << objUpperBound << endl;
    cout << "Upper bound - Lower bound: " << objUpperBound - objLowerBound << endl;
}

int Node::getVariableIndex(VariablePtr variable) const
{
    auto vars = constraints->getVariables();
    auto it = std::find(vars.begin(), vars.end(), variable);
    assert(it != vars.end());
    return it - vars.begin();
}

} // namespace BB

} // namespace CENSO
