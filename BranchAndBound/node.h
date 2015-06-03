/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef NODE_H
#define NODE_H

#include <list>
#include <unordered_map>

#include "Utils/definitions.h"
#include "OptimizationProblem/constraint.h"
#include "bbutils.h"

namespace CENSO
{

namespace BB
{

class Node
{
public:

    Node(ConstraintPtr constraints);

    Node(const Node &copy);
    Node& operator=(const Node &assign) = delete;

    virtual ~Node() {}

    // Inherit from parent node
    Node* inherit() const;

    // Check if root node
    bool isRootNode() const { return depth == 0; }

    // Convergence check
    bool hasConverged(double epsilon) const;

    // Local refinement at lower bound solution
    void localRefinement();

    /*
     * Getters
     */
    ConstraintPtr getConstraints() { return constraints; }
    std::vector<double> getLowerBoundSolution() const { return solLowerBound; }
    std::vector<double> getUpperBoundSolution() const { return solUpperBound; }
    double getLowerBound() const { return objLowerBound; }
    double getUpperBound() const { return objUpperBound; }
    int getNodeId() const { return nodeId; }
    int getDepth() const { return depth; }
    int getVariableIndex(VariablePtr variable) const; // From pointer to index
    std::vector<VariablePtr> getBranchingVariables() const { return findBranchingVariables(constraints); }
    std::unordered_map<int,double> getPseudoCostsContinuous() { return pCostsContinuous; }

    /*
     * Setters
     */
    void setDepth(int depth){ this->depth = depth; }
    void setParentBranchingVariable(int var) { parentBranchingVariable = var; }

    // Bound update functions
    void updateLowerBound(double objLowerBound, std::vector<double> solLowerBound);
    void updateUpperBound(double objUpperBound, std::vector<double> solUpperBound);

    void printNode () const;

private: 

    ConstraintPtr constraints;

    // Node ID, depth and branching variable
    static int nextId; // Static ID counter for nodes
    const int nodeId; // Unique ID of node
    int depth; // Node depth (depth 0 at root node)
    int parentBranchingVariable; // Variable that was branched on to create this node

    // Solution variables
    double objLowerBound, objUpperBound;
    std::vector<double> solLowerBound, solUpperBound;

    // Node branching history
    std::unordered_map<int,double> pCostsInteger; // Not in use now
    std::unordered_map<int,double> pCostsContinuous;

    // Return next node ID and increment counter
    static int getNextNodeId();
};

typedef std::shared_ptr<Node> NodePtr;



} // namespace BB

} // namespace CENSO

#endif // NODE_H
