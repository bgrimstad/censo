/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHANDBOUND_H
#define BRANCHANDBOUND_H

#include "Utils/definitions.h"
#include "Utils/timer.h"

#include "node.h"
#include "nodelist.h"
#include "OptimizationProblem/constraint.h"
#include "SolverInterface/solver.h"

#include "Branching/brancheradmin.h"
#include "BoundsTightening/fbbt.h"
#include "BoundsTightening/obbt.h"

namespace CENSO
{

namespace BB
{

class BranchAndBound : public Solver
{
public:

    BranchAndBound(ConstraintPtr constraints);

    BranchAndBound(const BranchAndBound &copy) = delete;
    BranchAndBound& operator=(const BranchAndBound &assign) = delete;

    ~BranchAndBound() {} // Smart pointers will delete themselves

    /*
     * Setters
     */
    void setEpsilon(double eps) { this->epsilon = eps; }
    void setFeasibilityTol(double feasibilityTol) { this->feasibilityTol = feasibilityTol; }
    void setSolveUpperBoundFlag(bool solveUpperBoundFlag) { this->solveUpperBoundFlag = solveUpperBoundFlag; }
    void setSolveUpperBoundDepth(int solveUpperBoundDepth) { this->solveUpperBoundDepth = solveUpperBoundDepth; }
    void setBoundsDeductionFlag(bool boundsDeductionFlag) { this->fbbtFlag = boundsDeductionFlag; }
    void setBoundsTighteningFlag(bool boundsTighteningFlag) { this->obbtFlag = boundsTighteningFlag; }
    void setBoundsTighteningDepth(bool boundsTighteningDepth) { this->obbtDepth = boundsTighteningDepth; }
    void setBoundsTighteningOccurence(bool boundsTighteningOccurence) { this->obbtOccurence = boundsTighteningOccurence; }
    void setPrintLevel(int printLevel) { this->printLevel = printLevel; }

protected:
    void initialize() override;
    SolverResult runOptimizer() override;

private:

    // Solves lower and upper bound of a node
    bool solveNodeLowerBound(NodePtr node);
    bool solveNodeUpperBound(NodePtr node);

    // Incumbent update
    void updateIncumbent(const NodePtr incumbentCandidate);

    // Update global lower bound (lowest bound of all nodes)
    void updateGlobalLowerBound();

    // Check if tree search can be terminated
    bool checkTerminationCriteria();

    // Bounds tightening
    bool doOBBT(NodePtr node) const;

    // Spatial measures of remaining search space
    double calculateSearchSpaceVolume() const; // Average interval length measure
    double calculateSearchSpaceVolume2() const; // Union of intervals
//    double calculateSearchSpaceVolume3() const; // Lebesque measure

    // Print function
    void printDebugInfo(std::string text, int level) const;

    // Node list
    NodeList nodeList;

    // Incumbent points to the node with the best upper bound found
    NodePtr incumbent;

    // Branchers
    BrancherAdmin brancherAdmin;

    // Bounds tighteners
    FBBT fbbt;
    OBBT obbt;

    // NLP solver
    // Possible improvement is to maintain only one solver and update it with new objectives and constraints.
    // Solver and interface classes must be updated to enable this.
    // SolverPtr nlpSolver;

    // Global (search tree) lower bound
    double globalLowerBound;

    // Internal variables
    unsigned int numIntegerVariables; // Number of binary/integer variables
    unsigned int numIterations; // Number of iterations - corresponds to number of nodes processed
    bool feasibleSolutionFound; // True if a feasible solution (upper bound) has been found
    unsigned int continuousBranchingDepth; // Depth where continuous branching/bisection starts

    // Algorithm parameters
    double epsilon; // epsilon convergence criteria: globalUpperBound - globalLowerBound < epsilon
    double epsilonrel; // Relative epsilon convergence criteria: |globalupper-globallower|/max(|globallower|, |globalupper|) < epsilonrel)
    unsigned int timeLimit; // Time limit in seconds.
    unsigned int maxIterations; // Maximum number of iterations - corresponds to number of branchings
    unsigned int maxInfeasibleIterations; // Maximum number of iterations without a feasible solution (upper bound)
    unsigned int maxDepth; // Maximum depth (stops the algorithm in case relaxation does not get tight resulting in indefinite branching)
    double feasibilityTol; // Important: must be forwarded to subproblem solvers!
    // Default response to solver error is to fathom node

    // Options for Heuristics
//    int rootNodeHeuristic; // Solve root node to optimality (Bonmin, etc) or feasibility (specialized flow heuristic, feasibility pump, user defined integers, or other construction heuristics)

    // Solve for upper bounds (heuristic)
    bool solveUpperBoundFlag;

    // Solve for upper bound only on certain depths:
    // may increase solution speed if upper bounds problems take long to solve
    unsigned int solveUpperBoundDepth; // 1: every node, 2: every other node, 3: every third node, etc

    // Bounds deduction options
    bool fbbtFlag;

    // Bounds tightening options
    bool obbtFlag;
    unsigned int obbtDepth; // Do bound tightening on nodes down to this depth
    unsigned int obbtOccurence; // On nodes deeper than boundsTighteningDepth do boundtightening every boundTighteningOccurence depth

    // Debug data
    unsigned int numLowerBoundsSolved;   // Lower bound problems solved
    unsigned int numUpperBoundsSolved;   // Upper bound problems solved
    unsigned int numInfeasibleChildren;  // Number of nodes that is feasible, but has two infeasible children
    unsigned int numFathoms;             // Number of fathomed/pruned nodes
    unsigned int numLeafNodes;           // Number of leaf nodes found

    /*
     * Print levels:
     * 0: prints simple summary
     * 1: prints detailed summary and issues during branching
     * 2: prints detailed summary and branching information
     * 3: prints detailed summary and branching information, including node details
     * 4: prints everything that could be of interest
     */
    int printLevel;

    // Timer
    Timer timer;
};

} // namespace BB

} // namespace CENSO

#endif // BRANCHANDBOUND_H
