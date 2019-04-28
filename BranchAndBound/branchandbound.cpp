/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "branchandbound.h"

#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"
#include "SolverInterface/solverbonmin.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <ctime>
#include "Utils/timer.h"
#include "Utils/metrics.h"
#include "Utils/eigen_utils.h"

namespace CENSO
{

namespace BB
{

using std::cout;
using std::endl;
using std::string;
using std::stringstream;

BranchAndBound::BranchAndBound(ConstraintPtr constraints)
    : Solver(constraints),

      // Internal variables
      nodeList(NodeList(BEST_FIRST)),
      brancherAdmin(BrancherAdmin()),
      fbbt(FBBT(0.5, 100)),
      obbt(OBBT(20, 2)),
      globalLowerBound(-INF),
      numIntegerVariables(0),
      numIterations(0),
      feasibleSolutionFound(false),

      // Settings for Branch-and-Bound
      epsilon(1e-6), // 1e-6
      epsilonrel(1e-6),
      timeLimit(1e6), // Seconds
      maxIterations(1e6), // 1e6
      maxInfeasibleIterations(1e6),
      maxDepth(1000), // Not in use
      feasibilityTol(1e-6), // Subproblem solvers must enforce the same value

      // Settings for heuristics
      solveUpperBoundFlag(true),
      solveUpperBoundDepth(10),
      fbbtFlag(true),
      obbtFlag(false), // Works best when there are many (connected) branching variables
      obbtDepth(1),
      obbtOccurence(1),

      // Internal counters
      numLowerBoundsSolved(0),
      numUpperBoundsSolved(0),
      numInfeasibleChildren(0),
      numFathoms(0),
      numLeafNodes(0),

      // Setting print level
      printLevel(1),
      timer(Timer())
{

    // Deduce bounds before checking?
//    if (boundsDeductionFlag)
//        constraints->reduceVariablesBounds();

    auto bvars = findBranchingVariables(constraints);

    // Check variables bounds
    for (unsigned int i = 0; i < bvars.size(); i++)
    {
        auto bvar = bvars.at(i);

        // Check that branching variable is bounded
        assert(bvar->getLowerBound() > -INF);
        assert(bvar->getUpperBound() < INF);

        if (bvar->getType() == VariableType::BINARY
            || bvar->getType() == VariableType::INTEGER)
        {
            // Correct bounds on binary variables
            if (bvar->getType() == VariableType::BINARY)
            {
                // Correct bounds to [0,1]
                if (bvar->getLowerBound() != 0 && bvar->getLowerBound() != 1)
                    bvar->setLowerBound(0);
                if (bvar->getUpperBound() != 0 && bvar->getUpperBound() != 1)
                    bvar->setUpperBound(1);
            }

            // Remove excess fractional values from bounds on integer variables
            // E.g. -0.5 <= x <= 2.5 is rounded to 0 <= x <= 2
            if (!isInteger(bvar->getLowerBound()))
            {
                cout << "Lower bound on variable " << i << " is rounded from " << bvar->getLowerBound();
                bvar->setLowerBound(std::ceil(bvar->getLowerBound()));
                cout << " to " << bvar->getLowerBound() << endl;
            }
            if (!isInteger(bvar->getUpperBound()))
            {
                cout << "Upper bound on variable " << i << " is rounded from " << bvar->getUpperBound();
                bvar->setUpperBound(std::floor(bvar->getUpperBound()));
                cout << " to " << bvar->getUpperBound() << endl;
            }
        }
    }

    // Count number of binary/integer variables and calculate depth where continuous branching starts
    // NOTE: Continuous branching may start higher in the tree due to bounds tightening!
    continuousBranchingDepth = 0;
    for (unsigned int i = 0; i < bvars.size(); i++)
    {
        auto bvar = bvars.at(i);

        if (bvar->getType() == VariableType::BINARY || bvar->getType() == VariableType::INTEGER)
        {
            numIntegerVariables++;

            if (bvar->getType() == VariableType::BINARY
                    && bvar->getLowerBound() != bvar->getUpperBound())
            {
                continuousBranchingDepth++;
            }

            if (bvar->getType() == VariableType::INTEGER
                    && bvar->getLowerBound() != bvar->getUpperBound())
            {
                // Calculate required number of bisections
                double x = std::log(bvar->getUpperBound() - bvar->getLowerBound()) / std::log(2);
                x = std::ceil(x);
                continuousBranchingDepth += (int)x;
            }
        }
    }

    // Check problem class
    ProblemClass problemClass = constraints->assessClass();
    if (problemClass == ProblemClass::LP)
    {
        printDebugInfo("Attempting to solve problem of class LP...", 1);
    }
    else if (problemClass == ProblemClass::CNLP)
    {
        printDebugInfo("Attempting to solve problem of class Convex NLP...", 1);
    }
    else if (problemClass == ProblemClass::NLP)
    {
        if (constraints->hasConvexRelaxation())
        {
            printDebugInfo("Attempting to solve problem of class Non-convex NLP with convex relaxations.", 1);
        }
        else
        {
            printDebugInfo("Attempting to solve problem of class Non-convex NLP without convex relaxations... The branch-and-bound is run as a heuristic and cannot guarantee a solution!", 1);
        }
    }
    else if (problemClass == ProblemClass::MILP)
    {
        printDebugInfo("Attempting to solve problem of class MILP...", 1);
    }
    else if (problemClass == ProblemClass::CMINLP)
    {
        printDebugInfo("Attempting to solve problem of class Convex MINLP...", 1);
    }
    else if (problemClass == ProblemClass::MINLP)
    {
        if (constraints->hasConvexRelaxation())
        {
            printDebugInfo("Attempting to solve problem of class Non-convex MINLP with convex relaxations...", 1);
        }
        else
        {
            printDebugInfo("Attempting to solve problem of class Non-convex MINLP without convex relaxations... The branch-and-bound is run as a heuristic and cannot guarantee a solution...", 1);
        }
    }
    else
    {
        printDebugInfo("CENSO cannot solve unknown problem class.", 1);
        exit(1);
    }
}

void BranchAndBound::initialize()
{
    // Create root node and put it in the node list (constraints are copied)
    NodePtr rootNode(new Node(constraints));
    nodeList.addNode(rootNode);

    // Initialize incumbent node - copy of root node
    NodePtr rootNode2(new Node(constraints));
    incumbent = rootNode2;

    // Finished initializing
    setInitialized(true);
}

/* Starts the Branch-and-Bound algorithm
 * Solves MILP and convex MINLP to global optimality.
 * Solves non-convex MINLP to global optimality (epsilon convergence)
 *  if valid convex relaxations exist for all non-convex constraints.
 * Minimizes the objective.
 *
 * Fathoming rules:
 * 1) Bounds tightening proved problem infeasible
 * 2) The lower bound problem is infeasible
 * 3) A node can be fathomed if its lower bound is greater than the global upper bound (the incumbent)
 * 4) A node can be fathomed if it is a leaf node, meaning that its upper and lower bound have converged
 */
SolverResult BranchAndBound::runOptimizer()
{
    // Check if initialized
    if(!getInitialized())
        initialize();

    //printDebugInfo("\nStarting Branch-and-Bound machinery!\n", 1);
    // Maybe print BB settings here

    std::ofstream myfile("results.csv");

    // Start timer
    timer.start();

    // Calculate search space volume
    double volume = calculateSearchSpaceVolume();

    while (!checkTerminationCriteria())
    {
        // Increment iteration counter
        numIterations++;

        // Calculate search space volume
        //volume = calculateSearchSpaceVolume();

        // Select node to process
        NodePtr currentNode = nodeList.selectNode();

        // Print solver info
        string ldel("| ");
        string cdel(" | ");
        string rdel(" |");

        // Default terminal size on linux is 80 chars
        if (printLevel > 0 && numIterations == 1)
        {
            // Print header
            cout << std::right;
            cout << ldel;
            cout << std::setw(8) << "Iter" << cdel;
            cout << std::setw(8) << "Time" << cdel; // Run time in seconds
            cout << std::setw(8) << "Gap" << cdel; // Absolute gap
            cout << std::setw(8) << "Upper" << cdel; // Upper bound
            cout << std::setw(8) << "Lower" << cdel; // Lower bound
            cout << std::setw(8) << "Depth" << cdel; // Node depth
            cout << std::setw(8) << "|L|" << cdel; // List size
            cout << std::setw(8) << "Volume" << rdel; // Remaining search space volume (%)
            cout << endl;

            int w = 0; w += 8*8; w += 2+7*3+2;
            cout << std::setfill('-') << std::setw(w) << "" << endl;
            cout << std::setfill(' ');
        }

        if (printLevel > 0 && (numIterations == 1 || numIterations % 100 == 0))
        {
            const int prec = cout.precision(); //
            cout << ldel;
            cout << std::setw(8) << std::setprecision(6) << numIterations << cdel;
            cout << std::setw(8) << std::setprecision(6) << timer.getSeconds() << cdel;
            cout << std::setw(8) << std::setprecision(6) << incumbent->getUpperBound() - globalLowerBound << cdel;
            cout << std::setw(8) << std::setprecision(6) << incumbent->getUpperBound() << cdel;
            cout << std::setw(8) << std::setprecision(6) << globalLowerBound << cdel;
            cout << std::setw(8) << std::setprecision(6) << currentNode->getDepth() << cdel;
            cout << std::setw(8) << std::setprecision(6) << nodeList.size() << cdel;
            cout << std::setw(8) << std::setprecision(6) << volume << rdel;
            cout << endl;
        }

        // Write to file
        if (myfile.is_open())
        {
            if (numIterations == 1)
            {
                myfile << "iter,time,gap,ub,lb,nodedepth,listsize,volume" << "\n";
            }

            double lb = std::max(globalLowerBound, -1e9);
            double ub = std::min(incumbent->getUpperBound(), 1e9);
            myfile << numIterations << ",";
            myfile << timer.getSeconds() << ",";
            myfile << ub - lb << ",";
            myfile << ub << ",";
            myfile << lb << ",";
            myfile << currentNode->getDepth() << ",";
            myfile << nodeList.size() << ",";
            myfile << volume << "\n";
            myfile.flush();
        }
        else cout << "Unable to open file";

        /*
         * sBB logic below
         */

        // Check if node can be fahomed (check is performed early to avoid unnecessary calculations)
        if (currentNode->getLowerBound() >= incumbent->getUpperBound())
        {
            numFathoms++;
            printDebugInfo("Node fathomed: lower bound is higher than incumbent.", 2);
            continue;
        }

        // Do feasibility-based bounds tightening
        {
            bool fbbtSuccess = true;
            if (fbbtFlag)
            {
                ConstraintSetPtr constraintSet = std::make_shared<ConstraintSet>();
                ConstraintPtr objCut = objectiveCut(currentNode->getConstraints(), incumbent->getUpperBound());
                constraintSet->add(currentNode->getConstraints());
                constraintSet->add(objCut);

                fbbtSuccess = fbbt.run(constraintSet);
            }

            // Check if node can be fathomed (TODO: test for robustness)
//            if (!fbbtSuccess)
//            {
//                numFathoms++;
//                printDebugInfo("Node fathomed: bounds deduction proved node infeasible.", 2);
//                continue;
//            }
        }

        // Perform optimality-based bounds tightening
        {
            bool obbtSuccess = true;
            if (doOBBT(currentNode))
            {
                ConstraintSetPtr constraintSet = std::make_shared<ConstraintSet>();
                ConstraintPtr objCut = objectiveCut(currentNode->getConstraints(), incumbent->getUpperBound());
                constraintSet->add(currentNode->getConstraints());
                constraintSet->add(objCut);

                obbtSuccess = obbt.run(constraintSet);
            }

            // Check if node can be fathomed (will be detected when solving lower bound problem)
//            if (!boundsTighteningSuccess)
//            {
//                numFathoms++;
//                printDebugInfo("Node fathomed: bounds tightening proved node infeasible.", 2);
//                continue;
//            }
        }

        // Solve lower bound problem
        bool solveSuccess = solveNodeLowerBound(currentNode);

        if (!solveSuccess)
        {
            numFathoms++;
            printDebugInfo("Node fathomed: lower bound problem infeasible.", 2);
            continue;
        }
        else if (currentNode->getLowerBound() >= incumbent->getUpperBound())
        {
            // Lower bound may have increased above incumbent
            numFathoms++;
            printDebugInfo("Node fathomed: lower bound is higher than incumbent.", 2);
            continue;
        }

        // Refine locally
        if (solveSuccess)
        {
            /*
             * NOTE: No knots inserted since global refinement
             * will have used all available knots.
             *
             * NOTE: This should not affect the solution time,
             * but it does!
             */
            currentNode->localRefinement();
        }

        // Check if lower bound solution is primal feasible
        std::vector<double> xLower = currentNode->getLowerBoundSolution();
        //Eigen::Map<DenseVector> xLower3(xLower.data(),xLower.size());
        DenseVector xLower2 = stdToEigenVec(xLower);

        if (isSolutionIntegerFeasible(currentNode->getConstraints()->getVariables(), xLower)
                && constraints->checkFeasibility(xLower2, feasibilityTol))
        {
            currentNode->updateUpperBound(currentNode->getLowerBound(), xLower);
        }
        else
        {
            // Optional primal heuristic: local search for upper bound
            solveNodeUpperBound(currentNode);
        }

        // Check if new incumbent has been found
        updateIncumbent(currentNode);

        // Check if node has converged; if not, branch
        if (currentNode->hasConverged(epsilon))
        {
            // NOTE: This will only happen when a new incumbent has converged
            // since nodes with a lower bound above the incumbent are discarded.
            numFathoms++;
            printDebugInfo("Node fathomed: node converged, no further branching required.", 2);
            continue;
        }
        else
        {
            // Branch
            NodeList children = brancherAdmin.branch(currentNode);
            if (!children.isEmpty())
            {
                nodeList.addNodes(children);
            }
            else
            {
                cout << "Could not branch on node!" << endl;
                // TODO: throw exception!
                exit(1);
            }
        }

    } // End While

    timer.stop();

    // Close file
    if (myfile.is_open())
    {
        double lb = std::max(globalLowerBound, -1e9);
        double ub = std::min(incumbent->getUpperBound(), 1e9);
        myfile << numIterations << ",";
        myfile << timer.getSeconds() << ",";
        myfile << ub - lb << ",";
        myfile << ub << ",";
        myfile << lb << ",";
        myfile << "-" << ",";
        myfile << nodeList.size() << ",";
        myfile << volume << "\n";
        myfile.close();
    }
    else cout << "Unable to open file";

    // Get final search space volume
    volume = calculateSearchSpaceVolume();

    // Print summary of exploration
    if (printLevel > 0)
    {
        string ldel("| ");
        string cdel(" | ");
        string rdel(" |");

        int w = 0; w += 8*8; w += 2+7*3+2;
        cout << std::setfill('-') << std::setw(w) << "" << endl;
        cout << std::setfill(' ');

        cout << ldel;
        cout << std::setw(8) << std::setprecision(6) << numIterations << cdel;
        cout << std::setw(8) << std::setprecision(6) << timer.getSeconds() << cdel;
        cout << std::setw(8) << std::setprecision(6) << incumbent->getUpperBound() - globalLowerBound << cdel;
        cout << std::setw(8) << std::setprecision(6) << incumbent->getUpperBound() << cdel;
        cout << std::setw(8) << std::setprecision(6) << globalLowerBound << cdel;
        cout << std::setw(8) << "-" << cdel;
        cout << std::setw(8) << std::setprecision(6) << nodeList.size() << cdel;
        cout << std::setw(8) << std::setprecision(6) << volume << rdel;
        cout << endl;

        cout << std::setfill('-') << std::setw(w) << "" << endl;
        cout << std::setfill(' ');
        cout << endl;
    }

    SolverResult result(SolverStatus::ERROR, incumbent->getUpperBound(), incumbent->getUpperBoundSolution());

    if (feasibleSolutionFound)
    {
        // NOTE: if the upper and lower bound of the incumbent has not converged, it means that
        // all child nodes of the incumbent node had worse upper bound solutions (which must be due to very small numerical differences in the solution)
        // A partial fix for this behaviour is to update the incumbent node if new upper bound <= incumbent

        if (assertNear(incumbent->getUpperBound(), globalLowerBound, epsilon, epsilonrel))
        {
            stringstream ss;
            ss << "CENSO terminated with certified optimal solution." << endl;
            printDebugInfo(ss.str(), 1);
            result.status = SolverStatus::OPTIMAL;
        }
        else
        {
            stringstream ss;
            ss << "CENSO terminated with feasible solution, but could not guarantee epsilon (" << epsilon << ") optimality." << endl;
            printDebugInfo(ss.str(), 1);
            result.status = SolverStatus::FEASIBLE;
        }

        if (printLevel > 1)
        {
            cout << endl;
            cout << "Optimal point x* = ";
            cout << incumbent->getUpperBoundSolution() << endl;
            cout << endl;
        }
    }
    else
    {
        cout << "CENSO unsuccessfully terminated without finding any feasible solutions." << endl;
        result.status = SolverStatus::INFEASIBLE;
    }

    return result;
}

/*
 * Solve lower bounding problem:
 * - For MILPs and convex MINLPs the integrality constraint is relaxed
 * - For non-convex MINLPs the integrality constraint is relaxed and convex relaxations must be available
 * - For non-convex NLPs convex relaxations must be available
 * - For convex NLPs and LPs no relaxation is done -> the problems are solved directly
 */
bool BranchAndBound::solveNodeLowerBound(NodePtr node)
{
    bool solveSuccessful = false;

    // Solve root node for possible lower bound
    ConstraintPtr nodeConstraints = node->getConstraints();

    // Get convex relaxation and add objective cut
    //ConstraintPtr nodeConstraintsRelaxed = nodeConstraints->getConvexRelaxation();
    ConstraintSetPtr nodeConstraintsRelaxed = std::make_shared<ConstraintSet>();
    nodeConstraintsRelaxed->add(nodeConstraints->getConvexRelaxation());
    if (node->getLowerBound() > -INF)
    {
        ConstraintPtr objCut = objectiveCut(nodeConstraints, node->getLowerBound(), false);
        //nodeConstraintsRelaxed->add(objCut); // NOTE: This results in more iterations in some cases...
    }

    // Set all variable types to continuous
    std::vector<VariableType> varTypes;
    {
        auto relVars = nodeConstraintsRelaxed->getVariables();
        for (unsigned int i = 0; i < relVars.size(); i++)
        {
            auto var = relVars.at(i);
            varTypes.push_back(var->getType());
            var->setType(VariableType::CONTINUOUS);
        }
    }

    // Select, initialize and run solver
    Solver* solver;

    if (nodeConstraintsRelaxed->isConstraintLinear())
    {
        solver = new SolverGurobi(nodeConstraintsRelaxed);
    }
    else
    {
        SolverIpopt* nlpSolver = new SolverIpopt(nodeConstraintsRelaxed);
        nlpSolver->setBoundRelaxFactor(0);
        solver = nlpSolver;
    }

    SolverResult lbReturn = solver->optimize();

    numLowerBoundsSolved++;

    // Reset variable types
    {
        auto relVars = nodeConstraintsRelaxed->getVariables();
        for (unsigned int i = 0; i < relVars.size(); i++)
        {
            auto var = relVars.at(i);
            var->setType(varTypes.at(i));
        }
    }

    if (lbReturn.status == SolverStatus::OPTIMAL)
    {
        solveSuccessful = true;
        std::vector<double> zopt_aux = lbReturn.primalVariables;
        //std::vector<double> zopt_pri(zopt_aux.begin(), zopt_aux.begin() + nodeConstraints->getNumVariables());
        std::vector<double> zopt_pri;

        // Extract values of original variables (excluding the auxiliary variables)
        auto relvars = nodeConstraintsRelaxed->getVariables();
        for (auto var : nodeConstraints->getVariables())
        {
            auto varit = std::find(relvars.begin(), relvars.end(), var);
            assert(varit != relvars.end());
            zopt_pri.push_back(zopt_aux.at(varit-relvars.begin()));
        }

        node->updateLowerBound(lbReturn.objectiveValue, zopt_pri);

        // Reduced cost bounds tightening
        // If using dual-adequate solver add optimality-based cuts of Ryoo and Sahinidis
        // Gurobi returns the reduced costs which can be used in the same manner as the range constraint duals
        bool dualSolver = solver->dualAdequate();
        double deltaBounds = incumbent->getUpperBound() - node->getLowerBound();

        if (dualSolver && deltaBounds < INF && deltaBounds > 0) // Negative deltaBounds for nodes to be fathomed
        {
            // Get dual variables
            std::vector<double> lowerBoundDuals, upperBoundDuals, constraintDuals;
            lowerBoundDuals = lbReturn.lowerBoundDualVariables;
            upperBoundDuals = lbReturn.upperBoundDualVariables;
            constraintDuals = lbReturn.constraintDualVariables;
            assert(lowerBoundDuals.size() == nodeConstraintsRelaxed->getNumVariables());
            assert(lowerBoundDuals.size() == upperBoundDuals.size());
            assert(constraintDuals.size() == nodeConstraintsRelaxed->getNumConstraints());

            // Update variable bounds (all variables)
            auto varsRelaxed = nodeConstraintsRelaxed->getVariables();
            for (unsigned int i = 0; i < nodeConstraintsRelaxed->getNumVariables(); i++)
            {
                auto var = varsRelaxed.at(i);

                if (var->getLowerBound() > -INF
                    && assertNear(zopt_aux.at(i), var->getLowerBound())
                    && lowerBoundDuals.at(i) > 0)
                {
                    double newub = var->getLowerBound() + deltaBounds/lowerBoundDuals.at(i);
                    var->updateUpperBound(newub);
                }

                if (var->getUpperBound() < INF
                    && assertNear(zopt_aux.at(i), var->getUpperBound())
                    && upperBoundDuals.at(i) > 0)
                {
                    double newlb = var->getUpperBound() - deltaBounds/upperBoundDuals.at(i);
                    var->updateLowerBound(newlb);
                }
            }

//            {
//                // Update constraint bounds
//                std::vector<double> clb, cub, nclb, ncub;
//                nodeConstraintsRelaxed->getCodomainBounds(clb, cub);
//                nclb = clb;
//                ncub = cub;

//                assert(constraintDuals.size() == clb.size());

//                for (unsigned int i = 0; i < clb.size(); i++)
//                {
//                    /*
//                     * Let lambda_i = constraintDuals.at(i) be the dual variable for
//                     * the convex constraint g_i(x) <= 0
//                     * if (lambda_i > 0) then a valid inequality is: g_i(x) >= -(U-L)/lambda_i
//                     */
//                    if (constraintDuals.at(i) > 0)
//                    {
//                        nclb.at(i) = std::max(clb.at(i), cub.at(i) - deltaBounds/constraintDuals.at(i));
//                    }
//                }

//                // Set new constraint bounds
//                nodeConstraintsRelaxed->setCodomainBounds(nclb, ncub);
//            }

        }
        // END of Reduced cost
    }

    delete solver;

    return solveSuccessful;
}

/*
 * Solves a node's upper bounding problem:
 * - For convex mixed-integer problems (MILPs and convex MINLPs) a heuristical search may be performed
 * - For non-convex MINLP problems a convex MINLPs solver (e.g. Bonmin) may be run to local optimality
 * - For non-convex NLPs a NLP solver (e.g. Ipopt) may be run to local optimality
 * - Not needed for convex NLP and LP problems (they are solved directly by the lower bound problem of the root node)
 */
bool BranchAndBound::solveNodeUpperBound(NodePtr node)
{
    bool solveSuccessful = false;
    bool mixedIntegerProblem = isIntegerProblem(node->getConstraints());

    // TODO: move this logic to a separate function
    if (!solveUpperBoundFlag)
        return solveSuccessful;

    // Do not solve for upper bound
    if (node->getDepth() > 0 && node->getDepth() % solveUpperBoundDepth != 0)
        return solveSuccessful;

    // Do not attempt to solve minlp problems
    if (mixedIntegerProblem && node->getDepth() != 0)
        return solveSuccessful;

//    if (node->getDepth() == 0)
//        mixedIntegerProblem = true; // Use MINLP solver to solve root node

    // Node constraints are cloned since bounds will be altered here
    ConstraintPtr nodeConstraints = node->getConstraints()->clone(true);

    // Lock integer variables to starting point values
    if (node->getDepth() == 0)
    {
//        // Root node heuristic
//        // Use initial guess to lock integer variables
//        std::vector<double> clb, cub;
//        nodeConstraints->getDomainBounds(clb, cub);
//        for (unsigned int i = 0; i < variableTypes.size(); i++)
//        {
//            if (variableTypes.at(i) == BINARY || variableTypes.at(i) == INTEGER)
//            {
//                if (clb.at(i) != cub.at(i))
//                {
//                    // Root node heuristic
//                    if (isInteger(z0.at(i)) && z0.at(i) <= cub.at(i) && z0.at(i) >= clb.at(i))
//                    {
//                        cub.at(i) = z0.at(i); // lock upper bound to z0
//                        clb.at(i) = z0.at(i); // lock lower bound to z0
//                    }
//                    else
//                    {
//                        int lock = randomInteger(clb.at(i), cub.at(i));
//                        if (isInteger(lock))
//                        {
//                            clb.at(i) = lock;
//                            cub.at(i) = lock;
//                        }
//                        else
//                        {
//                            // Else, lock at upper bound. Should not happen.
//                            clb.at(i) = cub.at(i);
//                        }
//                    }
//                }
//            }
//        }

        // Set bounds
        //nodeConstraints->setDomainBounds(clb, cub);

        // Integer variables locked
        //mixedIntegerProblem = false;
    }

    // If integer problem switch to Bonmin w/ multi-start heuristic
    if (mixedIntegerProblem)
    {
        return solveSuccessful;
//        if (node->getDepth() != 0)
//            return solveSuccessful;

        SolverBonmin solver(nodeConstraints);
        SolverResult result = solver.optimize();
        if (result.status == SolverStatus::OPTIMAL
            && isSolutionIntegerFeasible(nodeConstraints->getVariables(), result.primalVariables))
        {
            node->updateUpperBound(result.objectiveValue, result.primalVariables);
            solveSuccessful = true;
        }
    }
    else
    {
        SolverIpopt solver(nodeConstraints);
        // Important for B-spline constraints which cannot be evaluated outside domain.
        // However, it can slow convergence of optimality gap for B-spline (and non-B-spline) problems.
        solver.setBoundRelaxFactor(0);

        SolverResult result = solver.optimize();
        if (result.status == SolverStatus::OPTIMAL
            && isSolutionIntegerFeasible(nodeConstraints->getVariables(), result.primalVariables))
        {
            node->updateUpperBound(result.objectiveValue, result.primalVariables);
            solveSuccessful = true;
        }
    }

    numUpperBoundsSolved++;

    if (solveSuccessful)
    {
        printDebugInfo("Upper bound problem solved", 2);

        // Test feasibility here (tolerance decided by boundRelaxFactor)
        auto x = stdToEigenVec(node->getUpperBoundSolution());
        if (!node->getConstraints()->checkFeasibility(x, feasibilityTol))
        {
            cout << "Upper bound solution is not feasible!" << endl;
            exit(1);
        }
    }
    else
    {
        printDebugInfo("Found no feasible solution for upper bound problem", 2);        
    }

    return solveSuccessful;
}

void BranchAndBound::updateIncumbent(const NodePtr incumbentCandidate)
{
    // Check if candidate is new incumbent
    if (isSolutionIntegerFeasible(incumbentCandidate->getConstraints()->getVariables(), incumbentCandidate->getUpperBoundSolution())
            && incumbentCandidate->getUpperBound() < incumbent->getUpperBound())
    {
        // Print info
        stringstream ss;
        ss  << "Incumbent update: " << endl
            << "Old upper bound =\t" << incumbent->getUpperBound() << endl
            << "New upper bound =\t" << incumbentCandidate->getUpperBound()
            << endl;
        printDebugInfo(ss.str(), 2);

        // Update incumbent
        incumbent.reset();
        incumbent = incumbentCandidate;

        feasibleSolutionFound = true;
    }
}

void BranchAndBound::updateGlobalLowerBound()
{
    globalLowerBound = incumbent->getUpperBound(); // When node list is empty the lower bound will equal the incumbent
    for (NodeList::const_iterator it = nodeList.begin(); it != nodeList.end(); it++)
    {
        if (globalLowerBound > (*it)->getLowerBound())
        {
            globalLowerBound = (*it)->getLowerBound();
        }
    }
}

bool BranchAndBound::checkTerminationCriteria()
{
    /*
     * Check remaining nodes and update global lower bound.
     */
    updateGlobalLowerBound();

    if (nodeList.size() <= 0)
    {
        // Node list empty
        printDebugInfo("Terminating: empty node list!", 1);
        return true;
    }
    else if (feasibleSolutionFound && assertNear(incumbent->getUpperBound(), globalLowerBound, epsilon, epsilonrel))
    {
        printDebugInfo("Terminating: epsilon convergence!", 1);
        return true;
    }
    else if (numIterations >= maxIterations)
    {
        // Maximum iterations
        printDebugInfo("Terminating: maximum iterations!", 1);
        return true;
    }
    else if (!feasibleSolutionFound && numIterations >= maxInfeasibleIterations)
    {
        // Infeasible iterations
        printDebugInfo("Terminating: maximum infeasible iterations!", 1);
        return true;
    }
    else if (timer.getSeconds() > timeLimit)
    {
        // Time limit reached
        printDebugInfo("Terminating: Reached time limit!", 1);
        return true;
    }

    // Not finished
    return false;
}

bool BranchAndBound::doOBBT(NodePtr node) const
{
    if (obbtFlag)
    {
        int depth = node->getDepth();
        if (depth >= obbtDepth && (depth-obbtDepth) % obbtOccurence == 0)
        //if (depth >= boundsTighteningDepth && depth >= continuousBranchingDepth && (depth-continuousBranchingDepth-boundsTighteningDepth) % boundsTighteningOccurence == 0)
        {
            return true;
        }
    }

    return false;
}

//double BranchAndBound::calculateSearchSpaceVolume3() const
//{
//    double relativeVolume = 0; // Final search space volume

//    std::vector<double> olb,oub;
//    constraints->getDomainBounds(olb,oub);

//    // Calculate volume of remaining nodes
//    for(NodeList::const_iterator it = nodeList.begin(); it != nodeList.end(); it++)
//    {
//        ConstraintPtr cs = (*it)->getConstraints();
//        std::vector<double> nlb,nub;
//        cs->getDomainBounds(nlb,nub);

//        double nodeContinuousVolume = 1;
//        double nodeIntegerVolume = 1;

//        for (unsigned int i = 0; i < branchingVariables.size(); i++)
//        {
//            int bvar = branchingVariables.at(i);

//            double odelta = oub.at(bvar) - olb.at(bvar);
//            double ndelta = nub.at(bvar) - nlb.at(bvar);

//            // Check if variable is considered branchable
//            if (odelta == 0)
//                continue;

//            // Check if continuous branching variable
//            // Volumetric does not work for integer variables (points)
//            if (variableTypes.at(bvar) != VariableType::CONTINUOUS)
//            {
//                nodeIntegerVolume *= (ndelta+1)/(odelta+1);
//            }
//            else
//            {
//                // Relative continuous volume using the Lebesgue measure: |I_1|*|I_2|*...*|I_n|, for intervals I_i
//                nodeContinuousVolume *= (ndelta)/(odelta);
//            }
//        }

//        relativeVolume += nodeIntegerVolume*nodeContinuousVolume;
//    }

//    return 100*relativeVolume;
//}

double BranchAndBound::calculateSearchSpaceVolume2() const
{
    if (nodeList.isEmpty())
        return 0;

    // Calculate sum of variable intervals as union over all nodes
    double relCombinations = 1;
    double relInterval = 0;

    unsigned int numContBranchingVars = 0;

    auto vars = constraints->getVariables();
    auto bvars = findBranchingVariables(constraints);

    for (auto bvar : bvars)
    {
        // Get index
        auto it = std::find(vars.begin(), vars.end(), bvar);
        assert(it != vars.end());
        int bvari = it - vars.begin();

        // Compute interval
        double origInterval = bvar->getUpperBound() - bvar->getLowerBound();

        if (origInterval == 0)
            continue;

        double unionInterval = 0;

        // Get variable interval for all nodes
        std::vector<std::pair<double,double>> intervals;
        for (NodeList::const_iterator it = nodeList.begin(); it != nodeList.end(); it++)
        {
            ConstraintPtr cs = (*it)->getConstraints();
            auto varNode = cs->getVariables().at(bvari);
            intervals.push_back(std::pair<double,double>(varNode->getLowerBound(), varNode->getUpperBound()));
        }

        // Sort the intervals
        std::sort(intervals.begin(),intervals.end(),
                  [](const std::pair<double,double> &a, const std::pair<double,double> &b) -> bool
        {
            return a.first < b.first;
        });

        // Calculate union of intervals
        std::vector<std::pair<double,double>> intervalsMerged;
        double from = intervals.at(0).first;
        double to = intervals.at(0).second;

        for (unsigned int i = 0; i < intervals.size()-1; i++)
        {

            double lb = intervals.at(i+1).first;
            double ub = intervals.at(i+1).second;

            if (to < lb)
            {
                intervalsMerged.push_back(std::pair<double,double>(from,to));
                from = lb;
            }

            to = std::max(ub,to);
        }

        intervalsMerged.push_back(std::pair<double,double>(from,to));

        if (bvar->getType() == VariableType::CONTINUOUS)
        {
            // Sum intervals
            for (auto interval : intervalsMerged)
            {
                unionInterval += interval.second - interval.first;
            }

            relInterval += unionInterval/origInterval;
            numContBranchingVars++;
        }
        else
        {
            // Sum intervals
            for (auto interval : intervalsMerged)
            {
                unionInterval += interval.second - interval.first + 1;
            }
            relCombinations *= unionInterval/(origInterval + 1);
        }
    }
    //cout << "Rel combinations: " << relCombinations << endl;

    return 100*relCombinations*relInterval/numContBranchingVars;
}

double BranchAndBound::calculateSearchSpaceVolume() const
{
    if (nodeList.isEmpty())
        return 0;

    // Calculate sum of variable intervals as union over all nodes
    double totalInterval = 0;

    unsigned int numContBranchingVars = 0;

    auto vars = constraints->getVariables();
    auto bvars = findBranchingVariables(constraints);

    for (auto bvar : bvars)
    {
        auto it = std::find(vars.begin(), vars.end(), bvar);
        assert(it != vars.end());
        int bvari = it - vars.begin();

        double origInterval = bvar->getUpperBound() - bvar->getLowerBound();

        if (origInterval == 0)
            continue;

        double averageInterval = 0;

        // Get variable interval for all nodes
        for (NodeList::const_iterator it = nodeList.begin(); it != nodeList.end(); it++)
        {
            ConstraintPtr cs = (*it)->getConstraints();
            auto bvarNode = cs->getVariableAt(bvari);
            averageInterval += bvarNode->getUpperBound() - bvarNode->getLowerBound();
        }

        averageInterval /= nodeList.size();

        averageInterval /= origInterval;

        totalInterval += averageInterval;

        numContBranchingVars++;
    }

    return 100*totalInterval/numContBranchingVars;
}

void BranchAndBound::printDebugInfo(string text, int level) const
{
    if (level <= printLevel)
    {
        cout << text << endl;
    }
}

} // namespace BB

} // namespace CENSO
