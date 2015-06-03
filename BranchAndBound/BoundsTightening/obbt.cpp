/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "obbt.h"

#include "Utils/definitions.h"
#include "OptimizationProblem/constraintlinear.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"

#include "thread"
#include "mutex"
#include "Utils/timer.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

// Global variable for now
std::mutex boundTightenerMutex;

OBBT::OBBT(double threshold, unsigned int maxIterations)
    : BoundsTightener(threshold, maxIterations),
      doParallelComputing(false),
      doAggressiveBoundsTightening(false)
{
}

bool OBBT::doTightening(ConstraintPtr constraints)
{
    // Get convex relaxation
    ConstraintPtr convexConstraints = constraints->getConvexRelaxation();

    assert(convexConstraints != nullptr);
    assert(convexConstraints->isConstraintConvex());

    // Tighten bounds of all complicating variables
    std::vector<VariablePtr> variables;
    for (auto var : constraints->getComplicatingVariables())
    {
        if (assertNear(var->getUpperBound(), var->getLowerBound()))
            continue;

        variables.push_back(var);
    }

    // Check if there are any variables to tighten
    if (variables.size() < 1)
        return true;

    // Tighten bounds
    return tightenBoundsSequential(convexConstraints, variables);
//    bool success = true;

//    if (doParallelComputing)
//    {
//        tightenBoundsParallel(convexConstraints, variables);
//    }
//    else
//    {
//        success = tightenBoundsSequential(convexConstraints, variables);
//    }

//    return success;
}

bool OBBT::tightenBoundsSequential(ConstraintPtr cs, std::vector<VariablePtr> variables)
{
    for (auto var : variables)
    {
        // Tighten bounds of single variable
        if (!tightenVariableBound(cs, var))
            return false;
    }

    return true;
}

void OBBT::tightenBoundsParallel(ConstraintPtr cs, std::vector<VariablePtr> variables)
{
    assert(variables.size() > 0);

    std::vector<std::thread> threads;
    int numThreads = std::thread::hardware_concurrency();
    if (numThreads < 1) numThreads = 1;

    int varsPerThread = variables.size()/numThreads;

//    printVec(variables);
//    cout << "Variables: " << variables.size() << endl;
//    cout << "Threads: " << numThreads << endl;

    for (int i = 0; i < numThreads; i++)
    {
        int start = i*varsPerThread;
        int end = (i+1)*varsPerThread;
        if (i == numThreads - 1) end = variables.size();

        //cout << "Thread " << i << ": " << start << " to " << end << endl;

        std::vector<VariablePtr>::const_iterator varStart = variables.begin() + start;
        std::vector<VariablePtr>::const_iterator varEnd = variables.begin() + end;
        std::vector<VariablePtr> threadVariables(varStart, varEnd);
        //printVec(threadVariables);

        threads.push_back(std::thread(&OBBT::tightenVariableBounds, this, cs, threadVariables));
//        threads.push_back(thread(&BoundTightener::test, this, convexConstraints, z0_aug, ref(threadVariables)));
    }

    for (int i = 0; i < numThreads; i++)
    {
        threads.at(i).join();
    }
}

void OBBT::tightenVariableBounds(ConstraintPtr cs, std::vector<VariablePtr> variables)
{
    for (auto var : variables)
    {
        // Read variable bounds
        boundTightenerMutex.lock();
        tightenVariableBound(cs, var);
        boundTightenerMutex.unlock();
    }
}

bool OBBT::tightenVariableBound(ConstraintPtr cs, VariablePtr variable)
{
    assert(cs->hasVariable(variable));

    auto vars = cs->getVariables();

    // Store and set objective costs to zero
    std::vector<double> costs;

    for (auto var : vars)
    {
        costs.push_back(var->getCost());
    }

    // Set costs for lower bound problem
    for (auto var : vars)
    {
        if (var == variable)
            var->setCost(1); // Min. variable
        else
            var->setCost(0);
    }

    //SolverIpopt solver_min(cs);
    SolverGurobi solver_min(cs);
    SolverResult result_min = solver_min.optimize();

    // Set costs for upper bound problem
    for (auto var : vars)
    {
        if (var == variable)
            var->setCost(-1); // Max. variable
        else
            var->setCost(0);
    }

    //SolverIpopt solver_max(cs);
    SolverGurobi solver_max(cs);
    SolverResult result_max = solver_max.optimize();

    // Reset costs
    int counter = 0;
    for (auto var : vars)
    {
        var->setCost(costs.at(counter));
        counter++;
    }

    // Check for infeasibility
    if (result_min.status == SolverStatus::INFEASIBLE
        || result_max.status == SolverStatus::INFEASIBLE)
        return false;

    // Update lower bound
    if (result_min.status == SolverStatus::OPTIMAL)
    {
        if (!variable->updateLowerBound(result_min.objectiveValue))
        {
            // This should never happen!
            cout << "Min bound" << endl;
            cout << *variable << endl;
            cout << result_min << endl;
            return false;
        }
    }

    // Update upper bound
    if (result_max.status == SolverStatus::OPTIMAL)
    {
        if (!variable->updateUpperBound(-result_max.objectiveValue))
        {
            cout << std::setprecision(10) << -result_max.objectiveValue << endl;
            cout << std::setprecision(10) << variable->getLowerBound() << endl;
            cout << std::setprecision(10) << variable->getLowerBound() + result_max.objectiveValue << endl;
            // This should never happen!
            cout << "Max bound" << endl;
            cout << *variable << endl;
            cout << result_max << endl;
            return false;
        }
    }

    // No update
    return true;
}

//void FBBT::tightenBoundsAggressive(ConstraintPtr cs,
//                                             std::vector<double> z0,
//                                             std::vector<int> variables,
//                                             std::vector<double> &lb,
//                                             std::vector<double> &ub)
//{
//    cout << "Performing aggressive bound tightening..." << endl;

//    // z0 is assumed to be the lower bound solution

//    assert(variables.size() > 0);

//    for (unsigned int i = 0; i < variables.size(); i++)
//    {
//        int var = variables.at(i);

//        double varopt = z0.at(var);
//        double varlb = lb.at(var);
//        double varub = ub.at(var);
//        double dist = varub - varlb;
//        //varopt = (varub + varlb)/2.0;

//        if (dist < 0.1) continue;

//        cout << "varlb: " << varlb << endl;
//        cout << "varlb: " << varub << endl;
//        cout << "var: " << varopt << endl;

//        // Try to update lower bound by introducing a fake upper bound
//        double ratio = 2.0/3.0; // 0 < ratio < 1
//        double fakeub = varlb + ratio*(varopt - varlb);
//        if ((fakeub - varlb)/dist > 0.1)
//        {
//            ConstraintPtr cs_copy(cs->clone());
//            std::vector<double> lb_copy;
//            std::vector<double> ub_copy;
//            cs_copy->getDomainBounds(lb_copy,ub_copy);

//            ub_copy.at(var) = fakeub;

//            cs_copy->setDomainBounds(lb_copy, ub_copy);

//            // Solve for new lower bound on variable
//            DenseMatrix cmin; cmin.setZero(1,z0.size()); cmin(var) = 1;
//            ObjectivePtr minObj(new ObjectiveLinear(cmin));
//            //SolverIpopt ip_min(minObj, cs_copy, z0);
//            SolverGurobi solver_min(minObj, cs_copy, z0);
//            solver_min.initialize();
//            int status_min = solver_min.optimize();

//            if (status_min != 1)
//            {
//                cout << "Able to tighten lower bound by aggressive FBBT!" << endl;
//                cout << "old lb: " << lb.at(var) << endl;
//                cout << "new lb: " << fakeub << endl;
//                cout << "ub: " << ub.at(var) << endl;
//                lb.at(var) = fakeub;
//                varlb = fakeub;
//                dist = varub - varlb;
//            }
//        }

//        // Try to update upper bound by introducing a fake lower bound
//        double fakelb = varub - ratio*(varub - varopt);
//        if ((varub - fakelb)/dist > 0.1)
//        {
//            ConstraintPtr cs_copy(cs->clone());
//            std::vector<double> lb_copy;
//            std::vector<double> ub_copy;
//            cs_copy->getDomainBounds(lb_copy,ub_copy);

//            lb_copy.at(var) = fakelb;

//            cs_copy->setDomainBounds(lb_copy, ub_copy);

//            // Solve for new lower bound on variable
//            DenseMatrix cmin; cmin.setZero(1,z0.size()); cmin(var) = 1;
//            ObjectivePtr minObj(new ObjectiveLinear(cmin));
//            //SolverIpopt ip_min(minObj, cs_copy, z0);
//            SolverGurobi solver_min(minObj, cs_copy, z0);
//            solver_min.initialize();
//            int status_min = solver_min.optimize();

//            if (status_min != 1)
//            {
//                cout << "Able to tighten upper bound by aggressive FBBT!" << endl;
//                ub.at(var) = fakelb;

//            }
//        }
//    }
//}

} // namespace BB

} // namespace CENSO
