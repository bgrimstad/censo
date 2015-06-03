/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "asaadi1.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintpolynomial.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solverbonmin.h"
#include "BranchAndBound/branchandbound.h"

using std::cout;
using std::endl;

namespace CENSO
{

Asaadi1::Asaadi1()
    : zopt_known(std::vector<double>(0)),
      zopt_found(std::vector<double>(0)),
      fopt_known(-INF),
      fopt_found(INF)
{
    testName = "Asaadi1";

    zopt_known.clear();
    zopt_known.push_back(0);
    zopt_known.push_back(1);
    zopt_known.push_back(2);
    zopt_known.push_back(0);

    fopt_known = -38.0;
}

void Asaadi1::runProblem()
{
    int numVars = 5; // 4 + 1 epigraph variable

//    std::vector<VariableType> varTypes(numVars,VariableType::CONTINUOUS);

//    std::vector<VariableType> varTypes = {
//        VariableType::INTEGER,
//        VariableType::INTEGER,
//        VariableType::CONTINUOUS,
//        VariableType::INTEGER,
//        VariableType::CONTINUOUS
//    };

    std::vector<VariableType> varTypes = {
        VariableType::INTEGER,
        VariableType::INTEGER,
        VariableType::INTEGER,
        VariableType::INTEGER,
        VariableType::CONTINUOUS
    };

    std::vector<double> costs = {0, 0, 0, 0, 1};
    std::vector<double> lb = {0, 0, 0, 0, -INF};
    std::vector<double> lb_unbounded(numVars, -INF);
    std::vector<double> ub = {100, 100, 100, 100, INF};

    std::vector<VariablePtr> vars;

    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i), varTypes.at(i));
        vars.push_back(var);
    }

    std::vector<VariablePtr> xvars = {
        vars.at(0),
        vars.at(1),
        vars.at(2),
        vars.at(3)
    };

    // Constraints
    ConstraintSetPtr constraints = std::make_shared<ConstraintSet>();

    // Constraint due to epigraph form
    DenseVector c0(9);
    c0.setZero();
    DenseMatrix C0(9,5);
    C0.setZero();

    c0(0) = 1; C0(0,0) = 2;
    c0(1) = 1; C0(1,1) = 2;
    c0(2) = 2; C0(2,2) = 2;
    c0(3) = 1; C0(3,3) = 2;
    c0(4) = -5; C0(4,0) = 1;
    c0(5) = -5; C0(5,1) = 1;
    c0(6) = -21; C0(6,2) = 1;
    c0(7) = 7; C0(7,3) = 1;
    c0(8) = -1; C0(8,4) = 1; // -t (epigraph)

    ConstraintPtr cPoly0 = std::make_shared<ConstraintPolynomial>(vars, c0, C0, false);
    constraints->add(cPoly0);

    DenseVector c1(9);
    c1.setZero();
    DenseMatrix C1(9,4);
    C1.setZero();

    c1(0) = -1; C1(0,0) = 2;
    c1(1) = -1; C1(1,1) = 2;
    c1(2) = -1; C1(2,2) = 2;
    c1(3) = -1; C1(3,3) = 2;
    c1(4) = -1; C1(4,0) = 1;
    c1(5) = 1; C1(5,1) = 1;
    c1(6) = -1; C1(6,2) = 1;
    c1(7) = 1; C1(7,3) = 1;
    c1(8) = 8; // -8 <= p(x) <=> 0 <= 8 + p(x)

    ConstraintPtr cPoly1 = std::make_shared<ConstraintPolynomial>(xvars, -c1, C1, false);
    constraints->add(cPoly1);

    DenseVector c2(7);
    c2.setZero();
    DenseMatrix C2(7,4);
    C2.setZero();

    c2(0) = -1; C2(0,0) = 2;
    c2(1) = -2; C2(1,1) = 2;
    c2(2) = -1; C2(2,2) = 2;
    c2(3) = -2; C2(3,3) = 2;
    c2(4) = 1; C2(4,0) = 1;
    c2(5) = 1; C2(5,3) = 1;
    c2(6) = 10; // -10 <= p(x)

    ConstraintPtr cPoly2 = std::make_shared<ConstraintPolynomial>(xvars, -c2, C2, false);
    constraints->add(cPoly2);

    DenseVector c3(7);
    c3.setZero();
    DenseMatrix C3(7,4);
    C3.setZero();

    c3(0) = -2; C3(0,0) = 2;
    c3(1) = -1; C3(1,1) = 2;
    c3(2) = -1; C3(2,2) = 2;
    c3(3) = -2; C3(3,0) = 1;
    c3(4) = 1; C3(4,1) = 1;
    c3(5) = 1; C3(5,3) = 1;
    c3(6) = 5; // -5 <= p(x)

    ConstraintPtr cPoly3 = std::make_shared<ConstraintPolynomial>(xvars, -c3, C3, false);
    constraints->add(cPoly3);

    // Problem 1: all variables are continuous and unbounded (not bounded from below)
    // f* = -44, x* = (0, 1, 2, -1)

    // Problem 2: all variables are continuous and non-negative
    // f* = -40.96, x* = (0, 1.038, 2.227, 0)

    // Problem 3: variables x1, x2, and x4 are non-negative integers, x3 is non-negative and continous
    // f* = -40.957, x* = (0, 1, 2.236, 0)

    // Problem 4: all variables are non-negative integers
    // f* = -38.000, x* = (0, 1, 2, 0)

    cout << constraints->isConstraintConvex() << endl;

    if (!constraints->isConstraintConvex())
        cout << "NOOOOOOOOOOOOOT CONVEEEEEEEEEEEEEEEX!!!!!!!!!!" << endl;
    else
        cout << "CONVEEEEEEEEEEEEEEEX!!!!!!!!!!" << endl;
    constraints->writeToGAMS("asaadi.gms");

    //SolverBonmin solver(constraints);
    BB::BranchAndBound solver(constraints);
    SolverResult res = solver.optimize();

    cout << res << endl;
    cout << res.primalVariables << endl;

    zopt_found = res.primalVariables;
    fopt_found = res.objectiveValue;
}

bool Asaadi1::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
    {
        return true;
    }
    else
    {
        return false;
    }
}

} // namespace CENSO
