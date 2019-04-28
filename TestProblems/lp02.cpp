/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "lp02.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"

using std::cout;
using std::endl;

namespace CENSO
{

LP02::LP02()
    : fopt_known(7425),
      fopt_found(INF)
{
    testName = "Linear Program 02";
}

void LP02::runProblem()
{
    int numVars = 12;
    std::vector<VariablePtr> variables;
    std::vector<double> costs = {20, 40, 35, 120, 50, 60, 20, 70, 90, 35, 70, 40};
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), 0, INF);
        variables.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    {
        std::vector<VariablePtr> vars = {variables.at(0),
                                         variables.at(1),
                                         variables.at(2),
                                         variables.at(3)};

        DenseMatrix A(1,4); A << 1, 1, 1, 1;
        DenseVector b(1);   b << 100;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(4),
                                         variables.at(5),
                                         variables.at(6),
                                         variables.at(7)};

        DenseMatrix A(1,4); A << 1, 1, 1, 1;
        DenseVector b(1);   b << 75;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(8),
                                         variables.at(9),
                                         variables.at(10),
                                         variables.at(11)};

        DenseMatrix A(1,4); A << 1, 1, 1, 1;
        DenseVector b(1);   b << 90;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(0),
                                         variables.at(4),
                                         variables.at(8)};

        DenseMatrix A(1,3); A << -1, -1, -1;
        DenseVector b(1);   b << -30;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(1),
                                         variables.at(5),
                                         variables.at(9)};

        DenseMatrix A(1,3); A << -1, -1, -1;
        DenseVector b(1);   b << -75;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(2),
                                         variables.at(6),
                                         variables.at(10)};

        DenseMatrix A(1,3); A << -1, -1, -1;
        DenseVector b(1);   b << -90;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(3),
                                         variables.at(7),
                                         variables.at(11)};

        DenseMatrix A(1,3);  A << -1, -1, -1;
        DenseVector b(1);   b << -50;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    //cout << *cs << endl;

    SolverGurobi solver(cs);
    //SolverIpopt solver(cs);
    auto res = solver.optimize();
    cout << res << endl;

    fopt_found = res.objectiveValue;
}

bool LP02::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
