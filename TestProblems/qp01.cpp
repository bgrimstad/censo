/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "qp01.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintquadratic.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"

using std::cout;
using std::endl;

namespace CENSO
{

QP01::QP01()
    : fopt_known(-17),
      fopt_found(INF)
{
    testName = "Quadratic Program 01";
}

void QP01::runProblem()
{
    int numVars = 6;
    std::vector<VariablePtr> variables;
    std::vector<double> costs = {0, 0, 0, 0, 0, 1};
    std::vector<double> start = {1, 1, 0, 1, 0, 0};
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), 0, 1);
        var->setValue(start.at(i));
        variables.push_back(var);
    }
    variables.at(5)->setLowerBound(-INF);
    variables.at(5)->setUpperBound(INF);

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    {
        //std::vector<VariablePtr> vars = {variables.at()

        DenseMatrix A = -50*DenseMatrix::Identity(numVars,numVars);
        A(numVars-1, numVars-1) = 0;
        DenseVector b = DenseVector::Zero(numVars);
        b << 42, 44, 45, 47, 47.5, -1;
        double c = 0;

        auto quadcon = std::make_shared<ConstraintQuadratic>(variables, A, b, c, -INF, 0);

        cs->add(quadcon);
    }

    {
        std::vector<VariablePtr> vars = {variables.at(0),
                                         variables.at(1),
                                         variables.at(2),
                                         variables.at(3),
                                         variables.at(4)};

        DenseMatrix A(1,5);
        A << 20, 12, 11, 7, 4;
        DenseVector b = DenseVector::Zero(1);
        b << 40;

        auto lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);

        cs->add(lincon);
    }

    //SolverGurobi solver(cs);
    SolverIpopt solver(cs);
    auto res = solver.optimize();
    cout << res << endl;

    fopt_found = res.objectiveValue;
}

bool QP01::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
