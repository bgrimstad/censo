/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "lp01.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"

using std::cout;
using std::endl;

namespace CENSO
{

LP01::LP01()
    : fopt_known(10417.29129),
      fopt_found(INF)
{
    testName = "Linear Program 01";
}

void LP01::runProblem()
{
    int numVars = 6;
    std::vector<VariablePtr> variables;
    std::vector<double> costs = {67, 66, 66.3, 80, 78.5, 78.4};
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(-costs.at(i), 0, INF);
        variables.push_back(var);
    }

    DenseMatrix A = DenseMatrix::Zero(4,numVars);
    A << 0.8, 1.3, 0.2, 1.2, 1.7, 0.5,
         0.5, 0.2, 1.3, 0.7, 0.3, 1.5,
         0.4, 0.4, 0.4, 1.0, 1.0, 1.0,
         1.0, 1.05, 1.1, 0.8, 0.82, 0.84;
    DenseVector b = DenseVector::Zero(4);
    b << 140, 90, 120, 125;

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();
    ConstraintPtr lincon = std::make_shared<ConstraintLinear>(variables, A, b, false);
    ConstraintPtr lincon2 = std::make_shared<ConstraintLinear>(variables, A, b, false);

    cs->add(lincon);
    cs->add(lincon2);

    SolverGurobi solver(cs);
    //SolverIpopt solver(cs);
    auto res = solver.optimize();
    cout << res << endl;

    fopt_found = -res.objectiveValue;
}

bool LP01::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
