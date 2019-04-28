/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "lp03.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintquadratic.h"
#include "SolverInterface/solveripopt.h"
#include "SolverInterface/solvergurobi.h"

using std::cout;
using std::endl;

namespace CENSO
{

LP03::LP03()
    : fopt_known(143750),
      fopt_found(INF)
{
    testName = "Linear Program 03";
}

void LP03::runProblem()
{
    int numVars = 11;
    std::vector<VariablePtr> variables;
    std::vector<double> costs = {4, 0.5, -5, -5, -5, -5, -5, -5, -5, -2, -0.2};
    for (int i = 0; i < numVars; i++)
    {
        auto var = std::make_shared<Variable>(-costs.at(i), 0, INF);
        variables.push_back(var);
    }

    DenseMatrix A = DenseMatrix::Zero(5,numVars);
    A <<    1,  0,  -30,    -50,    -65,    -75,    -80,    -80,    -75,    0,  0,
            0,  1,  -10,    -17,    -22,    -26,    -29,    -31,    -32,    0,  0,
            0,  0,  0,      5,      10,     15,     20,     25,     30,     -1, 0,
            0,  0,  10,     10,     10,     10,     10,     10,     10,     0,  -1,
            0,  0,  1,      1,      1,      1,      1,      1,      1,      0,  0;
    DenseVector b = DenseVector::Zero(5);
    b << 0, 0, 0, 0, 500;

    auto cs = std::make_shared<ConstraintLinear>(variables, A, b, false);

    SolverGurobi solver(cs);
    //SolverIpopt solver(cs);
    auto res = solver.optimize();
    cout << res << endl;

    fopt_found = -res.objectiveValue;
}

bool LP03::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
