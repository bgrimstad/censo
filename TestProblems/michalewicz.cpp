/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "michalewicz.h"

#include "datatable.h"
#include "OptimizationProblem/constraintbspline.h"
#include "SolverInterface/solveripopt.h"
#include "BranchAndBound/branchandbound.h"

using std::cout;
using std::endl;

using Splinter::BSpline;
using Splinter::DataTable;
using Splinter::BSplineType;

namespace CENSO
{

Michalewicz::Michalewicz()
    : zopt_known(std::vector<double>(0)),
      zopt_found(std::vector<double>(0)),
      fopt_known(-INF),
      fopt_found(INF)
{
    testName = "Michalewicz";

    zopt_known.clear();
    zopt_known.push_back(2.2029);
    zopt_known.push_back(1.5708);
    // Include aux variable

    fopt_known = -1.8013;
}

DenseVector Michalewicz::michalewiczFunction(DenseVector x)
{
    assert(x.rows() == 2);
    double pi = atan(1)*4;
    DenseVector y; y.setZero(1);
    y(0) = -sin(x(0))*pow(sin(x(0)*x(0)/pi), 20) -sin(x(1))*pow(sin(2*x(1)*x(1)/pi), 20);
    return y;
}

void Michalewicz::runProblem()
{
    double pi = atan(1)*4;

    std::vector<VariablePtr> vars = {
        std::make_shared<Variable>(0, 0, pi),
        std::make_shared<Variable>(0, 0, pi),
        std::make_shared<Variable>(1)
    };

    // Set starting points
    vars.at(0)->setValue(1.0);
    vars.at(1)->setValue(1.0);

    DataTable data;

    double dx = 0.05;
    for (double x1 = 0; x1 <= pi; x1+=dx)
    {
        for (double x2 = 0; x2 <= pi; x2+=dx)
        {
            std::vector<double> x = {x1, x2};

            DenseVector xd(2); xd << x1, x2;
            DenseVector yd = michalewiczFunction(xd);

            data.addSample(x,yd(0));
        }
    }

    // Test accuracy of B-spline
//    DenseVector (*foo)(DenseVector);
//    foo = &michalewiczFunction;
//    BSpline* bs = new BSpline(*data, 3);
//    bool testpassed = bs->testBspline(foo);
//    if (testpassed)
//    {
//        cout << "B-spline is very accurate:)" << endl;
//    }
//    else
//    {
//        cout << "B-spline is NOT very accurate:(" << endl;
//    }

    BSpline bs(data, BSplineType::CUBIC_FREE);
    auto constraint = std::make_shared<ConstraintBSpline>(vars, bs, false);

    //SolverIpopt solver(constraint);
    BB::BranchAndBound solver(constraint);

    // Optimize
    SolverResult result = solver.optimize();

    cout << result << endl;

    fopt_found = result.objectiveValue;
    zopt_found = result.primalVariables;

    cout << zopt_found << endl;
}

bool Michalewicz::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
