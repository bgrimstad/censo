/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bilinearrelaxationtest.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintbspline.h"
#include "OptimizationProblem/constraintbilinear.h"
#include "BranchAndBound/branchandbound.h"
#include "SolverInterface/solveripopt.h"

using std::cout;
using std::endl;
using Splinter::BSpline;
using Splinter::DataTable;
using Splinter::BSplineType;

namespace CENSO
{

BilinearRelaxationTest::BilinearRelaxationTest()
    : zopt_known(std::vector<double>(0)),
      zopt_found(std::vector<double>(0)),
      fopt_known(-INF),
      fopt_found(INF)
{
    testName = "Bilinear relaxation test";

    zopt_known.clear();
    zopt_known.push_back(-1);
    zopt_known.push_back(2.5);
    // Include aux variable

    fopt_known = -2.5;
}

DenseVector BilinearRelaxationTest::bilinearFunction(DenseVector x)
{
    assert(x.rows() == 2);
    DenseVector y; y.setZero(1);
    y(0) = x(0)*x(1);
    return y;
}

void BilinearRelaxationTest::runProblem()
{

    std::vector<double> cost, lb, ub, z0;

    cost.push_back(0);
    cost.push_back(0);
    cost.push_back(1);

    lb.push_back(-1.);
    lb.push_back(-1.);
    lb.push_back(-INF);

    ub.push_back(1);
    ub.push_back(2.5);
    ub.push_back(INF);

    z0.push_back(1);
    z0.push_back(-0.5);
    z0.push_back(-0.5);

    std::vector<VariablePtr> vars;
    for (unsigned int i = 0; i < 3; i++)
    {
        auto var = std::make_shared<Variable>(cost.at(i), lb.at(i), ub.at(i));
        var->setValue(z0.at(i));
        vars.push_back(var);
    }

    // Constraints
    ConstraintSetPtr constraints = std::make_shared<ConstraintSet>();

    ConstraintPtr myBilinearConstraint = std::make_shared<ConstraintBilinear>(vars,1,0,0);
    constraints->add(myBilinearConstraint);

    DataTable data;

    double dx = 0.5;
    for (double x1 = lb.at(0); x1 <= ub.at(0); x1+=dx)
    {
        for (double x2 = lb.at(1); x2 <= ub.at(1); x2+=dx)
        {
            std::vector<double> x;
            x.push_back(x1);
            x.push_back(x2);

            DenseVector xd; xd.setZero(2);
            xd(0) = x1;
            xd(1) = x2;
            DenseVector yd = bilinearFunction(xd);

            data.addSample(x,yd(0));
        }
    }

    ConstraintSetPtr constraints2 = std::make_shared<ConstraintSet>();
    BSpline bs(data, BSplineType::CUBIC_FREE);
    ConstraintPtr cbspline = std::make_shared<ConstraintBSpline>(vars, bs, true);
    constraints2->add(cbspline);

    // Test accuracy of B-spline
//    DenseVector (*foo)(DenseVector);
//    foo = &bilinearFunction;
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

    BB::BranchAndBound solver(constraints);
    //SolverIpopt solver(constraints);
    SolverResult res = solver.optimize();

    cout << res << endl;

    zopt_found = res.primalVariables;
    fopt_found = res.objectiveValue;

}

bool BilinearRelaxationTest::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
