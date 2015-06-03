/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "sixhumpcamel.h"

#include "datatable.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintbspline.h"
#include "BranchAndBound/branchandbound.h"
#include "SolverInterface/solveripopt.h"

using std::cout;
using std::endl;
using Splinter::BSpline;
using Splinter::DataTable;
using Splinter::BSplineType;

namespace CENSO
{

SixHumpCamel::SixHumpCamel()
    : zopt_known(std::vector<double>(0)),
      zopt_found(std::vector<double>(0)),
      fopt_known(-INF),
      fopt_found(INF)
{
    testName = "Six-hump camel back";

    // One of two global optima
    zopt_known.clear();
    zopt_known.push_back(0.089842);
    zopt_known.push_back(-0.712656);
    // Include aux variable

    fopt_known = -1.031628453;
}

DenseVector SixHumpCamel::sixHumpCamelFunction(DenseVector x)
{
    assert(x.rows() == 2);
    DenseVector y; y.setZero(1);
    y(0) = (4 - 2.1*x(0)*x(0) + (1/3.)*x(0)*x(0)*x(0)*x(0))*x(0)*x(0) + x(0)*x(1) + (-4 + 4*x(1)*x(1))*x(1)*x(1);
    return y;
}

void SixHumpCamel::runProblem()
{

    std::vector<double> costs = {0, 0, 1};
    std::vector<double> lb = {-1.9, -1.1, -INF};
    std::vector<double> ub = {1.9, 1.1, INF};
    std::vector<double> z0 = {0, 0, 0};

    std::vector<VariablePtr> vars;
    for (int i = 0; i < 3; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    DataTable data;

    double dx = 0.05;
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
            DenseVector yd = sixHumpCamelFunction(xd);

            data.addSample(x,yd(0));
        }
    }


//    DenseVector (*foo)(DenseVector);
//    foo = &sixHumpCamelFunction;
//    Bspline* bs = new Bspline(*data, 3);
//    bool testpassed = bs->testBspline(foo);
//    if (testpassed)
//    {
//        cout << "B-spline is very accurate:)" << endl;
//    }
//    else
//    {
//        cout << "B-spline is NOT very accurate:(" << endl;
//    }

    //    DenseMatrix controlpoints = bs->getControlPolygonVertices(bs->getDomainLowerBound(), bs->getDomainUpperBound());
    //    DenseMatrix bsmin = controlpoints.rowwise().minCoeff();
    //    cout << "Printing min coeff: " << endl;
    //    cout << bsmin << endl;

    ConstraintSetPtr constraints2 = std::make_shared<ConstraintSet>();
    BSpline bs(data, BSplineType::CUBIC_FREE);
    ConstraintPtr cbspline = std::make_shared<ConstraintBSpline>(vars, bs, true);
    constraints2->add(cbspline);

    // Solve Problem 1
    SolverIpopt ip3(constraints2);
    SolverResult status3 = ip3.optimize();

    // The optimal value
    cout << "Solved using B-spline approximation!" << endl;
    cout << status3 << endl;
    cout << status3.primalVariables;

    // Attempt to solve to global optimality
    BB::BranchAndBound bnb(constraints2);
    SolverResult bnbres = bnb.optimize();

    zopt_found = bnbres.primalVariables;
    fopt_found = bnbres.objectiveValue;
}

bool SixHumpCamel::validateResult()
{
    // Test if problem is solved maybe
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;

    return false;
}

} // namespace CENSO
