/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop01.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintbspline.h"
#include "BranchAndBound/branchandbound.h"
#include "Utils/bsplinepoly.h"
#include "Utils/definitions.h"

using std::cout;
using std::endl;

namespace CENSO
{

// P01 from B-spline paper
POP01::POP01()
{
    testName = "P01";
    fopt_known = -5.50801;
    fopt_found = INF;
}

void POP01::runProblem()
{
    // Problem 1 (Nataraj)
    // cout << "\n\nSolving problem P01..." << endl;

    int dim = 4;

    // x1,x2,l1,l2
    std::vector<double> costs = {-1, -1, 0, 0};
    std::vector<double> lb = {0,0,-INF,-INF};
    std::vector<double> ub = {3,4,INF,INF};
    //std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // 2*x1^4 - 8*x1^3 + 8*x1^2 + 2

        std::vector<VariablePtr> cvars = {vars.at(0), vars.at(2)};
        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {4};

        DenseVector c(5);
        c.setZero();
        c(0) = 2;
        c(1) = 0;
        c(2) = 8;
        c(3) = -8;
        c(4) = 2;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // 4*x1^4 - 32*x1^3 + 88*x1^2 - 96*x1 + 36

        std::vector<VariablePtr> cvars = {vars.at(0), vars.at(3)};
        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {4};

        DenseVector c(5); c.setZero();
        c(0) = 36;
        c(1) = -96;
        c(2) = 88;
        c(3) = -32;
        c(4) = 4;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // x2 - l1 <= 0, x2 - l2 <= 0

        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(3)
        };

        DenseMatrix A(2,3);
        A.setZero();
        A(0,0) = 1; A(0,1) = -1;
        A(1,0) = 1; A(1,2) = -1;

        DenseVector b; b.setZero(2);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }    

    BB::BranchAndBound solver(cs);
    Timer timer;
    timer.start();
    SolverResult res = solver.optimize();
    timer.stop();
    cout << res << endl;
    cout << res.primalVariables << endl;
    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP01::validateResult()
{
    return std::abs(fopt_found - fopt_known) <= 1e-3;
}

} // namespace CENSO
