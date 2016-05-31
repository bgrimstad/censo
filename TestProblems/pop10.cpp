/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop10.h"
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

POP10::POP10()
{
    testName = "P10";
    fopt_known = -213;
    fopt_found = INF;
}

void POP10::runProblem()
{
    // Problem 10 (Nataraj), QP
    // Five 1-D B-splines
    cout << "\n\nSolving problem P10..." << endl;

    int dim = 6+5; // x1,..,x5,y,l1,...,l5

    // x1,x2,x3,x4,l1,l2,l3
    std::vector<double> costs = {0,0,0,0,0,-10,1,1,1,1,1};
    std::vector<double> lb = {0,0,0,0,0,0,-INF,-INF,-INF,-INF,-INF};
    std::vector<double> ub = {1,1,1,1,1,INF,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    std::vector<double> a = {-10.5,-7.5,-3.5,-2.5,-1.5};

    // Add one B-spline for each variable
    for (int i = 0; i < 5; i++)
    {
        std::vector<VariablePtr> cvars = {
            vars.at(i),
            vars.at(i+6)
        };

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        // Poly coeffs
        DenseVector c(3);
        c.setZero();
        c(0) = 0;
        c(1) = a.at(i);
        c(2) = -0.5; // -50 or -0.5 (Floudas' problem has -50)

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(4),
            vars.at(5)
        };

        DenseMatrix A = DenseMatrix::Zero(2,6);
        A(0,0) = 6; A(0,1) = 3; A(0,2) = 3; A(0,3) = 2; A(0,4) = 1;
        A(1,0) = 10; A(1,2) = 10; A(1,5) = 1;

        DenseVector b;
        b.setZero(2);
        b(0) = 6.5;
        b(1) = 20;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    SolverResult res = bnb.optimize();
    timer.stop();
    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << "Time: " << timer.getMicroSeconds() << " (us)" << endl;

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP10::validateResult()
{
    return std::abs(fopt_found - fopt_known) <= 1e-3;
}

} // namespace CENSO
