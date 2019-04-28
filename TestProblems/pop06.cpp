/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop06.h"
#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintbspline.h"
#include "BranchAndBound/branchandbound.h"
#include "Utils/bsplinepoly.h"
#include "Utils/bspline_wrapper.h"
#include "Utils/definitions.h"

using std::cout;
using std::endl;

namespace CENSO
{

POP06::POP06()
{
    testName = "P06";
    fopt_known = 6395.5078;
    fopt_found = INF;
}

void POP06::runProblem()
{
    // Problem 6 (Nataraj)
    // 5 2-D B-splines
    int dim = 4+2;

    // x1,x2,x3,x4,l1,l2
    std::vector<double> costs = {0, 0, 0, 0, 1, 0};
    std::vector<double> lb = {1,0.625,47.5,90,-INF,-INF};
    std::vector<double> ub = {1.1375,1,52.5,112,INF,INF};
    std::vector<double> z0(dim,0);

    // x1,x2,x3,x4,l1,l2,l3,l4,l5
//    std::vector<double> lb = {1,0.625,47.5,90,-INF,-INF,-INF,-INF,-INF};
//    std::vector<double> ub = {1.1375,1,52.5,112,INF,INF,INF,INF,INF};

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // obj = l1

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(4)
        };

        std::vector<double> thislb = {
            cvars.at(0)->getLowerBound(),
            cvars.at(1)->getLowerBound(),
            cvars.at(2)->getLowerBound(),
            cvars.at(3)->getLowerBound()
        };

        std::vector<double> thisub = {
            cvars.at(0)->getUpperBound(),
            cvars.at(1)->getUpperBound(),
            cvars.at(2)->getUpperBound(),
            cvars.at(3)->getUpperBound()
        };

        std::vector<unsigned int> deg = {2,1,2,1};

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(0) = 0.6224;
        c(1) = 1.7781;
        c(2) = 3.1661;
        c(3) = 19.84;

        // Poly exponents
        DenseMatrix E(4,4);
        E.setZero();
        E(0,2) = 1; E(0,3) = 1;
        E(1,1) = 1; E(1,2) = 2;
        E(2,0) = 2; E(2,3) = 1;
        E(3,0) = 2; E(3,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);

        //DenseMatrix cpoints = bs.getControlPoints();
        //cout << cpoints << endl;
    }

    { // noncon = l2

        std::vector<VariablePtr> cvars = {
            vars.at(2),
            vars.at(3),
            vars.at(5)
        };

        std::vector<double> thislb = {
            cvars.at(0)->getLowerBound(),
            cvars.at(1)->getLowerBound()
        };

        std::vector<double> thisub = {
            cvars.at(0)->getUpperBound(),
            cvars.at(1)->getUpperBound()
        };

        std::vector<unsigned int> deg = {3,1};

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = -1;
        c(1) = -(4.0/3.0);

        // Poly exponents
        DenseMatrix E(2,2);
        E.setZero();
        E(0,0) = 2; E(0,1) = 1;
        E(1,0) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints of auxiliary variables
        DenseMatrix A(4,6);
        A.setZero();
        A(0,0) = -1; A(0,2) = 0.0193;
        A(1,1) = -1; A(1,2) = 0.00954;
        A(2,5) = 1;
        A(3,3) = 1;

        DenseVector b;
        b.setZero(4);
        b(2) = -750.1728/3.14159265359;
        b(3) = 240;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    SolverResult res = bnb.optimize();
    timer.stop();
    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;
    cout << "Time: " << timer.getMicroSeconds() << " (us)" << endl;

    //1 0.625 47.5 90 6395.50783 -345958.333

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP06::validateResult()
{
    return std::abs(fopt_found - fopt_known) <= 1e-3;
}

} // namespace CENSO
