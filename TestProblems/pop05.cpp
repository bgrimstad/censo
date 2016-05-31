/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop05.h"
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

POP05::POP05()
{
    testName = "P05";
    fopt_known = 0;
    fopt_found = INF;
}

void POP05::runProblem()
{
    // Problem 5 (Nataraj), Himmelblau root finding problem
    // cout << "\n\nSolving problem P05..." << endl;

    int dim = 3+2;

    // x1,x2,x3,l1,l2
    std::vector<double> costs = {0, 0, 1, 0, 0};
    std::vector<double> lb = {-5,-5,-5,-INF,-INF};
    std::vector<double> ub = {5,5,5,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // f1(x1,x2) = 2*x1^2 + 4*x1*x2 - 42*x1 +4*x1^3 = l1

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(3)
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

        DenseVector c(8);
        c.setZero();
        c(2) = -42; // x1
        c(3) = 4; // x1*x2
        c(4) = 2; // x1^2
        c(6) = 4; // x1^3

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // f2(x1,x2) = 2*x1^2 + 4*x1*x2 - 26*x2 +4*x2^3 = l2

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(4)
        };

        std::vector<double> thislb = {
            cvars.at(0)->getLowerBound(),
            cvars.at(1)->getLowerBound()
        };

        std::vector<double> thisub = {
            cvars.at(0)->getUpperBound(),
            cvars.at(1)->getUpperBound()
        };

        std::vector<unsigned int> deg = {2,3};

        DenseVector c(12);
        c.setZero();
        c(1) = -26; // x2
        c(3) = 4; // x2^3
        c(5) = 4; // x1*x2
        c(8) = 2; // x1^2

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // Linear constraints

        DenseMatrix A(4,5);
        A.setZero();
        A(0,2) = -1; A(0,3) = 1;
        A(1,2) = -1; A(1,3) = -1;
        A(2,2) = -1; A(2,4) = 1;
        A(3,2) = -1; A(3,4) = -1;

        DenseVector b;
        b.setZero(4);
        b(0) = 14;
        b(1) = -14;
        b(2) = 22;
        b(3) = -22;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    Timer timer;
    timer.start();
    SolverResult res = bnb.optimize();
    timer.stop();
    cout << "Time: " << timer.getMilliSeconds() << " (ms)" << endl;

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP05::validateResult()
{
    return std::abs(fopt_found - fopt_known) <= 1e-3;
}

} // namespace CENSO
