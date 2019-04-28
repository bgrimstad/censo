/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop07.h"
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

POP07::POP07()
{
    testName = "P07";
    fopt_known = 1.0899;
    fopt_found = INF;
}

void POP07::runProblem()
{
    P07_1();
}

bool POP07::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

void POP07::P07_1()
{
    // Problem 7 (Nataraj)
    // Two 2-D B-splines
    cout << "\n\nSolving problem P07..." << endl;

    int dim = 4+2;

    // x1,x2,x3,x4,l1,l2
    std::vector<double> costs = {0,0,0,1,0,0};
    std::vector<double> lb = {0,0,0,0,-INF,-INF};
    std::vector<double> ub = {5,5,5,5,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // x1^4*x2^4 - x1^4 = l1

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(4)};

        std::vector<double> thislb = {
            cvars.at(0)->getLowerBound(),
            cvars.at(1)->getLowerBound()
        };

        std::vector<double> thisub = {
            cvars.at(0)->getUpperBound(),
            cvars.at(1)->getUpperBound()
        };

        std::vector<unsigned int> deg = {4,4};

        DenseVector c(5*5); c.setZero();
        c(20) = -1;
        c(24) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // x2^4*x3 = l1
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
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

        std::vector<unsigned int> deg = {4,1};

        DenseVector c(10); c.setZero();
        c(9) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints of auxiliary variables
        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(5)};

        DenseMatrix A(1,2);
        A(0,0) = 1; A(0,1) = -1;

        DenseVector b; b.setZero(1);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

        cs->add(lincon);
    }

    { // Linear constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3)
        };

        DenseMatrix A(6,4);
        A.setZero();
        A(0,0) = -1; A(0,3) = -0.25;
        A(1,0) = 1; A(1,3) = -0.25;
        A(2,1) = -1; A(2,3) = -0.2;
        A(3,1) = 1; A(3,3) = -0.2;
        A(4,2) = -1; A(4,3) = -0.2;
        A(5,2) = 1; A(5,3) = -0.2;

        DenseVector b; b.setZero(6);
        b(0) = -1.4;
        b(1) = 1.4;
        b(2) = -1.5;
        b(3) = 1.5;
        b(4) = -0.8;
        b(5) = 0.8;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

void POP07::P07_2()
{
    // Problem 7 (Nataraj)
    // One 3-D B-spline

    int dim = 4+1;

    // x1,x2,x3,x4,l1,l2
    std::vector<double> costs = {0,0,0,1,0};
    std::vector<double> lb = {0,0,0,0,-INF};
    std::vector<double> ub = {5,5,5,5,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // x1^4*x2^4 - x1^4 = l1
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(4)
        };

        std::vector<double> thislb = {
            cvars.at(0)->getLowerBound(),
            cvars.at(1)->getLowerBound(),
            cvars.at(2)->getLowerBound()
        };

        std::vector<double> thisub = {
            cvars.at(0)->getUpperBound(),
            cvars.at(1)->getUpperBound(),
            cvars.at(2)->getUpperBound()
        };

        std::vector<unsigned int> deg = {4,4,1};

        // Poly coeffs
        DenseVector c(3); c.setZero();
        c(0) = 1;
        c(1) = -1;
        c(2) = -1;

        // Poly exponents
        DenseMatrix E(3,3); E.setZero();
        E(0,0) = 4; E(0,1) = 4; E(0,2) = 0;
        E(1,0) = 4; E(1,1) = 0; E(1,2) = 0;
        E(2,0) = 0; E(2,1) = 4; E(2,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints of auxiliary variables

        std::vector<VariablePtr> cvars = {vars.at(4)};

        DenseMatrix A(1,1);
        A(0,0) = 1;

        DenseVector b; b.setZero(1);

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

        cs->add(lincon);
    }

    { // Linear constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3)
        };

        DenseMatrix A(6,4);
        A.setZero();
        A(0,0) = -1; A(0,3) = -0.25;
        A(1,0) = 1; A(1,3) = -0.25;
        A(2,1) = -1; A(2,3) = -0.2;
        A(3,1) = 1; A(3,3) = -0.2;
        A(4,2) = -1; A(4,3) = -0.2;
        A(5,2) = 1; A(5,3) = -0.2;

        DenseVector b; b.setZero(6);
        b(0) = -1.4;
        b(1) = 1.4;
        b(2) = -1.5;
        b(3) = 1.5;
        b(4) = -0.8;
        b(5) = 0.8;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

} // namespace CENSO
