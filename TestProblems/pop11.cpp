/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop11.h"
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

POP11::POP11()
{
    testName = "P11";
    fopt_known = -310;
    fopt_found = INF;
}

void POP11::runProblem()
{
    // Problem 11 (Nataraj) - Test problem 3, Ch. 3.3 (Floudas)
    // cout << "\n\nSolving problem P11..." << endl;

    int dim = 14;

    // x1,...,x6,l1,...,l8
    std::vector<double> costs = {0,0,0,0,0,0,-25,-1,-1,-1,-1,-1,0,0};
    std::vector<double> lb = {0,0,1,0,1,0,0,0,0,0,0,0,0,0};
    std::vector<double> ub = {5,5,5,6,5,10,INF,INF,INF,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // (x1-2)^2

        std::vector<VariablePtr> cvars = {vars.at(0),vars.at(6)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 4;
        c(1) = -4;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x2-2)^2

        std::vector<VariablePtr> cvars = {vars.at(1), vars.at(7)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 4;
        c(1) = -4;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x3-1)^2

        std::vector<VariablePtr> cvars = {vars.at(2), vars.at(8)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 1;
        c(1) = -2;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x4-4)^2

        std::vector<VariablePtr> cvars = {vars.at(3), vars.at(9)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 16;
        c(1) = -8;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x5-1)^2

        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(10)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 1;
        c(1) = -2;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x6-4)^2

        std::vector<VariablePtr> cvars = {vars.at(5), vars.at(11)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 16;
        c(1) = -8;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x3-3)^2
        std::vector<VariablePtr> cvars = {vars.at(2), vars.at(12)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 9;
        c(1) = -6;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x5-3)^2

        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(13)};

        std::vector<double> thislb = {cvars.at(0)->getLowerBound()};
        std::vector<double> thisub = {cvars.at(0)->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 9;
        c(1) = -6;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

//        BSpline bs(coeffs, knots, deg);
        BSpline bs = BSplineWrap::build_bspline(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // l7 + x4 >= 4

        std::vector<VariablePtr> cvars = {vars.at(12), vars.at(3)};

        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

        DenseVector b(1); b(0) = -4;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    { // l8 + x6 >= 4

        std::vector<VariablePtr> cvars = {vars.at(13), vars.at(5)};

        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

        DenseVector b(1); b(0) = -4;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    { // x1 - 3x2 <= 2, -x1 + x2 <= 2, x1 + x2 <= 6, x1 + x2 >= 2

        std::vector<VariablePtr> cvars = {vars.at(0), vars.at(1)};

        DenseMatrix A(4,2);
        A(0,0) = 1; A(0,1) = -3;
        A(1,0) = -1; A(1,1) = 1;
        A(2,0) = 1; A(2,1) = 1;
        A(3,0) = -1; A(3,1) = -1;

        DenseVector b(4);
        b(0) = 2;
        b(1) = 2;
        b(2) = 6;
        b(3) = -2;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP11::validateResult()
{
    return std::abs(fopt_found - fopt_known) <= 1e-3;
}

} // namespace CENSO
