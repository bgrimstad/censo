/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop08.h"
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

POP08::POP08()
{
    testName = "P08";
    fopt_known = 42.444;
    fopt_found = INF;
}

void POP08::runProblem()
{
    // Problem 8 (Nataraj)
    // Three 4-D B-splines
    // Runs fast when refinement is turned off!
    // A lot of iterations when using bounding box relaxation (as expected)
    // cout << "\n\nSolving problem P08..." << endl;

    int dim = 4+3;

    // x1,x2,x3,x4,l1,l2,l3
    std::vector<double> costs = {0,0,0,0,0,0,1};
    std::vector<double> lb = {3,2,0.125,0.25,-INF,-INF,-INF};
    std::vector<double> ub = {20,15,0.75,1.25,INF,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // f(x) = l3

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(6)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,1,1};

        // Poly coeffs
        double a = 27.264;
        DenseVector c(3); c.setZero();
        c(0) = 2*a;
        c(1) = a;
        c(2) = -2*a;

        // Poly exponents
        DenseMatrix E(3,4); E.setZero();
        E(0,0) = 0; E(0,1) = 1; E(0,2) = 0; E(0,3) = 1;
        E(1,0) = 1; E(1,1) = 0; E(1,2) = 1; E(1,3) = 0;
        E(2,0) = 0; E(2,1) = 0; E(2,2) = 1; E(2,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // I1(x) = l1

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(4)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3,1,3,1};

        // Poly coeffs
        DenseVector c(7); c.setZero();
        c(0) = 6;
        c(1) = -12;
        c(2) = 8;
        c(3) = 1;
        c(4) = -6;
        c(5) = 12;
        c(6) = -8;

        // Poly exponents
        DenseMatrix E(7,4); E.setZero();
        E(0,0) = 2; E(0,1) = 1; E(0,2) = 1; E(0,3) = 0;
        E(1,0) = 1; E(1,1) = 1; E(1,2) = 2; E(1,3) = 0;
        E(2,0) = 0; E(2,1) = 1; E(2,2) = 3; E(2,3) = 0;
        E(3,0) = 3; E(3,1) = 0; E(3,2) = 0; E(3,3) = 1;
        E(4,0) = 2; E(4,1) = 0; E(4,2) = 1; E(4,3) = 1;
        E(5,0) = 1; E(5,1) = 0; E(5,2) = 2; E(5,3) = 1;
        E(6,0) = 0; E(6,1) = 0; E(6,2) = 3; E(6,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // I2(x) = l2

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(5)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3,1,4,2};

        // Poly coeffs
        double a = -3.5;
        DenseVector c(12); c.setZero();
        c(0) = 6*a;
        c(1) = -12*a;
        c(2) = 8*a;
        c(3) = 1*a;
        c(4) = -6*a;
        c(5) = 12*a;
        c(6) = -8*a;

        c(7) = 1;
        c(8) = -1;
        c(9) = 1;
        c(10) = 1;
        c(11) = -2;

        // Poly exponents
        DenseMatrix E(12,4); E.setZero();
        E(0,0) = 2; E(0,1) = 1; E(0,2) = 2; E(0,3) = 0;
        E(1,0) = 1; E(1,1) = 1; E(1,2) = 3; E(1,3) = 0;
        E(2,0) = 0; E(2,1) = 1; E(2,2) = 4; E(2,3) = 0;
        E(3,0) = 3; E(3,1) = 0; E(3,2) = 1; E(3,3) = 1;
        E(4,0) = 2; E(4,1) = 0; E(4,2) = 2; E(4,3) = 1;
        E(5,0) = 1; E(5,1) = 0; E(5,2) = 3; E(5,3) = 1;
        E(6,0) = 0; E(6,1) = 0; E(6,2) = 4; E(6,3) = 1;

        E(7,0) = 1;  E(7,1) = 1;  E(7,2) = 0;  E(7,3) = 1;
        E(8,0) = 0;  E(8,1) = 1;  E(8,2) = 0;  E(8,3) = 2;
        E(9,0) = 2;  E(9,1) = 0;  E(9,2) = 1;  E(9,3) = 0;
        E(10,0) = 0; E(10,1) = 0; E(10,2) = 1; E(10,3) = 2;
        E(11,0) = 1; E(11,1) = 0; E(11,2) = 1; E(11,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints of auxiliary variables

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(4),
            vars.at(5)
        };

        DenseMatrix A = DenseMatrix::Zero(7,6);
        A(0,4) = -1;
        A(1,0) = 8; A(1,4) = -1;
        A(2,5) = 1;
        A(3,0) = 1; A(3,1) = -3;
        A(4,1) = 2; A(4,0) = -1;
        A(5,2) = 1; A(5,3) = -1.5;
        A(6,3) = 0.5; A(6,2) = -1;

        DenseVector b; b.setZero(7);
        b(0) = -61.01627586;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);

        cs->add(lincon);
    }

//    DenseVector x(7);
//    x(0) = 3.76117027;
//    x(1) = 2;
//    x(2) = 0.125;
//    x(3) = 0.25;
//    x(4) = 30.662118;
//    x(5) = -10.1180522;
//    x(6) = 38.3780683;

//    bool feas = cs->checkFeasibility(x, 1e-5);
//    if (feas) cout << "YAAAI!" << endl;
//    else cout << "NEEEY!" << endl;
//    exit(1);

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

bool POP08::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
