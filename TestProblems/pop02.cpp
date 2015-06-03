/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop02.h"
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

POP02::POP02()
{
    testName = "P02";
    fopt_known = -6961.815;
    fopt_found = INF;
}

void POP02::runProblem()
{
    // Problem 2 (Nataraj)
    // cout << "\n\nSolving problem P02..." << endl;

    int dim = 7;

    // x1,x2,l1,...,l5
    std::vector<double> costs = {0, 0, 1, 1, 0, 0, 0};
    std::vector<double> lb = {13,0,-INF,-INF,0,0,0};
    std::vector<double> ub = {100,100,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // (x1-10)^3
        VariablePtr var = vars.at(0);
        std::vector<double> thislb = {var->getLowerBound()};
        std::vector<double> thisub = {var->getUpperBound()};

        std::vector<unsigned int> deg = {3};

        DenseVector c(4); c.setZero();
        c(0) = -1000;
        c(1) = 300;
        c(2) = -30;
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        std::vector<VariablePtr> cvars = {var,vars.at(2)};
        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // (x2-20)^3
        VariablePtr var = vars.at(1);
        std::vector<double> thislb = {var->getLowerBound()};
        std::vector<double> thisub = {var->getUpperBound()};

        std::vector<unsigned int> deg = {3};

        DenseVector c(4); c.setZero();
        c(0) = -8000;
        c(1) = 1200;
        c(2) = -60;
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        std::vector<VariablePtr> cvars = {var,vars.at(3)};
        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // (x1-5)^2
        VariablePtr var = vars.at(0);
        std::vector<double> thislb = {var->getLowerBound()};
        std::vector<double> thisub = {var->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 25;
        c(1) = -10;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);


        std::vector<VariablePtr> cvars = {var,vars.at(4)};
        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // (x2-5)^2
        VariablePtr var = vars.at(1);
        std::vector<double> thislb = {var->getLowerBound()};
        std::vector<double> thisub = {var->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 25;
        c(1) = -10;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        std::vector<VariablePtr> cvars = {var, vars.at(5)};
        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // (x1-6)^2
        VariablePtr var = vars.at(0);
        std::vector<double> thislb = {var->getLowerBound()};
        std::vector<double> thisub = {var->getUpperBound()};

        std::vector<unsigned int> deg = {2};

        DenseVector c(3); c.setZero();
        c(0) = 36;
        c(1) = -12;
        c(2) = 1;

        DenseMatrix T = getTransformationMatrix(deg,thislb,thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        std::vector<VariablePtr> cvars = {var,vars.at(6)};
        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);
        cs->add(cbs);
    }

    { // -l3 - l4 <= -100

        DenseMatrix A(1,2); A(0,0) = -1; A(0,1) = -1;

        DenseVector b(1); b(0) = -100;

        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(5)};
        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
        cs->add(lincon);
    }

    { // l5 + l4 <= 82.81

        DenseMatrix A(1,2); A(0,0) = 1; A(0,1) = 1;

        DenseVector b(1); b(0) = 82.81;

        std::vector<VariablePtr> cvars = {vars.at(6), vars.at(5)};
        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, false);
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

bool POP02::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-2)
        return true;
    return false;
}

} // namespace CENSO
