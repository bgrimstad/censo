/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop12.h"
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

POP12::POP12()
{
    testName = "P12";
    fopt_known = -450;
    fopt_known = -4.5; // for scaled problem P12_3()
    fopt_found = INF;
}

void POP12::runProblem()
{
    //P12_1();
    //P12_2();
    P12_3(); // Scaled version which solves faster
}

bool POP12::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

void POP12::P12_1()
{
    // Problem 12 (Nataraj), Bilinear problem
    // Six 2-D B-splines
    // cout << "\n\nSolving problem P12..." << endl;

    int dim = 7+6; // x1,..,x7,l1,...,l6

    std::vector<double> costs = {0,0,0,-9,-15,1,-5,  6,16,15,6,16,15};
    std::vector<double> lb = {0,0,0,0,0,0,0,            0,0,0,0,0,0};
    std::vector<double> ub = {1,1,1,100,200,100,200,    100,100,100,200,200,200};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // x_i * x_4 = li
    for (int i = 0; i < 3; i++)
    {
        std::vector<VariablePtr> cvars = {
            vars.at(i),
            vars.at(3),
            vars.at(i+7)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1};

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x_i * x_5 = li
    for (int i = 0; i < 3; i++)
    {
        std::vector<VariablePtr> cvars = {
            vars.at(i),
            vars.at(4),
            vars.at(i+7+3)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1};

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints
        DenseMatrix A = DenseMatrix::Zero(5,dim);
        A(0,9) = 1; A(0,12) = 1;
        A(1,3) = 1; A(1,5) = 1;
        A(2,4) = 1; A(2,6) = 1;
        A(3,7) = 3; A(3,8) = 1; A(3,9) = 1; A(3,3) = -2.5; A(3,5) = -0.5;
        A(4,10) = 3; A(4,11) = 1; A(4,12) = 1; A(4,4) = -1.5; A(4,6) = 0.5;

        DenseVector b;
        b.setZero(5);
        b(0) = 50;
        b(1) = 100;
        b(2) = 200;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false); // All variables

        cs->add(lincon);
    }

    { // Linear equality constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2)
        };

        DenseMatrix A = DenseMatrix::Zero(1,3);
        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

        DenseVector b(1);
        b(0) = 1;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

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

void POP12::P12_2()
{
    // Problem 12 (Nataraj), Bilinear problem
    // cout << "\n\nSolving problem P12..." << endl;

    int dim = 7+5; // x1,..,x7,l1,...,l5

    std::vector<double> costs = {0,0,0,0,0,1,-5,        -1,-1,0,0,0};
    std::vector<double> lb = {0,0,0,0,0,0,0,            -INF,-INF,-INF,-INF,-INF};
    std::vector<double> ub = {1,1,1,100,200,100,200,    INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // x4*(9 - 6*x1 - 16*x2 - 15*x3) = l1
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(7)
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
        DenseVector c(4);
        c.setZero();
        c(0) = 9;
        c(1) = -6;
        c(2) = -16;
        c(3) = -15;

        // Poly exponents
        DenseMatrix E(4,4); E.setZero();
        E(0,3) = 1;
        E(1,0) = 1; E(1,3) = 1;
        E(2,1) = 1; E(2,3) = 1;
        E(3,2) = 1; E(3,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5*(15 - 6*x1 - 16*x2 - 15*x3) = l2
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(4),
            vars.at(8)
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
        DenseVector c(4); c.setZero();
        c(0) = 15;
        c(1) = -6;
        c(2) = -16;
        c(3) = -15;

        // Poly exponents
        DenseMatrix E(4,4); E.setZero();
        E(0,3) = 1;
        E(1,0) = 1; E(1,3) = 1;
        E(2,1) = 1; E(2,3) = 1;
        E(3,2) = 1; E(3,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x3*(x4 + x5) = l3
    {
        std::vector<VariablePtr> cvars = {
            vars.at(2),
            vars.at(3),
            vars.at(4),
            vars.at(9)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,1};

        // Poly coeffs
        DenseVector c(2); c.setZero();
        c(0) = 1;
        c(1) = 1;

        // Poly exponents
        DenseMatrix E(2,3); E.setZero();
        E(0,0) = 1; E(0,1) = 1;
        E(1,0) = 1; E(1,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x4*(3*x1 + x2 + x3 - 2.5) = l4
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(10)
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
        DenseVector c(4); c.setZero();
        c(0) = 3;
        c(1) = 1;
        c(2) = 1;
        c(3) = -2.5;

        // Poly exponents
        DenseMatrix E(4,4); E.setZero();
        E(0,0) = 1; E(0,3) = 1;
        E(1,1) = 1; E(1,3) = 1;
        E(2,2) = 1; E(2,3) = 1;
        E(3,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5*(3*x1 + x2 + x3 - 1.5) = l4
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(4),
            vars.at(11)
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
        DenseVector c(4); c.setZero();
        c(0) = 3;
        c(1) = 1;
        c(2) = 1;
        c(3) = -1.5;

        // Poly exponents
        DenseMatrix E(4,4); E.setZero();
        E(0,0) = 1; E(0,3) = 1;
        E(1,1) = 1; E(1,3) = 1;
        E(2,2) = 1; E(2,3) = 1;
        E(3,3) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints

        DenseMatrix A = DenseMatrix::Zero(5,dim);
        A(0,9) = 1;
        A(1,3) = 1; A(1,5) = 1;
        A(2,4) = 1; A(2,6) = 1;
        A(3,10) = 1; A(3,5) = -0.5;
        A(4,11) = 1; A(4,6) = 0.5;

        DenseVector b;
        b.setZero(5);
        b(0) = 50;
        b(1) = 100;
        b(2) = 200;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false); // All variables

        cs->add(lincon);
    }

    { // Linear equality constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2)
        };

        DenseMatrix A = DenseMatrix::Zero(1,3);
        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

        DenseVector b(1);
        b(0) = 1;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

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

void POP12::P12_3()
{
    // Problem 12 (Nataraj), Bilinear problem
    // Six 2-D B-splines
    cout << "\n\nSolving problem P12..." << endl;

    int dim = 7+6; // x1,..,x7,l1,...,l6

    std::vector<double> costs = {0,0,0,-9,-15*2,1,-5*2,6,16,15,6*2,16*2,15*2};
    std::vector<double> lb = {0,0,0,0,0,0,0,            0,0,0,0,0,0};
    std::vector<double> ub = {1,1,1,1,1,1,1,    1,1,1,1,1,1};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // x_i * x_4 = li
    for (int i = 0; i < 3; i++)
    {
        std::vector<VariablePtr> cvars = {
            vars.at(i),
            vars.at(3),
            vars.at(i+7)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1};

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x_i * x_5 = li
    for (int i = 0; i < 3; i++)
    {
        std::vector<VariablePtr> cvars = {
            vars.at(i),
            vars.at(4),
            vars.at(i+7+3)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1};

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(3) = 1;

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints
        DenseMatrix A = DenseMatrix::Zero(5,dim);
        A(0,9) = 0.5; A(0,12) = 1;
        A(1,3) = 1; A(1,5) = 1;
        A(2,4) = 1; A(2,6) = 1;
        A(3,7) = 3; A(3,8) = 1; A(3,9) = 1; A(3,3) = -2.5; A(3,5) = -0.5;
        A(4,10) = 3; A(4,11) = 1; A(4,12) = 1; A(4,4) = -1.5; A(4,6) = 0.5;

        DenseVector b;
        b.setZero(5);
        b(0) = 50/200.0;
        b(1) = 1;
        b(2) = 1;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false); // All variables

        cs->add(lincon);
    }

    { // Linear equality constraints

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2)
        };

        DenseMatrix A = DenseMatrix::Zero(1,3);
        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;

        DenseVector b(1);
        b(0) = 1;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(cvars, A, b, true);

        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

} // namespace CENSO
