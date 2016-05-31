/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop13.h"
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

POP13::POP13()
{
    testName = "P13";
    fopt_known = 2994.42447;
    fopt_found = INF;
}

void POP13::runProblem()
{
    P13_3();
}

bool POP13::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-2)
        return true;
    return false;
}

void POP13::P13_1()
{
    // Problem 13 (Nataraj)
    // cout << "\n\nSolving problem P13..." << endl;

    int dim = 7+8; // x1,..,x7,l1,...,l5

    std::vector<double> costs(dim,0.0); costs.at(7) = 1;
    std::vector<double> lb = {2.6,0.7,17,7.3,7.3,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
    std::vector<double> ub = {3.6,0.8,28,8.3,8.3,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    // Testing manual bounds tightening
    lb.at(0) = 3.5;
    ub.at(1) = 0.72;
    lb.at(4) = 7.4;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // obj = l1
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(4),
            vars.at(5),
            vars.at(6),
            vars.at(7)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2,1,1,3,3};

        double a1 = 0.7854;
        double a2 = 3.3333;
        double a3 = 14.9334;
        double a4 = -43.0934;
        double a5 = -1.508;
        double a6 = 7.477;
        double a7 = 0.7854;

        // Poly coeffs
        DenseVector c(9);
        c.setZero();
        c(0) = a1*a2;
        c(1) = a1*a3;
        c(2) = a1*a4;
        c(3) = a5;
        c(4) = a5;
        c(5) = a6;
        c(6) = a6;
        c(7) = a7;
        c(8) = a7;

        // Poly exponents
        DenseMatrix E(9,7);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
        E(2,0) = 1; E(2,1) = 2;
        E(3,0) = 1; E(3,5) = 2;
        E(4,0) = 1; E(4,6) = 2;
        E(5,5) = 3;
        E(6,6) = 3;
        E(7,3) = 1; E(7,5) = 2;
        E(8,4) = 1; E(8,6) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3 = l2
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(8)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,1};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3*x3 = l3
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(9)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x6^4*x3 - 1.93*x4^3 = l4
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(5),
            vars.at(10)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,3,4};

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = 1;
        c(1) = -1.93;

        // Poly exponents
        DenseMatrix E(2,4);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,3) = 4;
        E(1,2) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x7^4*x3 - 1.93*x5^3 = l5
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(4),
            vars.at(6),
            vars.at(11)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,3,4};

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = 1;
        c(1) = -1.93;

        // Poly exponents
        DenseMatrix E(2,4);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,3) = 4;
        E(1,2) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // g(x2,x3,x4,x6) = l6
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(3),
            vars.at(5),
            vars.at(12)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,2,6};

        // Poly coeffs
        DenseVector c(3);
        c.setZero();
        c(0) = 745*745/16.911e6;
        c(1) = 1;
        c(2) = -110*110/16.911e6;

        // Poly exponents
        DenseMatrix E(3,4);
        E.setZero();
        E(0,2) = 2;
        E(1,0) = 2; E(1,1) = 2;
        E(2,0) = 2; E(2,1) = 2; E(2,3) = 6;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // g(x2,x3,x5,x7) = l7
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(4),
            vars.at(6),
            vars.at(13)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,2,6};

        // Poly coeffs
        DenseVector c(3);
        c.setZero();
        c(0) = 745*745/157.50e6;
        c(1) = 1;
        c(2) = -85*85/157.50e6;

        // Poly exponents
        DenseMatrix E(3,4);
        E.setZero();
        E(0,2) = 2;
        E(1,0) = 2; E(1,1) = 2;
        E(2,0) = 2; E(2,1) = 2; E(2,3) = 6;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3 = l8
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(14)
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

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints
        DenseMatrix A = DenseMatrix::Zero(11,dim);
        A(0,8) = -1;
        A(1,9) = -1;
        A(2,10) = -1;
        A(3,11) = -1;
        A(4,12) = 1;
        A(5,13) = 1;
        A(6,14) = 1;
        A(7,0) = -1; A(7,1) = 5;
        A(8,0) = 1; A(8,1) = -12;
        A(9,3) = -1; A(9,5) = 1.5;
        A(10,4) = -1; A(10,6) = 1.1;

        DenseVector b;
        b.setZero(11);
        b(0) = -27;
        b(1) = -397.5;
        b(6) = 40;
        b(9) = -1.9;
        b(10) = -1.9;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false); // All variables

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

void POP13::P13_2()
{
    // Problem 13 (Nataraj)
    // cout << "\n\nSolving problem P13..." << endl;

    int dim = 7+15; // x1,..,x7,l1,...,l15

    std::vector<double> costs(dim, 0.0);
    costs.at(7) = 1; costs.at(8) = 1; costs.at(9) = 1; costs.at(10) = 1;
    std::vector<double> lb = {2.6,0.7,17,7.3,7.3,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
    std::vector<double> ub = {3.6,0.8,28,8.3,8.3,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    // Testing manual bounds tightening
    lb.at(0) = 3.5;
    ub.at(1) = 0.72;
    lb.at(4) = 7.4;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // obj = l1
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(7)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2};

        double a1 = 0.7854;
        double a2 = 3.3333;
        double a3 = 14.9334;
        double a4 = -43.0934;

        // Poly coeffs
        DenseVector c(3);
        c.setZero();
        c(0) = a1*a2;
        c(1) = a1*a3;
        c(2) = a1*a4;

        // Poly exponents
        DenseMatrix E(3,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
        E(2,0) = 1; E(2,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l2
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(5),
            vars.at(6),
            vars.at(8)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,3,3};

        double a1 = -1.508;
        double a2 = 7.477;

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(0) = a1;
        c(1) = a1;
        c(2) = a2;
        c(3) = a2;

        // Poly exponents
        DenseMatrix E(4,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;
        E(1,0) = 1; E(1,2) = 2;
        E(2,1) = 3;
        E(3,2) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l3
    {
        std::vector<VariablePtr> cvars = {
            vars.at(3),
            vars.at(5),
            vars.at(9)
        };

        std::vector<double> thislb;
        std::vector<double> thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2};

        double a1 = 0.7854;

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = a1;

        // Poly exponents
        DenseMatrix E(1,2);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l4
    {
        std::vector<VariablePtr> cvars = {
            vars.at(4),
            vars.at(6),
            vars.at(10)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2};

        double a1 = 0.7854;

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = a1;

        // Poly exponents
        DenseMatrix E(1,2);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3 = l5
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(11)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,1};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3*x3 = l6
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(12)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3*x6^4 = l7
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(5),
            vars.at(13)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,4};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x4^3 = l8
    {
        std::vector<VariablePtr> cvars = {
            vars.at(3),
            vars.at(14)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3*x7^4 = l9
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(6),
            vars.at(15)};

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,4};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5^3 = l10
    {
        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(16)};

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x4^2 = l11
    {
        std::vector<VariablePtr> cvars = {vars.at(3), vars.at(17)};

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // l12
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(5),
            vars.at(18)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,6};

        double a1 = 12100/(745.0*745.0);
        double a2 = -16.91e6/(745.0*745.0);

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = a1;
        c(1) = a2;

        // Poly exponents
        DenseMatrix E(2,3);
        E.setZero();
        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
        E(1,0) = 2; E(1,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5^2 = l13
    {
        std::vector<VariablePtr> cvars = {vars.at(4), vars.at(19)};

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // l14
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(6),
            vars.at(20)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,6};

        double a1 = 7225/(745.0*745.0);
        double a2 = -157.5e6/(745.0*745.0);

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = a1;
        c(1) = a2;

        // Poly exponents
        DenseMatrix E(2,3);
        E.setZero();
        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
        E(1,0) = 2; E(1,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3 = l15
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(21)
        };

        std::vector<double> thislb, thisub;
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

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints
        DenseMatrix A = DenseMatrix::Zero(11,dim);
        A(0,11) = -1;
        A(1,12) = -1;
        A(2,13) = -1; A(2,14) = 1.93;
        A(3,15) = -1; A(3,16) = 1.93;
        A(4,17) = 1; A(4,18) = -1;
        A(5,19) = 1; A(5,20) = -1;
        A(6,21) = 1;
        A(7,0) = -1; A(7,1) = 5;
        A(8,0) = 1; A(8,1) = -12;
        A(9,3) = -1; A(9,5) = 1.5;
        A(10,4) = -1; A(10,6) = 1.1;

        DenseVector b;
        b.setZero(11);
        b(0) = -27;
        b(1) = -397.5;
        b(6) = 40;
        b(9) = -1.9;
        b(10) = -1.9;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);// All variables

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

void POP13::P13_3()
{
    // Problem 13 (Nataraj)
    // Exploiting additive separability and scaling (faster, scaling is very important here)
    // cout << "\n\nSolving problem P13..." << endl;

    int dim = 7+15; // x1,..,x7,l1,...,l15

    double x3s = 28.0;
    double x4s = 8.3;
    double x5s = 8.3;
    std::vector<double> costs(dim, 0.0);
    costs.at(7) = 1; costs.at(8) = 1; costs.at(9) = 1; costs.at(10) = 1;
    std::vector<double> lb = {2.6,0.7,17/x3s,7.3/x4s,7.3/x5s,2.9,5,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF,-INF};
    std::vector<double> ub = {3.6,0.8,28/x3s,8.3/x4s,8.3/x5s,3.9,5.5,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF,INF};
    std::vector<double> z0(dim,0);

    // Testing manual bounds tightening
    lb.at(0) = 3.5;
    ub.at(1) = 0.72;
    lb.at(4) = 7.4/x5s;

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    // obj = l1
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(7)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2};

        double a1 = 0.7854;
        double a2 = 3.3333;
        double a3 = 14.9334;
        double a4 = -43.0934;

        // Poly coeffs
        DenseVector c(3);
        c.setZero();
        c(0) = a1*a2*x3s*x3s;
        c(1) = a1*a3*x3s;
        c(2) = a1*a4;

        // Poly exponents
        DenseMatrix E(3,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;
        E(1,0) = 1; E(1,1) = 2; E(1,2) = 1;
        E(2,0) = 1; E(2,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l2
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(5),
            vars.at(6),
            vars.at(8)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,3,3};

        double a1 = -1.508;
        double a2 = 7.477;

        // Poly coeffs
        DenseVector c(4);
        c.setZero();
        c(0) = a1;
        c(1) = a1;
        c(2) = a2;
        c(3) = a2;

        // Poly exponents
        DenseMatrix E(4,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;
        E(1,0) = 1; E(1,2) = 2;
        E(2,1) = 3;
        E(3,2) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l3
    {
        std::vector<VariablePtr> cvars = {
            vars.at(3),
            vars.at(5),
            vars.at(9)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2};

        double a1 = 0.7854;

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = a1*x4s;

        // Poly exponents
        DenseMatrix E(1,2);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // obj = l4
    {
        std::vector<VariablePtr> cvars = {
            vars.at(4),
            vars.at(6),
            vars.at(10)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2};

        double a1 = 0.7854*x5s;

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = a1;

        // Poly exponents
        DenseMatrix E(1,2);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3 = l5
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(11)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,1};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 1;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x1*x2*x2*x3*x3 = l6
    {
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(12)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,2,2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 2; E(0,2) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3*x6^4 = l7
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(5),
            vars.at(13)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,4};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x4^3 = l8
    {
        std::vector<VariablePtr> cvars = {
            vars.at(3),
            vars.at(14)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1*x4s*x4s*x4s;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3*x7^4 = l9
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(6),
            vars.at(15)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {1,1,4};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1;

        // Poly exponents
        DenseMatrix E(1,3);
        E.setZero();
        E(0,0) = 1; E(0,1) = 1; E(0,2) = 4;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5^3 = l10
    {
        std::vector<VariablePtr> cvars = {
            vars.at(4),
            vars.at(16)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {3};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1*x5s*x5s;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 3;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x4^2 = l11
    {
        std::vector<VariablePtr> cvars = {
            vars.at(3),
            vars.at(17)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1*x4s*x4s;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // l12
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(5),
            vars.at(18)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,6};

        double a1 = 12100/(745.0*745.0)*x3s*x3s;
        double a2 = -16.91e6/(745.0*745.0)*x3s*x3s;

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = a1;
        c(1) = a2;

        // Poly exponents
        DenseMatrix E(2,3);
        E.setZero();
        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
        E(1,0) = 2; E(1,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x5^2 = l13
    {
        std::vector<VariablePtr> cvars = {
            vars.at(4),
            vars.at(19)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2};

        // Poly coeffs
        DenseVector c(1);
        c.setZero();
        c(0) = 1*x5s*x5s;

        // Poly exponents
        DenseMatrix E(1,1);
        E.setZero();
        E(0,0) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // l14
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(6),
            vars.at(20)
        };

        std::vector<double> thislb, thisub;
        for (unsigned int i = 0; i < cvars.size()-1; i++)
        {
            thislb.push_back(cvars.at(i)->getLowerBound());
            thisub.push_back(cvars.at(i)->getUpperBound());
        }

        std::vector<unsigned int> deg = {2,2,6};

        double a1 = 7225/(745.0*745.0)*x3s*x3s;
        double a2 = -157.5e6/(745.0*745.0)*x3s*x3s;

        // Poly coeffs
        DenseVector c(2);
        c.setZero();
        c(0) = a1;
        c(1) = a2;

        // Poly exponents
        DenseMatrix E(2,3);
        E.setZero();
        E(0,0) = 2; E(0,1) = 2; E(0,2) = 6;
        E(1,0) = 2; E(1,1) = 2;

        DenseMatrix coeffs = getBSplineBasisCoefficients(c, E, thislb, thisub);

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    // x2*x3 = l15
    {
        std::vector<VariablePtr> cvars = {
            vars.at(1),
            vars.at(2),
            vars.at(21)
        };

        std::vector<double> thislb, thisub;
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

        BSpline bs(coeffs, knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear inequality constraints
        DenseMatrix A = DenseMatrix::Zero(11,dim);
        A(0,11) = -1;
        A(1,12) = -1;
        A(2,13) = -1*x3s; A(2,14) = 1.93;
        A(3,15) = -1*x3s; A(3,16) = 1.93;
        A(4,17) = 1; A(4,18) = -1;
        A(5,19) = 1; A(5,20) = -1;
        A(6,21) = 1;
        A(7,0) = -1; A(7,1) = 5;
        A(8,0) = 1; A(8,1) = -12;
        A(9,3) = -1*x4s; A(9,5) = 1.5;
        A(10,4) = -1*x5s; A(10,6) = 1.1;

        DenseVector b;
        b.setZero(11);
        b(0) = -27/x3s;
        b(1) = -397.5/(x3s*x3s);
        b(6) = 40/x3s;
        b(9) = -1.9;
        b(10) = -1.9;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false); // All variables

        //std::vector<int> vars = {0,1,2,3,4,5};
        cs->add(lincon);
    }

    BB::BranchAndBound bnb(cs);
    SolverResult res = bnb.optimize();

    if (res.status == SolverStatus::OPTIMAL)
        fopt_found = res.objectiveValue;
}

} // namespace CENSO
