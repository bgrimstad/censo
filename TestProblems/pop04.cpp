/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "pop04.h"
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

POP04::POP04()
{
    testName = "P04";
    fopt_known = -4;
    fopt_found = INF;
}

void POP04::runProblem()
{
    // Problem 4 (Nataraj)
    // cout << "\n\nSolving problem P04..." << endl;

    int dim = 3+1;

    // x1,x2,x3,l1, l1 >= 0
    std::vector<double> costs = {-2, 1, -1, 0};
    std::vector<double> lb = {0,0,0,0};
    std::vector<double> ub = {2,10,3,INF};
    std::vector<double> z0(dim,0);

    std::vector<VariablePtr> vars;
    for (int i = 0; i < dim; i++)
    {
        auto var = std::make_shared<Variable>(costs.at(i), lb.at(i), ub.at(i));
        vars.push_back(var);
    }

    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    { // Quadratic constraint = l1
        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2),
            vars.at(3)
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

        std::vector<unsigned int> deg = {2,2,2};

        DenseVector c(27);
        c.setZero();
        c(0) = 24; // Constant
        c(1) = -13; // x3
        c(2) = 2; // x3^2
        c(3) = 9; // x2
        c(4) = -2; // x2*x3
        c(6) = 2; // x2^2
        c(9) = -20; // x1
        c(10) = 4; // x1*x3
        c(12) = -4; // x1*x2
        c(18) = 4; // x1^2

        DenseMatrix T = getTransformationMatrix(deg, thislb, thisub);
        DenseMatrix coeffs = T*c;

        std::vector< std::vector<double> > knots = getRegularKnotVectors(deg, thislb, thisub);

        BSpline bs(coeffs.transpose(), knots, deg);

        ConstraintPtr cbs = std::make_shared<ConstraintBSpline>(cvars, bs, true);

        cs->add(cbs);
    }

    { // Linear constraints

        DenseMatrix A(2,3);
        A.setZero();
        A(0,0) = 1; A(0,1) = 1; A(0,2) = 1;
        A(1,1) = 3; A(1,2) = 1;

        DenseVector b;
        b.setZero(2);
        b(0) = 4;
        b(1) = 6;

        std::vector<VariablePtr> cvars = {
            vars.at(0),
            vars.at(1),
            vars.at(2)
        };
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

bool POP04::validateResult()
{
    if (std::abs(fopt_found - fopt_known) <= 1e-3)
        return true;
    return false;
}

} // namespace CENSO
