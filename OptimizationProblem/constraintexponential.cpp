/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintexponential.h"

namespace CENSO
{

ConstraintExponential::ConstraintExponential(const ConstraintExponential &copy, bool deep)
    : Constraint(copy, deep), a(copy.a), b(copy.b)
{
}

ConstraintExponential::ConstraintExponential(std::vector<VariablePtr> variables, double a, double b, double lb, double ub)
    : Constraint(variables), a(a), b(b)
{
    assert(variables.size() == 2);
    numConstraints = 1;

    this->lb.push_back(lb);
    this->ub.push_back(ub);

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = false;

    // Sub-level sets of -exp(x) are actually convex, even though -exp(x) is strictly concave!
    // NOTE: Consider changing convexity check to include this exception

    // Depends on coefficient a
    if (a > 0 && lb != ub && lb <= -INF && ub < INF)
    {
        // Convex and sublevel set = convex
        constraintConvex = true;
    }
    else if (a < 0 && lb != ub && lb > -INF && ub >= INF)
    {
        // Concave and superlevel set = convex
        constraintConvex = true;
    }

    nnzJacobian = 2;
    nnzHessian  = 1;

    checkConstraintSanity();
}

DenseVector ConstraintExponential::eval(const DenseVector &x) const
{
    DenseVector y = DenseVector::Zero(numConstraints);
    y(0) = a*exp(b*x(0)) - x(1);
    return y;
}

DenseVector ConstraintExponential::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);
    dx(0) = a*b*exp(b*x(0));
    dx(1) = -1;
    return dx;
}

DenseVector ConstraintExponential::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);
    ddx(0) = a*b*b*exp(b*x(0));
    return ddx;
}

void ConstraintExponential::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    iRow.push_back(0);  jCol.push_back(0);
    iRow.push_back(0);  jCol.push_back(1);
}

void ConstraintExponential::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    eqnr.push_back(0);
    iRow.push_back(0);
    jCol.push_back(0);
}

} // namespace CENSO
