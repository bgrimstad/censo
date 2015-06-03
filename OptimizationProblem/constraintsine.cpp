/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintsine.h"
#include <cmath>

namespace CENSO
{

ConstraintSine::ConstraintSine(const ConstraintSine &copy, bool deep)
    : Constraint(copy, deep), a(copy.a), b(copy.b)
{
}

ConstraintSine::ConstraintSine(std::vector<VariablePtr> variables, double a, double b, bool equality)
    : Constraint(variables), a(a), b(b)
{
    assert(variables.size() == 1);
    numConstraints = 1;

    if (equality)
        lb.push_back(0);
    else
        lb.push_back(-INF);

    ub.push_back(0);

    nnzJacobian = 2;
    nnzHessian = 1;

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = false;

    constraintName = "Constraint Sine";

    checkConstraintSanity();
}

ConstraintSine::~ConstraintSine()
{
}

DenseVector ConstraintSine::eval(const DenseVector &x) const
{
    DenseVector y = DenseVector::Zero(numConstraints);
    y(0) = a*sin(b*x(0)) - x(1);
    return y;
}

DenseVector ConstraintSine::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);

    dx(0) = a*b*std::cos(b*x(0));
    dx(1) = -1;

    return dx;
}

DenseVector ConstraintSine::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);
    ddx(0) = -a*b*b*sin(b*x(0));
    return ddx;
}

void ConstraintSine::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    iRow.push_back(0);  jCol.push_back(0);
    iRow.push_back(0);  jCol.push_back(1);
}

void ConstraintSine::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    eqnr.push_back(0);
    iRow.push_back(0);
    jCol.push_back(0);
}

} // namespace CENSO
