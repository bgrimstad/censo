/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTQUADRATIC_H
#define CONSTRAINTQUADRATIC_H

#include "constraint.h"

namespace CENSO
{

/*
 * Quadratic constraints of the type:
 * lb <= x'*A*x + b'*x + c <= ub
 * where A is a square matrix of size nx*nx (A does not have to be symmetric),
 * b is a column vector of size nx, and c is a double
 */
class ConstraintQuadratic : public Constraint
{
public:
    ConstraintQuadratic(std::vector<VariablePtr> variables, DenseMatrix A, DenseMatrix b, double c, double lb, double ub);
    ConstraintQuadratic(const ConstraintQuadratic &copy, bool deep = false);
    ConstraintQuadratic& operator = (const ConstraintQuadratic &assign) = delete;

    ~ConstraintQuadratic() override {}

    // Clone function - uses copy constructor
    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintQuadratic>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int>& iRow, std::vector<int>& jCol) override;

    void structureHessian(std::vector<int>& eqnr,std::vector<int>& iRow, std::vector<int>& jCol) override;

    ConstraintPtr getConvexRelaxation() override;

private:
    // f(x,y) = x'*A*x + b'*x + c - y
    DenseMatrix A;
    DenseMatrix b;
    double c;

    // Hessian matrix H = A.transpose() + A
    DenseMatrix H;

};

} // namespace CENSO

#endif // CONSTRAINTQUADRATIC_H
