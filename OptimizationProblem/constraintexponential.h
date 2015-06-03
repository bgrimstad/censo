/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTEXPONENTIAL_H
#define CONSTRAINTEXPONENTIAL_H

#include "constraint.h"

namespace CENSO
{

/*
 * Exponential constraint:
 * lb <= a*exp(b*x(0)) - x(1) <= ub
 * where a and b are constants
 */
class ConstraintExponential : public Constraint
{
public:
    ConstraintExponential(std::vector<VariablePtr> variables, double a, double b, double lb, double ub);
    ConstraintExponential(const ConstraintExponential &copy, bool deep = false);
    ConstraintExponential& operator = (const ConstraintExponential &assign) = delete;

    ~ConstraintExponential() override {}

    // Clone function - uses copy constructor
    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintExponential>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

    void structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol) override;

private:
    double a, b;
};

} // namespace CENSO

#endif // CONSTRAINTEXPONENTIAL_H
