/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTLINEAR_H
#define CONSTRAINTLINEAR_H

#include "constraint.h"

namespace CENSO
{

/**********************************
 * Linear constraints on the form:
 * lb <= Ax <= ub = b
 **********************************/
class ConstraintLinear : public Constraint
{
public:
    ConstraintLinear(std::vector<VariablePtr> variables, DenseMatrix &A, DenseVector &b, bool equality);
    //ConstraintLinear(const ConstraintLinear &copy) = default;
    ConstraintLinear(const ConstraintLinear &copy, bool deep = false);
    ConstraintLinear& operator=(ConstraintLinear const& assign) = delete;

    ~ConstraintLinear() override {}

    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintLinear>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override { DenseVector ddx; ddx.resize(0); return ddx; }

    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

    void structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol) override {}

    bool reduceVariableRanges() const override;

    void writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const override;

private:
    DenseMatrix A;
    bool equality;

    bool intervalAnalysis() const;

};

} // namespace CENSO

#endif // CONSTRAINTLINEAR_H
