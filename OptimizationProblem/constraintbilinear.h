/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTBILINEAR_H
#define CONSTRAINTBILINEAR_H

#include "constraint.h"

namespace CENSO
{

/**********************************
 * Bilinear constraint:
 * lb <= a*x(0)*x(1) - x(2) <= ub
 * where a is a constant
 **********************************/

class ConstraintBilinear : public Constraint
{
public:
    ConstraintBilinear(std::vector<VariablePtr> variables, double a, double lb, double ub);
    ConstraintBilinear(const ConstraintBilinear &copy, bool deep = false);
    ConstraintBilinear& operator=(const ConstraintBilinear &assign) = delete;
    ~ConstraintBilinear() override {}

    // Clone function - uses copy constructor
    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintBilinear>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int>& iRow, std::vector<int>& jCol) override;

    void structureHessian(std::vector<int>& eqnr,std::vector<int>& iRow, std::vector<int>& jCol) override;

    std::vector<VariablePtr> getComplicatingVariables() const override;

    bool reduceVariableRanges() const override;

    ConstraintPtr getConvexRelaxation() override;

private:
    double a;

};

} // namespace CENSO

#endif // CONSTRAINTBILINEAR_H
