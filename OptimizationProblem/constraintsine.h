/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTSINE_H
#define CONSTRAINTSINE_H

#include "Utils/definitions.h"
#include "constraint.h"

namespace CENSO
{

class ConstraintSine : public Constraint
{
public:
    ConstraintSine(std::vector<VariablePtr> variables, double a, double b, bool equality);
    ConstraintSine(const ConstraintSine &copy, bool deep = false); // Default copy constructor is ok
    ConstraintSine& operator = ( ConstraintSine const& assign) = delete;

    ~ConstraintSine() override;

    // Clone function - uses copy constructor
    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintSine>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const;

    DenseVector evalJacobian(const DenseVector &x) const;

    DenseVector evalHessian(const DenseVector &x) const;

    virtual void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol);

    virtual void structureHessian(std::vector<int> &eqnr,std::vector<int> &iRow, std::vector<int> &jCol);

private:
    double a;
    double b;

};

} // namespace CENSO

#endif // CONSTRAINTSINE_H
