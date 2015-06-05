/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTAPPROXIMANT_H
#define CONSTRAINTAPPROXIMANT_H

//#include "approximant.h"
//#include "constraint.h"

//using SPLINTER::Approximant;

//namespace CENSO
//{

///*
// * In the current implementation a approximant constraint
// * is implemented as: f(x) - y = 0, where f(x) is the approximant
// */
//class ConstraintApproximant : public Constraint
//{
//public:
//    ConstraintApproximant(std::vector<VariablePtr> variables, const Approximant &approx, bool equality);
//    ConstraintApproximant(const ConstraintApproximant &copy, bool deep = false);
//    ConstraintApproximant& operator=(const ConstraintApproximant &assign) = delete;
//    ~ConstraintApproximant() override {}

//    ConstraintPtr clone(bool deep = true) const override
//    {
//        return std::make_shared<ConstraintApproximant>(*this, deep);
//    }

//    void init(bool equality);

//    DenseVector eval(const DenseVector &x) const override;

//    DenseVector evalJacobian(const DenseVector &x) const override;

//    DenseVector evalHessian(const DenseVector &x) const override;

//    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

//    void structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol) override;

//    //bool reduceVariableRanges() const override;

//    std::vector<VariablePtr> getComplicatingVariables() const override;

//    //ConstraintPtr getConvexRelaxation() override;

//private:

//    // B-spline object
//    Appoximant approx;

//    // Indicate if equality constraint
//    bool equality;
//};

//} // namespace CENSO

#endif // CONSTRAINTAPPROXIMANT_H
