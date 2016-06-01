/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTBSPLINE_H
#define CONSTRAINTBSPLINE_H

#include "constraint.h"
#include "bsplinebuilder.h"

using SPLINTER::BSpline;

namespace CENSO
{

/*
 * In the current implementation a B-spline constraint
 * is implemented as: f(x) - y = 0
 */
class ConstraintBSpline : public Constraint
{
public:
    ConstraintBSpline(std::vector<VariablePtr> variables, const BSpline &bspline, bool equality);
    ConstraintBSpline(const ConstraintBSpline &copy, bool deep = false);
    ConstraintBSpline& operator=(const ConstraintBSpline &assign) = delete;
    ~ConstraintBSpline() override {}

    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintBSpline>(*this, deep);
    }

    void init(bool equality);

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

    void structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol) override;

    bool reduceVariableRanges() const override;

    std::vector<VariablePtr> getComplicatingVariables() const override;

    ConstraintPtr getConvexRelaxation() override;

    void localRefinement(const DenseVector &x) override;

private:

    // B-spline object
    // NOTE: control point matrix is assumed to be R^(m x (n + 1)),
    // with m = number of control points and n = number of variables (x)
    BSpline bspline;

    // Indicate if equality constraint
    bool equality;

    // Related to convex relaxation
    int maxNumAuxiliaryVariables;

    // Related to variable bounds check
    std::vector<Variable> storedVariables;

    // Relaxation is the hypercube defined by the extremal values in each dimension.
    // This is a poor convex relaxation of most B-splines and should be avoided.
    ConstraintPtr computeRelaxationHyperrectangle();

    // Relaxation is the convex hull of the control points.
    // This is the tightest possible convex relaxation of the B-spline
    // obtained from the control points alone.
    // The method requires one auxiliary variables per control point.
    ConstraintPtr computeRelaxationConvexHull();

    // Relaxation is the convex hull of the control points.
    // Uses QuickHull (with Akl-Toussaint heuristic) to retrieve
    // a subset S of the control points C that fulfills:
    // S <= C and conv(S) = conv(C)
    //void computeRelaxationQuickHull();

    // Auxililary function that can be removed
    // bool testConvexHullMembership(const DenseMatrix &C, const DenseVector &p) const;

    // Reduce domain of B-spline
    void reduceBSplineDomain();

    // Control point bounds deduction
    bool controlPointBoundsDeduction() const;

    // Store old variables
    void storeVariables();

    // Check if variable bounds have changed
    bool haveVariableBoundsChanged() const;

    // Central difference (for testing only)
    DenseVector centralDifference(const DenseVector &x) const;
};

} // namespace CENSO

#endif // CONSTRAINTBSPLINE_H
