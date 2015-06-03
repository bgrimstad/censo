/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTPOLYNOMIAL_H
#define CONSTRAINTPOLYNOMIAL_H

#include "constraint.h"

namespace CENSO
{

/* Constraint polynomial implements a constraint on the form:
 * lb <= f(x) <= ub,
 * The polynomial f(x), with x in R^n, has m terms on the form
 * f(x) = c(0)*x(0)^E(0,0)*x(1)^E(0,1)*...*x(n-1)^E(0,n-1)
 *       +c(1)*x(0)^E(1,0)*x(1)^E(1,1)*...*x(n-1)^E(1,n-1)
 *       +...
 *       +c(m-1)*x(0)^E(m-1,0)*x(1)^E(m-1,1)*...*x(n-1)^E(m-1,n-1)
 * where c in R^m is a vector with coefficients for each of the m monomials,
 * and E in N^(mxn) is a matrix with the exponents of each variable in each of the m terms,
 * e.g. the first row of E defines the first term with variable exponents E(0,0) to E(0,n-1).
 *
 * The polynomial is thus completely defined by lb, ub, c, and E.
 *
 * NOTE: Only positive integral exponents are supported! The pow() function is a bit picky:
 * "If base is negative and exponent is not an integral value, or if base is zero and exponent is negative,
 * a domain error occurs, setting the global variable errno to the value EDOM."
 */
class ConstraintPolynomial : public Constraint
{
public:
    ConstraintPolynomial(std::vector<VariablePtr> variables, DenseVector c, DenseMatrix E, bool equality);
    ConstraintPolynomial(const ConstraintPolynomial &copy, bool deep = false);
    ConstraintPolynomial& operator=(const ConstraintPolynomial &assign) = delete;

    ~ConstraintPolynomial() override {}

    // Clone function - uses copy constructor
    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintPolynomial>(*this, deep);
    }

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

    void structureHessian(std::vector<int> &eqnr,std::vector<int> &iRow, std::vector<int> &jCol) override;

    void writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const override;

    /*
     * Getter specific for polynomial constraints
     */
    DenseVector getCoefficients() const
    {
        return c;
    }

private:
    DenseVector c;
    DenseMatrix E;
    bool equality;

    bool isTermConvex(std::vector<int> e) const;
    bool isTermLinear(std::vector<int> e) const;
    bool isTermConstant(std::vector<int> e) const;
};

} // namespace CENSO

#endif // CONSTRAINTPOLYNOMIAL_H
