/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintpolynomial.h"
#include "Eigen/Eigenvalues"

using std::cout;
using std::endl;

namespace CENSO
{

ConstraintPolynomial::ConstraintPolynomial(const ConstraintPolynomial &copy, bool deep)
    : Constraint(copy, deep), c(copy.c), E(copy.E), equality(copy.equality)
{
}

// Implementes polynomial p(x) <= 0 or p(x) = 0
ConstraintPolynomial::ConstraintPolynomial(std::vector<VariablePtr> variables, DenseVector c, DenseMatrix E, bool equality)
    : Constraint(variables), c(c), E(E), equality(equality)
{
    assert(c.rows() >= 1); // At least one monomial (assumes that c is a column vector)
    assert(E.cols() >= 1); // At least one variable
    assert(c.rows() == E.rows());

    assert(E.cols() == (int)variables.size());
    numConstraints = 1;

    this->ub = {0};
    if (equality)
        this->lb = {0};
    else
        this->lb = {-INF};

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = false;

    // Check if polynomial has non-negative integer exponents
    for (int i = 0; i < E.rows(); i++)
    {
        // If coefficient is zero, continue
        if (c(i) == 0) continue;

        // First column not considered, coefficients do not affect linearity
        for (int j = 0; j < E.cols(); j++)
        {
            assert(E(i,j) == floor(E(i,j)));
            assert(E(i,j) >= 0);
        }
    }

    // Check if polynomial is linear
    // Does not check for cancelling terms
    bool linearityCheck = true;
    for (int i = 0; i < E.rows(); i++)
    {
        // If coefficient is zero, continue
        if (c(i) == 0) continue;

        // power of monomial/term
        double termPower = 0;

        // First column not considered, coefficients do not affect linearity
        for (int j = 0; j < E.cols(); j++)
        {
            if (E(i,j) != 0 && E(i,j) != 1)
            {
                linearityCheck = false;
                break;
            }
            else
            {
                termPower += E(i,j);
            }
        }
        if (termPower != 0 && termPower != 1)
        {
            linearityCheck = false;
            break;
        }
        if (!linearityCheck) break;
    }
    if (linearityCheck) constraintLinear = true;

    // Check if polynomial is (globally) convex
    // Local convexity, using domain bounds, is not investigated
    if (constraintLinear)
    {
        constraintConvex = true;
    }
    else
    {
        // Check for inequality
        if (!equality)
        {
            // Every term must be convex
            bool convexityCheck = true;
            for (int i = 0; i < E.rows(); i++)
            {
                if (c(i) == 0) continue;

                std::vector<int> e;
                for (int j = 0; j < E.cols(); j++) e.push_back(E(i,j));

                if (isTermLinear(e) || isTermConstant(e)) continue;

                // Check if coefficient is positive and term is convex
                if (!(c(i) > 0 && isTermConvex(e)))
                {
                    convexityCheck = false;
                    break;
                }
            }
            if (convexityCheck) constraintConvex = true;
        }
        else
        {
            // Equality + nonlinearity = non-convex constraint
        }
    }

    // Calculate number of non-zeros in Jacobian and Hessian
    nnzJacobian = variables.size();
    nnzHessian = 0; // Lower left triangle only

    // Hessian is symmetric so fill out lower left triangle only
    for (unsigned int row = 0; row < variables.size(); row++)
    {
        for (unsigned int col = 0; col <= row; col++)
        {
            nnzHessian++;
        }
    }

    // If the Hessian is a quadratic and symmetric matrix
    // the dense number of non-zeros can be calculated as:
//    int nnzHessian2 = 0;
//    for (int i = 1; i <= dimensionDomain; i++) nnzHessian2 += i;

    constraintName = "Polynomial";

    // Check sanity of constraint
    checkConstraintSanity();
}

DenseVector ConstraintPolynomial::eval(const DenseVector &x) const
{
    assert(x.size() == (int)variables.size());

    DenseVector y = DenseVector::Zero(numConstraints);

    // Calculate polynomial
    y(0) = 0;

    // One term at the time
    for (int i = 0; i < E.rows(); i++)
    {
        // Calculate term and add it to y
        double monomial = c(i);
        if (monomial == 0) continue;

        for (int j = 0; j < E.cols(); j++)
        {
            if (E(i,j) == 0) continue;
            monomial = monomial*pow(x(j),E(i,j));
        }
        y(0) += monomial;
    }

    return y;
}

DenseVector ConstraintPolynomial::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);

    for (unsigned int k = 0; k < variables.size(); k++)
    {
        dx(k) = 0;

        // For each term (monomial) in the polynomial
        for (unsigned int i = 0; i < E.rows(); i++)
        {
            // Calculate term and add it to dx
            double monomial = c(i);
            if (monomial == 0) continue;

            // For each variable in the term
            for (unsigned int j = 0; j < E.cols(); j++)
            {
                // Differentiate
                if (j == k)
                {
                    if (E(i,j) < 1)
                    {
                        // Variable differentiates to zero
                        monomial = 0;
                        break;
                    }
                    else
                    {
                        monomial = monomial*E(i,j)*pow(x(j),E(i,j)-1);
                    }
                }
                else
                {
                    monomial = monomial*pow(x(j),E(i,j));
                }
            }
            dx(k) += monomial;
        }
    }

    return dx;
}


DenseVector ConstraintPolynomial::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);

    // Calculate Hessian
    // Hessian is symmetric so fill out lower left triangle only
    int idx = 0;
    for (unsigned int row = 0; row < variables.size(); row++)
    {
        for (unsigned int col = 0; col <= row; col++)
        {
            // Calculate Hessian element
            ddx(idx) = 0;

            // For each term (monomial) in the polynomial
            for (unsigned int i = 0; i < E.rows(); i++)
            {
                // Calculate term and add it to ddx
                double monomial = c(i);
                if (monomial == 0) continue;

                // For each variable in term
                for (unsigned int j = 0; j < E.cols(); j++)
                {
                    // Differentiate
                    if (j == row && row == col)
                    {
                        // Diagonal - differentiate twice
                        if (E(i,j) < 2)
                        {
                            // Variable differentiates to zero
                            monomial = 0;
                            break;
                        }
                        else
                        {
                            // Two differentiations
                            // This pow() may throw errors if E contains non-integer values, be careful
                            monomial = monomial*E(i,j)*(E(i,j)-1)*pow(x(j),E(i,j)-2);
                        }
                    }
                    else if (j == row || j == col)
                    {
                        if (E(i,j) == 0)
                        {
                            // Variable differentiates to zero
                            monomial = 0;
                            break;
                        }
                        else
                        {
                            // One differantiation
                            monomial = monomial*E(i,j)*pow(x(j),E(i,j)-1);
                        }
                    }
                    else
                    {
                        monomial = monomial*pow(x(j),E(i,j));
                    }
                }

                // Add term to ddx
                ddx(idx) += monomial;
            }

            idx++;
        }
    }

    return ddx;
}

void ConstraintPolynomial::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    for (unsigned int i = 0; i < variables.size(); i++)
    {
        iRow.push_back(0); jCol.push_back(i);
    }
}

void ConstraintPolynomial::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    // Calculate Hessian
    // Hessian is symmetric so fill out lower left triangle only
    for (unsigned int row = 0; row < variables.size(); row++)
    {
        for (unsigned int col = 0; col <= row; col++)
        {
            eqnr.push_back(0);
            iRow.push_back(row);
            jCol.push_back(col);
        }
    }
}

// Currently this function only asserts global convexity by checking that epxonents are 0, 1, or even.
// It does not consider domain bounds to check local convexity. Function tree expantion and interval analysis
// would be powerful tools to use for this and the generation of convex relaxations of general polynomials.
bool ConstraintPolynomial::isTermConvex(std::vector<int> e) const
{
    if (isTermLinear(e) || isTermConstant(e)) return true;

    // Check number of variables in monomial
    int variables = 0;
    int first = -1; // Position of the first exponent != 0 in e
    for (unsigned int i = 0; i < e.size(); i++)
    {
        if (e.at(i) != 0)
        {
            variables++;
            if (first < 0) first = i;
        }
    }

    if (variables == 0)
    {
        // Nothing to do here, monomial is a constant
        return true;
    }
    else if (variables == 1)
    {
        // Univariate
        // Check if exponent is even or 1 (disregards domain, e.g. that x^3 is convex for x >= 0)
        if (e.at(first) % 2 == 0 || e.at(first) == 1) return true;
    }
    else
    {
        // Multivariate
        return false;
    }

    return false;
}

bool ConstraintPolynomial::isTermLinear(std::vector<int> e) const
{
    // Check number of variables in term
    int variables = 0;
    for (unsigned int i = 0; i < e.size(); i++)
    {
        if (e.at(i) != 0)
        {
            variables++;
            if (e.at(i) != 1) return false;
        }
        if (variables > 1) return false;
    }
    return true;
}

bool ConstraintPolynomial::isTermConstant(std::vector<int> e) const
{
    // Check number of variables in term
    for (unsigned int i = 0; i < e.size(); i++)
    {
        if (e.at(i) != 0) return false;
    }
    return true;
}

void ConstraintPolynomial::writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const
{

    os << "e" << std::to_string(start) << ".. ";

    for (int i = 0; i < c.size(); ++i)
    {
        auto ci = c(i);

        if (assertNear(ci, 0.0))
        {
            continue;
        }
        else if (ci > 0)
        {
            os << " + " << std::to_string(ci);
        }
        else
        {
            os << " - " << std::to_string(std::abs(ci));
        }

        for (int j = 0; j < E.cols(); ++j)
        {
            for (int k = 0; k < E(i,j); ++k)
                os << "*" << variables.at(j)->getName();
        }
    }

    if (equality) os << " =E= ";
    else os << " =L= ";

    os << ub.at(0) << ";" << endl;
}

} // namespace CENSO
