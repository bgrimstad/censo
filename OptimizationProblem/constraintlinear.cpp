/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintlinear.h"

using std::cout;
using std::endl;

namespace CENSO
{

ConstraintLinear::ConstraintLinear(const ConstraintLinear &copy, bool deep)
    : Constraint(copy, deep), A(copy.A), equality(copy.equality)
{
}

/*
 * TODO: Sanitize A
 * For example, rows with only zeros cause trouble for Ipopt
 */
ConstraintLinear::ConstraintLinear(std::vector<VariablePtr> variables, DenseMatrix &A, DenseVector &b, bool equality)
    : Constraint(variables), A(A), equality(equality)
{
    assert(A.rows() == b.rows());
    assert(A.cols() == (int)variables.size());
    numConstraints = A.rows();

    for (unsigned int i = 0; i < numConstraints; i++)
    {
        if (equality)
            lb.push_back(b(i));
        else
            lb.push_back(-INF);

        ub.push_back(b(i));
    }

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = true;
    constraintConvex = true;
    convexRelaxationAvailable = true;

    nnzJacobian = 0;
    nnzHessian = 0;
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            if (A(i,j) != 0)
            {
                nnzJacobian++;
            }
        }
    }

    constraintName = "Linear";

    checkConstraintSanity();
}

DenseVector ConstraintLinear::eval(const DenseVector &x) const
{
    DenseVector y = A*x;
    return y;
}

DenseVector ConstraintLinear::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);

    unsigned int k = 0;
    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            if (A(i,j) != 0)
            {
                dx(k) = A(i,j);
                k++;
            }
        }
    }
    assert(k == nnzJacobian);
    return dx;
}

void ConstraintLinear::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    iRow.clear();
    jCol.clear();

    for (int i = 0; i < A.rows(); i++)
    {
        for (int j = 0; j < A.cols(); j++)
        {
            if (A(i,j) != 0)
            {
                iRow.push_back(i);
                jCol.push_back(j);
            }
        }
    }

    assert(iRow.size() == jCol.size());
    assert(iRow.size() == nnzJacobian);
}

bool ConstraintLinear::reduceVariableRanges() const
{
    return intervalAnalysis();
}

/*
 * Interval analysis of a linear constraint (Gauss-Seidel method)
 * Ax <= b or Ax = b, where lb <= x <= ub
 *
 * NOTE: this is the same as "Poor man's linear programming"
 */
bool ConstraintLinear::intervalAnalysis() const
{
    // For each variable k
    for (unsigned int k = 0; k < A.cols(); k++)
    {
        // Check bounds on variables i != k
        bool infBounds = false;
        for (unsigned int i = 0; i < A.cols(); i++)
        {
            if (i != k && (variables.at(i)->getLowerBound() <= -INF || variables.at(i)->getUpperBound() >= INF))
            {
                infBounds = true;
                break;
            }
        }

        if (infBounds)
            continue;

        // Do interval analysis on each row
        for (unsigned int i = 0; i < A.rows(); i++)
        {
            if (A(i,k) == 0)
                continue;

            // For row i calculate range of s = sum_j (Aij*xj) - Aik*xk - b,
            // so that s <= -aik*xk, and imin <= s <= imax

            double imin = -ub.at(i); // -b(i)
            double imax = -ub.at(i); // -b(i)

            for (unsigned int j = 0; j < A.cols(); j++)
            {
                if (j == k)
                    continue;

                if (A(i,j) > 0)
                {
                    imin += A(i,j)*variables.at(j)->getLowerBound();
                    imax += A(i,j)*variables.at(j)->getUpperBound();
                }
                else if (A(i,j) < 0)
                {
                    imin += A(i,j)*variables.at(j)->getUpperBound();
                    imax += A(i,j)*variables.at(j)->getLowerBound();
                }
                else // A(i,j) = 0
                {
                    continue;
                }
            }

            // Update variable bounds
            if (!equality)
            {
                // Inequality constraint Ax <= b
                double ib = -imin/A(i,k);

                if (A(i,k) > 0)
                {
                    // Calculate new upper bound on x(k)
                    if (!variables.at(k)->updateUpperBound(ib))
                        return false;
                }
                else // A(i,k) < 0
                {
                    // Calculate new lower bound on x(k)
                    if (!variables.at(k)->updateLowerBound(ib))
                        return false;
                }
            }
            else
            {
                // Equality constraint Ax = b
                double ilb = -INF;
                double iub = INF;

                if (A(i,k) > 0)
                {
                    ilb = -imax/A(i,k);
                    iub = -imin/A(i,k);
                }
                else // A(i,k) < 0
                {
                    ilb = -imin/A(i,k);
                    iub = -imax/A(i,k);
                }

                if (!variables.at(k)->updateBounds(ilb, iub))
                    return false;
            }
        }
    }

    return true;
}

void ConstraintLinear::writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const
{
    unsigned int counter = start;

    for (int i = 0; i < A.rows(); ++i)
    {
        os << "e" << std::to_string(counter) << ".. ";

        for (int j = 0; j < A.cols(); ++j)
        {
            double aij = A(i,j);

            if (assertNear(aij, 0.0))
            {
                continue;
            }
            else if (aij >= 0)
            {
                os << " + " << std::to_string(aij) << "*";
            }
            else
            {
                os << " - " << std::to_string(std::abs(aij)) << "*";
            }
            os << variables.at(j)->getName();
        }

        if (equality) os << " =E= ";
        else os << " =L= ";

        os << ub.at(i) << ";" << endl;

        counter++;
    }
}

} // namespace CENSO
