/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintquadratic.h"
#include "Eigen/Eigenvalues"

using std::cout;
using std::endl;

namespace CENSO
{

ConstraintQuadratic::ConstraintQuadratic(const ConstraintQuadratic &copy, bool deep)
    : Constraint(copy, deep), A(copy.A), b(copy.b), c(copy.c), H(copy.H)
{
}

ConstraintQuadratic::ConstraintQuadratic(std::vector<VariablePtr> variables, DenseMatrix A, DenseMatrix b, double c, double lb, double ub)
    : Constraint(variables), A(A), b(b), c(c)
{
    assert(A.cols() == (int)variables.size());
    assert(A.rows() == b.rows());
    assert(b.cols() == 1);

    numConstraints = 1;

    this->lb.push_back(lb);
    this->ub.push_back(ub);

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = true;

    nnzJacobian = A.rows();
    nnzHessian = 0;

    constraintName = "Constraint Quadratic";

    //    // Check for parameters for NaN
    //    for (int i = 0; i < A.rows(); i++)
    //    {
    //        for (int j = 0; j < A.cols(); j++)
    //        {
    //            bool nanA = false;
    //            if (A(i,j) != A(i,j)) nanA = true;
    //            assert(nanA == false);
    //        }
    //        bool nanb = false;
    //        if (b(i) != b(i)) nanb = true;
    //        assert(nanb == false);
    //    }

    // Calculate and store Hessian
    H = A.transpose() + A;

    // H is symmetric so fill out lower left triangle only
    for (int row = 0; row < H.rows(); row++)
    {
        for (int col = 0; col <= row; col++)
        {
            if (H(row,col) != 0)
            {
                nnzHessian++;
            }
        }
    }

    // Check convexity using Hessian
    Eigen::EigenSolver<DenseMatrix> es(H);
    DenseVector eigs = es.eigenvalues().real();
    double minEigVal = 0;
    for (int i = 0; i < eigs.rows(); i++)
    {
        if (eigs(i) < minEigVal) minEigVal = eigs(i);
    }

    if (minEigVal >= 0 && lb <= -INF)
    {
        constraintConvex = true;
    }
    // Note could also check that max. eigen value <= 0 and ub = INF

    checkConstraintSanity();
}

DenseVector ConstraintQuadratic::eval(const DenseVector &x) const
{
    // Convert x to matrix so that it can be transposed
    DenseMatrix z(x);
    DenseMatrix zt = z.transpose();
    DenseMatrix bt = b.transpose();
    DenseVector cvec(1); cvec(0) = c;

    DenseVector y = zt*A*z + bt*z + cvec;
    return y;
}

DenseVector ConstraintQuadratic::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = H*x + b; // H = A.transpose() + A
    return dx;
}


DenseVector ConstraintQuadratic::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);

    // H is symmetric so fill out lower left triangle only
    int idx = 0;
    for (int row = 0; row < H.rows(); row++)
    {
        for (int col = 0; col <= row; col++)
        {
            if (H(row,col) != 0)
            {
                ddx(idx++) = H(row,col);
            }
        }
    }

    return ddx;
}

void ConstraintQuadratic::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    for (int i = 0; i < A.rows(); i++)
    {
        iRow.push_back(0); jCol.push_back(i);
    }
}

void ConstraintQuadratic::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    // Calculate Hessian
    DenseMatrix H = A.transpose() + A;

    // H is symmetric so fill out lower left triangle only
    for (int row = 0; row < A.rows(); row++)
    {
        for (int col = 0; col <= row; col++)
        {
            if (H(row,col) != 0)
            {
                eqnr.push_back(0);
                iRow.push_back(row);
                jCol.push_back(col);
            }
        }
    }
}


ConstraintPtr ConstraintQuadratic::getConvexRelaxation()
{
    if (!constraintConvex)
    {
        // Temp fix
        bool isQuadConvex = false;
        assert(isQuadConvex);
        return this->clone();
        // Temp fix end

        // TODO: this code has not been used or tested in a long while!
        Eigen::EigenSolver<DenseMatrix> es(H);
        DenseVector eigs = es.eigenvalues().real();

        double max_alpha = 0;
        for (int i = 0; i < eigs.rows(); i++)
        {
            cout <<"eig nr " << i << " = " << eigs(i) << endl;
            max_alpha = std::max(max_alpha, -0.5 * eigs(i));
        }

        if (max_alpha != 0)
        {
            cout << "Max alpha = " << max_alpha <<endl;

            std::vector<double> alpha;
            for (unsigned int i = 0; i < variables.size(); i++)
            {
                alpha.push_back(max_alpha);
            }

            cout << "Sour about that negative hessian." <<endl;

            //ConstraintPtr quad_org(this);
            //ConstraintPtr quad_org(this->clone());
            //return new ConstraintDecoratorQuadraticConvexRelaxation(quad_org,alpha);
        }
    }

    // Constraint is convex
    return this->clone(false);
}

} // namespace CENSO
