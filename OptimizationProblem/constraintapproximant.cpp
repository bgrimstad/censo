/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintapproximant.h"
#include "constraintlinear.h"
#include "constraintset.h"

using std::cout;
using std::endl;

namespace CENSO
{

//ConstraintApproximant::ConstraintApproximant(std::vector<VariablePtr> variables, const Approximant &approximant, bool equality)
//    : Constraint(variables),
//      approximant(approximant),
//      equality(true) // Ineq not implemented!
//{
//    init(equality);
//}

//ConstraintApproximant::ConstraintApproximant(const ConstraintApproximant &copy, bool deep)
//    : Constraint(copy, deep),
//      approximant(copy.approximant),
//      equality(copy.equality)
//{
//}

//void ConstraintApproximant::init(bool equality)
//{
//    numConstraints = 1;
//    assert(variables.size() == approximant.getNumVariables() + numConstraints);

//    // Update variable bounds
//    auto varlb = approximant.getDomainLowerBound();
//    auto varub = approximant.getDomainUpperBound();

//    for (unsigned int i = 0; i < approximant.getNumVariables(); i++)
//    {
//        variables.at(i)->updateBounds(varlb.at(i), varub.at(i));
//    }

//    for (unsigned int i = 0; i < numConstraints; i++)
//    {
//        if (equality)
//        {
//            lb.push_back(0);
//        }
//        else
//        {
//            lb.push_back(-INF);
//        }

//        ub.push_back(0);

//        constraintTypes.push_back(ConstraintType::NONLINEAR_NONCONVEX);
//    }

//    nnzJacobian = approximant.getNumVariables() + 1;
//    nnzHessian = 0;

//    jacobianCalculated = true;
//    hessianCalculated = false;
//    constraintLinear = false;
//    constraintConvex = false;
//    convexRelaxationAvailable = true;

//    for (unsigned int row = 0; row < approximant.getNumVariables(); row++)
//        for (unsigned int col = 0; col <= row; col++)
//            nnzHessian++;

//    constraintName = "Approximant";

//    // Update bounds
//    reduceVariableRanges();

//    checkConstraintSanity();
//}

//DenseVector ConstraintApproximant::eval(const DenseVector &x) const
//{
//    // Only x-variables (B-spline input variables) are adjusted to bounds
//    DenseVector xy = x;
//    DenseVector xa = adjustToDomainBounds(x);

//    DenseVector xx = xa.block(0,0,variables.size()-numConstraints,1);
//    DenseVector yy = xy.block(variables.size()-numConstraints,0,numConstraints,1);

//    double by = approximant.eval(xx);

//    DenseVector y = DenseVector::Zero(numConstraints);
//    y(0) = by - yy(0);

//    //    adjustToDomainBounds(x);
//    //    y.resize(numConstraints);

//    //    VecD xx = x.block(0,0,variables.size()-numConstraints,1);
//    //    VecD yy = x.block(variables.size()-numConstraints,0,numConstraints,1);

//    //    VecD by = bspline->evaluate(xx);

//    //    y = by - yy;

//    return y;
//}

//DenseVector ConstraintApproximant::evalJacobian(const DenseVector &x) const
//{
//    DenseVector xa = adjustToDomainBounds(x);
//    DenseVector dx = DenseVector::Zero(nnzJacobian);

//    //return centralDifference(xa);

//    // Get x-values
//    DenseVector xx = xa.block(0,0,bspline.getNumVariables(),1);

//    // Evaluate Jacobian
//    DenseMatrix jac = approximant.evalJacobian(xx);

//    // Derivatives on inputs x
//    int k = 0;
//    for (int i = 0; i < jac.rows(); i++)
//    {
//        for (int j = 0; j < jac.cols(); j++)
//        {
//            dx(k++) = jac(i,j);
//        }
//    }

//    // Derivatives on outputs y
//    for (unsigned int i = 0; i < numConstraints; i++)
//    {
//        dx(k++) = -1;
//    }

//    return dx;
//}

//DenseVector ConstraintApproximant::evalHessian(const DenseVector &x) const
//{
//    DenseVector xa = adjustToDomainBounds(x);
//    DenseVector ddx = DenseVector::Zero(nnzHessian);

//    // Get x-values
//    DenseVector xx = xa.block(0,0,bspline.getNumVariables(),1);

//    // Calculate Hessian
//    DenseMatrix H = approximant.evalHessian(xx);

//    // H is symmetric so fill out lower left triangle only
//    int idx = 0;
//    for (int row = 0; row < H.rows(); row++)
//    {
//        for (int col = 0; col <= row; col++)
//        {
//            //if (H(row,col) != 0)
//            //{
//                ddx(idx++) = H(row,col);
//            //}
//        }
//    }

//    return ddx;
//}

//void ConstraintApproximant::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
//{
//    for (unsigned int i = 0; i < numConstraints; i++)
//    {
//        for (unsigned int j = 0; j < variables.size(); j++)
//        {
//            iRow.push_back(i);
//            jCol.push_back(j);
//        }
//    }
//}

//void ConstraintApproximant::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
//{
//    // H is symmetric so fill out lower left triangle only
//    // Also, not neccessary to fill out for outputs y (=0)
//    for (unsigned int row = 0; row < approximant.getNumVariables(); row++)
//    {
//        for (unsigned int col = 0; col <= row; col++)
//        {
//            eqnr.push_back(0);
//            iRow.push_back(row);
//            jCol.push_back(col);
//        }
//    }
//}

//std::vector<VariablePtr> ConstraintBSpline::getComplicatingVariables() const
//{
//    std::vector<VariablePtr> vars;

//    // Only the variables participating in the B-spline are complicating
//    for (unsigned int i = 0; i < variables.size()-numConstraints; i++)
//        vars.push_back(variables.at(i));

//    return vars;
//}

} // namespace CENSO
