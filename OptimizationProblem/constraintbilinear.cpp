/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintbilinear.h"
#include "constraintlinear.h"

using std::cout;
using std::endl;

namespace CENSO
{

ConstraintBilinear::ConstraintBilinear(const ConstraintBilinear &copy, bool deep)
    : Constraint(copy, deep), a(copy.a)
{
}

ConstraintBilinear::ConstraintBilinear(std::vector<VariablePtr> variables, double a, double lb, double ub)
    : Constraint(variables), a(a)
{
    assert(variables.size() == 3);
    numConstraints = 1;

    this->lb.push_back(lb);
    this->ub.push_back(ub);

    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = true;

    nnzJacobian = 3;
    nnzHessian  = 1;

    constraintName = "Constraint Bilinear";

    checkConstraintSanity();
}

DenseVector ConstraintBilinear::eval(const DenseVector &x) const
{
    DenseVector y = DenseVector::Zero(numConstraints);
    y(0) = a*x(0)*x(1) - x(2);
    return y;
}

DenseVector ConstraintBilinear::evalJacobian(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);
    dx(0) = ( a*x(1) );
    dx(1) = ( a*x(0) );
    dx(2) = ( -1 );
    return dx;
}

DenseVector ConstraintBilinear::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);
    ddx(0) = ( a ); // Hessian is symmetric
    //ddx(1) = ( a );
    return ddx;
}

void ConstraintBilinear::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    iRow.push_back(0); jCol.push_back(0);
    iRow.push_back(0); jCol.push_back(1);
    iRow.push_back(0); jCol.push_back(2);
}

void ConstraintBilinear::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    eqnr.push_back(0); iRow.push_back(1); jCol.push_back(0);
    //eqnr.push_back(0); iRow.push_back(0); jCol.push_back(1);

}

std::vector<VariablePtr> ConstraintBilinear::getComplicatingVariables() const
{
    std::vector<VariablePtr> vars;
    vars.push_back(variables.at(0));
    vars.push_back(variables.at(1));
    return vars;
}

// TODO: Only implemented for constraints a*x0*x1 = x2 (lb = ub = 0)
bool ConstraintBilinear::reduceVariableRanges() const
{
    // Calculate bounds on x2 from bounds on x0 and x1 using interval analysis
    double lb0 = variables.at(0)->getLowerBound();
    double lb1 = variables.at(1)->getLowerBound();
    double ub0 = variables.at(0)->getUpperBound();
    double ub1 = variables.at(1)->getUpperBound();
    std::vector<double> endpoints = {lb0*lb1, lb0*ub1, ub0*lb1, ub0*ub1};
    double newlb = *std::min_element(endpoints.begin(), endpoints.end());
    double newub = *std::max_element(endpoints.begin(), endpoints.end());

    if (a >= 0)
    {
        newlb = newlb*a;
        newub = newub*a;
    }
    else
    {
        // Negative a, flip interval
        double temp = newlb;
        newlb = newub*a;
        newub = temp*a;
    }

    return variables.at(2)->updateBounds(newlb, newub);
}

ConstraintPtr ConstraintBilinear::getConvexRelaxation()
{
    /* Return convex hull of bilinear constraint
     * Ax <= b, where
     *
     * A =
     * | a*xmin(1),     a*xmin(0),      -1 |
     * | a*xmax(1),     a*xmax(0),      -1 |
     * | -a*xmin(1),    -a*xmax(0),      1 |
     * | -a*xmax(1),    -a*xmin(0),      1 |
     *
     * b = [ub + a*xmin(0)*xmin(1)
     *      ub + a*xmax(0)*xmax(1)
     *      -lb - a*xmax(0)*xmin(1)
     *      -lb - a*xmin(0)*xmax(1)]
     */

    double lb0 = variables.at(0)->getLowerBound();
    double ub0 = variables.at(0)->getUpperBound();

    double lb1 = variables.at(1)->getLowerBound();
    double ub1 = variables.at(1)->getUpperBound();

    if (lb0 <= -INF
        || lb1 <= -INF
        || ub0 >= INF
        || ub1 >= INF)
    {
        cout << "Cannot relax bilinear constraint with unbounded domain!" << endl;
        exit(1);
    }

    DenseMatrix Ac(4,3); Ac.setZero();
    Ac(0,0) = a*lb1;     Ac(0,1) = a*lb0;     Ac(0,2) = -1;
    Ac(1,0) = a*ub1;     Ac(1,1) = a*ub0;     Ac(1,2) = -1;
    Ac(2,0) = -a*lb1;    Ac(2,1) = -a*ub0;    Ac(2,2) = 1;
    Ac(3,0) = -a*ub1;    Ac(3,1) = -a*lb0;    Ac(3,2) = 1;

    DenseVector bc(4); bc.setZero();
    bc(0) = ub.at(0) + a*lb0*lb1;
    bc(1) = ub.at(0) + a*ub0*ub1;
    bc(2) = -lb.at(0) - a*ub0*lb1;
    bc(3) = -lb.at(0) - a*lb0*ub1;

    return std::make_shared<ConstraintLinear>(variables, Ac, bc, false);
}

} // namespace CENSO
