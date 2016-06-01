/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintbspline.h"
#include "constraintlinear.h"
#include "constraintset.h"
//#include "SolverInterface/solveripopt.h"
#include "unsupported/Eigen/KroneckerProduct"

using std::cout;
using std::endl;

namespace CENSO
{

ConstraintBSpline::ConstraintBSpline(std::vector<VariablePtr> variables, const BSpline &bspline, bool equality)
    : Constraint(variables),
      bspline(bspline),
      equality(true), // Ineq not implemented!
      maxNumAuxiliaryVariables(1e9) //1e6
{
    init(equality);
}

ConstraintBSpline::ConstraintBSpline(const ConstraintBSpline &copy, bool deep)
    : Constraint(copy, deep),
      bspline(copy.bspline),
      equality(copy.equality),
      maxNumAuxiliaryVariables(copy.maxNumAuxiliaryVariables),
      storedVariables(copy.storedVariables)
{
}

void ConstraintBSpline::init(bool equality)
{
    numConstraints = 1;
    assert(variables.size() == bspline.getNumVariables() + numConstraints);

    // Update variable bounds
    auto varlb = bspline.getDomainLowerBound();
    auto varub = bspline.getDomainUpperBound();

    for (unsigned int i = 0; i < bspline.getNumVariables(); i++)
    {
        variables.at(i)->updateBounds(varlb.at(i), varub.at(i));
    }

    // Set bound on y
//    for (unsigned int i = bspline.getNumVariables(); i < variables.size(); i++)
//    {
//        variables.at(i)->updateBounds(-INF, INF);
//    }

    for (unsigned int i = 0; i < numConstraints; i++)
    {
        if (equality)
        {
            lb.push_back(0);
        }
        else
        {
            lb.push_back(-INF);
        }

        ub.push_back(0);

        constraintTypes.push_back(ConstraintType::NONLINEAR_NONCONVEX);
    }

    nnzJacobian = bspline.getNumVariables() + 1;
    nnzHessian = 0;

    jacobianCalculated = true;
    hessianCalculated = false;
    constraintLinear = false;
    constraintConvex = false;
    convexRelaxationAvailable = true;

    for (unsigned int row = 0; row < bspline.getNumVariables(); row++)
        for (unsigned int col = 0; col <= row; col++)
            nnzHessian++;

    constraintName = "B-spline";

    // Update bounds
    reduceVariableRanges();

    checkConstraintSanity();
}

DenseVector ConstraintBSpline::eval(const DenseVector &x) const
{
    // Only x-variables (B-spline input variables) are adjusted to bounds
    DenseVector xy = x;
    DenseVector xa = adjustToDomainBounds(x);

    DenseVector xx = xa.block(0, 0, variables.size()-numConstraints, 1);
    DenseVector yy = xy.block(variables.size()-numConstraints, 0, numConstraints, 1);

    double by = bspline.eval(xx);

    DenseVector y = DenseVector::Zero(numConstraints);
    y(0) = by - yy(0);

    //    adjustToDomainBounds(x);
    //    y.resize(numConstraints);

    //    VecD xx = x.block(0,0,variables.size()-numConstraints,1);
    //    VecD yy = x.block(variables.size()-numConstraints,0,numConstraints,1);

    //    VecD by = bspline->evaluate(xx);

    //    y = by - yy;

    return y;
}

DenseVector ConstraintBSpline::evalJacobian(const DenseVector &x) const
{
    DenseVector xa = adjustToDomainBounds(x);
    DenseVector dx = DenseVector::Zero(nnzJacobian);

    //return centralDifference(xa);

    // Get x-values
    DenseVector xx = xa.block(0, 0, bspline.getNumVariables(), 1);

    // Evaluate B-spline Jacobian
    DenseMatrix jac = bspline.evalJacobian(xx);

    // Derivatives on inputs x
    int k = 0;
    for (int i = 0; i < jac.rows(); i++)
    {
        for (int j = 0; j < jac.cols(); j++)
        {
            dx(k++) = jac(i,j);
        }
    }

    // Derivatives on outputs y
    for (unsigned int i = 0; i < numConstraints; i++)
    {
        dx(k++) = -1;
    }

    return dx;
}

DenseVector ConstraintBSpline::evalHessian(const DenseVector &x) const
{
    DenseVector xa = adjustToDomainBounds(x);
    DenseVector ddx = DenseVector::Zero(nnzHessian);

    // Get x-values
    DenseVector xx = xa.block(0, 0, bspline.getNumVariables(), 1);

    // Calculate Hessian
    DenseMatrix H = bspline.evalHessian(xx);

    // H is symmetric so fill out lower left triangle only
    int idx = 0;
    for (int row = 0; row < H.rows(); row++)
    {
        for (int col = 0; col <= row; col++)
        {
            //if (H(row,col) != 0)
            //{
                ddx(idx++) = H(row,col);
            //}
        }
    }

    return ddx;
}

void ConstraintBSpline::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{
    for (unsigned int i = 0; i < numConstraints; i++)
    {
        for (unsigned int j = 0; j < variables.size(); j++)
        {
            iRow.push_back(i);
            jCol.push_back(j);
        }
    }
}

void ConstraintBSpline::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    // H is symmetric so fill out lower left triangle only
    // Also, not neccessary to fill out for outputs y (=0)
    for (unsigned int row = 0; row < bspline.getNumVariables(); row++)
    {
        for (unsigned int col = 0; col <= row; col++)
        {
            eqnr.push_back(0);
            iRow.push_back(row);
            jCol.push_back(col);
        }
    }
}

bool ConstraintBSpline::reduceVariableRanges() const
{
    // Update bound of y in f(x) - y = 0
    auto cp = bspline.getControlPoints();
    DenseVector minControlPoints = cp.colwise().minCoeff();
    DenseVector maxControlPoints = cp.colwise().maxCoeff();

    for (unsigned int i = variables.size()-numConstraints; i < variables.size(); i++)
    {
        assert(minControlPoints(i) <= maxControlPoints(i));

        // Fix for B-splines with collapsed bounds
        double xlb = variables.at(i)->getLowerBound();
        double xub = variables.at(i)->getUpperBound();

        if (assertNear(xlb, xub))
            continue;

        double newlb = std::max(xlb, minControlPoints(i));
        double newub = std::min(xub, maxControlPoints(i));

        // Detect constraint infeasibility or update bounds
        if (!variables.at(i)->updateBounds(newlb, newub))
            return false;
    }

    // Compute and check variable bounds
    return controlPointBoundsDeduction();
}

std::vector<VariablePtr> ConstraintBSpline::getComplicatingVariables() const
{
    std::vector<VariablePtr> vars;

    // Only the variables participating in the B-spline are complicating
    for (unsigned int i = 0; i < variables.size()-numConstraints; i++)
        vars.push_back(variables.at(i));

    return vars;
}

ConstraintPtr ConstraintBSpline::getConvexRelaxation()
{
    if (constraintConvex)
        return this->clone(false);

    if (haveVariableBoundsChanged())
    {
        storeVariables();

        reduceBSplineDomain();
    }

    // Compute B-spline relaxation
    if (bspline.getNumControlPoints() > maxNumAuxiliaryVariables)
        return computeRelaxationHyperrectangle();

    return computeRelaxationConvexHull();
}

void ConstraintBSpline::reduceBSplineDomain()
{
    std::vector<double> varlb, varub;
    for (unsigned int i = 0; i < bspline.getNumVariables(); i++)
    {
        varlb.push_back(variables.at(i)->getLowerBound());
        varub.push_back(variables.at(i)->getUpperBound());
    }

    std::vector<double> bslb = bspline.getDomainLowerBound();
    std::vector<double> bsub = bspline.getDomainUpperBound();

    // Hack for fixed input variables
    for (unsigned int i = 0; i < varlb.size(); i++)
    {
        if (assertNear(varlb.at(i), varub.at(i)))
        {
            /*
             * NOTE: Expand B-spline domain to avoid knot multiplicity
             * (the B-spline cannot have an empty domain).
             *
             * This is especially important when doing integer optimization,
             * where bounds often are collapsed (lb = ub).
             *
             * The bound threshold should be a very small number, but not
             * too small! It must allow a at least 1000 distinct knot values
             * between the lower and upper bound!
             */
            double boundThreshold = 100000*std::numeric_limits<double>::epsilon();
            varlb.at(i) = std::max(bslb.at(i), varlb.at(i)-boundThreshold);
            varub.at(i) = std::min(bsub.at(i), varub.at(i)+boundThreshold);
        }
    }

    // Reduce domain of B-spline
    bspline.reduceSupport(varlb, varub);

    // Refinement for low dimensional B-spline
    //if (bspline.getNumVariables() <= 2)
        bspline.globalKnotRefinement();
}

void ConstraintBSpline::localRefinement(const DenseVector &x)
{
    DenseVector xadj = adjustToDomainBounds(x);
    bspline.localKnotRefinement(xadj);
}

ConstraintPtr ConstraintBSpline::computeRelaxationHyperrectangle()
{
    /*
     * Hyperrectangle model:
     *
     * X =   [x' y']'
     *
     * A =   [-I ]   b = [-lb ]
     *       [ I ]       [ ub ]
     *
     * AX <= b
     */

    int dim = bspline.getNumVariables() + 1;
    assert(dim == (int)variables.size());

    auto cp = bspline.getControlPoints();
    DenseVector minControlPoints = cp.colwise().minCoeff();
    DenseVector maxControlPoints = cp.colwise().maxCoeff();

    DenseMatrix Idim;
    Idim.setIdentity(dim,dim);

    DenseMatrix A(2*dim, dim);
    A.block(0,0, dim, dim) = - Idim;
    A.block(dim,0, dim, dim) = Idim;

    DenseVector b(2*dim);

    b.block(0, 0, dim, 1) = - minControlPoints;
    b.block(dim, 0, dim, 1) = maxControlPoints;

    ConstraintPtr relaxedConstraint = std::make_shared<ConstraintLinear>(variables, A ,b, false);
    relaxedConstraint->setName("B-spline hypercube relaxation (linear)");

    return relaxedConstraint;
}

ConstraintPtr ConstraintBSpline::computeRelaxationConvexHull()
{
    /*
     * Convex combination model:
     *
     * X = [x' y' auxiliary variables]'
     *
     * A = [ I -C ]   b = [ 0 ]
     *     [ 0  1 ]       [ 1 ]
     *
     * AX = b, [lb' 0]' <= X <= [ub' inf]'
     */

    // Consider renaming here (control points matrix was transposed in old implementation)
    int rowsC = bspline.getNumVariables() + 1;
    int colsC = bspline.getNumControlPoints();
    int rowsA = rowsC + 1;
    int colsA = rowsC + colsC;

    DenseMatrix I;
    I.setIdentity(rowsC,rowsC);

    DenseMatrix zeros;
    zeros.setZero(1,rowsC);

    DenseMatrix ones;
    ones.setOnes(1,colsC);

    DenseMatrix A(rowsA, colsA);
    A.block(0,     0,     rowsC, rowsC) = I;
    A.block(0,     rowsC, rowsC, colsC) = -bspline.getControlPoints().transpose();
    A.block(rowsC, 0,     1,     rowsC) = zeros;
    A.block(rowsC, rowsC, 1,     colsC) = ones;

    zeros.setZero(rowsC,1);
    ones.setOnes(1,1);

    DenseVector b(rowsA);
    b.block(0,     0, rowsC, 1) = zeros;
    b.block(rowsC, 0, 1,     1) = ones;

    auto auxVariables = variables;

    // Number of auxiliary variables equals number of control points
    for (int i = 0; i < bspline.getNumControlPoints(); i++)
        auxVariables.push_back(std::make_shared<Variable>(0, 0, 1));

    ConstraintPtr relaxedConstraint = std::make_shared<ConstraintLinear>(auxVariables, A ,b, true);
    relaxedConstraint->setName("B-spline convex hull relaxation");

    return relaxedConstraint;
}

//bool ConstraintBSpline::testConvexHullMembership(const DenseMatrix& C, const DenseVector& p) const
//{
//    //  Test if point p is in convex hull of the points in C
//    //
//    //  A = [ C ]   b = [ p ]
//    //      [ 1 ]       [ 1 ]
//    //
//    //  p is a member if Ax = b, s.t. x >= 0, has a solution.

//    // Create matrix with extremal control points and a row of ones
//    int rowsC = C.rows();
//    int colsC = C.cols();

//    assert(p.rows() == rowsC);

//    DenseMatrix ones;
//    ones.setOnes(1, colsC);

//    DenseMatrix A(rowsC+1, colsC);
//    A.block(0, 0,  rowsC,  colsC) = C;
//    A.block(rowsC, 0,  1,  colsC) = ones;

//    DenseVector b(rowsC+1);
//    b.block(0,0,rowsC,1) = p;
//    b(rowsC) = 1;

//    // Objective (feasibility problem)
//    DenseMatrix c; c.setOnes(1, colsC);
//    ObjectivePtr obj(new ObjectiveLinear(c));

//    // x >= 0
//    std::vector<double> lb(colsC,0);
//    std::vector<double> ub(colsC,1); // Could also be set to +inf
//    std::vector<double> z0(colsC,0);

//    ConstraintSetPtr cHull(new ConstraintComposite(colsC, lb, ub));
//    cHull->add(new ConstraintLinear(A, b, true));

//    // Determine membership of p in conv(C)
//    ProblemPtr prob(new Problem(obj,cHull));
//    SolverIpopt ipopt(prob, z0);
//    ipopt.setBoundRelaxFactor(1e-8); // Important since problem is very close to infeasibility
//    SolverResult result = ipopt.optimize();

//    return (result.status == SolverStatus::OPTIMAL);
//}

bool ConstraintBSpline::controlPointBoundsDeduction() const
{
    // Get variable bounds
    auto xlb = bspline.getDomainLowerBound();
    auto xub = bspline.getDomainUpperBound();

    // Use these instead?
//    for (unsigned int i = 0; i < bspline.getNumVariables(); i++)
//    {
//        xlb.at(i) = variables.at(i)->getLowerBound();
//        xub.at(i) = variables.at(i)->getUpperBound();
//    }

    double lowerBound = variables.back()->getLowerBound(); // f(x) = y > lowerBound
    double upperBound = variables.back()->getUpperBound(); // f(x) = y < upperBound

    // Get knot vectors and basis degrees
    auto knotVectors = bspline.getKnotVectors();
    auto basisDegrees = bspline.getBasisDegrees();

    // Compute n value for each variable
    // Total number of control points is ns(0)*...*ns(d-1)
    std::vector<unsigned int> numBasisFunctions = bspline.getNumBasisFunctionsPerVariable();

    // Get matrix of coefficients
//    DenseMatrix cps = bspline.getControlPoints().transpose();
//    DenseMatrix coeffs = cps.block(bspline.getNumVariables(), 0, 1, cps.cols());
    DenseMatrix coeffs = bspline.getCoefficients().transpose();

    for (unsigned int d = 0; d < bspline.getNumVariables(); d++)
    {
        if (assertNear(xlb.at(d), xub.at(d)))
            continue;

        auto n = numBasisFunctions.at(d);
        auto p = basisDegrees.at(d);
        std::vector<double> knots = knotVectors.at(d);
        assert(knots.size() == n+p+1);

        // Tighten lower bound
        unsigned int i = 1;
        for (; i <= n; i++)
        {
            // Knot interval of interest: [t_0, t_i]

            // Selection matrix
            DenseMatrix S = DenseMatrix::Ones(1,1);

            for (unsigned int d2 = 0; d2 < bspline.getNumVariables(); d2++)
            {
                DenseMatrix temp(S);

                DenseMatrix Sd_full = DenseMatrix::Identity(numBasisFunctions.at(d2), numBasisFunctions.at(d2));
                DenseMatrix Sd(Sd_full);
                if (d == d2)
                    Sd = Sd_full.block(0,0,n,i);

                S = kroneckerProduct(temp, Sd);
            }

            // Control points that have support in [t_0, t_i]
            DenseMatrix selc = coeffs*S;
            DenseVector minCP = selc.rowwise().minCoeff();
            DenseVector maxCP = selc.rowwise().maxCoeff();
            double minv = minCP(0);
            double maxv = maxCP(0);

            // Investigate feasibility
            if (minv > upperBound || maxv < lowerBound)
                continue; // infeasible
            else
                break; // feasible
        }

        // New valid lower bound on x(d) is knots(i-1)
        if (i > 1)
        {
            if (!variables.at(d)->updateLowerBound(knots.at(i-1)))
                return false;
        }

        // Tighten upper bound
        i = 1;
        for (; i <= n; i++)
        {
            // Knot interval of interest: [t_{n+p-i}, t_{n+p}]

            // Selection matrix
            DenseMatrix S = DenseMatrix::Ones(1,1);

            for (unsigned int d2 = 0; d2 < bspline.getNumVariables(); d2++)
            {
                DenseMatrix temp(S);

                DenseMatrix Sd_full = DenseMatrix::Identity(numBasisFunctions.at(d2), numBasisFunctions.at(d2));
                DenseMatrix Sd(Sd_full);
                if (d == d2)
                    Sd = Sd_full.block(0,n-i,n,i);

                S = kroneckerProduct(temp, Sd);
            }

            // Control points that have support in [t_{n+p-i}, t_{n+p}]
            DenseMatrix selc = coeffs*S;
            DenseVector minCP = selc.rowwise().minCoeff();
            DenseVector maxCP = selc.rowwise().maxCoeff();
            double minv = minCP(0);
            double maxv = maxCP(0);

            // Investigate feasibility
            if (minv > upperBound || maxv < lowerBound)
                continue; // infeasible
            else
                break; // feasible
        }

        // New valid lower bound on x(d) is knots(n+p-(i-1))
        if (i > 1)
        {
            if (!variables.at(d)->updateUpperBound(knots.at(n+p-(i-1))))
                return false;
            // NOTE: the upper bound seems to not be tight! can we use knots.at(n+p-i)?
        }

    }

    return true;
}

void ConstraintBSpline::storeVariables()
{
    storedVariables.clear();
    for (const auto var : variables)
        storedVariables.push_back(*var);
}

bool ConstraintBSpline::haveVariableBoundsChanged() const
{
    if (variables.size() != storedVariables.size())
        return true;

    for (unsigned int i = 0; i < variables.size()-1; i++)
    {
        auto var = *variables.at(i);
        auto svar = storedVariables.at(i);

        if (var.getUpperBound() != svar.getUpperBound()
                || var.getLowerBound() != svar.getLowerBound())
            return true;
    }

    return false;
}

DenseVector ConstraintBSpline::centralDifference(const DenseVector &x) const
{
    DenseVector dx = DenseVector::Zero(nnzJacobian);

    double h = 1e-6; // perturbation step size

    int k = 0;
    for (unsigned int i = 0; i < getNumVariables(); i++)
    {
        double hForward = 0.5*h;
        DenseVector xForward(x);
        if (xForward(i) + hForward > variables.at(i)->getUpperBound())
        {
            hForward = 0;
        }
        else
        {
            xForward(i) = xForward(i) + hForward;
        }

        double hBackward = 0.5*h;
        DenseVector xBackward(x);
        if (xBackward(i) - hBackward < variables.at(i)->getLowerBound())
        {
            hBackward = 0;
        }
        else
        {
            xBackward(i) = xBackward(i) - hBackward;
        }

        DenseVector yForward = eval(xForward);
        DenseVector yBackward = eval(xBackward);

        for (unsigned int j = 0; j < numConstraints; ++j)
        {
            dx(k++) = (yForward(j) - yBackward(j))/(hBackward + hForward);
        }
    }

    return dx;
}

} // namespace CENSO
