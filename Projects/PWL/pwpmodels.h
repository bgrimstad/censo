/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef PWPMODELS_H
#define PWPMODELS_H

#include "Utils/definitions.h"
#include "Utils/bsplinepoly.h"
#include "Projects/PWL/polytope.h"
#include "OptimizationProblem/constraint.h"
#include "OptimizationProblem/constraintpolynomial.h"
#include "bspline.h"
#include "unsupported/Eigen/KroneckerProduct"

using namespace CENSO;

template <unsigned int N>
struct PiecewisePolynomial
{
    PiecewisePolynomial(std::vector<std::shared_ptr<ConstraintPolynomial>> polynomials_, std::vector<Polytope<N>> polytopes_)
        : polynomials(polynomials_),
          polytopes(polytopes_)
    {
        assert(polynomials.size() == polytopes.size());
    }

    // Polynomials in polytopes
    std::vector<std::shared_ptr<ConstraintPolynomial>> polynomials;
    std::vector<Polytope<N>> polytopes;
};

SparseMatrix getSelectionMatrix(unsigned int k, unsigned int degree, unsigned int knotIntervals);

/*
 * Decompose B-spline into set of polytopes (domain) and polynomials (power basis).
 * NOTE: Using dummy variables in the polynomial constraints.
 */
template <unsigned int N>
PiecewisePolynomial<N> decomposeBSplinePower(SPLINTER::BSpline bspline)
{
    assert(bspline.getNumVariables() == N);

    // Decompose B-spline
    bspline.decomposeToBezierForm();
    auto knots = bspline.getKnotVectors();

    // Extract unique knots
    std::vector<std::vector<double>> uniqueKnots = getUniqueKnots(knots);

    // Find bounds on x
    std::vector<double> xlb, xub;
    for (auto knotvec : uniqueKnots)
    {
        xlb.push_back(knotvec.front());
        xub.push_back(knotvec.back());
    }

    // Get B-spline coefficients
    DenseMatrix cpoints = bspline.getControlPoints();
    DenseMatrix coeffs = cpoints.block(bspline.getNumVariables(),0,1,cpoints.cols());

    // Compute number of knot spans
    std::vector<unsigned int> numKnotSpans;
    for (unsigned int i = 0; i < uniqueKnots.size(); ++i)
        numKnotSpans.push_back(uniqueKnots[i].size()-1);

    // Index each hyperrectangle using a permutation matrix
    DenseMatrix Mhr = getPermutationMatrix(numKnotSpans);

    // Get power matrix for polynomials
    std::vector<unsigned int> degrees = bspline.getBasisDegrees();
    DenseMatrix Mpow = getPowersMatrix(degrees);

    // Compute number of variables
    // Order of variables: [x z x_P z_P lambda_{P,v} y_P]
    int npol = Mhr.rows(); // Number of polytopes

    // For each hyperrectangle compute corresponding polynomial
    std::vector<std::shared_ptr<ConstraintPolynomial>> polynomials;
    std::vector<Polytope<N>> polytopes;

    for (unsigned int i = 0; i < npol; ++i)
    {
        // Get hyperrectangle bounds
        std::vector<double> lb, ub;

        for (unsigned int j = 0; j < Mhr.cols(); ++j)
        {
            unsigned int k = Mhr(i,j);
            lb.push_back(uniqueKnots[j][k]);
            ub.push_back(uniqueKnots[j][k+1]);
        }

        // Compute basis function selection matrix for hyperrectangle i
        SparseMatrix S(1,1);
        S.insert(0,0) = 1;

        for (unsigned int j = 0; j < Mhr.cols(); ++j)
        {
            SparseMatrix temp = S;

            unsigned int k = Mhr(i,j);

            SparseMatrix Sk = getSelectionMatrix(k, degrees.at(j), numKnotSpans.at(j));

            S = Eigen::kroneckerProduct(temp, Sk);
        }

        // Multiply with B-spline coefficients to get relevant coefficients
        assert(coeffs.cols() == S.rows());
        DenseMatrix Ci = coeffs*S;
        assert(Ci.rows() == 1);

        // Multiply with transformation matrix to get power basis coefficients
        DenseMatrix Cit = Ci.transpose();
        DenseMatrix Li = getPowerBasisCoefficients(Cit, degrees, lb, ub);

        // Want power basis on [0,1]
        {
            // Calculate bounds on t variable (in [0,1])
            std::vector<double> tlb, tub;

            for (int k = 0; k < N; ++k)
            {
                double dx = xub.at(k) - xlb.at(k);
                assert(dx > 0);
                tlb.push_back((lb.at(k)-xlb.at(k))/dx);
                tub.push_back((ub.at(k)-xlb.at(k))/dx);
            }

            // From [lb,ub] to [0,1]
            DenseMatrix R1 = getReparameterizationMatrixND(degrees, lb, ub);

            // From [0,1] to [tlb,tub]
            DenseMatrix R2 = getReparameterizationMatrixND(degrees, tlb, tub);

            DenseMatrix Li2 = R2.inverse()*(R1*Li);
            Li = Li2;
        }

        // Clean up Li
        for (unsigned int j = 0; j < Li.rows(); ++j)
            if (assertNear(Li(j,0), 0.0))
                Li(j,0) = 0;

        //cout << Li.transpose() << endl;

        // Should remove zero monomials from Li and Mpow!
//        cout << Li << endl;
//        cout << "--" << endl;
//        cout << Mpow << endl;
//        cout << "--" << endl;

        // Create polynomial for hyperrectangle
        // Need to add one term for z variable
        DenseMatrix MpowExt = DenseMatrix::Zero(Mpow.rows()+1, Mpow.cols()+1);
        MpowExt.block(0, 0, Mpow.rows(), Mpow.cols()) = Mpow;
        MpowExt(Mpow.rows(), Mpow.cols()) = 1; // z

        DenseMatrix LiExt = DenseMatrix::Zero(Li.rows()+1, Li.cols());
        LiExt.block(0, 0, Li.rows(), Li.cols()) = Li;
        LiExt(Li.rows(), 0) = -1; // -z

        // Dummy variables
        Variables pVars;
        for (unsigned int j = 0; j < N; ++j)
            pVars.push_back(std::make_shared<Variable>(0)); // x_P
        pVars.push_back(std::make_shared<Variable>(0)); // z_P

        auto poly = std::make_shared<ConstraintPolynomial>(pVars, LiExt, MpowExt, false);
        polynomials.push_back(poly);

        // Create polytope for hyperrectangle
        Polytope<N> pol = createPolytopeFromHyperRectangle<N>(lb, ub);
        polytopes.push_back(pol);
    }

    PiecewisePolynomial<N> pwp(polynomials, polytopes);

    return pwp;
}

/*
 * DCC model with additional x_P variables
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDCC(Variables vars, SPLINTER::BSpline bspline);

/*
 * DCC model
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDCC2(Variables vars, SPLINTER::BSpline bspline);

/*
 * DLog model with additional x_P variables
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog(Variables vars, SPLINTER::BSpline bspline);

/*
 * DLog model
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog2(Variables vars, SPLINTER::BSpline bspline);

/*
 * DLog model basis functions
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog3(Variables vars, SPLINTER::BSpline bspline);

/*
 * DLog model basis functions and McCormick
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog4(Variables vars, SPLINTER::BSpline bspline);

/*
 * DLog model basis functions and McCormick
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog5(Variables vars, SPLINTER::BSpline bspline, bool equality = false);

/*
 * DLog model basis functions and McCormick
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog6(Variables vars, SPLINTER::BSpline bspline, bool equality = false);

/*
 * DLog model using basis functions and McCormick, and assuming a rectangular grid
 */
template <unsigned int N>
ConstraintPtr constraintPiecewisePolynomialDLog7(Variables vars, SPLINTER::BSpline bspline, bool equality = false);

void testSplineDecomposition();

unsigned int polyTest(int numSamples = 10);

void michaTest();

void pumpSynthesis2(unsigned int grid);

#endif // PWPMODELS_H
