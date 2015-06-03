/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bsplinepoly.h"
#include "math.h"
#include "set"
#include "unsupported/Eigen/KroneckerProduct"

#include "iostream"
using std::cout;
using std::endl;

namespace CENSO
{

/*
 * Calculate coefficients of B-spline representing a multivariate polynomial
 *
 * The polynomial f(x), with x in R^n, has m terms on the form
 * f(x) = c(0)*x(0)^E(0,0)*x(1)^E(0,1)*...*x(n-1)^E(0,n-1)
 *       +c(1)*x(0)^E(1,0)*x(1)^E(1,1)*...*x(n-1)^E(1,n-1)
 *       +...
 *       +c(m-1)*x(0)^E(m-1,0)*x(1)^E(m-1,1)*...*x(n-1)^E(m-1,n-1)
 * where c in R^m is a vector with coefficients for each of the m terms,
 * and E in N^(mxn) is a matrix with the exponents of each variable in each of the m terms,
 * e.g. the first row of E defines the first term with variable exponents E(0,0) to E(0,n-1).
 *
 * Note: E must be a matrix of nonnegative integers
 */
DenseMatrix getBSplineBasisCoefficients(DenseVector c, DenseMatrix E, std::vector<double> lb, std::vector<double> ub)
{
    unsigned int dim = E.cols();
    unsigned int terms = E.rows();
    assert(dim >= 1); // At least one variable
    assert(terms >= 1); // At least one term (assumes that c is a column vector)
    assert(terms == c.rows());
    assert(dim == lb.size());
    assert(dim == ub.size());

    // Get highest power of each variable
    DenseVector powers = E.colwise().maxCoeff();

    // Store in std vector
    std::vector<unsigned int> powers2;
    for (unsigned int i = 0; i < powers.size(); ++i)
        powers2.push_back(powers(i));

    // Calculate tensor product transformation matrix T
    DenseMatrix T = getTransformationMatrix(powers2, lb, ub);

    // Compute power basis coefficients (lambda vector)
    SparseMatrix L(T.cols(),1);
    L.setZero();

    for (unsigned int i = 0; i < terms; i++)
    {
        SparseMatrix Li(1,1);
        Li.insert(0,0) = 1;

        for (unsigned int j = 0; j < dim; j++)
        {
            int e = E(i,j);
            SparseVector li(powers(j)+1);
            li.reserve(1);
            li.insert(e) = 1;

            SparseMatrix temp = Li;
            Li = kroneckerProduct(temp, li);
        }

        L += c(i)*Li;
    }

    // Compute B-spline coefficients
    DenseMatrix C = T*L;

    return C;
}

/*
 * Compute power basis coefficients from B-spline basis coefficients
 */
DenseMatrix getPowerBasisCoefficients(DenseVector c, std::vector<unsigned int> degrees, std::vector<double> lb, std::vector<double> ub)
{
    // Calculate tensor product transformation matrix T*lambda = c
    DenseMatrix T = getTransformationMatrix(degrees, lb, ub);

    // Compute power basis coefficients from B-spline basis coefficients
    DenseMatrix L = T.colPivHouseholderQr().solve(c);

    return L;
}

/*
 * Transformation matrix taking a power basis to a B-spline basis
 */
DenseMatrix getTransformationMatrix(std::vector<unsigned int> degrees, std::vector<double> lb, std::vector<double> ub)
{
    unsigned int dim = degrees.size();
    assert(dim >= 1); // At least one variable
    assert(dim == lb.size());
    assert(dim == ub.size());

    // Calculate tensor product transformation matrix T
    SparseMatrix T(1,1);
    T.insert(0,0) = 1;

    for (unsigned int i = 0; i < dim; i++)
    {
        SparseMatrix temp(T);
        unsigned int deg = degrees.at(i);
        DenseMatrix Mi = getBSplineToPowerBasisMatrix1D(deg);
        DenseMatrix Ri = getReparameterizationMatrix1D(deg, lb.at(i), ub.at(i));
        DenseMatrix Ti = Mi.colPivHouseholderQr().solve(Ri);

        T = kroneckerProduct(temp, Ti);
    }
    return T;
}

/*
 * Calculates the matrix that transforms a B-spline basis to a power basis
 * Input: p (degree of basis)
 * Output: M (power matrix)
 */
DenseMatrix getBSplineToPowerBasisMatrix1D(unsigned int p)
{
    assert(p > 0); // m = n+p+1

    // M is a lower triangular matrix of size (p+1)x(p+1)
    DenseMatrix M; M.setZero(p+1,p+1);

    for (unsigned int j = 0; j <= p; j++) // cols
    {
        for (unsigned int i = j; i <= p; i++) // rows
        {
            M(i,j) = pow(-1,i-j)*binomialCoeff(p,j)*binomialCoeff(p-j,i-j);
        }
    }
    return M;
}

/*
 * Calculates the reparameterization matrix that changes the domain from [0,1] to [a,b]
 */
DenseMatrix getReparameterizationMatrixND(std::vector<unsigned int> degrees, std::vector<double> lb, std::vector<double> ub)
{
    unsigned int dim = degrees.size();
    assert(dim >= 1); // At least one variable
    assert(dim == lb.size());
    assert(dim == ub.size());

    // Calculate tensor product transformation matrix T
    SparseMatrix R(1,1);
    R.insert(0,0) = 1;

    for (unsigned int i = 0; i < dim; i++)
    {
        SparseMatrix temp(R);
        unsigned int deg = degrees.at(i);
        DenseMatrix Ri = getReparameterizationMatrix1D(deg, lb.at(i), ub.at(i));
        R = kroneckerProduct(temp, Ri);
    }
    return R;
}

/*
 * Calculates the reparameterization matrix that changes the domain of a basis
 * Input: p (degree of basis), a (left bound/first knot), b (right bound/last knot)
 * Output: Rinv (reparameterization matrix that changes the domain from [0,1] to [a,b])
 */
DenseMatrix getReparameterizationMatrix1D(int p, double a, double b)
{
    assert(p > 0 && a < b);

    // R is a upper triangular matrix of size (p+1)x(p+1)
    DenseMatrix R; R.setZero(p+1,p+1);
    for (int i = 0; i <= p; i++) // rows
    {
        for (int j = i; j <= p; j++) // cols
        {
            R(i,j) = binomialCoeff(j,i)*pow(b-a,i)*pow(a,j-i);
        }
    }
    return R;
}

// Returns the binomial coefficient
double binomialCoeff(int n, int k)
{
    assert(n >= 0 && k >= 0 && n >= k);
    std::vector<double> c; c.push_back(1);
    for (double i = 0; i < k; ++i)
        c.push_back((c.back()*(n-i)) / (i+1));
    return c.back();
}

// Get regular knot sequence for one spline piece
std::vector< std::vector<double> > getRegularKnotVectors(std::vector<unsigned int> deg, std::vector<double> lb, std::vector<double> ub)
{
    unsigned int dim = deg.size();
    assert(dim > 0);
    assert(dim == lb.size());
    assert(dim == ub.size());

    std::vector< std::vector<double> > knots;
    for (unsigned int i = 0; i < dim; i++)
    {
        std::vector<double> ks;
        for (unsigned int j = 0; j < deg.at(i)+1; j++)
        {
            ks.push_back(lb.at(i));
        }
        for (unsigned int j = 0; j < deg.at(i)+1; j++)
        {
            ks.push_back(ub.at(i));
        }

        knots.push_back(ks);
    }
    return knots;
}

// Extract unique knots
std::vector<std::vector<double> > getUniqueKnots(std::vector<std::vector<double> > knots)
{
    // Extract unique knots
    std::vector<std::vector<double>> uniqueKnots;
    for (unsigned int i = 0; i < knots.size(); ++i)
    {
        auto knotsi = knots.at(i);

        // Remove duplicates
        std::set<double> s(knotsi.begin(), knotsi.end());
        std::vector<double> knotsu;
        knotsu.assign(s.begin(), s.end());
        uniqueKnots.push_back(knotsu);
    }

    return uniqueKnots;
}

/*
 * Compute power matrix E(m,n) given polynomial degrees.
 */
DenseMatrix getPowersMatrix(std::vector<unsigned int> degrees)
{
    assert(degrees.size() > 0);

    // Number of powers = degree+1
    auto powers = degrees;
    for (auto &pow : powers)
        ++pow;

    return getPermutationMatrix(powers);
}

/*
 * Compute matrix P(m,n) with all possible (m) permutations of
 * a vector x of size n which can have values x(i) in [0,...,num(i)-1]
 */
DenseMatrix getPermutationMatrix(std::vector<unsigned int> num)
{
    unsigned int n = num.size();
    unsigned int m = 1;

    for (unsigned int val : num)
        m *= val;

    DenseMatrix Perm(m, n);

    for (unsigned int i = 0; i < n; ++i)
    {
        DenseVector vec(1); vec(0) = 1;

        for (unsigned int j = 0; j < n; ++j)
        {
            DenseVector temp = vec;
            DenseVector vecj = DenseVector::Ones(num.at(j));

            if (i == j)
            {
                // Values [0,...,num(i)-1]
                for (unsigned int k = 0; k < vecj.size(); ++k)
                    vecj(k) = k;
            }

            vec = Eigen::kroneckerProduct(temp, vecj);
        }

        assert(vec.size() == m);

        for (unsigned int k = 0; k < m; ++k)
        {
            Perm(k,i) = vec(k);
        }
    }

    return Perm;
}

} // namespace CENSO
