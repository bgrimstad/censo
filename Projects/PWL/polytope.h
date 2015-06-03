/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef POLYTOPE_H
#define POLYTOPE_H

#include <stdio.h>
#include <algorithm>
#include "point.h"
#include "Utils/linearsolvers.h"

using namespace CENSO;

/*
 * An N-polytope is defined by a set of N-dimensional points
 * called vertices. The minimum number of points required
 * for construction is N+1, fin which case the polytope is a simplex.
 *
 * For example:
 * A 1-polytope is an edge defined by two 1D points
 * A 2-polytope is a polygon defined by 2D points
 * A 3-polytope is a polyhedron defined by 3D points
 * Etc.
 */
template <unsigned int N>
struct Polytope
{
    Polytope(std::vector<Point<N>> vertices_)
        : vertices(vertices_)//, qHull(QHullInterface())
    {
    }

    //Polytope(Polytope const& copy) = delete;
    Polytope& operator = (Polytope const& assign) = delete;

    std::vector<Point<N>> vertices;
//    QHullInterface qHull;

    /*
     * Check if point is a vertex of the polytope
     */
    bool isVertex(const Point<N> &p) const
    {
        if (std::find(vertices.begin(), vertices.end(), p) == vertices.end())
            return false;
        return true;
    }

    /*
     * Check if point is in polytope
     */

    DenseMatrix computeInterpolationCoefficients() const;

    /*
     * Compute facets of polytope.
     * The facets are represented by
     * half-spaces A with offset b.
     */
//    void computeFacets(DenseMatrix &A, DenseVector &b)
//    {
//        if (!qHull.isInitialized())
//            computeConvexHull();
//        assert(qHull.extractFacets(A, b));
//    }

    /*
     * Compute convex hull of polytope vertices
     * NOTE: does not include y-values
     */
//    void computeConvexHull()
//    {
//        DenseMatrix P = makePointMatrix(vertices, false);
//        qHull.computeConvexHull(P);
//    }

    /*
     * Operator overloading
     */
    template <unsigned int T>
    friend std::ostream &operator<<(std::ostream &outputStream, const Polytope<T> &polytope);
};

/*
 * Compute interpolation coefficients
 * For vertices at (x -> y) solve the linear system
 * with rows [x^T 1] [m^T c]^T = y
 * for m in R^N and c in R so that
 * m^T*x + c = y for all x in the polytope.
 */
template <unsigned int N>
DenseMatrix Polytope<N>::computeInterpolationCoefficients() const
{
    /*
     * Triangulation required for unique solution
     * Polytope must be a simplex with N+1 vertices
     */
    assert(vertices.size() == N+1);

    DenseMatrix A = DenseMatrix::Zero(vertices.size(), N+1);
    DenseMatrix b = DenseMatrix::Zero(vertices.size(), 1);
    int i = 0;
    for (const auto &v : vertices)
    {
        for (unsigned int j = 0; j < N; j++)
            A(i,j) = v.x.at(j);
        A(i,N) = 1;
        b(i) = v.y;
        i++;
    }

    /*
     * Solve A [m c] = b
     */
    DenseMatrix mc;
    DenseQR s;
    if (!s.solve(A, b, mc))
    {
        cout << A << endl;
        cout << b << endl;
        cout << "Could not solve for interpolation coefficients!" << endl;
        exit(1);
    }

    // Last element of mc is c
    return mc;
}

/*
 * Specialization for the special case of 1D polytopes
 */
//template <>
//inline void Polytope<1>::computeFacets(DenseMatrix &A, DenseVector &b)
//{
//    assert(vertices.size() == 2);
//    double xl = vertices.at(0).x.front();
//    double xu = vertices.at(1).x.front();

//    A.resize(2,1);
//    A(0,0) = -1;
//    A(1,0) = 1;
//    b.resize(2,1);
//    b(0,0) = -xl;
//    b(1,0) = xu;
//}

/*
 * Specialization for the special case of 1D polytopes
 */
//template <>
//inline void Polytope<1>::computeConvexHull()
//{
//}

template <unsigned int T>
std::ostream &operator<<(std::ostream &outputStream, const Polytope<T> &polytope)
{
    outputStream << T << "D polytope vertices:" << endl;
    for (const auto &vertex : polytope.vertices)
        outputStream << vertex << endl;

    return outputStream;
}

/*
 * For addYCol = false
 * Makes a matrix P with m rows and N cols,
 * where m is the number of points.
 * For addYCol = true
 * Makes a matrix P with m rows and N+1 cols, where
 * the last column of P contain the y values at the
 * x points.
 */
template <unsigned int N>
DenseMatrix makePointMatrix(const std::vector<Point<N>> &points, bool addYCol = true)
{
    DenseMatrix P;
    if (addYCol)
        P = DenseMatrix::Zero(points.size(), N+1);
    else
        P = DenseMatrix::Zero(points.size(), N);

    int i = 0;
    for (const auto &p : points)
    {
        for (unsigned int j = 0; j < N; j++)
            P(i,j) = p.x.at(j);
        if (addYCol)
            P(i,N) = p.y;
        i++;
    }
    return P;
}

/*
 * Makes a vector of Point<N> from
 * a matrix P with m rows and N+1 columns,
 * where m is the number of points.
 * The last column of P is the y values
 * stored at a point x.
 */
template <unsigned int N>
std::vector<Point<N>> makePointVector(const DenseMatrix &P)
{
    assert(P.cols() == N+1);
    std::vector<Point<N>> points;

    for (int i = 0; i < P.rows(); i++)
    {
        std::array<double, N> x;
        for (unsigned int j = 0; j < N; j++)
            x.at(j) = P(i,j);

        double y = P(i,N);
        points.push_back(Point<N>(x,y));
    }
    return points;
}

/*
 * Builder functions
 */
template <unsigned int N>
Polytope<N> createPolytopeFromHyperRectangle(std::vector<double> lb, std::vector<double> ub)
{
    assert(lb.size() == N);
    assert(ub.size() == N);

    /*
     * Check if digit l of integer b is 1,
     * where l is nonnegative (zero-indexed).
     */
    auto isDigitOne = [](int b, int l)
    {
        assert(l >= 0);
        int mask = 1 << l;
        int res = (b & mask) >> l;
        return (res == 1);
    };

    // Number of vertices of a hyperrectangle is 2^N
    std::vector<Point<N>> vertices;

    // Permutate using N bits (0: lb, 1: ub)
    for (unsigned int i = 0; i < std::pow(2, N); ++i)
    {
        std::array<double, N> x;

        for (unsigned int j = 0; j < N; ++j)
        {
            if (isDigitOne(i, j))
                x.at(j) = ub.at(j);
            else
                x.at(j) = lb.at(j);
        }

        Point<N> p(x);
        vertices.push_back(p);
    }

    Polytope<N> poly(vertices);

    return poly;
}

typedef Polytope<1> Polytope1D;
typedef Polytope<2> Polytope2D;
typedef Polytope<3> Polytope3D;

#endif // POLYTOPE_H
