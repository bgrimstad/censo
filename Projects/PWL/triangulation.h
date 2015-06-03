/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef DELAUNAY_H
#define DELAUNAY_H

#include "point.h"
#include "polytope.h"
//#include "qhullinterface.h"

/*
 * Implementation of Delaunay triangulation for
 * n-dimensional point sets (not working).
 */
template <unsigned int N>
std::vector<Polytope<N> > delaunay(const std::vector<Point<N>> &points)
{
//    // Project points to the R(N+1) paraboloid
//    // Note: the point value y is used.
//    std::vector<Point<N>> projpoints;
//    for (const auto &p : points)
//    {
//        auto x = p.x;
//        double norm = 0;
//        for (unsigned int i = 0; i < N; i++)
//        {
//            norm += x.at(i)*x.at(i);
//        }
//        projpoints.push_back(Point<N>(x, norm));
//    }

//    // Compute convex hull of the R(N+1) paraboloid
//    QHullInterface qHull;
//    DenseMatrix P = makePointMatrix(projpoints);
//    qHull.computeConvexHull(P);

//    qHull.printQhull();

//    DenseMatrix vertices;
//    qHull.extractVertices(vertices);
//    cout << "Vertices: " << endl;
//    cout << vertices << endl;

//    qHull.printQhull();

//    // Extract facets of convex hull
//    DenseMatrix A;
//    DenseVector b;
//    qHull.extractFacets(A, b);
//    cout << "A" << endl;
//    cout << A << endl;
//    cout << "b" << endl;
//    cout << b << endl;

//    // Find point with the maximum (N+1)th coordinate
//    DenseVector p0 = DenseVector::Zero(N,1);
//    double w0 = 0;
//    for (const auto &p : projpoints)
//    {
//        double w = p.y;
//        if (w > w0)
//        {
//            for (unsigned int i = 0; i < N; i++)
//                p0(i) = p.x.at(i);
//            w0 = w;
//        }
//    }

//    /*
//     * Find the point where the plane tangent to the point (p0,w0) on the
//     * paraboloid crosses the w axis. This is the point that can see the entire
//     * lower hull.
//     */
//    double w_optimal = w0 - 2*p0.squaredNorm();

//    /*
//     * Subtract 1000 times the absolute value of w_optimal to ensure that the
//     * point where the tangent plane crosses the w axis will see all points on
//     * the lower hull. This avoids numerical roundoff errors.
//     */
//    w_optimal -= 1000*std::abs(w_optimal);

//    // Set the point where the tangent plane crosses the w axis
//    DenseVector p = DenseVector::Zero(N+1);
//    p(N) = w_optimal;

//    // Find all faces that are visible from this point
//    // Visible if A*p + b > 0
//    DenseVector visible = A*p+b;
//    cout << "Visible" << endl;
//    cout << visible << endl;

//    cout << "Trying to extract vertices of facets: " << endl;
//    //qHull.extractFacetVertices(A);

//    qHull.printQhull();

//    // Construct the polytopes of the convex hull

//    std::vector<Polytope<N>> polytopes;
//    return polytopes;
}

/*
 * Triangulation
 * Creates a triangulation of a point set.
 * The triangulation is represented as a set of polytopes,
 * that is, f : point set -> family of polytopes.
 * 1D: sorted list
 * nD: Delaunay triangulation
 */
template <unsigned int N>
std::vector<Polytope<N> > triangulate(std::vector<Point<N>> &points)
{
    return delaunay(points);
}

/*
 * Triangulation in 1D is simply a sorting of the points.
 *
 * NOTE: fully specialized so it needs to be inline.
 */
template <>
inline std::vector<Polytope<1> > triangulate<1>(std::vector<Point<1> > &points)
{
    // Sort segments
    std::sort(points.begin(), points.end());

    // Create list of polytopes (line segments in 1D)
    std::vector<Polytope1D> polytopes;
    for (unsigned int i = 0; i < points.size()-1; i++)
    {
        std::vector<Point1D> pts;
        pts.push_back(points.at(i));
        pts.push_back(points.at(i+1));
        polytopes.push_back(Polytope1D(pts));
    }
    return polytopes;
}

/*
 * Triangulation in 2D assuming regular grid.
 *
 * NOTE: fully specialized so it needs to be inline.
 */
template <>
inline std::vector<Polytope<2> > triangulate<2>(std::vector<Point<2> > &points)
{
    // Sort segments
    std::sort(points.begin(), points.end());
    std::vector<double> x1, x2;
    for (const auto &p : points)
    {
        x1.push_back(p.x.front());
        x2.push_back(p.x.back());
    }

    std::sort(x1.begin(), x1.end());
    auto it1 = std::unique(x1.begin(), x1.end());
    x1.resize(std::distance(x1.begin(), it1));

    std::sort(x2.begin(), x2.end());
    auto it2 = std::unique(x2.begin(), x2.end());
    x2.resize(std::distance(x2.begin(), it2));

    // Create list of polytopes (line segments in 2D)
    std::vector<Polytope2D> polytopes;
    for (unsigned int i = 0; i < x1.size()-1; i++)
    {
        for (unsigned int j = 0; j < x2.size()-1; j++)
        {
            // Add two triangles per rectangle
            std::vector<Point<2>> pl1, pl2;

            // Point in first triangle
            Point2D p11({x1.at(i), x2.at(j)});
            auto p11o = std::find(points.begin(), points.end(), p11);
            assert(p11o != points.end());
            pl1.push_back(*p11o);

            Point2D p12({x1.at(i), x2.at(j+1)});
            auto p12o = std::find(points.begin(), points.end(), p12);
            assert(p12o != points.end());
            pl1.push_back(*p12o);

            Point2D p13({x1.at(i+1), x2.at(j)});
            auto p13o = std::find(points.begin(), points.end(), p13);
            assert(p13o != points.end());
            pl1.push_back(*p13o);

            // Point in second triangle
            Point2D p21({x1.at(i+1), x2.at(j)});
            auto p21o = std::find(points.begin(), points.end(), p21);
            assert(p21o != points.end());
            pl2.push_back(*p21o);

            Point2D p22({x1.at(i), x2.at(j+1)});
            auto p22o = std::find(points.begin(), points.end(), p22);
            assert(p22o != points.end());
            pl2.push_back(*p22o);

            Point2D p23({x1.at(i+1), x2.at(j+1)});
            auto p23o = std::find(points.begin(), points.end(), p23);
            assert(p23o != points.end());
            pl2.push_back(*p23o);

            // Add the two polytopes
            // y-values of points are preserved
            polytopes.push_back(pl1);
            polytopes.push_back(pl2);
        }
    }
    return polytopes;
}

#endif // DELAUNAY_H
