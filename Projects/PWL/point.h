/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef POINT_H
#define POINT_H

#include "stdio.h"

template <unsigned int N>
struct Point
{
    // Initialize with origin
    Point() : x({}), y(0) {}

    Point(std::array<double, N> x_) : x(x_), y(0) {}

    Point(std::array<double, N> x_, double y_) : x(x_), y(y_) {}

    // Point coordinate
    std::array<double, N> x;

    // Point value
    double y;

    /*
     * Operator overloading
     */
    bool operator==(const Point<N> &rhs) const;
    bool operator>(const Point<N> &rhs) const;
    bool operator<(const Point<N> &rhs) const;

    template <unsigned int T>
    friend std::ostream &operator<<(std::ostream &outputStream, const Point<T> &point);
};

/*
 * Points are equal if
 * x-coordinates are equal
 */
template <unsigned int N>
bool Point<N>::operator==(const Point<N> &rhs) const
{
    for (unsigned int i = 0; i < N; i++)
        if (x.at(i) != rhs.x.at(i))
            return false;
    return true;
}

template <unsigned int N>
bool Point<N>::operator>(const Point<N> &rhs) const
{
    for (unsigned int i = 0; i < N; i++)
    {
        if (x.at(i) > rhs.x.at(i))
            return true;
        else if (x.at(i) < rhs.x.at(i))
            return false;
    }
    return false;
}


template <unsigned int N>
bool Point<N>::operator<(const Point<N> &rhs) const
{
    for (unsigned int i = 0; i < N; i++)
    {
        if (x.at(i) < rhs.x.at(i))
            return true;
        else if (x.at(i) > rhs.x.at(i))
            return false;
    }
    return false;
}

template <unsigned int T>
std::ostream &operator<<(std::ostream &outputStream, const Point<T> &point)
{
    outputStream << "(";

    bool firstLoop = true;
    for (const auto &xi : point.x)
    {
        if (!firstLoop)
            outputStream << ", ";

        outputStream << xi;
        firstLoop = false;
    }

    outputStream << ")";

    return outputStream;
}

typedef Point<1> Point1D;
typedef Point<2> Point2D;
typedef Point<3> Point3D;

#endif // POINT_H
