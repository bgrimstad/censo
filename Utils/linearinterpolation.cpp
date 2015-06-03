/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "linearinterpolation.h"
#include <iostream>
#include <cassert>

using std::cout;
using std::endl;

namespace CENSO
{

/* Linear interpolation on line from (x1,y1) to (x2,y2)
 */
double linearInterpolation(double x, double x1, double x2, double y1, double y2)
{
    // Check that x is within x-interval [x1,x2]
    assert(x1 <= x && x <= x2);

    // Interpolate
    if (x1 == x2)
    {
        cout << "linearInterpolation: x1 = x2, returning y1." << endl;
        return y1;
    }
    else return y1 + (x - x1)/(x2 - x1)*(y2 - y1);
}

/* Bilinear interpolation (on regular grid):
 * A function f(x,y) is interpolated at the point (x,y)
 * using the function values at the neighbouring grid points
 * where x1 <= x <= x2 and y1 <= y <= y2
 * and f11 = f(x1,y1), f21 = f(x2,y1), f12 = f(x1,y2), f22 = f(x2,y2)
 */
double bilinearInterpolation(double x, double y, double x1, double x2, double y1, double y2, double f11, double f21, double f12, double f22)
{
    // Check that (x,y) is within the square
    assert(x1 <= x); assert(x <= x2);
    assert(y1 <= y); assert(y <= y2);
    assert(x1 < x2); assert(y1 < y2);

    // Check grid
    assert(x1 < x2);
    assert(y1 < y2);

    // Interpolate in x-direction
    double fr1 = f11*(x2-x)/(x2-x1) + f21*(x-x1)/(x2-x1);
    double fr2 = f12*(x2-x)/(x2-x1) + f22*(x-x1)/(x2-x1);

    // Interpolate in y-direction
    return fr1*(y2-y)/(y2-y1) + fr2*(y-y1)/(y2-y1);
}

} // namespace CENSO
