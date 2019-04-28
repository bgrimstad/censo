/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bspline_wrapper.h"
#include "Utils/eigen_utils.h"

namespace CENSO {

BSpline BSplineWrap::build_bspline(DenseMatrix coefficients, std::vector<std::vector<double>> knots,
        std::vector<unsigned int> degrees)
{
    auto cp = eigMatToStdVecVec(coefficients);
    return BSpline(cp, knots, degrees);
}

BSpline BSplineWrap::fit_bspline(DataTable data, unsigned int degree)
{
    return BSpline::Builder(data.getDimX(), data.getDimY()).degree(degree).fit(data);
}

} // namespace CENSO