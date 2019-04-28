/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CENSO_BSPLINE_WRAPPER_H
#define CENSO_BSPLINE_WRAPPER_H

#include "definitions.h"
#include "bspline_builder.h"


using SPLINTER::BSpline;
using SPLINTER::DataTable;

namespace CENSO {

class BSplineWrap {
public:
    static BSpline build_bspline(DenseMatrix coefficients, std::vector<std::vector<double>> knots,
                                 std::vector<unsigned int> degrees);

    static BSpline fit_bspline(DataTable data, unsigned int degree = 3);
};

} // namespace CENSO

#endif //CENSO_BSPLINE_WRAPPER_H
