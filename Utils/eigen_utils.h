/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CENSO_EIGEN_UTILS_H
#define CENSO_EIGEN_UTILS_H

#include "definitions.h"

namespace CENSO {

// Converters
std::vector<double> eigenToStdVec(const DenseVector vec);
DenseVector stdToEigenVec(const std::vector<double> vec);
std::vector<std::vector<double>> eigMatToStdVecVec(const DenseMatrix &mat);
DenseMatrix stdVecVecToEigMat(const std::vector<std::vector<double>> &vec);

} // namespace CENSO

#endif //CENSO_EIGEN_UTILS_H
