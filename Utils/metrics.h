/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef METRICS_H
#define METRICS_H

#include "Utils/definitions.h"

namespace CENSO
{

/*
 * Various metrics acting on
 * vectors and matrices.
 */

// Relative difference between variable bounds
double relativeChangeIntervals(std::vector<double> lb1, std::vector<double> ub1, std::vector<double> lb2, std::vector<double> ub2);

} // namespace CENSO

#endif // METRICS_H
