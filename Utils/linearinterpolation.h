/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef LINEARINTERPOLATION_H
#define LINEARINTERPOLATION_H

namespace CENSO
{

double linearInterpolation(double x, double x1, double x2, double y1, double y2);
double bilinearInterpolation(double x, double y, double x1, double x2, double y1, double y2, double f11, double f21, double f12, double f22);

} // namespace CENSO

#endif // LINEARINTERPOLATION_H
