/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "metrics.h"

namespace CENSO
{

/*
 * Measures the relative difference in the sum of variable intervals
 * (not counting infinite intervals)
 */
double relativeChangeIntervals(std::vector<double> lb1, std::vector<double> ub1, std::vector<double> lb2, std::vector<double> ub2)
{
    assert(lb1.size() == ub1.size());
    assert(lb2.size() == ub2.size());
    assert(lb1.size() == lb2.size());

    double diff1 = 0;
    double diff2 = 0;

    for (unsigned int i = 0; i < lb1.size(); i++)
    {
        /* NOTE: not counting infinite bounds.
         * This may be a problem when an inifite bound has become finite.
         * An alternative measure is the average relative difference.
         */
        if (ub1.at(i) >= INF
            || lb1.at(i) <= -INF
            || ub2.at(i) >= INF
            || lb2.at(i) <= -INF)
            continue;

        diff1 += ub1.at(i) - lb1.at(i);
        diff2 += ub2.at(i) - lb2.at(i);
    }

    if (diff1 > 0)
        return (1 - diff2/diff1);

    return 0;
}

} // namespace CENSO
