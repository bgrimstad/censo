/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "Utils/definitions.h"

namespace CENSO
{

bool isInteger(double value)
{
    return assertNear(value, std::floor(value));
    //return (value == std::floor(value)); // Fast, to int precision
    // return (std::abs(value - round(value)) < 1e-9); // Slow, specified precision
}

// Returns a random integer in the range [min, max]
int randomInteger(int min, int max)
{
    if (!randomSeedSet)
    {
        std::srand(time(NULL)); // Random seed
        randomSeedSet = true;
    }
    return std::rand() % (max - min + 1) + min;
}

std::vector<double> linspace(double start, double stop, unsigned int points)
{
    std::vector<double> ret;
    double dx = 0;
    if (points > 1)
        dx = (stop - start)/(points-1);
    for (unsigned int i = 0; i < points; ++i)
        ret.push_back(start + i*dx);
    return ret;
}

} // namespace CENSO
