/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef FBBT_H
#define FBBT_H

#include "boundstightener.h"
#include "OptimizationProblem/constraint.h"

#include "mutex"

namespace CENSO
{

namespace BB
{

class FBBT : public BoundsTightener
{

public:
    FBBT(double threshold, unsigned int maxIterations)
        : BoundsTightener(threshold, maxIterations)
    {}

    virtual ~FBBT() {}

private:
    bool doTightening(ConstraintPtr constraints) override;
};

} // namespace BB

} // namespace CENSO

#endif // FBBT_H
