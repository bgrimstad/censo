/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "fbbt.h"
#include "Utils/metrics.h"

namespace CENSO
{

namespace BB
{

bool FBBT::doTightening(ConstraintPtr constraints)
{
    // Reduce variable ranges
    return constraints->reduceVariableRanges();
}

} // namespace BB

} // namespace CENSO
