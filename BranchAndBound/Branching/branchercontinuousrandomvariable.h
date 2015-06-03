/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHERCONTINUOUSRANDOMVARIABLE_H
#define BRANCHERCONTINUOUSRANDOMVARIABLE_H

#include "branchercontinuous.h"

namespace CENSO
{

namespace BB
{

class BrancherContinuousRandomVariable : public BrancherContinuous
{
public:
    BrancherContinuousRandomVariable() : BrancherContinuous() {}
    ~BrancherContinuousRandomVariable() {}

    virtual bool selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const;
};

} // namespace BB

} // namespace CENSO

#endif // BRANCHERCONTINUOUSRANDOMVARIABLE_H
