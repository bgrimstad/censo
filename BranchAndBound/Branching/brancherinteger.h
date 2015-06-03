/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHERINTEGER_H
#define BRANCHERINTEGER_H

#include "brancher.h"

namespace CENSO
{

namespace BB
{

// Branching rules for integer variables
enum BranchingRuleInteger
{
    MOST_FRACTIONAL, // Branch on the most infeasible integer/binary variable (default)
    PSEUDOCOST, // Not implemented.
    STRONG_BRANCHING, // Not implemented. Strong braching is a very popular approach.
    USER_SPECIFIED // Not implemented, order in branching_variables determines branching priority
};

class BrancherInteger : public Brancher
{
public:
    BrancherInteger() : Brancher() {}
    virtual ~BrancherInteger() {}

    virtual NodeList branchOnVariable(const NodePtr node, VariablePtr branchingVariable) const;

protected:
    virtual bool validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const;
};

} // namespace BB

} // namespace CENSO

#endif // BRANCHERINTEGER_H
