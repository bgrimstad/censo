/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHERADMIN_H
#define BRANCHERADMIN_H

#include "brancherinteger.h"
#include "branchercontinuous.h"

#include "../nodelist.h"

namespace CENSO
{

namespace BB
{

// Branching priorities (not in use)
enum class BranchingPriority
{
    BINARY_INTEGER_CONTINUOUS, // Branching in that order. Currently the only one that is implemented (default)
    RANDOM_MIX // Basically removes prioritization
};

// Add option for different division schemes
// Division into disjunct sub-problems: default is dichotomy (2), but polytomic (n) could be implemented

class BrancherAdmin
{
public:
    BrancherAdmin();
    BrancherAdmin(BranchingRuleInteger branchingRuleInteger,
                  BranchingRuleContinuous branchingRuleContinuous);
    //BrancherAdmin(BrancherAdmin const& copy) = delete;
    //BrancherAdmin& operator = (BrancherAdmin const& assign) = delete;

    ~BrancherAdmin() {}

    NodeList branch(const NodePtr node) const;

private:

    BranchingRuleInteger branchingRuleInteger;
    BranchingRuleContinuous branchingRuleContinuous;

    BrancherPtr brancherInteger, brancherContinuous;

    // Factory method for integer branchers
    Brancher* createIntegerBrancher(BranchingRuleInteger branchingRule);

    // Factory method for continuous branchers
    Brancher* createContinuousBrancher(BranchingRuleContinuous branchingRule);

};

} // namespace BB

} // namespace CENSO

#endif // BRANCHERADMIN_H
