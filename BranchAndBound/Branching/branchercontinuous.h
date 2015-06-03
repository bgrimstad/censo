/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHERCONTINUOUS_H
#define BRANCHERCONTINUOUS_H

#include "brancher.h"

namespace CENSO
{

namespace BB
{

// Branching rules for continuous variables
// These rules are only used in global optimization of non-convex problems
// NOTE: could add rule that uses circular buffer (mostly for testing)
enum BranchingRuleContinuous
{
    MOST_PROMISING, // Default. Selects branching variable based on branching costs.
    LONGEST_INTERVAL, // Selects branching variable with bound intervals (argmax(ub-lb)). Note: dependent on variable scaling!
    RANDOM_VARIABLE // Selects a random branching variable. Note: gives stochastic run time!
};

class BrancherContinuous : public Brancher
{
public:
    BrancherContinuous();
    virtual ~BrancherContinuous() {}

    virtual NodeList branchOnVariable(const NodePtr node, VariablePtr branchingVariable) const;

protected:
    virtual bool validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const;

    double branchingThreshold; // Branch until variable bounds meet: ub-lb <= branchingThreshold
};

} // namespace BB

} // namespace CENSO

#endif // BRANCHERCONTINUOUS_H
