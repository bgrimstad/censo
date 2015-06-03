/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "brancher.h"

namespace CENSO
{

namespace BB
{

bool Brancher::validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const
{
    // Check if variable is a branching variable
    auto branchingVariables = node->getBranchingVariables();
    auto it = std::find(branchingVariables.begin(), branchingVariables.end(), branchingVariable);
    if (it == branchingVariables.end())
    {
        std::cout << "Brancher: Invalid branching variable! Variable is not in branching variables list!" << std::endl;
        return false;
    }

    // Looks OK
    return true;
}

} // namespace BB

} // namespace CENSO
