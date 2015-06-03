/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BRANCHER_H
#define BRANCHER_H

#include "../node.h"
#include "../nodelist.h"

namespace CENSO
{

namespace BB
{

// Virtual class for branchers
// A brancher is state-less and all information should be stored in the node object
class Brancher
{
public:
    Brancher() = default;
    virtual ~Brancher() {}

    // Selects next branching variable for node. Returns true if a branching variable was found.
    virtual bool selectBranchingVariable(const NodePtr node, VariablePtr &branchingVariable) const = 0;

    // Branches on variable and return a list of children nodes.
    virtual NodeList branchOnVariable(const NodePtr node, VariablePtr branchingVariable) const = 0;

protected:

    virtual bool validateBranchingVariable(const NodePtr node, VariablePtr branchingVariable) const;

};

typedef std::shared_ptr<Brancher> BrancherPtr;

} // namespace BB

} // namespace CENSO

#endif // BRANCHER_H
