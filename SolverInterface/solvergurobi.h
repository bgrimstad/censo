/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SOLVERGUROBI_H
#define SOLVERGUROBI_H

#include "solver.h"
#include "gurobi_c++.h"

namespace CENSO
{

class SolverGurobi : public Solver
{
public:
    SolverGurobi(ConstraintPtr constraints);

    SolverGurobi(const SolverGurobi &copy) = delete;
    SolverGurobi& operator=(const SolverGurobi &assign) = delete;
    ~SolverGurobi() override {}

protected:
    void initialize() override;
    SolverResult runOptimizer() override;
};

} // namespace CENSO

#endif // SOLVERGUROBI_H
