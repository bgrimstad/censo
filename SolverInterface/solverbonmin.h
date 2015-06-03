/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef OPTIMIZERBONMIN_H
#define OPTIMIZERBONMIN_H

#include "solver.h"
#include "interfacebonmin.h"
#include "BonCbc.hpp"
#include "BonBonminSetup.hpp"

namespace CENSO
{

class SolverBonmin : public Solver
{
public:
    SolverBonmin(ConstraintPtr constraints)
        : Solver(constraints)
    {}

    SolverBonmin(const SolverBonmin &copy) = delete;
    SolverBonmin& operator=(const SolverBonmin &assign) = delete;
    ~SolverBonmin() override {}

protected:
    void initialize() override;
    CENSO::SolverResult runOptimizer() override;

private:
    SmartPtr<InterfaceBonmin> bonminInterface;
    Bonmin::BonminSetup bonmin;
};

} // namespace CENSO

#endif // OPTIMIZERBONMIN_H
