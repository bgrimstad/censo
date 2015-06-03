/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef IPOPTOPTIMIZER_H
#define IPOPTOPTIMIZER_H

#include "solver.h"
#include "interfaceipopt.h"
#include "IpIpoptApplication.hpp"
//#include "IpSolveStatistics.hpp"

using Ipopt::SmartPtr;
using Ipopt::IpoptApplication;

namespace CENSO
{

class SolverIpopt : public Solver
{
public:

    /*
     * Default constructor
     */
    SolverIpopt(ConstraintPtr constraints);

    SolverIpopt(const SolverIpopt &copy) = delete;
    SolverIpopt& operator=(const SolverIpopt &assign) = delete;
    ~SolverIpopt() override {}

    void setBoundRelaxFactor(double relaxFactor)
    {
        if (relaxFactor < 0) relaxFactor = 0;
        boundRelaxFactor = relaxFactor;
    }

protected:
    void initialize() override;
    SolverResult runOptimizer() override;

private:
    SmartPtr<InterfaceIpopt> ipoptInterface;
    SmartPtr<IpoptApplication> ipoptApplication;

    double boundRelaxFactor;
};

} // namespace CENSO

#endif // IPOPTOPTIMIZER_H
