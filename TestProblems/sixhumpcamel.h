/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef SIXHUMPCAMEL_H
#define SIXHUMPCAMEL_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Testing using Six hump camel back function (this function is NOT very nice to Ipopt):
 * Global optima at x* = (0.089842, -0.712656) and x* = (-0.089842, 0.712656), with f(x*) = -1.031628453
 */
class SixHumpCamel : public TestProblem
{
public:
    SixHumpCamel();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:
    static double sixHumpCamelFunction(DenseVector x); // Static allows for function pointers

    std::vector<double> zopt_known; // The known optimal point of the problem
    std::vector<double> zopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // SIXHUMPCAMEL_H
