/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef MICHALEWICZ_H
#define MICHALEWICZ_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Testing using Michalewics (this function is very non-linear):
 * Global optimum f(x*) = -1.8013, at x* = (2.2029, 1.5708)
 */
class Michalewicz : public TestProblem
{
public:
    Michalewicz();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:
    static double michalewiczFunction(DenseVector x, unsigned int m); // Static allows for function pointers

    std::vector<double> zopt_known; // The known optimal point of the problem
    std::vector<double> zopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // MICHALEWICZ_H
