/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BILINEARRELAXATIONTEST_H
#define BILINEARRELAXATIONTEST_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Test of convex relaxation of the bilinear constraint: xy = z
 * The problem is simply to minimize z in the domain x = [-1, 1], y = [-1, 2].
 * Two problems are created:
 * 1) non-convex problem with the bilinear constraint
 *    the starting point of (1,-0.5,-0,5) should trick the optimizer
 *    to find the local optimum (1,-1,-1)
 * 2) convex problem with relaxation of the bilinear constraint
 *    the optimizer should find the global optimum at (-1,2,-2)
 */
class BilinearRelaxationTest : public TestProblem
{
public:
    BilinearRelaxationTest();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:
    static DenseVector bilinearFunction(DenseVector x); // Static allows for function pointers

    std::vector<double> zopt_known; // The known optimal point of the problem
    std::vector<double> zopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // BILINEARRELAXATIONTEST_H
