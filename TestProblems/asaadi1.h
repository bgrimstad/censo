/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef ASAADI1_H
#define ASAADI1_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Asaadi1 is a convex MINLP problem from Leyffer's PhD thesis
 *
 * The problem is solved three times:
 * 1) All integer variables relaxed (Ipopt)
 * 2) One integer variable relaxed (Bonmin and BnB)
 * 3) Only integer variables (Bonmin and BnB)
 */
class Asaadi1 : public TestProblem
{
public:
    Asaadi1();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:

    std::vector<double> zopt_known; // The known optimal point of the problem
    std::vector<double> zopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // ASAADI1_H
