/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef POLYNOMIALPROBLEM04_H
#define POLYNOMIALPROBLEM04_H

#include "testproblem.h"

namespace CENSO
{

class POP04 : public TestProblem
{
public:
    POP04();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:
    std::vector<double> xopt_known; // The known optimal point of the problem
    std::vector<double> xopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // POLYNOMIALPROBLEM04_H
