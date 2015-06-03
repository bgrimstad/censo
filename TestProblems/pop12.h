/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef POLYNOMIALPROBLEM12_H
#define POLYNOMIALPROBLEM12_H

#include "testproblem.h"

namespace CENSO
{

class POP12 : public TestProblem
{
public:
    POP12();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:
    std::vector<double> xopt_known; // The known optimal point of the problem
    std::vector<double> xopt_found; // The optimal point found by the optimizer

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point

    void P12_1();
    void P12_2(); // Used in JOGO paper
    void P12_3();
};

} // namespace CENSO

#endif // POLYNOMIALPROBLEM12_H
