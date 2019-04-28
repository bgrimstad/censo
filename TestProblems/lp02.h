/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef LP02_H
#define LP02_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Transportation problem from GAMS
 * Optimal solution: 7,425
 */
class LP02 : public TestProblem
{
public:
    LP02();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO

#endif // LP02_H
