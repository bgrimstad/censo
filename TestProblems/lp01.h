/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef LP01_H
#define LP01_H

#include "testproblem.h"

namespace CENSO
{

/*
 * Chair maker problem from GAMS
 * Optimal solution: 10,417.29
 */
class LP01 : public TestProblem
{
public:
    LP01();

protected:
    virtual void runProblem();

    virtual bool validateResult();

private:

    double fopt_known; // The objective function value at the known optimal point
    double fopt_found; // The objective function value at the found optimal point
};

} // namespace CENSO


#endif // LP01_H
