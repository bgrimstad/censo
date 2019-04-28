/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef TESTPROBLEM_H
#define TESTPROBLEM_H

#include <string>
#include "Utils/timer.h"
#include "SolverInterface/solveripopt.h"

namespace CENSO
{

/*
 * Class for general test problems (need not be an optimization problems)
 * - runs test problems
 * - times run
 * - validates result
 * - displays summary
 */
class TestProblem
{
public:
    TestProblem();
    TestProblem(TestProblem const& copy) = delete;
    TestProblem& operator = (TestProblem const& assign) = delete;

    virtual ~TestProblem() {}

    void run();

protected:

    std::string testName;

    // Timer
    Timer timer;

    // Solving problem
    virtual void runProblem() = 0;

    virtual bool validateResult() = 0;
};

} // namespace CENSO

#endif // TESTPROBLEM_H
