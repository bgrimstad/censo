/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef LINEARINTERVALANALYSISTEST_H
#define LINEARINTERVALANALYSISTEST_H

#include "testproblem.h"

namespace CENSO
{

class LinearIntervalAnalysisTest : public TestProblem
{
public:
    LinearIntervalAnalysisTest()
        : testResult(false)
    {
        testName = "Linear Interval Analysis";
    }

protected:
    virtual void runProblem();

    virtual bool validateResult()
    {
        return testResult;
    }

private:
    bool testResult;
};

} // namespace CENSO

#endif // LINEARINTERVALANALYSISTEST_H
