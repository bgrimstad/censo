/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "linearintervalanalysistest.h"
#include "OptimizationProblem/constraintlinear.h"

using std::cout;
using std::endl;

namespace CENSO
{

void LinearIntervalAnalysisTest::runProblem()
{
    // NOTE: should test with INF bounds!
    std::vector<VariablePtr> vars = {
        std::make_shared<Variable>(1, 0, 1),
        std::make_shared<Variable>(1, 0, 2),
        std::make_shared<Variable>(1, 0, 3)
    };

    {
        DenseMatrix A(2,3);
        A << 1, -0.1, 0,
             0, 0, 1;
        DenseVector b(2); b << 0.5, 0.1;

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);
        testResult = lincon->reduceVariableRanges();
    }

    for (auto var : vars)
        cout << *var << endl;
}

} // namespace CENSO
