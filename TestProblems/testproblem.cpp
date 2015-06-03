/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "testproblem.h"
#include <iostream>

using std::cout;
using std::endl;

namespace CENSO
{

TestProblem::TestProblem()
    : testName(" "),
      timer(Timer())
{
}

void TestProblem::run()
{
    cout << "Running " << testName << " problem..." << endl;

//    This code can be used to mute solveProblem() console output
//    std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
//    std::ofstream   fout("/dev/null");
//    std::cout.rdbuf(fout.rdbuf()); // redirect 'cout' to a 'fout'
//    // solveProblem()...
//    std::cout.rdbuf(cout_sbuf); // restore the original stream buffer

    timer.start();
    runProblem();
    timer.stop();

    bool success = validateResult();

    if (success)
    {
        // Print info
        cout << testName << " problem successfully solved in " << timer.getMilliSeconds() << " (ms)" << endl;
    }
    else
    {
        cout << "Could not solve " << testName << " problem!" << endl;
    }
    cout << endl;

    assert(success);
}

} // namespace CENSO
