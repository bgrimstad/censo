/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "Utils/definitions.h"
#include "TestProblems/alltestproblems.h"

int main(int argc, char *argv[])
{
    CENSO::Timer timer;
    timer.start();

    CENSO::POP01 p01; p01.run();
    CENSO::POP02 p02; p02.run();
    CENSO::POP03 p03; p03.run();
    CENSO::POP04 p04; p04.run();
    CENSO::POP05 p05; p05.run();
    CENSO::POP06 p06; p06.run();
    CENSO::POP07 p07; p07.run();
    CENSO::POP08 p08; p08.run();
    CENSO::POP09 p09; p09.run();
    CENSO::POP10 p10; p10.run();
    CENSO::POP11 p11; p11.run();
    CENSO::POP12 p12; p12.run();
    CENSO::POP13 p13; p13.run();

    timer.stop();
    std::cout << "Total time: " << timer.getMilliSeconds() << " sec" << std::endl;

    return 0;
}


