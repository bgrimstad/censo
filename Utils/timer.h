/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef TIMER_H
#define TIMER_H

#include "chrono"

namespace CENSO
{

class Timer
{
public:
    Timer()
        : tStart(std::chrono::steady_clock::now()),
          tDur(std::chrono::duration<double>(0)),
          running(false)
    {
    }

    void start()
    {
        if (!running)
        {
            tStart = std::chrono::steady_clock::now();
            running = true;
        }
    }

    void stop()
    {
        if (running)
        {
            auto tNow = std::chrono::steady_clock::now();
            auto dur = (tNow-tStart);
            tDur += dur;
            running = false;
        }
    }

    void reset()
    {
        tDur = std::chrono::duration<double>(0);
        running = false;
    }

    unsigned int getSeconds()
    {
        return std::chrono::duration_cast<std::chrono::seconds> (getDuration()).count();
    }

    unsigned int getMilliSeconds()
    {
        return std::chrono::duration_cast<std::chrono::milliseconds> (getDuration()).count();
    }

    unsigned int getMicroSeconds()
    {
        return std::chrono::duration_cast<std::chrono::microseconds> (getDuration()).count();
    }

private:
    std::chrono::steady_clock::time_point tStart;
    std::chrono::duration<double> tDur;

    bool running;

    std::chrono::duration<double> getDuration()
    {
        if (running)
        {
            auto tNow = std::chrono::steady_clock::now();
            auto dur = (tNow-tStart);
            return (tDur + dur);
        }
        else
        {
            return tDur;
        }
    }

};

} // namespace CENSO

#endif // TIMER_H
