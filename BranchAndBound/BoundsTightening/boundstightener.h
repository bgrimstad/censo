/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef BOUNDSTIGHTENER_H
#define BOUNDSTIGHTENER_H

#include "OptimizationProblem/constraint.h"
#include "Utils/metrics.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

/*
 * Interface for bounds tightening algorithms
 */
class BoundsTightener
{
public:
    BoundsTightener()
        : BoundsTightener(0.1, 1)
    {}

    BoundsTightener(double threshold, unsigned int maxIterations)
        : threshold(threshold),
          maxIterations(maxIterations)
    {}

    virtual ~BoundsTightener() {}

    bool run(ConstraintPtr constraints)
    {
        // Retrieve old variable bounds
        std::vector<double> olb, oub, nlb, nub;

        auto vars = constraints->getVariables();

        for (auto var : vars)
        {
            olb.push_back(var->getLowerBound());
            oub.push_back(var->getUpperBound());
        }

        // Run sequential tightening
        bool success = true;
        double reduction = 0; // 0-100 %
        unsigned int iter = 0;

        do
        {
            success = doTightening(constraints);

            reduction = 0;

            // Retrieve new bounds
            if (success)
            {
                nlb.clear();
                nub.clear();

                for (auto var : vars)
                {
                    nlb.push_back(var->getLowerBound());
                    nub.push_back(var->getUpperBound());
                }

                // Compute reduction
                reduction = 100*relativeChangeIntervals(olb, oub, nlb, nub);

                // Save bounds
                olb = nlb;
                oub = nub;

//                cout << "BT reduction: " << reduction << "%" << endl;
            }

            assert(reduction >= 0 && reduction <= 100);
            ++iter;
        }
        while (success && reduction > threshold && iter < maxIterations);

        return success;
    }

    double getThreshold() const
    {
        return threshold;
    }

    unsigned int getMaxIterations() const
    {
        return maxIterations;
    }

    void setThreshold(double threshold)
    {
        assert(threshold >= 0 && threshold <= 100);
        this->threshold = threshold;
    }

    void setMaxIterations(unsigned int maxIterations)
    {
        this->maxIterations = maxIterations;
    }

private:
    double threshold; // Stop when improvement is < threshold (0,1)
    unsigned int maxIterations; // Stop after maxIterations

    virtual bool doTightening(ConstraintPtr constraints) = 0;
};

} // namespace BB

} // namespace CENSO

#endif // BOUNDSTIGHTENER_H
