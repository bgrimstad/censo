/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Utils/definitions.h"
#include "OptimizationProblem/constraint.h"

namespace CENSO
{

// Return statuses of optimizer
enum class SolverStatus
{
    OPTIMAL,        // Success
    FEASIBLE,
    INFEASIBLE,
    UNBOUNDED,      // Continuous unbounded
    LIMIT_EXCEEDED,
    ERROR
};

// Print solver status (inline to avoid multiple definitions)
inline std::ostream& operator<<(std::ostream& os, const SolverStatus status)
{
    switch (status)
    {
        case SolverStatus::OPTIMAL: os << "OPTIMAL"; break;
        case SolverStatus::FEASIBLE: os << "FEASIBLE"; break;
        case SolverStatus::INFEASIBLE: os << "INFEASIBLE"; break;
        case SolverStatus::UNBOUNDED: os << "UNBOUNDED"; break;
        case SolverStatus::LIMIT_EXCEEDED: os << "LIMIT_EXCEEDED"; break;
        case SolverStatus::ERROR: os << "ERROR"; break;
        default: assert(false); break; // Throw exception (unreachable)
    }
    return os;
}

/*
 * Struct for solver results
 */
struct SolverResult
{
    // Default results is error
    SolverResult()
        : status(SolverStatus::ERROR),
          objectiveValue(INF)
    {
    }

    SolverResult(SolverStatus status,
                 double objectiveValue,
                 std::vector<double> primalVariables)
        : status(status),
          objectiveValue(objectiveValue),
          primalVariables(primalVariables)
    {
    }

    // Printing solver result
    friend std::ostream& operator<<(std::ostream& os, const SolverResult result)
    {
        os << "Status:          " << result.status << std::endl;
        os << "Objective value: " << std::setprecision(10) << result.objectiveValue << std::endl;
        return os;
    }

    SolverStatus status;
    double objectiveValue; // Optimal objective function value
    std::vector<double> primalVariables; // Optimal point. Length n.
    std::vector<double> lowerBoundDualVariables; // Dual variables for variable lower bounds. Length n.
    std::vector<double> upperBoundDualVariables; // Dual variables for variable upper bounds. Length n.
    std::vector<double> constraintDualVariables; // Dual variables for constraints. Length m.
};

/*
 * Solver provides an interface to all solvers.
 * The Solver class and all its derived classes
 * should be in the CENSO namespace. However, each
 * solver class should have an interface object
 * which is kept outside of the CENSO namespace.
 * The reason for this is to eliminate any overlap
 * between CENSO and external libraries/codes.
 */
class Solver
{
public:

    Solver(ConstraintPtr constraints)
        : constraints(constraints),
          isInitialized(false),
          isDualAdequate(false)
    {
    }

    Solver(const Solver &copy) = delete;
    Solver& operator=(const Solver &assign) = delete;

    virtual ~Solver() {}

    SolverResult optimize()
    {
        if (!isInitialized)
        {
            initialize(); // This can fail - catch exceptions!
            isInitialized = true;
        }

        SolverResult result = runOptimizer();

        if (!verifySolverResult(result))
            result.status = SolverStatus::ERROR;

        //assert(result.status != SolverStatus::ERROR);

        return result;
    }

    bool dualAdequate() const { return isDualAdequate; }

protected:

    virtual void initialize() = 0;
    virtual SolverResult runOptimizer() = 0;

    bool verifySolverResult(const SolverResult &result) const
    {
        unsigned int n = constraints->getNumVariables();
        unsigned int m = constraints->getNumConstraints();
        if (n != result.primalVariables.size())
            return false;
//        if (n != result.lowerBoundDualVariables)
//            return false;
//        if (n != result.upperBoundDualVariables)
//            return false;
//        if (m != result.constraintDualVariables)
//            return false;

        return true;
    }

    void setInitialized(bool init) { isInitialized = init; }
    bool getInitialized() const { return isInitialized; }

    ConstraintPtr constraints; // Solver should not alter this object (shared pointer)
    bool isInitialized;
    bool isDualAdequate;

    // TODO: add flag to indicate if solver is dual-adequate (produces dual variables)

};

//typedef std::shared_ptr<Solver> OptimizerPtr;

} // namespace CENSO

#endif // OPTIMIZER_H
