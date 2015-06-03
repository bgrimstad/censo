/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "bbutils.h"
#include "OptimizationProblem/constraintlinear.h"

namespace CENSO
{

namespace BB
{

bool isIntegerProblem(ConstraintPtr constraints)
{

    auto vars = constraints->getVariables();

    for (unsigned int i = 0; i < vars.size(); i++)
    {
        if (vars.at(i)->getType() == VariableType::BINARY
            || vars.at(i)->getType() == VariableType::INTEGER)
        {
            if (vars.at(i)->getLowerBound() != vars.at(i)->getUpperBound())
                return true;
        }
    }

    return false;
}

bool isSolutionIntegerFeasible(std::vector<VariablePtr> variables, std::vector<double> z)
{
    assert(z.size() == variables.size());

    for (unsigned int i = 0; i < z.size(); i++)
    {
        if (variables.at(i)->getType() == VariableType::BINARY
            || variables.at(i)->getType() == VariableType::INTEGER)
        {
            if (!isInteger(z.at(i)))
                return false;
        }
    }
    return true;
}

std::vector<VariablePtr> findBranchingVariables(ConstraintPtr constraints)
{
    std::vector<VariablePtr> bvars = constraints->getComplicatingVariables();

    auto vars = constraints->getVariables();
    for (auto var : vars)
    {
        if (var->getType() == VariableType::BINARY
            || var->getType() == VariableType::INTEGER)
        {
            // Add to bv if not already there
            if (std::find(bvars.begin(), bvars.end(), var) == bvars.end())
            {
                bvars.push_back(var);
            }
        }
    }

    return bvars;
}

ConstraintPtr objectiveCut(ConstraintPtr constraints, double objectiveBound, bool upperCut)
{
    auto allVars = constraints->getVariables();

    std::vector<VariablePtr> vars;

    for (auto var : allVars)
    {
        if (var->getCost() != 0)
            vars.push_back(var);
    }

    DenseMatrix A = DenseMatrix::Zero(1, vars.size());

    int i = 0;
    for (auto var : vars)
    {
        if (upperCut)
            A(0,i++) = var->getCost();
        else
            A(0,i++) = -var->getCost();
    }

    DenseVector b(1);
    if (upperCut)
        b(0) = objectiveBound;
    else
        b(0) = -objectiveBound;

    return std::make_shared<ConstraintLinear>(vars, A, b, false);
}

} // namespace BB

} // namespace CENSO
