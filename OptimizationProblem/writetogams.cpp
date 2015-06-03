/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraint.h"
#include <fstream>

using std::endl;

namespace CENSO
{

/*
 * The following code is inspired by the
 * writeGAMS.cpp file in the Couenne project
 * (https://projects.coin-or.org/Couenne).
 */
void Constraint::writeToGAMS(const std::string &fname) const
{
    std::ofstream f(fname.c_str());

    f << std::setprecision(10);

    const int nline = 20;

    // Header
    f << "* Problem written by CENSO (https://github.com/bgrimstad/censo)" << endl
      << "* Problem name: " << constraintName << endl
      << "* " << getNumVariables() << " variables, " << getNumConstraints() << " constraints" << endl << endl;

    // Give names to variables
    for (unsigned int i = 0; i < getNumVariables(); ++i)
    {
        std::string name = "x" + std::to_string(i);
        auto var = getVariableAt(i);
        var->setName(name);
    }

    // Variables
    f << "variables ";

    for (unsigned int i = 0, j = 0; i < getNumVariables(); ++i)
    {
        if (j != 0) f << ",";
        if (++j % nline == 0) f << endl << "          ";
        auto var = getVariableAt(i);
        f << var->getName();
    }

    f << ",objvar;" << endl << endl;

    bool have_nonnegative = false, have_binary = false, have_integer = false;

    for (auto var : getVariables())
    {
        // Check if there are nonnegative variables
        if (var->getLowerBound() >= 0
                && assertNear(var->getLowerBound(), 0.0)
                && var->getType() == VariableType::CONTINUOUS)
        {
            have_nonnegative = true;
        }

        // Check if there are binary variables
        if (var->getType() == VariableType::BINARY)
        {
            have_binary = true;
        }

        // Check if there are integer variables
        if (var->getType() == VariableType::INTEGER)
        {
            have_integer = true;
        }
    }

    // Write nonnegative variables
    if (have_nonnegative)
    {
        f << "positive variables ";

        for (unsigned int i = 0, j = 0; i < getNumVariables(); ++i)
        {
            auto var = getVariableAt(i);

            // Only write continuous variables with lower bound equal to zero
            if (assertNear(var->getLowerBound(), 0.0)
                    && var->getType() == VariableType::CONTINUOUS)
            {
                if (j != 0) f << ",";
                if (++j % nline == 0) f << endl << "                   ";
                auto var = getVariableAt(i);
                f << var->getName();
            }
        }

        f << ";" << endl << endl;
    }

    // Write binary variables
    if (have_binary)
    {
        f << "binary variables ";

        for (unsigned int i = 0, j = 0; i < getNumVariables(); ++i)
        {
            auto var = getVariableAt(i);

            // Only write continuous variables with lower bound equal to zero
            if (var->getType() == VariableType::BINARY)
            {
                if (j != 0) f << ",";
                if (++j % nline == 0) f << endl << "                 ";
                auto var = getVariableAt(i);
                f << var->getName();
            }
        }

        f << ";" << endl << endl;
    }

    // Write integer variables
    if (have_integer)
    {
        f << "integer variables ";

        for (unsigned int i = 0, j = 0; i < getNumVariables(); ++i)
        {
            auto var = getVariableAt(i);

            // Only write continuous variables with lower bound equal to zero
            if (var->getType() == VariableType::INTEGER)
            {
                if (j != 0) f << ",";
                if (++j % nline == 0) f << endl << "                  ";
                auto var = getVariableAt(i);
                f << var->getName();
            }
        }

        f << ";" << endl << endl;
    }

    // Equations
    f << "equations ";

    for (unsigned int i = 0; i < getNumConstraints(); ++i)
    {
        if (i != 0) f << ",";
        if ((i+1) % nline == 0) f << std::endl << "          ";
        f << "e" << std::to_string(i);
    }

    if (getNumConstraints() > 0) f << ",";
    f << "objdef;" << endl << endl;

    writeConstraintEquationsToGAMS(f, 0);

    // Objective
    f << "objdef.. - objvar";
    for (auto var : getVariables())
    {
        if (!assertNear(var->getCost(), 0.0))
        {
            double cost = var->getCost();
            if (cost > 0) f << " + " << std::to_string(cost) << "*";
            else f << " - " << std::to_string(std::abs(cost)) << "*";
            f << var->getName();
        }
    }
    f << " =E= 0;" << endl << endl;

    // Initial levels and variable bounds
    for (auto var : getVariables())
    {
        f << var->getName();
        f << ".l = " << var->getValue() << ";" << endl;

        if (var->getLowerBound() > -INF)
        {
            f << var->getName();
            f << ".lo = " << var->getLowerBound() << ";" << endl;
        }

        if (var->getUpperBound() < INF)
        {
            f << var->getName();
            f << ".up = " << var->getUpperBound() << ";" << endl;
        }
    }

    f << endl << endl;

    // Model
    f << "Model m / all /;" << endl
      << "m.limrow=0; m.limcol=0;" << endl
      << "Solve m using ";

    auto pClass = assessClass();

    if (pClass == ProblemClass::MINLP)
        f << "MINLP";
    else if (pClass == ProblemClass::MILP)
        f << "MIP";
    else if (pClass == ProblemClass::NLP)
        f << "NLP";
    else if (pClass == ProblemClass::LP)
        f << "LP";
    else
        f << "MINLP";

    f << " minimizing objvar;" << endl;

    f.close ();
}

} // namespace CENSO
