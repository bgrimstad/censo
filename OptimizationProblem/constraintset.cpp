/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraintset.h"

using std::cout;
using std::endl;

namespace CENSO
{

/*
 * By default, and opposed to a single constraint, constraint composites are assumed to:
 * - have an analytical gradient and hessian
 * - be linear
 * - be convex
 * - have a convex relaxation available
 *
 * These properties must be checked and updated when
 * adding new constraints to the composite.
 */
ConstraintSet::ConstraintSet()
{
    jacobianCalculated = true;
    hessianCalculated = true;
    constraintLinear = true;
    constraintConvex = true;
    convexRelaxationAvailable = true;

    checkConstraintSanity();
}

ConstraintSet::ConstraintSet(const ConstraintSet &copy, bool deep)
    : Constraint(copy, deep)
{
    for (unsigned int i = 0; i < copy.constraints.size(); i++)
    {
        constraints.push_back(copy.constraints.at(i)->clone(deep));
    }

    if (deep)
    {
        // Update variables
        for (unsigned int i = 0; i < copy.constraints.size(); i++)
        {
            std::vector<int> conVars = copy.getConstraintVariableIndices(i);

            std::vector<VariablePtr> newConVars;
            for (unsigned int j = 0; j < conVars.size(); j++)
                newConVars.push_back(variables.at(conVars.at(j)));

            constraints.at(i)->setVariables(newConVars);
        }
    }
}

DenseVector ConstraintSet::eval(const DenseVector &x) const
{
    setVariableValues(x);

    DenseVector y = DenseVector::Zero(numConstraints);

    int k = 0;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        DenseVector yi = constraints.at(i)->eval();

        for (unsigned int j = 0; j < yi.size(); j++)
        {
            y(k++) = yi(j);
        }
    }

    return y;
}

DenseVector ConstraintSet::evalJacobian(const DenseVector &x) const
{
    setVariableValues(x);

    DenseVector dx = DenseVector::Zero(nnzJacobian);

    int k = 0;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        DenseVector di = constraints.at(i)->evalJacobian();

        for (unsigned int j = 0; j < di.size(); j++)
        {
            dx(k++) = di(j);
        }
    }

    return dx;
}

DenseVector ConstraintSet::evalHessian(const DenseVector &x) const
{
    DenseVector ddx = DenseVector::Zero(nnzHessian);

    int k = 0;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        DenseVector ddi = constraints.at(i)->evalHessian();

        for (unsigned int j = 0; j < ddi.size(); j++)
        {
            ddx(k++) = ddi(j);
        }
    }

    return ddx;
}

void ConstraintSet::structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol)
{

    unsigned int globalConstraintIndex = 0;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        std::vector<int> iRowI;
        std::vector<int> jColI;

        constraints.at(i)->structureJacobian(iRowI, jColI);
        assert(iRowI.size() == jColI.size());

        std::vector<int> variableIndices = getConstraintVariableIndices(i);

        for (unsigned int j = 0; j < jColI.size(); j++)
        {
            iRow.push_back(globalConstraintIndex + iRowI.at(j));
            jCol.push_back(variableIndices.at(jColI.at(j)));
        }

        globalConstraintIndex += constraints.at(i)->getNumConstraints();
    }

    assert(iRow.size() == nnzJacobian);
    assert(iRow.size() == jCol.size());
    assert(globalConstraintIndex == numConstraints);
}

void ConstraintSet::structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol)
{
    unsigned int globalConstraintIndex = 0;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        std::vector<int> eqNrI;
        std::vector<int> iRowI;
        std::vector<int> jColI;

        constraints.at(i)->structureHessian(eqNrI,iRowI,jColI);
        assert(eqNrI.size() == iRowI.size());
        assert(iRowI.size() == jColI.size());

        std::vector<int> variableIndices = getConstraintVariableIndices(i);

        for (unsigned int j = 0; j < jColI.size(); j++)
        {
            eqnr.push_back(globalConstraintIndex + eqNrI.at(j));
            iRow.push_back(variableIndices.at(iRowI.at(j)));
            jCol.push_back(variableIndices.at(jColI.at(j)));
        }

        globalConstraintIndex += constraints.at(i)->getNumConstraints();
    }

    assert(eqnr.size() == iRow.size());
    assert(jCol.size() == iRow.size());
    assert(nnzHessian == iRow.size());
    assert(globalConstraintIndex == numConstraints);
}

bool ConstraintSet::checkFeasibility(const DenseVector &x, double tol) const
{
    setVariableValues(x);

    if (tol < 0) tol = 0;

    for (const auto &con : constraints)
    {
        if (!con->checkFeasibility(tol))
            return false;
    }
    return true;
}

void ConstraintSet::localRefinement(const DenseVector &x)
{
    setVariableValues(x);

    for (const auto &con : constraints)
        con->localRefinement();
}

void ConstraintSet::add(ConstraintPtr constraint)
{
    // Add new variables
    addVariables(constraint->getVariables());

    // Copy the constraint (no copy)
    ConstraintPtr cnew = constraint;

    // Store constraint object
    constraints.push_back(cnew);
    numConstraints += cnew->getNumConstraints();
    nnzJacobian += cnew->getNumNonZerosJacobian();
    nnzHessian += cnew->getNumNonZerosHessian();

    // Set constraint bounds
    std::vector<double> lbF,ubF;
    cnew->getConstraintBounds(lbF, ubF);
    for (unsigned int i = 0; i < lbF.size(); i++)
    {
        lb.push_back(lbF.at(i));
        ub.push_back(ubF.at(i));
    }

    // Check available calculations
    // Different from the constraint, the composite assumes that everything is true!
    if (!cnew->isJacobianCalculated())
    {
        jacobianCalculated = false;
    }

    if (!cnew->isHessianCalculated())
    {
        hessianCalculated = false;
    }

    if (!cnew->isConstraintLinear())
    {
        constraintLinear = false;
    }

    if (!cnew->isConstraintConvex())
    {
        constraintConvex = false;
    }

    if (!cnew->isConstraintLinear() && !cnew->isConstraintConvex() && !cnew->hasConvexRelaxation())
    {
        convexRelaxationAvailable = false;
    }

    // Check constraint
    checkConstraintSanity();
}

/*
 * Returns a set of constraints in the original variables and
 * auxiliary variables (shallow copy). The ordering of the
 * original variables will change when aux variables are added.
 */
std::shared_ptr<Constraint> ConstraintSet::getConvexRelaxation()
{
    // Create new constraint composite for the convex relaxation
    ConstraintSetPtr relConstraintSet = std::make_shared<ConstraintSet>();

    // Fill the new constraint composite
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        // Create convex relaxation of all non-convex constraints, and clone convex constraints
        ConstraintPtr relConstraint = constraints.at(i)->getConvexRelaxation();
        relConstraintSet->add(relConstraint);
    }

    assert(relConstraintSet->isConstraintConvex());

    return relConstraintSet;
}

std::vector<VariablePtr> ConstraintSet::getComplicatingVariables() const
{
    std::vector<VariablePtr> vars;

    if (constraintLinear || constraintConvex)
        return vars;

    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        std::vector<VariablePtr> varsi = constraints.at(i)->getComplicatingVariables();

        for (unsigned int j = 0; j < varsi.size(); j++)
        {
            // Check if variable has been added, if not, add it
            auto var = varsi.at(j);
            if (std::find(vars.begin(), vars.end(), var) == vars.end())
            {
                vars.push_back(var);
            }
        }
    }

    // Needed?
    //std::sort(vars.begin(), vars.end());

    return vars;
}

bool ConstraintSet::reduceVariableRanges() const
{
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        if (!constraints.at(i)->reduceVariableRanges())
            return false;
    }

    return true;
}

void ConstraintSet::addVariables(const std::vector<VariablePtr> &variables)
{
    for (unsigned int i = 0; i < variables.size(); i++)
    {
        if (!hasVariable(variables.at(i)))
            this->variables.push_back(variables.at(i));
    }
}

bool ConstraintSet::isVariableSubset(const std::vector<VariablePtr> &variables) const
{
    if (this->variables.size() < variables.size())
        return false;

    for (unsigned int i = 0; i < variables.size(); i++)
    {
        if (!hasVariable(variables.at(i)))
            return false;
    }

    return true;
}

bool ConstraintSet::hasVariable(const VariablePtr &variable) const
{
    auto res = std::find(variables.begin(), variables.end(), variable);
    if (res == variables.end())
        return false;
    return true;
}

/*
 * Return the global indices of local constraint variables
 */
std::vector<int> ConstraintSet::getConstraintVariableIndices(int constraintIndex) const
{
    std::vector<int> variableIndices;

    for (const auto &var : constraints.at(constraintIndex)->getVariables())
    {
        auto it = std::find(variables.begin(), variables.end(), var);
        assert(it != variables.end()); // Exception
        variableIndices.push_back(it - variables.begin());
    }

    assert(variableIndices.size() == constraints.at(constraintIndex)->getNumVariables());

    return variableIndices;
}

std::ostream& ConstraintSet::output(std::ostream &os) const
{
    os << "Constraint Composite" << endl;
    Constraint::output(os);

    os << "Contains the following " << constraints.size() << " constraint objects:"<< endl;
    for (unsigned int i = 0; i < constraints.size(); i++)
    {
        os << "Constraint object nr " << i << ":" << endl;

        //std::vector<int> xi = constraintsDomainIndex.at(i);
        //os << "x: \t\t";
        //printVector(xi);

        constraints.at(i)->output(os);

        os << endl;
    }

    return os;
}

void ConstraintSet::writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const
{
    unsigned int counter = start;

    for (auto con : constraints)
    {
        con->writeConstraintEquationsToGAMS(os, counter);
        counter += con->getNumConstraints();
    }
}

} // namespace CENSO
