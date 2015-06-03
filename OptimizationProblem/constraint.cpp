/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "constraint.h"
#include <cassert>

using std::cout;
using std::endl;

namespace CENSO
{

Constraint::Constraint(const Constraint &copy, bool deep)
    : variables(copy.variables),
      lb(copy.lb),
      ub(copy.ub),
      constraintTypes(copy.constraintTypes),
      numConstraints(copy.numConstraints),
      nnzJacobian(copy.nnzJacobian),
      nnzHessian(copy.nnzHessian),
      jacobianCalculated(copy.jacobianCalculated),
      hessianCalculated(copy.hessianCalculated),
      constraintLinear(copy.constraintLinear),
      constraintConvex(copy.constraintConvex),
      convexRelaxationAvailable(copy.convexRelaxationAvailable),
      constraintName(copy.constraintName)
{
    if (deep)
    {
        copyVariables();
    }
}

void Constraint::copyVariables()
{
    std::vector<VariablePtr> vars;
    for (unsigned int i = 0; i < variables.size(); i++)
    {
        auto var = std::make_shared<Variable>(*variables.at(i));
        vars.push_back(var);
    }

    variables = vars;
}

DenseVector Constraint::getVariableValues() const
{
    DenseVector x = DenseVector::Zero(variables.size());
    int i = 0;
    for (const auto &var : variables)
        x(i++) = var->getValue();
    return x;
}

void Constraint::setVariableValues(const DenseVector &x) const
{
    assert(x.size() == (int)variables.size());
    int i = 0;
    for (const auto &var : variables)
        var->setValue(x(i++));
}

DenseVector Constraint::eval() const
{
    DenseVector x = getVariableValues();
    return eval(x);
}

DenseVector Constraint::evalJacobian() const
{
    DenseVector x = getVariableValues();
    return evalJacobian(x);
}

DenseVector Constraint::evalHessian() const
{
    DenseVector x = getVariableValues();
    return evalHessian(x);
}

std::vector<ConstraintType> Constraint::getConstraintTypes()
{
    if (constraintTypes.size() == 0)
    {
        if (constraintLinear)
        {
            for (unsigned int i = 0; i < numConstraints; i++ )
                constraintTypes.push_back(ConstraintType::LINEAR);
        }
        else
        {
            if (constraintConvex)
            {
                for (unsigned int i = 0; i < numConstraints; i++ )
                    constraintTypes.push_back(ConstraintType::NONLINEAR_CONVEX);
            }
            else
            {
                for (unsigned int i = 0; i < numConstraints; i++ )
                    constraintTypes.push_back(ConstraintType::NONLINEAR_NONCONVEX);
            }
        }
    }

    assert(constraintTypes.size() == numConstraints);

    return constraintTypes;
}

bool Constraint::checkFeasibility(const DenseVector &x, double tol) const
{
    assert(x.size() == (int)variables.size());

    if (tol < 0) tol = 0;

    for (unsigned int i = 0; i < variables.size(); i++)
    {
        if (!variables.at(i)->isValueFeasible(x(i), tol))
            return false;
    }

    DenseVector y = eval(x);

    for (unsigned int i = 0; i < numConstraints; i++)
    {
        if (y(i) < lb.at(i) - tol || y(i) > ub.at(i) + tol)
            return false;
    }

    return true;
}

bool Constraint::checkFeasibility(double tol) const
{
    DenseVector x = getVariableValues();
    return checkFeasibility(x, tol);
}

void Constraint::localRefinement(const DenseVector &x)
{
    // Do nothing
}

void Constraint::localRefinement()
{
    DenseVector x = getVariableValues();
    return localRefinement(x);
}

void Constraint::getConstraintBounds(std::vector<double> &lb, std::vector<double> &ub) const
{
    lb = this->lb;
    ub = this->ub;
}

bool Constraint::reduceVariableRanges() const
{
    return true;
}

VariablePtr Constraint::getVariableAt(int index) const
{
    return variables.at(index);
}

std::vector<VariablePtr> Constraint::getVariables() const
{
    return variables;
}

std::vector<VariablePtr> Constraint::getComplicatingVariables() const
{
    if (!constraintConvex)
        return variables;

    std::vector<VariablePtr> vars;
    return vars;
}

/*
 * Return a convex relaxation of the constraint with the following properties:
 * - It is convex
 * - It will become tight when all variable bound intervals go to zero
 */
ConstraintPtr Constraint::getConvexRelaxation()
{
    if (!constraintLinear && !constraintConvex)
    {
        cout << "Constraint::getConvexRelaxation(): Error! Tried to return non-convex constraint." << endl;
        exit(1);
    }
    return this->clone(false); // Shallow copy (using the same variables)
}

// This function causes discontinuities.
// Should be used only when absolutely necessary!
DenseVector Constraint::adjustToDomainBounds(const DenseVector &x) const
{
    assert(x.size() == (int)variables.size());

    DenseVector xadj(x);

    for (unsigned int i = 0; i < variables.size(); i++)
    {
        if (xadj(i) > variables.at(i)->getUpperBound())
        {
            xadj(i) = variables.at(i)->getUpperBound();
        }
        else if(xadj(i) < variables.at(i)->getLowerBound())
        {
            xadj(i) = variables.at(i)->getLowerBound();
        }
    }

    return xadj;
}

void Constraint::checkConstraintSanity() const
{
    // Check codomain bounds
    assert(lb.size() == numConstraints);
    assert(ub.size() == numConstraints);

    for (unsigned int i = 0; i < numConstraints; i++)
    {
        assert(lb.at(i) <= ub.at(i));
    }

    // Check domain bounds
    for (const auto &var : variables)
    {
        // Check variable if needed
    }
}

ProblemClass Constraint::assessClass() const
{
    int numIntVars = 0;
    for (const auto &var : variables)
    {
        VariableType type = var->getType();
        if (type == VariableType::INTEGER || type == VariableType::BINARY)
            numIntVars++;
    }

    if (numIntVars == 0)
    {
        if (isConstraintLinear())
        {
            return ProblemClass::LP;
        }
        else if (isConstraintConvex())
        {
            return ProblemClass::CNLP;
        }
        else
        {
            return ProblemClass::NLP;
        }
    }
    else
    {
        if (isConstraintLinear())
        {
            return ProblemClass::MILP;
        }
        else if (isConstraintConvex())
        {
            return ProblemClass::CMINLP;
        }
        else
        {
            return ProblemClass::MINLP;
        }
    }
}

std::ostream& Constraint::output(std::ostream &os) const
{
    os << "f: \t\tR^" << variables.size() << " |--> R^"<< numConstraints << endl;

    if (constraintName.compare(" ") != 0)
    {
        os << "Name:    \t" << constraintName  << endl;
    }

    if (jacobianCalculated)
    {
        os << "Gradient:\tcalculated with  " << nnzJacobian << " non-zero elements" << endl;
    }
    else
    {
        os << "Gradient:\tnot available "  << endl;
    }

    if (hessianCalculated)
    {
        os << "Hessian:\tcalculated with  " << nnzHessian << " non-zero elements" << endl;
    }
    else
    {
        os << "Hessian:\tnot available "  << endl;
    }

//    os << "Domain lb: \t";
//    os << domainLowerBound;
//    os << "Domain ub: \t";
//    os << domainUpperBound;

    return os;
}

std::ostream& operator<<(std::ostream &os, const Constraint &cs)
{
    cs.output(os);
    return os;
}

void Constraint::writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const
{
}

} // namespace CENSO
