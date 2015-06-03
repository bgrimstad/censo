/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <cassert>
#include <iostream>
#include <iomanip>
#include <string>
#include <memory>

#include "Utils/definitions.h"
#include "variable.h"

namespace CENSO
{

enum class ConstraintType
{
    LINEAR,
    NONLINEAR_CONVEX,
    NONLINEAR_NONCONVEX
};

enum class ProblemClass
{
    LP,     // Linear Programming
    NLP,    // Non-Linear Programming
    CNLP,   // Convex Non-Linear Programming
    MILP,   // Mixed-Integer Linear Programming
    MINLP,  // Mixed-Integer Non-Linear Programming
    CMINLP  // Convex Mixed-Integer Non-Linear Programming
};

/*
 * Class for generic constraints
 */
class Constraint
{
public:

    Constraint(const Constraint &copy, bool deep = false);
    Constraint& operator=(const Constraint &assign) = delete;
    virtual ~Constraint() {}

    /*
     * Deep copy
     */
    virtual std::shared_ptr<Constraint> clone(bool deep = true) const = 0;

    virtual DenseVector eval(const DenseVector &x) const = 0;

    DenseVector eval() const;

    virtual DenseVector evalJacobian(const DenseVector &x) const = 0;

    DenseVector evalJacobian() const;

    virtual DenseVector evalHessian(const DenseVector &x) const = 0;

    DenseVector evalHessian() const;

    virtual void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) = 0;

    virtual void structureHessian(std::vector<int> &eqnr, std::vector<int> &iRow, std::vector<int> &jCol) = 0;

    virtual bool reduceVariableRanges() const; // Used by domain reduction techniques

    virtual bool checkFeasibility(const DenseVector &x, double tol = 0) const;

    virtual bool checkFeasibility(double tol = 0) const;

    virtual void localRefinement(const DenseVector &x);

    virtual void localRefinement();

    /*
     * Getters
     */
    VariablePtr getVariableAt(int index) const;
    std::vector<VariablePtr> getVariables() const;
    void getConstraintBounds(std::vector<double> &lb, std::vector<double> &ub) const;
    unsigned int getNumVariables() const { return variables.size(); }
    unsigned int getNumConstraints() const { return numConstraints; }
    unsigned int getNumNonZerosJacobian() const { return nnzJacobian; }
    unsigned int getNumNonZerosHessian() const { return nnzHessian; }
    std::vector<ConstraintType> getConstraintTypes();

    /*
     * Returns variables that participate nonlinearly in nonconvex constraints (complicating variables)
     */
    virtual std::vector<VariablePtr> getComplicatingVariables() const;

    /*
     * Returns a convex relaxation of the constraint
     * NOTE: the convex relaxation is in the original
     * variables.
     */
    virtual std::shared_ptr<Constraint> getConvexRelaxation();

    /*
     * Setters
     */
    void setVariables(std::vector<VariablePtr> variables)
    {
        assert(this->variables.size() == variables.size());
        this->variables = variables;
    }

    void setConstraintConvex (bool isConstraintConvex)
    {
        constraintConvex = isConstraintConvex;
    }

    void setName(std::string name)
    {
        constraintName = name;
    }

    /*
     * Checkers
     */

    bool isJacobianCalculated() const { return jacobianCalculated; }
    bool isHessianCalculated() const { return hessianCalculated; }
    bool isConstraintLinear() const { return constraintLinear; }
    bool isConstraintConvex() const { return constraintConvex; }
    bool hasConvexRelaxation() const { return convexRelaxationAvailable; }
    bool hasVariable(VariablePtr variable) const
    {
        return (std::find(variables.begin(), variables.end(), variable) != variables.end());
    }

    /*
     * Assert problem class
     */
    ProblemClass assessClass() const;

    /*
     * Printing
     */
    virtual std::ostream& output(std::ostream &os) const;
    friend std::ostream& operator<<(std::ostream &os, const Constraint &cs);

    void writeToGAMS(const std::string &fname) const;
    virtual void writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const;

protected:

    /*
     * Default constraint:
     * - do not have an analytical gradient or hessian
     * - are non-linear
     * - are non-convex
     * - do not have a convex relaxation available
     *
     * If any of these assumptions are wrong they must be
     * overloaded in the derived constraint class.
     */
    Constraint()
        : numConstraints(0),
          nnzJacobian(0),
          nnzHessian(0),
          jacobianCalculated(false),
          hessianCalculated(false),
          constraintLinear(false),
          constraintConvex(false),
          convexRelaxationAvailable(false),
          constraintName("")
    {}

    Constraint(std::vector<VariablePtr> variables)
        : Constraint()
    {
        this->variables = variables;
    }

    /*
     * Member variables
     */

    std::vector<VariablePtr> variables;

    std::vector<double> lb;
    std::vector<double> ub;

    std::vector<ConstraintType> constraintTypes;

    unsigned int numConstraints;
    unsigned int nnzJacobian;
    unsigned int nnzHessian;

    bool jacobianCalculated;
    bool hessianCalculated;
    bool constraintLinear;
    bool constraintConvex;
    bool convexRelaxationAvailable;

    std::string constraintName;

    /*
     * Helper functions
     */
    DenseVector getVariableValues() const;

    void setVariableValues(const DenseVector &x) const;

    DenseVector adjustToDomainBounds(const DenseVector &x) const;

    void checkConstraintSanity() const;

    void copyVariables();
};

typedef std::shared_ptr<Constraint> ConstraintPtr;

} // namespace CENSO

#endif // CONSTRAINT_H
