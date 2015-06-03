/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef CONSTRAINTSET_H
#define CONSTRAINTSET_H

#include <vector>
#include <cassert>

#include "Utils/definitions.h"
#include "constraint.h"

namespace CENSO
{

class ConstraintSet : public Constraint
{
public:
    // Constructors
    ConstraintSet();
    ConstraintSet(std::vector<VariablePtr> variables)
        : ConstraintSet()
    {
        this->variables = variables;
    }
    ConstraintSet(const ConstraintSet &copy, bool deep = false);
    ConstraintSet& operator=(const ConstraintSet &assign) = delete;

    ~ConstraintSet() override {}

    ConstraintPtr clone(bool deep = true) const override
    {
        return std::make_shared<ConstraintSet>(*this, deep);
    }

    // Must be made visible since they are overridden below
    using Constraint::eval;
    using Constraint::evalJacobian;
    using Constraint::evalHessian;
    using Constraint::checkFeasibility;

    DenseVector eval(const DenseVector &x) const override;

    DenseVector evalJacobian(const DenseVector &x) const override;

    DenseVector evalHessian(const DenseVector &x) const override;

    void structureJacobian(std::vector<int> &iRow, std::vector<int> &jCol) override;

    void structureHessian(std::vector<int> &eqnr,std::vector<int> &iRow, std::vector<int> &jCol) override;

    bool reduceVariableRanges() const override;

    bool checkFeasibility(const DenseVector &x, double tol = 1e-12) const override;

    void localRefinement(const DenseVector &x) override;

    // Prev. used for testing branching on variable infeasibility
    //std::vector<double> variableInfeasibility(DenseVector &x) const;

    // Add constraint to composite
    void add(ConstraintPtr constraint);

    /*
     * Getters
     */
    std::shared_ptr<Constraint> getConvexRelaxation() override;

    std::vector<VariablePtr> getComplicatingVariables() const override;

    unsigned int getNumConstraintObjects() const { return constraints.size(); }

    // Print function
    std::ostream &output(std::ostream &os) const override;

    void writeConstraintEquationsToGAMS(std::ostream &os, unsigned int start) const override;

private:

    /*
     * Member variables
     */
    std::vector<ConstraintPtr> constraints;

    /*
     * Private helper functions
     */
    void addVariables(const std::vector<VariablePtr> &variables);

    bool isVariableSubset(const std::vector<VariablePtr> &variables) const;

    bool hasVariable(const VariablePtr &variable) const;

    std::vector<int> getConstraintVariableIndices(int constraintIndex) const;
};

typedef std::shared_ptr<ConstraintSet> ConstraintSetPtr;

} // namespace CENSO

#endif // CONSTRAINTSET_H
