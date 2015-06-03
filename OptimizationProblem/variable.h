/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef VARIABLE_H
#define VARIABLE_H

namespace CENSO
{

enum class VariableType
{
    CONTINUOUS,
    BINARY,
    INTEGER
};

inline std::ostream& operator<<(std::ostream& os, const VariableType type)
{
    switch (type)
    {
        case VariableType::CONTINUOUS: os << "CONTINUOUS"; break;
        case VariableType::BINARY: os << "BINARY"; break;
        case VariableType::INTEGER: os << "INTEGER"; break;
        default: assert(false); break; // Throw exception (unreachable)
    }
    return os;
}

class Variable
{
public:

    /*
     * Default constructor
     */
    Variable(double cost, double lb, double ub, VariableType type)
        : value(0),
          cost(cost),
          lb(lb),
          ub(ub),
          type(type),
          name("x")
    {
        assert(lb <= ub);
    }

    /*
     * Assume continuous variable
     */
    Variable(double cost, double lb, double ub)
        : Variable(cost, lb, ub, VariableType::CONTINUOUS)
    {}

    /*
     * Assume unbounded and continuous variable
     */
    Variable(double cost)
        : Variable(cost, -INF, INF, VariableType::CONTINUOUS)
    {}

    /*
     * Assume unbounded and continuous variable with zero cost
     */
    Variable()
        : Variable(0, -INF, INF, VariableType::CONTINUOUS)
    {}

    virtual ~Variable() {}

    /*
     * Getters
     */
    double getValue() const
    {
        return value;
    }

    double getCost() const
    {
        return cost;
    }

    double getLowerBound() const
    {
        return lb;
    }

    double getUpperBound() const
    {
        return ub;
    }

    VariableType getType() const
    {
        return type;
    }

    std::string getName() const
    {
        return name;
    }

    /*
     * Setters
     */
    void setValue(double value)
    {
        this->value = value;
    }

    void setCost(double cost)
    {
        this->cost = cost;
    }

    // TODO: ceil if integer variable
    void setLowerBound(double lb)
    {
        this->lb = lb;
    }

    // TODO: floor if integer variable
    void setUpperBound(double ub)
    {
        this->ub = ub;
    }

    void setType(VariableType type)
    {
        this->type = type;
    }

    void setName(std::string name)
    {
        this->name = name;
    }

    bool updateLowerBound(double lb)
    {
        if (type == VariableType::BINARY || type == VariableType::INTEGER)
            return updateLowerBoundInteger(lb);

        if (lb > this->ub && !assertNear(lb, this->ub))
            return false;
        else if (lb > this->ub && assertNear(lb, this->ub))
            this->lb = this->ub;
        else if (lb > this->lb && !assertNear(lb, this->lb))
            this->lb = lb;

        return true;
    }

    bool updateUpperBound(double ub)
    {
        if (type == VariableType::BINARY || type == VariableType::INTEGER)
            return updateUpperBoundInteger(ub);

        if (ub < this->lb && !assertNear(ub, this->lb))
            return false;
        else if (ub < this->lb && assertNear(ub, this->lb))
            this->ub = lb;
        else if (ub < this->ub && !assertNear(ub, this->ub))
            this->ub = ub;

        return true;
    }

    /*
     * This function should be used to update
     * variable bounds. Checks correctness and
     * does not allow expansion of range.
     */
    bool updateBounds(double lb, double ub)
    {
        if (lb > ub) return false;

        return updateLowerBound(lb) && updateUpperBound(ub);
    }

    /*
     * Check if value is within bounds
     */
    bool isValueFeasible(double value, double tol = 0) const
    {
        if (value >= lb - tol && value <= ub + tol)
            return true;
        return false;
    }

//    bool isValueFeasible(double tol = 0) const
//    {
//        return isValueFeasible(value, tol);
//    }

    /*
     * Check if integer value
     */
    bool isValueInteger() const
    {
        if (isInteger(value))
            return true;
        return false;
    }

    friend std::ostream& operator<<(std::ostream &os, const Variable &var);

private:
    double value; // Can be used to store optimal value and starting point.
    double cost; // Objective function cost (minimization)
    double lb; // Lower bound
    double ub; // Upper bound
    VariableType type; // Variable type
    std::string name; // Variable name

    bool updateLowerBoundInteger(double lb)
    {
        if (assertNear(lb, this->lb))
            return true;

        if (assertNear(lb, this->ub))
        {
            this->lb = this->ub;
            return true;
        }

        double clb = std::ceil(lb);

        // Check for infeasibility
        if (clb > this->ub)
            return false;

        // Update bound
        if (clb > this->lb)
            this->lb = clb;

        return true;
    }

    bool updateUpperBoundInteger(double ub)
    {
        if (assertNear(ub, this->ub))
            return true;

        if (assertNear(ub, this->lb))
        {
            this->ub = this->lb;
            return true;
        }

        double fub = std::floor(ub);

        // Check for infeasibility
        if (fub < this->lb)
            return false;

        // Update bound
        if (fub < this->ub)
            this->ub = fub;

        return true;
    }
};

typedef std::shared_ptr<Variable> VariablePtr;
typedef std::vector<VariablePtr> Variables;

inline std::ostream& operator<<(std::ostream &os, const Variable &var)
{
    os << "Cost:        " << var.cost   << std::endl;
    os << "Value:       " << var.value  << std::endl;
    os << "Lower bound: " << var.lb     << std::endl;
    os << "Upper bound: " << var.ub     << std::endl;
    os << "Type:        " << var.type   << std::endl;

    return os;
}

} // namespace CENSO


#endif // VARIABLE_H
