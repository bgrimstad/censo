/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef FLOWNETWORKDATA_H
#define FLOWNETWORKDATA_H

#include "Utils/definitions.h"
#include "OptimizationProblem/constraint.h"
#include "OptimizationProblem/constraintset.h"

using std::cout;
using std::endl;
using CENSO::ConstraintPtr;
using CENSO::ConstraintSetPtr;

using namespace CENSO;

namespace FlowNetwork
{

// NOTE: units are assumed consistent in the network
enum class NodeVariable
{
    PRES        // Pressure
};

enum class EdgeVariable
{
    QOIL,       // Flowrate of oil
    QGAS,       // Flowrate of gas
    QWAT,       // Flowrate of water
    QLIQ,       // Flowrate of liquid
    WCT,        // Water-cut
    GOR,        // Gas-oil ratio
    PRES_DIFF,  // Differential pressure
    TEMP,       // Temperature upstream
    TEMP_DIFF,  // Temperature difference
    TEMP_DS,    // Temperature downstream
    ENTH,       // Enthalpy
    ENTH_DIFF,  // Enthalpy difference
    DISC,       // Discrete edge variable
    VEL_MIX     // Mixed in-situ velocity
};

// Some variables maps
typedef std::map<NodeVariable,VariablePtr> node_variable_map;
typedef std::map<EdgeVariable,VariablePtr> edge_variable_map;

template<class VarType>
class FlowNetworkData
{
public:
    FlowNetworkData()
    {
    }

    int getNumVariables() const
    {
        return this->variables.size();
    }

    bool hasVariable(VarType var) const
    {
        if (variables.count(var) > 0)
            return true;
        return false;
    }

    VariablePtr getVariable(VarType var) const
    {
        assert(hasVariable(var));
        return variables.at(var);
    }

    void addConstraint(ConstraintPtr constraint)
    {
        // Check that all constraint variables exist
        auto vars = constraint->getVariables();

        for (auto var : vars)
        {
            bool exists = false;
            for (auto existingVar : variables)
            {
                if (var == existingVar.second)
                {
                    exists = true;
                    break;
                }
            }

            assert(exists);
        }

        constraints.push_back(constraint);
    }

    bool setVariableBounds(VarType var, double lb, double ub, bool intersection = true)
    {
        if (!hasVariable(var))
            return false;

        if (intersection)
        {
            return variables.at(var)->updateBounds(lb, ub);
        }

        variables.at(var)->setLowerBound(lb);
        variables.at(var)->setUpperBound(ub);
        return true;
    }

    void addConstraints(ConstraintSetPtr &cs) const
    {
        for (auto con : constraints)
        {
            cs->add(con);
        }
    }

    bool addLocalVariable(VarType type, VariablePtr variable)
    {
        auto ret = variables.emplace(type, variable);
        return ret.second;
    }

protected:

    // Variable map (maps local to global variable indices)
    std::map<VarType, VariablePtr> variables;

    // Constraints
    std::vector<ConstraintPtr> constraints;
};

/*
 * Object holding node data
 * - associated constraints
 * - associated variables
 */
class NodeData : public FlowNetworkData<NodeVariable>
{
public:
    NodeData() : FlowNetworkData()
    {
        auto pres = std::make_shared<Variable>(0, 0, 1e3);
        variables   = {{NodeVariable::PRES, pres}};
    }
};

/*
 * Object holding edge data
 * - associated constraints
 * - associated variables
 * - if it is discrete
 */
class EdgeData : public FlowNetworkData<EdgeVariable>
{
public:
    EdgeData()
        : FlowNetworkData(),
          is_discrete(false)
    {
        auto qoil = std::make_shared<Variable>(0, 0, 1e3);
        auto qgas = std::make_shared<Variable>(0, 0, 1e3);
        auto qwat = std::make_shared<Variable>(0, 0, 1e3);
        auto pdiff = std::make_shared<Variable>(0, 0, 1e3);

        variables   = {{EdgeVariable::QOIL, qoil},
                       {EdgeVariable::QGAS, qgas},
                       {EdgeVariable::QWAT, qwat},
                       {EdgeVariable::PRES_DIFF, pdiff}};
    }

    bool isDiscrete() const
    {
        return is_discrete;
    }

protected:
    bool is_discrete;
};

/*
 * Discrete edges
 * have a discrete variable
 */
class DiscreteEdgeData : public EdgeData
{
public:
    DiscreteEdgeData()
        : EdgeData()
    {
        is_discrete = true;
        auto disc = std::make_shared<Variable>(0, 0, 1, VariableType::BINARY);
        variables.emplace(EdgeVariable::DISC, disc);
    }
};

} // namespace FlowNetwork

#endif // FLOWNETWORKDATA_H
