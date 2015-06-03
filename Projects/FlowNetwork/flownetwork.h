/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef FLOWNETWORK_H
#define FLOWNETWORK_H

#include "Utils/definitions.h"
#include "graph.h"
#include "flownetworkdata.h"

#include "OptimizationProblem/constraintset.h"
#include "OptimizationProblem/constraintlinear.h"
#include "OptimizationProblem/constraintbspline.h"

using namespace CENSO;
using Splinter::BSpline;
using Splinter::DataTable;

namespace FlowNetwork
{

/* Implementation of a flow network based on a bigraph
 * - three-phase (oil, gas, water) flow specific, but could be generalized to any number of flow components
 * - Currently gas lift can only be included by adding a node and an edge with gas flow only.
 */
class FlowNetworkProblem
{
    typedef unsigned int g_id;

public:
    FlowNetworkProblem() : energy_network(false) {}
    FlowNetworkProblem(bool energy) : energy_network(energy) {}

    bool addNode(g_id id, NodeData data);
    bool addEdge(g_id from, g_id to, EdgeData data);

    void addConstraint(ConstraintPtr constraint);
    bool addConstraintWCT(g_id from, double WCT);
    bool addConstraintWCT(g_id from, g_id to, double WCT);
    bool addConstraintGOR(g_id from, double GOR);
    bool addConstraintGOR(g_id from, g_id to, double GOR);
    bool addConstraintWPC(g_id from, const DataTable &samples);
    bool addConstraintWPC(g_id from, g_id to, const DataTable &samples);
    bool addConstraintWPC2(g_id from, g_id to, const DataTable &pwhToQoil, const DataTable &pwhToPbh);
    bool addConstraintTemperature(g_id from, g_id to, const DataTable &pwhToTwh);
    bool addConstraintFlowlineVLP(g_id from, g_id to, const BSpline &bspline);
    bool addConstraintFlowlineVLPTemp(g_id from, g_id to, const DataTable &samples);
    bool addConstraintJumper(g_id from, g_id to, double dp_a, double dp_b, double dt_a, double dt_b);
    bool addConstraintCapacity(g_id node, EdgeVariable variable, double capacity);
    bool addConstraintVelocity(g_id from, g_id to, const DataTable &samples, double limit);

    bool setNodeVariableBounds(g_id node, NodeVariable var, double lb, double ub, bool intersection = true);
    bool setEdgeVariableBounds(g_id from, g_id to, EdgeVariable var, double lb, double ub, bool intersection = true);

    ConstraintSetPtr buildFlowProblem() const;

    // Getters
    unsigned int getNumNodes() const { return graph.nodes.size(); }
    unsigned int getNumEdges() const { return graph.edges.size(); }

    VariablePtr getVariable(g_id node, NodeVariable variable) const;
    VariablePtr getVariable(g_id from, g_id to, EdgeVariable variable) const;

    void printSolution(std::vector<double> solution, ConstraintPtr cs) const;

    void printGraph() const
    {
        std::cout << graph << std::endl;
    }

private:
    Graph<NodeData,EdgeData> graph;

    bool energy_network;

    std::vector<ConstraintPtr> constraints;

    // Problem builders
    void setObjectiveCosts() const;
    void massConservation(ConstraintSetPtr &cs) const;
    void momentumConservation(ConstraintSetPtr &cs) const;
    void energyConservation(ConstraintSetPtr &cs) const;
    void routingConstraints(ConstraintSetPtr &cs) const;
    void localConstraints(ConstraintSetPtr &cs) const;
    void globalConstraints(ConstraintSetPtr &cs) const;

    bool checkNetwork() const;
};

} // namespace FlowNetwork

#endif // FLOWNETWORK_H
