/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "flownetwork.h"
#include "Utils/definitions.h"

#include "OptimizationProblem/constraintbspline.h"
#include "OptimizationProblem/constraintbilinear.h"
#include "bspline.h"

using std::cout;
using std::endl;
using Splinter::BSplineType;

namespace FlowNetwork
{

bool FlowNetworkProblem::addNode(g_id id, NodeData data)
{
    return graph.add_node(id, data);
}

bool FlowNetworkProblem::addEdge(g_id from, g_id to, EdgeData data)
{
    if (energy_network)
    {
        // Add enery related variables to edge
        auto temp = std::make_shared<Variable>(0, 0, 1e3);
        auto temp_diff = std::make_shared<Variable>(0, 0, 1e3);
        auto enth = std::make_shared<Variable>(0, 0, 1e3);
        auto enth_diff = std::make_shared<Variable>(0, 0, 1e3);

        assert(data.addLocalVariable(EdgeVariable::TEMP, temp));
        assert(data.addLocalVariable(EdgeVariable::TEMP_DIFF, temp_diff));
        assert(data.addLocalVariable(EdgeVariable::ENTH, enth));
        assert(data.addLocalVariable(EdgeVariable::ENTH_DIFF, enth_diff));
    }

    return graph.add_edge(from, to, data);
}

// Add a constraint on global variables
void FlowNetworkProblem::addConstraint(ConstraintPtr constraint)
{
    constraints.push_back(constraint);
}

/*
 * Adds WCT (%) constraint to all edges leaving node from
 */
bool FlowNetworkProblem::addConstraintWCT(g_id from, double WCT)
{
    assert(graph.node_exists(from));

    auto node = graph.nodes.at(from);

    for (const auto &edge : node->edges_out)
    {
        if (!addConstraintWCT(node->id, edge->to->id, WCT))
            return false;
    }

    return true;
}

/*
 * Adds WCT (%) constraint to edge
 */
bool FlowNetworkProblem::addConstraintWCT(g_id from, g_id to, double WCT)
{
    assert(WCT >= 0);

    if (!graph.edge_exists(from,to))
        return false;

    auto edge = graph.get_edge(from,to);

    if (edge == nullptr)
        return false;

    VariablePtr qoil = edge->data.getVariable(EdgeVariable::QOIL);
    VariablePtr qwat = edge->data.getVariable(EdgeVariable::QWAT);

    std::vector<VariablePtr> vars;
    vars.push_back(qoil);
    vars.push_back(qwat);

    DenseMatrix A; A.setZero(1,2);
    A(0,0) = WCT/(100 - WCT); A(0,1) = -1;
    DenseVector b; b.setZero(1);

    ConstraintPtr con = std::make_shared<ConstraintLinear>(vars, A, b, true);
    std::ostringstream cname;
    cname << "WCT, edge (" << from << "," << to << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    return true;
}

/*
 * Adds GOR (mmscf/mSTB) constraint to all edges leaving node from
 */
bool FlowNetworkProblem::addConstraintGOR(g_id from, double GOR)
{
    assert(graph.node_exists(from));

    auto node = graph.nodes.at(from);

    for (const auto &edge : node->edges_out)
    {
        if (!addConstraintGOR(node->id, edge->to->id, GOR))
            return false;
    }

    return true;
}

/*
 * Adds GOR (mmscf/mSTB) constraint to edge
 */
bool FlowNetworkProblem::addConstraintGOR(g_id from, g_id to, double GOR)
{
    assert(GOR >= 0);

    auto edge = graph.get_edge(from, to);

    if (edge == nullptr)
    {
        cout << "Edge (" << from << "," << to << ") does not exist!" << endl;
        exit(1);
    }

    VariablePtr qoil = edge->data.getVariable(EdgeVariable::QOIL);
    VariablePtr qgas = edge->data.getVariable(EdgeVariable::QGAS);

    std::vector<VariablePtr> vars;
    vars.push_back(qoil);
    vars.push_back(qgas);

    DenseMatrix A; A.setZero(1,2);
    A(0,0) = GOR; A(0,1) = -1;
    DenseVector b; b.setZero(1);

    ConstraintPtr con = std::make_shared<ConstraintLinear>(vars, A, b, true);
    std::ostringstream cname;
    cname << "GOR, edge (" << from << "," << to << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    return true;
}

/*
 * Adds WPC constraint, qoil = f(pwh), to all edges leaving node from
 */
bool FlowNetworkProblem::addConstraintWPC(g_id node, const DataTable &samples)
{
    if (!graph.node_exists(node))
        return false;

    auto n = graph.nodes.at(node);

    // Add to entering edges unless the node is a source node
//    auto edges = n->edges_in;
//    if(n->is_source_node())
//        edges = n->edges_out;

    for (const auto &edge : n->edges_out)
    {
        if (!addConstraintWPC(edge->from->id, edge->to->id, samples))
            return false;
    }
    return true;
}

/*
 * Adds WPC constraint, qoil = f(pwh), to edge (from,to)
 */
bool FlowNetworkProblem::addConstraintWPC(g_id from, g_id to, const DataTable &samples)
{
    assert(samples.isGridComplete());
    assert(samples.getNumVariables() == 1);

    assert(graph.node_exists(from) && graph.node_exists(to));
    assert(graph.edge_exists(from,to));
    auto node = graph.nodes.at(from);
    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    VariablePtr pi = node->data.getVariable(NodeVariable::PRES);
    VariablePtr qoil = edge->data.getVariable(EdgeVariable::QOIL);

    std::vector<VariablePtr> vars;
    vars.push_back(pi);
    vars.push_back(qoil);

    BSpline bs(samples, BSplineType::CUBIC_FREE); // LINEAR
    ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bs, true);
    std::ostringstream cname;
    cname << "WPC, edge (" << edge->from->id << "," << edge->to->id << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    return true;
}

/*
 * Adds WPC constraint, qoil = f(pwh), to edge (from,to)
 */
bool FlowNetworkProblem::addConstraintWPC2(g_id from, g_id to, const DataTable &pwhToQoil, const DataTable &pwhToPbh)
{
    assert(pwhToQoil.isGridComplete());
    assert(pwhToQoil.getNumVariables() == 1);
    assert(pwhToPbh.isGridComplete());
    assert(pwhToPbh.getNumVariables() == 1);

    assert(graph.node_exists(from) && graph.node_exists(to));
    assert(graph.edge_exists(from,to));
    auto node_from = graph.nodes.at(from);
    auto node_to = graph.nodes.at(to);
    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    VariablePtr p_from = node_from->data.getVariable(NodeVariable::PRES);
    VariablePtr p_to = node_to->data.getVariable(NodeVariable::PRES);
    VariablePtr qoil = edge->data.getVariable(EdgeVariable::QOIL);

    // Add qoil = f(pwh) constraint
    std::vector<VariablePtr> vars_rate;
    vars_rate.push_back(p_to);
    vars_rate.push_back(qoil);

    BSpline bs_rate(pwhToQoil, BSplineType::CUBIC_FREE);
    ConstraintPtr con_rate = std::make_shared<ConstraintBSpline>(vars_rate, bs_rate, true);
    std::ostringstream cname_rate;
    cname_rate << "WPC rate, edge (" << from << "," << to << ")";
    con_rate->setName(cname_rate.str());

    constraints.push_back(con_rate);

    // Add pbh = f(pwh) constraint
    std::vector<VariablePtr> vars_pres;
    vars_pres.push_back(p_to);
    vars_pres.push_back(p_from);

    BSpline bs_pres(pwhToPbh, BSplineType::CUBIC_FREE);
    ConstraintPtr con_pres = std::make_shared<ConstraintBSpline>(vars_pres, bs_pres, true);
    std::ostringstream cname_pres;
    cname_pres << "WPC pres, edge (" << from << "," << to << ")";
    con_pres->setName(cname_pres.str());

    constraints.push_back(con_pres);

    return true;
}

/*
 * WPC wellhead temperature
 * twh = f(pwh)
 */
bool FlowNetworkProblem::addConstraintTemperature(g_id from, g_id to, const DataTable &pwhToTwh)
{
    assert(energy_network);

    if (!graph.edge_exists(from, to))
        return false;

    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    // Add variable for downstream temperature
    if (!edge->data.hasVariable(EdgeVariable::TEMP_DS))
        assert(edge->data.addLocalVariable(EdgeVariable::TEMP_DS, std::make_shared<Variable>(0, 0, 1e3)));

    VariablePtr temp = edge->data.getVariable(EdgeVariable::TEMP);
    VariablePtr temp_diff = edge->data.getVariable(EdgeVariable::TEMP_DIFF);
    VariablePtr temp_ds = edge->data.getVariable(EdgeVariable::TEMP_DS);
    VariablePtr pres_to = edge->to->data.getVariable(NodeVariable::PRES);

    // Add relation t = t_diff + t_ds
    {
        std::vector<VariablePtr> vars;
        vars.push_back(temp);
        vars.push_back(temp_diff);
        vars.push_back(temp_ds);

        DenseMatrix Am; Am.setZero(1,3);
        Am(0,0) = 1; Am(0,1) = -1; Am(0,2) = -1;
        DenseVector bm(1); bm(0) = 0;
        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, Am, bm, true);

        std::ostringstream cname;
        cname << "Linear temperature constraint, edge (" << from << "," << to << ")";
        lincon->setName(cname.str());

        constraints.push_back(lincon);
    }

    // Add twh = f(pwh) constraint
    {
        assert(pwhToTwh.isGridComplete());
        assert(pwhToTwh.getNumVariables() == 1);

        std::vector<VariablePtr> vars;
        vars.push_back(pres_to);
        vars.push_back(temp_ds);

        BSpline bs(pwhToTwh, BSplineType::CUBIC_FREE);
        ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bs, true);
        std::ostringstream cname;
        cname << "WPC temp, edge (" << from << "," << to << ")";
        con->setName(cname.str());

        constraints.push_back(con);
    }

    return true;
}

bool FlowNetworkProblem::addConstraintFlowlineVLP(g_id from, g_id to, const BSpline &bspline)
{
    assert(bspline.getNumVariables() == 4);

    if (!graph.edge_exists(from, to))
        return false;

    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    // Augment edge with new local variables
    auto qliq = std::make_shared<Variable>(0, 0, INF);
    auto gor = std::make_shared<Variable>(0, 0, INF);
    auto wct = std::make_shared<Variable>(0, 0, 100);

    assert(edge->data.addLocalVariable(EdgeVariable::QLIQ, qliq));
    assert(edge->data.addLocalVariable(EdgeVariable::GOR, gor));
    assert(edge->data.addLocalVariable(EdgeVariable::WCT, wct));


    auto from_pres = edge->from->data.getVariable(NodeVariable::PRES);
    auto to_pres = edge->to->data.getVariable(NodeVariable::PRES);
    auto diff_pres = edge->data.getVariable(EdgeVariable::PRES_DIFF);
    auto qoil = edge->data.getVariable(EdgeVariable::QOIL);
    auto qgas = edge->data.getVariable(EdgeVariable::QGAS);
    auto qwat = edge->data.getVariable(EdgeVariable::QWAT);

    // This order is fixed
    // NOTE: change from to_pres to diff_pres
    //std::vector<int> vars = {qoil, qgas, qwat, from_pres, to_pres};
    std::vector<VariablePtr> vars = {qliq, from_pres, gor, wct, to_pres};

    ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bspline, true);
    std::ostringstream cname;
    cname << "Pressure correlation, edge (" << edge->from->id << "," << edge->to->id << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    // Add local constraints for pipeline LIQ, GOR, and WC
    std::vector<VariablePtr> liqVars = {qliq, qoil, qwat};
    DenseMatrix A1; A1.setZero(1,3);
    A1(0,0) = 1; A1(0,1) = -1; A1(0,2) = -1;
    DenseVector b1; b1.setZero(1);
    ConstraintPtr c1 = std::make_shared<ConstraintLinear>(liqVars, A1, b1, true);
    edge->data.addConstraint(c1);


    std::vector<VariablePtr> vGOR = {qoil, gor, qgas};
    ConstraintPtr bl_gor = std::make_shared<ConstraintBilinear>(vGOR, 1, 0, 0);
    edge->data.addConstraint(bl_gor);

    std::vector<VariablePtr> vWCT = {qliq, wct, qwat};
    ConstraintPtr bl_wc = std::make_shared<ConstraintBilinear>(vWCT, 1/100.0, 0, 0); // Divided by 100 since WC is in %
    edge->data.addConstraint(bl_wc);

    return true;
}

bool FlowNetworkProblem::addConstraintFlowlineVLPTemp(g_id from, g_id to, const DataTable &samples)
{
    assert(energy_network);
    assert(samples.isGridComplete());
    assert(samples.getNumVariables() == 5);

    if(!graph.edge_exists(from, to))
        return false;

    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    // Check edge variables
    assert(edge->data.hasVariable(EdgeVariable::QLIQ));
    assert(edge->data.hasVariable(EdgeVariable::GOR));
    assert(edge->data.hasVariable(EdgeVariable::WCT));

    auto from_pres = edge->from->data.getVariable(NodeVariable::PRES);
    auto qliq = edge->data.getVariable(EdgeVariable::QLIQ);
    auto gor = edge->data.getVariable(EdgeVariable::GOR);
    auto wct = edge->data.getVariable(EdgeVariable::WCT);
    auto temp = edge->data.getVariable(EdgeVariable::TEMP);
    auto temp_diff = edge->data.getVariable(EdgeVariable::TEMP_DIFF);

    // This order is fixed
    std::vector<VariablePtr> vars = {qliq, from_pres, gor, wct, temp, temp_diff};

    BSpline bs(samples, BSplineType::LINEAR); // CUBIC_FREE
    ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bs, true);
    std::ostringstream cname;
    cname << "Pressure correlation, edge (" << edge->from->id << "," << edge->to->id << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    return true;
}

/*
 * Constraint on jumper pressure and temperature drop
 */
bool FlowNetworkProblem::addConstraintJumper(g_id from, g_id to, double dp_a, double dp_b, double dt_a, double dt_b)
{
    if (!graph.edge_exists(from, to))
        return false;

    auto e = graph.get_edge(from, to);

    if (e == nullptr)
        return false;

    auto qoil = e->data.getVariable(EdgeVariable::QOIL);
    auto pres_diff = e->data.getVariable(EdgeVariable::PRES_DIFF);
    auto temp_diff = e->data.getVariable(EdgeVariable::TEMP_DIFF);

    {
        std::vector<VariablePtr> vars;
        vars.push_back(pres_diff);
        vars.push_back(qoil);

        DenseMatrix A; A.setOnes(1,2);
        A(0,0) = 1; A(0,1) = -dp_a;
        DenseVector b(1); b(0) = dp_b; // b cannot be negative since q=0 gives dp = -b < 0 during bound deduction
        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, true);

        std::ostringstream cname;
        cname << "Jumper pressure diff, edge (" << from << ", " << to << ")";
        lincon->setName(cname.str());

        constraints.push_back(lincon);
    }

    {
        std::vector<VariablePtr> vars;
        vars.push_back(temp_diff);
        vars.push_back(qoil);

        DenseMatrix A; A.setOnes(1,2);
        A(0,0) = 1; A(0,1) = -dt_a;
        DenseVector b(1); b(0) = dt_b;
        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, true);

        std::ostringstream cname;
        cname << "Jumper temperature diff, edge (" << from << ", " << to << ")";
        lincon->setName(cname.str());

//        constraints.push_back(lincon);
//        constraint_variables.push_back(vars);
    }

    return true;
}

/*
 * Constraint on flow rate through node
 */
bool FlowNetworkProblem::addConstraintCapacity(g_id node, EdgeVariable variable, double capacity)
{
    if (!graph.node_exists(node))
        return false;

    // Only rates are allowed
    if (variable != EdgeVariable::QOIL
    && variable != EdgeVariable::QGAS
    && variable != EdgeVariable::QWAT)
        return false;

    auto n = graph.nodes.at(node);

    auto edges = n->edges_out;
    if (n->is_sink_node())
        edges = n->edges_in;

    std::vector<VariablePtr> vars;
    for (const auto &e : edges)
    {
        auto ret = e->data.getVariable(variable);
        vars.push_back(ret);
    }

    DenseMatrix A; A.setOnes(1,vars.size());
    DenseVector b(1); b(0) = capacity;
    ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);

    std::ostringstream cname;
    cname << "Capacity constraint, node (" << node << ")";
    lincon->setName(cname.str());

    constraints.push_back(lincon);

    return true;
}

bool FlowNetworkProblem::addConstraintVelocity(g_id from, g_id to, const DataTable &samples, double limit)
{
    assert(energy_network);
    assert(samples.isGridComplete());
    assert(samples.getNumVariables() == 5);

    if (!graph.edge_exists(from, to))
        return false;

    auto edge = graph.get_edge(from, to);
    assert(edge != nullptr);

    // Augment edge with new local variables
    if (!edge->data.hasVariable(EdgeVariable::QLIQ))
        assert(edge->data.addLocalVariable(EdgeVariable::QLIQ, std::make_shared<Variable>(0, 0, INF)));
    if (!edge->data.hasVariable(EdgeVariable::GOR))
        assert(edge->data.addLocalVariable(EdgeVariable::GOR, std::make_shared<Variable>(0, 0, INF)));
    if (!edge->data.hasVariable(EdgeVariable::WCT))
        assert(edge->data.addLocalVariable(EdgeVariable::WCT, std::make_shared<Variable>(0, 0, 100)));
    if (!edge->data.hasVariable(EdgeVariable::VEL_MIX))
        assert(edge->data.addLocalVariable(EdgeVariable::VEL_MIX, std::make_shared<Variable>(0, 0, limit)));
    if (!edge->data.hasVariable(EdgeVariable::TEMP_DS))
        assert(edge->data.addLocalVariable(EdgeVariable::TEMP_DS, std::make_shared<Variable>(0, 0, 1e3)));

    auto from_pres = edge->from->data.getVariable(NodeVariable::PRES);
    auto to_pres = edge->to->data.getVariable(NodeVariable::PRES);

    auto from_temp = edge->data.getVariable(EdgeVariable::TEMP);
    auto temp_diff = edge->data.getVariable(EdgeVariable::TEMP_DIFF);
    auto to_temp = edge->data.getVariable(EdgeVariable::TEMP_DS);
    auto qliq = edge->data.getVariable(EdgeVariable::QLIQ);
    auto gor = edge->data.getVariable(EdgeVariable::GOR);
    auto wct = edge->data.getVariable(EdgeVariable::WCT);
    auto vel = edge->data.getVariable(EdgeVariable::VEL_MIX);

    // Constraint variables
    std::vector<VariablePtr> vars = {to_pres, to_temp, qliq, gor, wct, vel};

    BSpline bs(samples, BSplineType::LINEAR); // CUBIC_FREE
    ConstraintPtr con = std::make_shared<ConstraintBSpline>(vars, bs, true);
    std::ostringstream cname;
    cname << "Velocity constraint, edge (" << edge->from->id << "," << edge->to->id << ")";
    con->setName(cname.str());

    constraints.push_back(con);

    // Add local temperature constraints (tus - tdiff = tds)
    std::vector<VariablePtr> vars_temp = {from_temp, temp_diff, to_temp};
    DenseMatrix A1; A1.setZero(1,3);
    A1(0,0) = 1; A1(0,1) = -1; A1(0,2) = -1;
    DenseVector b1; b1.setZero(1);
    ConstraintPtr c1 = std::make_shared<ConstraintLinear>(vars_temp, A1, b1, true);

    edge->data.addConstraint(c1);

    return true;
}

bool FlowNetworkProblem::setNodeVariableBounds(g_id node, NodeVariable var, double lb, double ub, bool intersection)
{
    if (!graph.node_exists(node))
        return false;

    auto n = graph.nodes.at(node);

    if (!n->data.setVariableBounds(var, lb, ub, intersection))
        return false;

    return true;
}

bool FlowNetworkProblem::setEdgeVariableBounds(g_id from, g_id to, EdgeVariable var, double lb, double ub, bool intersection)
{
    if (!graph.edge_exists(from, to))
        return false;

    auto edge = graph.get_edge(from, to);

    if (edge == nullptr)
        return false;

    if (!edge->data.setVariableBounds(var, lb, ub, intersection))
        return false;

    return true;
}

// Passing shared_ptr by reference here since they will be assigned
ConstraintSetPtr FlowNetworkProblem::buildFlowProblem() const
{
    assert(checkNetwork());

    // NOTE: check for auxiliary variables
    ConstraintSetPtr cs = std::make_shared<ConstraintSet>();

    massConservation(cs);
    momentumConservation(cs);
    if (energy_network)
        energyConservation(cs);
    routingConstraints(cs);
    localConstraints(cs);
    globalConstraints(cs);

    setObjectiveCosts();

    return cs;
}

void FlowNetworkProblem::setObjectiveCosts() const
{
    // Maximize sum of source oil rates
    for (const auto& n : graph.nodes)
    {
        auto node = n.second;
        if (node->is_source_node())
        {
            for (const auto& e : node->edges_out)
            {
                auto var = e->data.getVariable(EdgeVariable::QOIL);
                var->setCost(-1);
            }
        }
    }
}

void FlowNetworkProblem::massConservation(ConstraintSetPtr &cs) const
{
    for (const auto &n : graph.nodes)
    {
        auto node = n.second;

        if (node->is_source_node() || node->is_sink_node())
            continue;

        int num_edges = node->edges_in.size()+node->edges_out.size();

        if (num_edges == 0)
            continue;

        int ai = 0;
        DenseMatrix A; A.setZero(1,num_edges);
        DenseVector b; b.setZero(1);

        std::vector<VariablePtr> qoil_vars;
        std::vector<VariablePtr> qgas_vars;
        std::vector<VariablePtr> qwat_vars;

        // Flow rates into node
        for (const auto &ein : node->edges_in)
        {
            A(0,ai++) = 1;

            auto qoil = ein->data.getVariable(EdgeVariable::QOIL);
            qoil_vars.push_back(qoil);

            auto qgas = ein->data.getVariable(EdgeVariable::QGAS);
            qgas_vars.push_back(qgas);

            auto qwat = ein->data.getVariable(EdgeVariable::QWAT);
            qwat_vars.push_back(qwat);
        }

        // Flow rates out of node
        for (const auto &eout : node->edges_out)
        {
            A(0,ai++) = -1;
            auto qoil = eout->data.getVariable(EdgeVariable::QOIL);
            qoil_vars.push_back(qoil);

            auto qgas = eout->data.getVariable(EdgeVariable::QGAS);
            qgas_vars.push_back(qgas);

            auto qwat = eout->data.getVariable(EdgeVariable::QWAT);
            qwat_vars.push_back(qwat);
        }

        {
            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(qoil_vars, A, b, true);

            std::ostringstream cname_oil;
            cname_oil << "Mass balance oil, node " << node->id;
            lincon->setName(cname_oil.str());
            cs->add(lincon);
        }

        {
            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(qgas_vars, A, b, true);

            std::ostringstream cname_gas;
            cname_gas << "Mass balance gas, node " << node->id;
            lincon->setName(cname_gas.str());
            cs->add(lincon);
        }

        {
            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(qwat_vars, A, b, true);

            std::ostringstream cname_wat;
            cname_wat << "Mass balance water, node " << node->id;
            lincon->setName(cname_wat.str());
            cs->add(lincon);
        }
    }
}

void FlowNetworkProblem::momentumConservation(ConstraintSetPtr &cs) const
{
    for (const auto &e : graph.edges)
    {
        auto edge = e.second;

        auto from_pres = edge->from->data.getVariable(NodeVariable::PRES);
        auto to_pres = edge->to->data.getVariable(NodeVariable::PRES);
        auto diff_pres = edge->data.getVariable(EdgeVariable::PRES_DIFF);

        // Add diff pressure relation
        if (edge->data.isDiscrete())
        {
            auto disc = edge->data.getVariable(EdgeVariable::DISC);

            /*
             * Big-M constraints:
             * -M(1-ye) <= pi - pj - Dpe <= M(1-ye),
             * where M = (pi_U - pi_L) + (pj_U - pj_L)
             */
            double M = 1000;

            DenseMatrix A; A.setZero(2,4);
            A(0,0) = 1; A(0,1) = -1; A(0,2) = -1; A(0,3) = M;
            A(1,0) = -1; A(1,1) = 1; A(1,2) = 1; A(1,3) = M;

            DenseVector b; b.setZero(2);
            b(0) = M; b(1) = M;

            std::vector<VariablePtr> vars = {from_pres, to_pres, diff_pres, disc};
            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, false);

            std::ostringstream cname;
            cname << "Pressure drop, edge (" << edge->from->id << "," << edge->to->id << ")";
            lincon->setName(cname.str());
            cs->add(lincon);
        }
        else
        {
            // to_pres - from_pres - diff_pres = 0
            DenseMatrix A; A.setZero(1,3);
            A(0,0) = 1; A(0,1) = -1; A(0,2) = -1;
            DenseVector b; b.setZero(1);

            std::vector<VariablePtr> vars = {from_pres, to_pres, diff_pres};
            ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, true);

            std::ostringstream cname;
            cname << "Pressure drop, edge (" << edge->from->id << "," << edge->to->id << ")";
            lincon->setName(cname.str());
            cs->add(lincon);
        }
    }

}

void FlowNetworkProblem::energyConservation(ConstraintSetPtr &cs) const
{
    assert(energy_network);

    // Energy balance
    for (const auto &n : graph.nodes)
    {
        auto node = n.second;

        if (node->is_source_node() || node->is_sink_node())
            continue;

        auto num_edges_in = node->edges_in.size();
        auto num_edges_out = node->edges_out.size();

        DenseMatrix A; A.setZero(1, 2*num_edges_in+num_edges_out);
        DenseVector b; b.setZero(1);
        std::vector<VariablePtr> vars;

        int count = 0;
        for (const auto &e_in : node->edges_in)
        {
            A(0,count++) = 1;
            A(0,count++) = -1;

            auto enth = e_in->data.getVariable(EdgeVariable::ENTH);
            vars.push_back(enth);

            auto enth_diff = e_in->data.getVariable(EdgeVariable::ENTH_DIFF);
            vars.push_back(enth_diff);
        }

        for (const auto &e_out : node->edges_out)
        {
            A(0,count++) = -1;

            auto enth = e_out->data.getVariable(EdgeVariable::ENTH);
            vars.push_back(enth);
        }

        ConstraintPtr lincon = std::make_shared<ConstraintLinear>(vars, A, b, true);

        std::ostringstream cname;
        cname << "Energy balance, node (" << n.first << ")";
        lincon->setName(cname.str());
        cs->add(lincon);
    }

    // Add enthalpy equations
    for (const auto &e : graph.edges)
    {
        auto edge = e.second;

        // Heat capacities
        double coil = 1/1000.0; // Scaled
        double cgas = 1/1000.0; // Scaled
        double cwat = 1/1000.0; // Scaled

        auto enth = edge->data.getVariable(EdgeVariable::ENTH);
        auto enth_diff = edge->data.getVariable(EdgeVariable::ENTH_DIFF);
        auto temp = edge->data.getVariable(EdgeVariable::TEMP);
        auto temp_diff = edge->data.getVariable(EdgeVariable::TEMP_DIFF);
        auto qoil = edge->data.getVariable(EdgeVariable::QOIL);
        auto qgas = edge->data.getVariable(EdgeVariable::QGAS);
        auto qwat = edge->data.getVariable(EdgeVariable::QWAT);

        // NOTE: bounds deduction is only implemented for bilinear equality constraints
        // Thus, if the bounds are relaxed CENSO may not converge!
        {
            std::vector<VariablePtr> vars = {temp, qoil, enth};
            ConstraintPtr blin = std::make_shared<ConstraintBilinear>(vars, coil, 0, 0.0001); // Prev. relaxed by 0.0001
            cs->add(blin);
        }

        {
            std::vector<VariablePtr> vars = {temp_diff, qoil, enth_diff};
            ConstraintPtr blin = std::make_shared<ConstraintBilinear>(vars, coil, 0, 0.0001); // Prev. relaxed by 0.0001
            cs->add(blin);
        }

    }
}

void FlowNetworkProblem::routingConstraints(ConstraintSetPtr &cs) const
{
    for(const auto& e : graph.edges)
    {
        auto edge = e.second;

        if (edge->data.isDiscrete())
        {
            std::vector<EdgeVariable> qvars;
            qvars.push_back(EdgeVariable::QOIL);
            qvars.push_back(EdgeVariable::QGAS);
            qvars.push_back(EdgeVariable::QWAT);

            auto bvar = edge->data.getVariable(EdgeVariable::DISC);

            for (const auto& qvar : qvars)
            {
                auto var = edge->data.getVariable(qvar);

                std::vector<VariablePtr> vars = {var, bvar};

                // Collapse of upper bound
                DenseMatrix A1; A1.setZero(1,2);
                A1(0,0) = 1; A1(0,1) = -var->getUpperBound();
                DenseVector b1; b1.setZero(1);
                ConstraintPtr c1 = std::make_shared<ConstraintLinear>(vars, A1, b1, false);
                std::ostringstream cname1;
                cname1 << "Flow restriction, edge (" << edge->from->id << "," << edge->to->id << ")";
                c1->setName(cname1.str());
                cs->add(c1);

                // Collapse of lower bound
                DenseMatrix A2; A2.setZero(1,2);
                A2(0,0) = -1; A2(0,1) = var->getLowerBound();
                DenseVector b2; b2.setZero(1);
                ConstraintPtr c2 = std::make_shared<ConstraintLinear>(vars, A2, b2, false);
                std::ostringstream cname2;
                cname2 << "Flow restriction, edge (" << edge->from->id << "," << edge->to->id << ")";
                c2->setName(cname2.str());
                cs->add(c2);
            }
        }
    }

    // Loop through nodes and add routing constraints if it has several outgoing edges
    for (const auto& n : graph.nodes)
    {
        auto node = n.second;

        std::vector<VariablePtr> disc_vars;

        for (const auto& e : node->edges_out)
        {
            if (e->data.isDiscrete())
            {
               auto ret = e->data.getVariable(EdgeVariable::DISC);
               disc_vars.push_back(ret);
            }
        }

        if (disc_vars.size() > 0)
        {
            // Add routing constraint
            DenseMatrix A;
            A.setOnes(1,disc_vars.size());
            DenseVector b(1); b(0) = 1;
            ConstraintPtr con = std::make_shared<ConstraintLinear>(disc_vars, A, b, false);
            std::ostringstream cname;
            cname << "Routing constraint, node " << node->id;
            con->setName(cname.str());
            cs->add(con);
        }
    }
}

void FlowNetworkProblem::localConstraints(ConstraintSetPtr &cs) const
{
    for(const auto &n : graph.nodes)
    {
        auto node = n.second;
        node->data.addConstraints(cs);
    }

    for(const auto &e : graph.edges)
    {
        auto edge = e.second;
        edge->data.addConstraints(cs);
    }
}

void FlowNetworkProblem::globalConstraints(ConstraintSetPtr &cs) const
{
    for (auto con : constraints)
        cs->add(con);
}

VariablePtr FlowNetworkProblem::getVariable(g_id node, NodeVariable variable) const
{
    assert(graph.node_exists(node));
    return graph.nodes.at(node)->data.getVariable(variable);
}

VariablePtr FlowNetworkProblem::getVariable(g_id from, g_id to, EdgeVariable variable) const
{
    assert(graph.edge_exists(from, to));

    auto edge = graph.get_edge(from, to);

    return edge->data.getVariable(variable);
}

void FlowNetworkProblem::printSolution(std::vector<double> solution, ConstraintPtr cs) const
{
    assert(solution.size() == cs->getNumVariables());

    // Assuming that solution is mapped directly to constraint variables
    std::vector<VariablePtr> vars = cs->getVariables();

    cout << endl << "EDGE FLOW RATES" << endl;
    for (const auto &e : graph.edges)
    {
        auto edge = e.second;

        auto qoil = edge->data.getVariable(EdgeVariable::QOIL);
        auto qgas = edge->data.getVariable(EdgeVariable::QGAS);
        auto qwat = edge->data.getVariable(EdgeVariable::QWAT);

        auto qoil_it = std::find(vars.begin(), vars.end(), qoil);
        auto qoil_i = qoil_it - vars.begin();

        auto qgas_it = std::find(vars.begin(), vars.end(), qgas);
        auto qgas_i = qgas_it - vars.begin();

        auto qwat_it = std::find(vars.begin(), vars.end(), qwat);
        auto qwat_i = qwat_it - vars.begin();

        cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
        cout << "(o,g,w,l) = (" << solution.at(qoil_i) << ", " << solution.at(qgas_i) << ", " << solution.at(qwat_i)
             << ", " << solution.at(qoil_i) + solution.at(qwat_i) << ")" << endl;
    }

    cout << endl << "EDGE ROUTING" << endl;
    for (const auto &e : graph.edges)
    {
        auto edge = e.second;

        if (!edge->data.isDiscrete())
            continue;

        auto disc = edge->data.getVariable(EdgeVariable::DISC);
        auto disc_it = std::find(vars.begin(), vars.end(), disc);
        auto disc_i = disc_it - vars.begin();

        cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
        cout << "(routing) = (" << solution.at(disc_i) << ")" << endl;
    }

    cout << endl << "EDGE DIFFERENTIAL PRESSURES" << endl;
    for (const auto &e : graph.edges)
    {
        auto edge = e.second;

        auto diff_pres = edge->data.getVariable(EdgeVariable::PRES_DIFF);
        auto dp_it = std::find(vars.begin(), vars.end(), diff_pres);
        auto dp = solution.at(dp_it - vars.begin());

        cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
        cout << dp << " (" << diff_pres->getLowerBound() << " < " << diff_pres->getUpperBound() << ")" << endl;
    }

    cout << endl << "NODE PRESSURES" << endl;
    for (const auto &n : graph.nodes)
    {
        auto node = n.second;
        auto pres = node->data.getVariable(NodeVariable::PRES);
        auto pres_it = std::find(vars.begin(), vars.end(), pres);
        auto pres_sol = solution.at(pres_it - vars.begin());

        cout << "Node (" << node->id << "): \t";
        cout << pres_sol << " (" << pres->getLowerBound() << " < " << pres->getUpperBound() << ")" << endl;
    }

    if (energy_network)
    {
        cout << endl << "EDGE TEMPERATURES" << endl;
        for (const auto &e : graph.edges)
        {
            auto edge = e.second;

            auto temp = edge->data.getVariable(EdgeVariable::TEMP);
            auto temp_it = std::find(vars.begin(), vars.end(), temp);
            auto temp_sol = solution.at(temp_it - vars.begin());

            auto temp_diff = edge->data.getVariable(EdgeVariable::TEMP_DIFF);
            auto td_it = std::find(vars.begin(), vars.end(), temp_diff);
            auto td_sol = solution.at(td_it - vars.begin());

            cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
            cout << temp_sol << " - " << td_sol << " = " << temp_sol - td_sol << endl;
        }
    }

    if (energy_network)
    {
        cout << endl << "EDGE ENTHALPIES" << endl;
        for (const auto &e : graph.edges)
        {
            auto edge = e.second;

            auto enth = edge->data.getVariable(EdgeVariable::ENTH);
            auto enth_it = std::find(vars.begin(), vars.end(), enth);
            auto enth_sol = solution.at(enth_it - vars.begin());

            auto enth_diff = edge->data.getVariable(EdgeVariable::ENTH_DIFF);
            auto ed_it = std::find(vars.begin(), vars.end(), enth_diff);
            auto ed_sol = solution.at(ed_it - vars.begin());

            cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
            cout << enth_sol << " - " << ed_sol << " = " << enth_sol - ed_sol << endl;
        }
    }

    if(energy_network)
    {
        cout << endl << "EDGE VELOCITIES" << endl;
        for (const auto &e : graph.edges)
        {
            auto edge = e.second;

            if (!edge->data.hasVariable(EdgeVariable::VEL_MIX))
                continue;

            auto vel = edge->data.getVariable(EdgeVariable::VEL_MIX);
            auto vel_it = std::find(vars.begin(), vars.end(), vel);
            auto vel_sol = solution.at(vel_it - vars.begin());

            cout << "Edge (" << edge->from->id << "," << edge->to->id << "): \t";
            cout << vel_sol << endl;
        }
    }

}

bool FlowNetworkProblem::checkNetwork() const
{
    bool sourceFound = false;
    bool sinkFound = false;
    for (auto n : graph.nodes)
    {
        if(n.second->is_source_node())
            sourceFound = true;
        if(n.second->is_sink_node())
            sinkFound = true;
        if(sourceFound && sinkFound)
            break;
    }
    if (!(sourceFound && sinkFound))
    {
        cout << "Network does not have a source and sink node!" << endl;
        return false;
    }

    return true;
}

} // namespace FlowNetwork
