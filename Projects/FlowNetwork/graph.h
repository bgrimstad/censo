/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <vector>
#include <algorithm>
#include <memory>
#include <unordered_map>

using std::cout;
using std::endl;

/*
 * Template class implementation of a bidirectional graph (bigraph)
 */
template<class NodeData = double, class EdgeData = double>
class Graph
{
    typedef unsigned int g_id;

    // Internal counter for edge IDs
    g_id next_edge_id;
    g_id get_next_edge_id() { return next_edge_id++; }

    // Define node and edge structs
    struct Node;
    struct Edge;
    typedef std::shared_ptr<Node> node_ptr;
    typedef std::shared_ptr<Edge> edge_ptr;

    struct Node
    {
        g_id id;
        std::vector<edge_ptr> edges_in, edges_out;
        NodeData data;

        Node(g_id id, NodeData data) : id(id), data(data) {}

        bool has_link_to(g_id other_node) const
        {
            for(const auto& e : edges_out)
                if(e->to->id == other_node)
                    return true;
            return false;
        }

        bool has_link_from(g_id other_node) const
        {
            for(const auto& e : edges_out)
                if(e->to->id == other_node)
                    return true;
            return false;
        }

        void remove_edge_to(g_id other_node)
        {
            edges_out.erase(std::remove_if(edges_out.begin(), edges_out.end(), [&other_node](edge_ptr e) { return e->to->id == other_node; }), edges_out.end());
        }

        void remove_edge_from(g_id other_node)
        {
            edges_in.erase(std::remove_if(edges_in.begin(), edges_in.end(), [&other_node](edge_ptr e) { return e->from->id == other_node; }), edges_in.end());
        }

        bool is_source_node() const
        {
            if(edges_in.empty())
                return true;
            return false;
        }

        bool is_sink_node() const
        {
            if(edges_out.empty())
                return true;
            return false;
        }
    };

    struct Edge
    {
        g_id id;
        node_ptr from, to;
        EdgeData data;

        Edge(g_id id, node_ptr from, node_ptr to, EdgeData data)
            : id(id), from(from), to(to), data(data) {}
    };

    // Print function
    friend std::ostream& operator<<(std::ostream& os, const Graph<NodeData,EdgeData> &graph)
    {
        for(const auto& n : graph.nodes)
        {
            os << "Node " << n.first << "-> ";
            for(const auto e : n.second->edges_out)
                os << e->to->id << ", ";
                //os << e->to->id << " (" << e->data << "), ";
            os << endl;
        }
        return os;
    }

public:
    Graph() { next_edge_id = 0; }

    // Map of nodes and edges (using id as key)
    std::unordered_map<g_id, node_ptr> nodes;
    std::unordered_map<g_id, edge_ptr> edges;
    // TODO: consider encapsulating the nodes and edges

    bool add_node(g_id id)
    {
        // Assumes that NodeData() exists
        return add_node(id, NodeData());
    }

    bool add_node(g_id id, NodeData data)
    {
        if(node_exists(id))
        {
            cout << "Node " << id << " already exist!" << endl;
            return false;
        }

        node_ptr n(new Node(id,data));
        nodes.emplace(id, n);
        return true;
    }

    bool remove_node(g_id id)
    {
        if(!node_exists(id))
        {
            cout << "Cannot remove nonexistent node!" << endl;
            return false;
        }

        // Remove all edges to and from node
        for(auto& n : nodes)
        {
            //cout << n.first << endl;
            if(edge_exists(n.first, id))
                remove_edge(n.first, id);
            if(edge_exists(id, n.first))
                remove_edge(id, n.first);
        }

        // Remove node
        nodes.erase(id);

        return true;
    }

    bool add_edge(g_id from, g_id to)
    {
        // Assumes that EdgeData() exists
        return add_edge(from, to, EdgeData());
    }

    bool add_edge(g_id from, g_id to, EdgeData data)
    {
        if(edge_exists(from, to))
        {
            cout << "Edge already exist!" << endl;
            return false;
        }
        else if(!node_exists(from) || !node_exists(to))
        {
            cout << "Nonexistent node(-s). Cannot add edge!" << endl;
            return false;
        }

        edge_ptr e(new Edge(get_next_edge_id(), nodes.at(from), nodes.at(to), data));
        edges.emplace(e->id,e);
        nodes.at(from)->edges_out.push_back(e);
        nodes.at(to)->edges_in.push_back(e);
        return true;
    }

    bool remove_edge(g_id from, g_id to)
    {
        if(!edge_exists(from, to))
        {
            cout << "Cannot remove nonexistent edge!" << endl;
            return false;
        }

        // TODO: remove this code duplication
        auto edge_it = std::find_if(edges.begin(), edges.end(), [from,to](std::pair<g_id, edge_ptr> e){ return (e.second->from->id == from && e.second->to->id == to); });
        if (edge_it != edges.end()) edges.erase(edge_it->first);
        nodes.at(from)->remove_edge_to(to);
        nodes.at(to)->remove_edge_from(from);
        return true;
    }

    bool node_exists(g_id node) const
    {
        return !(nodes.find(node) == nodes.end());
    }

    bool edge_exists(g_id from, g_id to) const
    {
        return std::find_if(edges.begin(), edges.end(), [from,to](std::pair<g_id, edge_ptr> e){ return (e.second->from->id == from && e.second->to->id == to); }) != edges.end();
    }

    edge_ptr get_edge(g_id from, g_id to) const
    {
        auto edge = std::find_if(edges.begin(), edges.end(), [from,to](std::pair<g_id, edge_ptr> e){ return (e.second->from->id == from && e.second->to->id == to); });
        if(edge == edges.end())
        {
            return nullptr;
        }
        else
        {
            return (*edge).second;
        }
    }

};

#endif // GRAPH_H
