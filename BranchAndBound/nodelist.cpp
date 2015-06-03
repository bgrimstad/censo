/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#include "nodelist.h"

using std::cout;
using std::endl;

namespace CENSO
{

namespace BB
{

NodeList::NodeList()
    : nodeSelectionStrategy(BEST_FIRST)
{
}

NodeList::NodeList(NodeSelectionStrategy strategy)
    : nodeSelectionStrategy(strategy)
{
}

void NodeList::addNode(NodePtr node)
{
    nodeList.push_back(node);
}

void NodeList::addNodes(NodeList nodes)
{
    for (iterator node = nodes.begin(); node != nodes.end(); node++)
    {
        addNode(*node);
    }
}

int NodeList::size() const
{
    return nodeList.size();
}

bool NodeList::isEmpty() const
{
    return (nodeList.size() <= 0);
}

// Selects from the list of active nodes which node to process next
NodePtr NodeList::selectNode()
{
    // Initialize next node object
    NodePtr nextNode(nullptr);

    if (isEmpty())
    {
        cout << "Trying to select node from empty list. Returning nullptr!" << endl;
        return nextNode;
    }

    // Remember node list size
    unsigned int nodeListSize = nodeList.size();

    // Select node to process
    switch (nodeSelectionStrategy)
    {
    case DEPTH_FIRST:
        // Depth first search strategy (LIFO)
        nextNode = nodeList.back(); // Select node from list
        nodeList.pop_back(); // Remove node from list
        break;

    case BREADTH_FIRST:
        // Breadth first search strategy (FIFO)
        nextNode = nodeList.front(); // Select node from list
        nodeList.pop_front(); // Remove node from list
        break;

    case RANDOM_NODE:
        { // Local scope to allow local variables in switch
            std::list<NodePtr>::iterator next = nodeList.begin();
            int random = CENSO::randomInteger(0,nodeList.size()-1);
            for (int i = 0; i < random; i++) {next++;}
            nextNode = *next;
            nodeList.erase(next); // Remove node from list
        } // End local scope
        break;

    default: // BEST_FIRST
        { // Local scope to allow local variables in switch
            std::list<NodePtr>::iterator best = nodeList.begin();
            for (auto it = nodeList.begin(); it != nodeList.end(); it++)
            {
                // cout << "NODE LB: " << (*it)->getLowerBound() << endl;
                if ((*it)->getLowerBound() <= (*best)->getLowerBound()) best = it;
            }
            nextNode = *best;
            nodeList.erase(best); // Remove node from list
        } // End of local scope
        break;
    }

    // Make sure that exactly one node has been removed
    assert(nodeList.size() == nodeListSize - 1);

    // Check for nullptr
    if (nextNode == nullptr)
    {
        cout << "Could not select node. Returning nullptr!" << endl;
    }

    // Return node
    return nextNode;
}

} // namespace BB

} // namespace CENSO
