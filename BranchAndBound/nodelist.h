/*
 * This file is part of the CENSO optimization framework.
 * Copyright (C) 2012 Bjarne Grimstad (bjarne.grimstad@gmail.com).
 *
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
*/

#ifndef NODELIST_H
#define NODELIST_H

#include "node.h"

#include <list>

namespace CENSO
{

namespace BB
{

// Tree search strategies
enum NodeSelectionStrategy
{
    BEST_FIRST, // Default
    DEPTH_FIRST,
    BREADTH_FIRST,
    RANDOM_NODE
};

class NodeList
{
public:
    NodeList();
    NodeList(NodeSelectionStrategy strategy);
    //NodeList(NodeList const& copy) = delete;
    //NodeList& operator = (NodeList const& assign) = delete;
    virtual ~NodeList() {} // All nodes are shared ptrs

    /* Node selection
     * To guarantee global convergence the branch-and-bound algorithm should select
     * the node with the worst bound (lowest lower bound) occasionally.
     *
     * It is important to restrict the number of active nodes to keep within an acceptable memory load.
     * This can be achieved by switching between different search strategies. Typically, BREADTH_FIRST
     * will maintain a large node list, while DEPTH_FIRST maintains a minimal node list. BEST_FIRST will
     * normally maintain more nodes than DEPTH_FIRST and fewer than BREADTH_FIRST.
     */
    NodePtr selectNode();

    void addNode(NodePtr node);

    void addNodes(NodeList nodes);

    int size() const;

    bool isEmpty() const;

    // Forward iterator definition of nodeList
    typedef std::list<NodePtr>::iterator iterator;
    typedef std::list<NodePtr>::const_iterator const_iterator;

    iterator begin() { return nodeList.begin(); }
    const_iterator begin() const { return nodeList.begin(); }
    iterator end() { return nodeList.end(); }
    const_iterator end() const { return nodeList.end(); }

private:

    // Node list (linked-list) - change to unordered map
    std::list<NodePtr> nodeList;

    NodeSelectionStrategy nodeSelectionStrategy;
};

typedef std::shared_ptr<NodeList> NodeListPtr;

} // namespace BB

} // namespace CENSO

#endif // NODELIST_H
