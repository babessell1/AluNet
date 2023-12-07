#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <Rcpp.h>
#include "GraphUtils.h"


class Community {
private:
    // properties
    int communityIndex; 
    std::vector<int> nodeIndices;  // node indices
    
public:
    // getters
    int getCommunityIndex() const {
        return communityIndex;
    }

    std::vector<int> getNodeIndices() const {
        return nodeIndices;
    }

    // setters
    void setCommunityIndex(int communityIndex) {
        communityIndex = communityIndex;
    }

    void setNodeIndices(std::vector<int> nodeIndices) {
        nodeIndices = nodeIndices;
    }
    
    // methods
    Community(const std::vector<int>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(const Graph& G) const; // sum weights of all edges in the community and count number of nodes
    size_t size() const; // number of nodes in the community
    int countPossibleEdges(const Graph& G) const; // count possible edges in the community
    double getClusterWeight(const Graph&G) const; // sum weights of all nodes in the community
    bool hasNode(int node_index) const; // check if the community has a node
    bool hasEdge(int node_index, int neighbor_index, const Graph& G) const; // check if the community has an edge
    bool isEmpty() const; // check if the community is empty
    void addNode(int node_index); // add a node to the community
    void removeNode(int node_index); // remove a node from the community
};

#endif