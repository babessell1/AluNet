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

    void addNodeIndex(int nodeIndex) {
        nodeIndices.push_back(nodeIndex);
    }

    void removeNodeIndex(int nodeIndex) {
        nodeIndices.erase(std::remove(nodeIndices.begin(), nodeIndices.end(), nodeIndex), nodeIndices.end());
    }
    
    // methods
    Community(const std::vector<int>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(const Graph& G) const; // sum weights of all edges in the community and count number of nodes
    size_t size() const; // number of nodes in the community
    int countPossibleEdges(const Graph& G) const; // count possible edges in the community
    double getClusterWeight(const Graph&G) const; // sum weights of all nodes in the community

};

#endif