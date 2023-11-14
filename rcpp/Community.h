#ifndef COMMUNITY_H
#define COMMUNITY_H

#include <Rcpp.h>
#include "GraphUtils.h"


class Community {
public:
    // properties
    int communityIndex; 
    std::vector<int> nodeIndices;  // node indices
    
    // methods
    Community(const std::vector<int>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(const Graph& G); // sum weights of all edges in the community and count number of nodes
    size_t size(); // number of nodes in the community
    int countPossibleEdges(const Graph& G); // count possible edges in the community
    double getClusterWeight(const Graph&G) const; // sum weights of all nodes in the community

};

#endif