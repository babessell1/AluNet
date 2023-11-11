#ifndef LEIDEN_H
#define LEIDEN_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    // properties
    size_t number_of_nodes;
    int communityIndex; 
    std::vector<std::vector<int>> nodeIndices;  // node indices nested in their subset group
    
    
    // methods
    Community(const std::vector<std::vector<int>>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
};

class Partition {
public:
    // properties
    std::vector<int> communityIndices;
    std::unordered_map<int, Community> communityIndexMap;
    std::vector<Community> communities;

    // methods
    Partition(const std::vector<Community>& communities);
    std::vector<int> getCommunityIndices();
    void flattenPartition();
    void updateCommunityMembership(int nodeIndex, int newCommunityIndex);
    //size_t number_of_nodes();
    void addCommunity(const Community& newCommunity);
    //void removeCommunity(int communityIndex);
    //void updateCommunityMembership(int nodeIndex, int newCommunityIndex);
    //std::vector<int> getCommunityIndices();
    //int getCommunityIndex(int nodeIndex);
    inline size_t number_of_nodes_in_community(size_t communityIndex) {
        if (communityIndex >= communities.size()) {
            Rcpp::stop("Invalid community index: " + std::to_string(communityIndex));
        }
        return communities[communityIndex].size();
    }
    double quality(double resolution_parameter);
};

class Optimizer {
public:
    // properties
    Graph G;
    Partition P;

    // methods
    Optimizer(const Graph& G, const Partition& P);
    void optimize();
    //moveNodesFast();
    //Partition refinePartition() const;
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    //double constantPotts(double gamma);
};

Partition initializePartition();

#endif
