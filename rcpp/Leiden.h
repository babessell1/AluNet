#ifndef LEIDEN_H
#define LEIDEN_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    // properties
    int communityIndex; 
    std::vector<std::vector<int>> nodeIndices;  // node indices nested in their subset group
    
    // methods
    Community(const std::vector<std::vector<int>>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(Graph& G); // sum of weights of all edges in the community
    size_t size(); // number of nodes in the community
};

class Partition {
public:
    // properties
    std::vector<int> communityIndices;
    std::unordered_map<int, Community> communityIndexMap;
    std::vector<Community> communities;

    // methods
    Partition(const std::vector<Community>& communities);
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembership(int nodeIndex, int newCommunityIndex);
    //size_t number_of_nodes();
    void addCommunity(const Community& newCommunity);
    //void removeCommunity(int communityIndex);
    //void updateCommunityMembership(int nodeIndex, int newCommunityIndex);
    //std::vector<int> getCommunityIndices();
    //int getCommunityIndex(int nodeIndex);
    //double quality(double resolution_parameter);
    //void renumber_communities();
    //void renumber_communities(std::vector<size_t> fixed_nodes, std::vector<size_t> fixed_membership);
};

class Optimizer {
public:
    // properties
    Graph& G;
    Partition& P;

    // methods
    Optimizer(Graph& G, Partition& P);
    void optimize();
    //moveNodesFast();
    //Partition refinePartition() const;
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    //double constantPotts(double gamma);
};

Partition initializePartition();

#endif
