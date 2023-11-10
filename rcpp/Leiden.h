#ifndef LEIDEN_H
#define LEIDEN_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    // properties
    size_t number_of_nodes;
    std::vector<int> nodeIndices;  // Store node indices associated with the community
    int communityIndex;  // Index for the community
    
    // methods
    Community(const std::vector<int>& nodes, int index);
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
    //size_t number_of_nodes();
    //void addCommunity(const Community& newCommunity);
    //void removeCommunity(int communityIndex);
    //void updateCommunityMembership(int nodeIndex, int newCommunityIndex);
    //std::vector<int> getCommunityIndices();
    //int getCommunityIndex(int nodeIndex);
};

class Optimizer {
public:
    // properties
    Graph G;
    Partition P;

    // methods
    Optimizer(const Graph& G, const Partition& P);
    //Partition optimize();
    //Partition moveNodesFast();
    //Partition refinePartition();
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    //double constantPotts(double gamma);
};

Partition initializePartition();

#endif