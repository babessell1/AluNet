#ifndef LEIDEN_H
#define LEIDEN_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    // properties
    int communityIndex; 
    std::vector<int> nodeIndices;  // node indices nested in their subset group
    //int nodeCount;  // number of nodes in the community
    
    // methods
    Community(const std::vector<int>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(Graph& G); // sum weights of all edges in the community and count number of nodes
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
size_t get_community_of_vertex(size_t vertex);
double diff_move(size_t vertex, size_t new_community);
    //void removeCommunity(int communityIndex);
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
    double gamma;

    // methods
    Optimizer(Graph& G, Partition& P, double gamma);
    void optimize();
    //bool moveNodesFast();
    //Partition refinePartition() const;
    //Partition mergeNodesSubset(const Community& subset);
    //Graph aggregateGraph(const Graph& G, const Partition& P);
    //double constantPotts(double gamma);
};

Partition initializePartition();

#endif
