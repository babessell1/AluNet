#ifndef LEIDEN_H
#define LEIDEN_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Community {
public:
    // properties
    int communityIndex; 
    std::vector<int> nodeIndices;  // node indices
    //int nodeCount;  // number of nodes in the community
    
    // methods
    Community(const std::vector<int>& nodes, int index); // remember to convert input to vector in a vector when constructing!!
    double aggregateWeights(Graph& G); // sum weights of all edges in the community and count number of nodes
    size_t size(); // number of nodes in the community
};

class Partition {
public:
    // properties
    std::unordered_map<int, Community> communityIndexMap;
    std::unordered_map<int, int> nodeCommunityMap;

    // methods
    Partition(const std::vector<Community>& communities);
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembershipSearch(int node_index, int new_community_index);
    void updateCommunityMembership(int node_index, int old_community_index, int new_community_index);
    void purgeEmptyCommunities();
    //size_t number_of_nodes(); // could define a inline function
    void addCommunity(const Community& newCommunity);
    size_t get_community_of_vertex(size_t vertex);
    double diff_move(size_t vertex, size_t new_community);
    //void removeCommunity(int communityIndex); // should define this unction
    //std::vector<int> getCommunityIndices(); // should define this unction
    //int getCommunityIndex(int node_index); // should define this unction
    //double quality(double resolution_parameter); // should define this unction
    //void renumber_communities(); // should define this unction
    //void renumber_communities(std::vector<size_t> fixed_nodes, std::vector<size_t> fixed_membership); // should define this unction
    //size_t membership(size_t vertex); function that returns vertex is in which community
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
