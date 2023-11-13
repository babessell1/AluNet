#ifndef PARTITION_H
#define PARTITION_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Partition {
public:
    // properties
    std::unordered_map<int, Community> communityIndexMap;
    std::unordered_map<int, int> nodeCommunityMap;
    double quality;

    // methods
    Partition(const std::vector<Community>& communities);
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembershipSearch(int node_index, int new_community_index);
    void updateCommunityMembership(int node_index, int old_community_index, int new_community_index);
    int getNodeCommunity(int node_index);
    void purgeEmptyCommunities(bool leave_one);
    //size_t number_of_nodes(); // could define a inline function
    void addCommunity(const Community& newCommunity);
    std::vector<double> getPartitionWeights(const Graph& G) const;
    double calcQuality(double gamma, const Graph& G) const;
    //void removeCommunity(int communityIndex); // should define this unction
    //std::vector<int> getCommunityIndices(); // should define this unction
    //int getCommunityIndex(int node_index); // should define this unction
    //double quality(double resolution_parameter); // should define this unction
    //void renumber_communities(); // should define this unction
    //void renumber_communities(std::vector<size_t> fixed_nodes, std::vector<size_t> fixed_membership); // should define this unction
    //size_t membership(size_t vertex); function that returns vertex is in which community
};

#endif