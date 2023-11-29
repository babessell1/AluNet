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
    Partition(const std::vector<Community>& communities, const Graph& G);
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembershipSearch(int node_index, int new_community_index);
    void updateCommunityMembership(int node_index, int old_community_index, int new_community_index);
    int getNodeCommunity(int node_index);
    void purgeEmptyCommunities(bool renumber);
    void addCommunity(const Community& newCommunity);
    std::unordered_map<int, double> getPartitionWeights(const Graph& G) const;
    double calcQuality(double gamma, const Graph& G) const;
    void makeSingleton(const Graph& G);
    bool inSingleton(int node_index) const;
};

#endif
