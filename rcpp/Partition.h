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
/*
bool isSingleton(const Partition& partition, int node) {
    int communityIndex = partition.nodeCommunityMap.at(node);
    const auto& community = partition.communityIndexMap.at(communityIndex);

    return community.nodeIndices.size() == 1;
}
*/
    Partition(const std::vector<Community>& communities);
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembershipSearch(int node_index, int new_community_index);
    void updateCommunityMembership(int node_index, int old_community_index, int new_community_index);
    int getNodeCommunity(int node_index);
    void purgeEmptyCommunities(bool leave_one);
    void addCommunity(const Community& newCommunity);
    std::vector<double> getPartitionWeights(const Graph& G) const;
    double calcQuality(double gamma, const Graph& G) const;
};

#endif
