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
bool Partition::isSingleton(int node) const {
    int communityIndex = nodeCommunityMap.at(node); // Get the community index of the node
    const auto& community = communityIndexMap.at(communityIndex); // Get the Community object

    // Check if the size of the community is 1, indicating a singleton
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
