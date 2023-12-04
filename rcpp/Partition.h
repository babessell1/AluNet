#ifndef PARTITION_H
#define PARTITION_H

#include <Rcpp.h>
#include "GraphUtils.h"

class Partition {
private:
    // attributes
    std::unordered_map<int, Community> communityIndexMap;
    std::unordered_map<int, int> nodeCommunityMap;  // node_idx, community_idx
    double quality;

public:
    // getters
    std::unordered_map<int, Community> getCommunityIndexMap() const {
        return communityIndexMap;
    }
    
    std::unordered_map<int, int> getNodeCommunityMap() const {
        return nodeCommunityMap;
    }

    double getQuality() const {
        return quality;
    }

    // setters
    void setCommunityIndexMap(std::unordered_map<int, Community> communityIndexMap) {
        communityIndexMap = communityIndexMap;
    }
    
    void setNodeCommunityMap(std::unordered_map<int, int> nodeCommunityMap) {
        nodeCommunityMap = nodeCommunityMap;
    }

    void setNodeCommunity(int node_index, int community_index) {
        nodeCommunityMap[node_index] = community_index;
    }

    void setQuality(double quality) {
        quality = quality;
    }

    // methods
    Partition(const std::vector<Community>& communities, const Graph& G);
    Partition();
    std::vector<int> getCommunityIndices() const;
    void flattenPartition();
    void updateCommunityMembershipSearch(int node_index, int new_community_index);
    void updateCommunityMembership(int node_index, int old_community_index, int new_community_index);
    int getNodeCommunity(int node_index) const;
    void purgeEmptyCommunities(bool renumber);
    void addCommunity(const Community& newCommunity);
    std::unordered_map<int, double> getPartitionWeights(const Graph& G) const;
    double calcQuality(double gamma, const Graph& G, bool quality) const;
    void makeSingleton(const Graph& G);
    bool inSingleton(int node_index) const;
};

#endif
