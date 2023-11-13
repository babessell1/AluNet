#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"
#include "Partition.h"

/*
################################################################################
############################ PARTITION CLASS ###################################
*/
// Construct a partition based on a set of communities
Partition::Partition(const std::vector<Community>& communities) {
    for (const auto& community : communities) {
        // add the community to the community map
        communityIndexMap.insert({community.communityIndex, community});
        // for each node in the community
        for (int node_index : community.nodeIndices) {
            // add the node to the node map
            nodeCommunityMap.insert({node_index, community.communityIndex});
        }
    }
    // initialize the quality
    quality = 0.0;
}
    
// get indices of communities in the partition (different from node indices!)
std::vector<int> Partition::getCommunityIndices() const {
    // get community indices
    // from communityIndexMap
    std::vector<int> indices;
    for (const auto& entry : communityIndexMap) {
        indices.push_back(entry.first);
    }

    return indices;
}

// flatten the partition, last step of the optimization step in the paper
 void Partition::flattenPartition() {
    // for each community index in the partition
    for (int communityIndex : getCommunityIndices()) {
        // get the community
        auto comm = communityIndexMap.find(communityIndex);
        // flatten the subsets
        std::vector<int> flat_set;
        for (int nidx : comm->second.nodeIndices) {
            flat_set.push_back(nidx);
        }
        // set the flat set as the community in the map
        // store as a vector of vectors!
        comm->second.nodeIndices = flat_set;
    }
}

void Partition::addCommunity(const Community& newCommunity) {
    // Add the new community to the communities vector
    communityIndexMap.insert({newCommunity.communityIndex, newCommunity});
}

// Function to get the community of a given vertex
int Partition::getNodeCommunity(int n_idx) {
    // get the community of the node
    auto node_it = nodeCommunityMap.find(n_idx);
    if (node_it != nodeCommunityMap.end()) {
        return node_it->second;
    } else {
        Rcpp::stop("Node not found in the node community map: " + std::to_string(n_idx));
    }
  }

std::vector<double> Partition::getPartitionWeights(const Graph& G) const {
    // for each cluster, get the sum of weights of all nodes in the community
    std::vector<double> cluster_weights;
    for (int c_idx : getCommunityIndices()) {
        cluster_weights.push_back(communityIndexMap.at(c_idx).getClusterWeight(G));
    }
    return cluster_weights;
}

// move a node from one community to another
void Partition::updateCommunityMembershipSearch(int node_index, int new_community_index) {
    // Find the current community of the node and remove the node from the community
    bool found = false;
    for (auto& entry : communityIndexMap) {
        // if the node is in the community, (throw error if not found)
        auto node_it = std::find(entry.second.nodeIndices.begin(), entry.second.nodeIndices.end(), node_index);
        if (node_it != entry.second.nodeIndices.end()) {
            // remove the node from the community
            entry.second.nodeIndices.erase(node_it);
            found = true;
        }
    }
    if (!found) {
        Rcpp::stop("Node not found in any community: " + std::to_string(node_index));
    }

    // Add the node to the new community
    communityIndexMap.at(new_community_index).nodeIndices.push_back(node_index);
}

void Partition::updateCommunityMembership(int node_index, int old_community_index, int new_community_index) {
    // Find the current community of the node and remove the node from the community
    auto old_comm = communityIndexMap.find(old_community_index);
    if (old_comm != communityIndexMap.end()) {
        // if the node is in the community, remove it
        auto node_it = std::find(old_comm->second.nodeIndices.begin(), old_comm->second.nodeIndices.end(), node_index);
        if (node_it != old_comm->second.nodeIndices.end()) {
            old_comm->second.nodeIndices.erase(node_it);
        } else {
            Rcpp::stop("Node not found in the old community: " + std::to_string(node_index));
        }
    } else {
        Rcpp::stop("Old community not found: " + std::to_string(old_community_index));
    }

    // Add the node to the new community
    communityIndexMap.at(new_community_index).nodeIndices.push_back(node_index);
}

// ################################
// [[test it]]

// purge empty communities
void Partition::purgeEmptyCommunities(bool leave_one) {
    // Remove empty communities
    for (auto it = communityIndexMap.begin(); it != communityIndexMap.end();) {
        if (it->second.size() == 0) {
            it = communityIndexMap.erase(it);
        } else {
            ++it;
        }
    }
    // if leave_one, create one empty community as an option to move nodes to
    if (leave_one) {
        // get the index of the last community
        int last_index = communityIndexMap.begin()->first;
        // create an empty community
        Community empty_community({}, last_index + 1);
        // add the empty community to the partition
        communityIndexMap.insert({last_index + 1, empty_community});
    }

    // renumber the communities to be consecutive
    int new_index = 0;
    for (auto& entry : communityIndexMap) {
        entry.second.communityIndex = new_index;
        // find index of the community to be updated in the nodeCommunityMap and update it
        auto node_it = nodeCommunityMap.find(entry.first);
        if (node_it != nodeCommunityMap.end()) {
            node_it->second = new_index;
        } else {
            Rcpp::stop("Node not found in the node community map: " + std::to_string(entry.first));
        }
        new_index++;
    }
}

double Partition::calcQuality(double gamma, const Graph& G) const {
   /*
    Calculate the quality of the partition 1 / (2 * m) * sum(d(c[i], c[j]) * (a[i][j] - gamma * n[i] * n[j])),
    where:
        a[i][j] is the weight of the edge between nodes i and j
        d(c[i], c[j]) is 1 if c[i] == c[j] and 0 otherwise
        gamma is the resolution parameter
        n[i] is the weight of node i
        m is the total edge weight
        sum is over all pairs of nodes i and j in the graph
    */

    double quality = 0.0;
    // for each node in the graph
    for (int node_idx : G.nodes) {
        // get the neighbors of the node
        for (int neigh_idx : G.getNeighbors(node_idx)) {
            // get the weight of the edge between the node and its neighbor
            double edge_weight = G.edgeWeights.at(node_idx).at(neigh_idx);
            // if the neighbor is in the same community as the node
            if (nodeCommunityMap.at(node_idx) == nodeCommunityMap.at(neigh_idx)) {
                // add the edge weight to the quality
                quality += edge_weight;
            }
        }
    }

    // our graph has no self links, so we do not need to handle them

    std::vector<double> p_weights = getPartitionWeights(G);
    // for each community in the partition, get the weight of the community and remove it from the quality
    for (int c_idx : getCommunityIndices()) {
        // get the weight of the community
        double c_weight = p_weights.at(c_idx);
        // subtract the weight of the community from the quality
        quality -= gamma * c_weight * c_weight;
    }

    quality /= 2 * G.totalEdgeWeight;

    return quality;   
}
