#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"
#include "Partition.h"

/*
################################################################################
############################ PARTITION CLASS ###################################
*/
// Construct a partition based on a set of communities
Partition::Partition(const std::vector<Community>& communities, const Graph& G) {
    for (const auto& community : communities) {
        // add the community to the community map
        communityIndexMap.insert({community.communityIndex, community});
        // for each node in the community
        for (int node_index : community.nodeIndices) {
            // add the node to the node map
            nodeCommunityMap.insert({node_index, community.communityIndex});
            // add the first community assignment to the node
            communityAssignments.insert({G.getNodeName(node_index), {community.communityIndex}});
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

std::unordered_map<int, double> Partition::getPartitionWeights(const Graph& G) const {
    // for each cluster, get the sum of weights of all nodes in the community
    std::unordered_map<int, double> cluster_weights;
    for (int c_idx : getCommunityIndices()) {
        // get the weight of the community
        double c_weight = communityIndexMap.at(c_idx).getClusterWeight(G);
        // add the weight of the community to the map
        cluster_weights.insert({c_idx, c_weight});
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

    // update the node community map
    nodeCommunityMap.at(node_index) = new_community_index;
}

void Partition::updateCommunityMembership(int node_index, int old_community_index, int new_community_index) {
    // Find the current community of the node and remove the node from the community
    auto old_comm = communityIndexMap.find(old_community_index);
    if (old_comm != communityIndexMap.end()) {
        // if the node is in the community, remove it
        auto node_it = std::find(old_comm->second.nodeIndices.begin(), old_comm->second.nodeIndices.end(), node_index);
        if (node_it != old_comm->second.nodeIndices.end()) {
            old_comm->second.nodeIndices.erase(node_it);
        }/* else {
            Rcpp::stop("Node not found in the old community: " + std::to_string(node_index));
        }*/ // I do not think that it is necessary but we can add the erase the old community to delete the old community
    } else {
        Rcpp::stop("Old community not found: " + std::to_string(old_community_index));
    }

    // Add the node to the new community
    communityIndexMap.at(new_community_index).nodeIndices.push_back(node_index);

    // update the node community map
    nodeCommunityMap.at(node_index) = new_community_index;
     // check and handle empty old community
    if (old_comm->second.nodeIndices.empty()) {
        // Remove or flag the empty community as needed
        communityIndexMap.erase(old_comm);
    }
}

// purge empty communities
void Partition::purgeEmptyCommunities(bool renumber) {
    // Remove empty communities
    for (auto it = communityIndexMap.begin(); it != communityIndexMap.end();) {
        // causing segfaults
        // print the community index and size
        if (it->second.size() == 0) {
            it = communityIndexMap.erase(it);
        } else {
            ++it;
        }
    }

    if (renumber) {
        // renumber the communities to be consecutive
        int new_index = 0;
        std::unordered_map<int, Community> new_community_index_map;
        std::unordered_map<int, int> old_index_to_new_index; // map old community index to new community index
        for (auto& entry : communityIndexMap) {
            // map the old index to the new index (old : new)
            old_index_to_new_index[entry.first] = new_index;
            new_index++;
        }
        // use the map to update community indices in the communityIndexMap
        for (auto& entry : communityIndexMap) {
            // get the community
            auto comm = entry.second;
            // update the community index
            comm.communityIndex = old_index_to_new_index[entry.first];
            // add the community to the new map
            new_community_index_map.insert({comm.communityIndex, comm});
        }
        // update the community index map
        communityIndexMap = new_community_index_map;
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
    // our graph has no self links, so we do not need to handle them for now
    std::unordered_map<int, double> p_weights = getPartitionWeights(G);

    // for each community in the partition, get the weight of the community and remove it from the quality
    for (int c_idx : getCommunityIndices()) {

        // get the weight of the community
        double c_weight = p_weights.at(c_idx);

        // subtract the weight of the community from the quality
        quality -= gamma * c_weight * c_weight;

    }

    // print total edge weight
    //Rcpp::Rcout << "Total edge weight: " << G.totalEdgeWeight << std::endl;

    quality /= 2 * G.totalEdgeWeight;

    // print the quality
    //Rcpp::Rcout << "Quality: " << quality << std::endl;

    return quality;   
}

void Partition::updateCommunityAssignments(const Graph& G) {
    for (auto& entry : nodeCommunityMap) {
        std::string node_name = G.getNodeName(entry.first);
        int community_index = entry.second;

        auto it = communityAssignments.find(node_name);
        if (it == communityAssignments.end()) {
            // If not found, insert a new entry with an empty vector
            it = communityAssignments.insert({node_name, {}}).first;
            // Optionally reserve space if you have an estimate of the size
            // it->second.reserve(estimated_size);
        }
        it->second.push_back(community_index);
    }
}
void Partition::makeSingleton(const Graph& G) {
    // initialize the community index map
    std::unordered_map<int, Community> single_community_index_map;
    std::unordered_map<int, int> single_node_community_map;
    // for each node in the graph
    for (int node_idx : G.nodes) {
        // create a community for the node
        Community comm({node_idx}, node_idx);
        // add the community to the community index map
        single_community_index_map.insert({node_idx, comm});
        // add the node to the node community map
        single_node_community_map.insert({node_idx, node_idx});
    }
    // update the community index map
    communityIndexMap = single_community_index_map;
    // update the node community map
    nodeCommunityMap = single_node_community_map;
}

bool Partition::inSingleton(const int node_idx) const {
    // if node is the only member of its community, return true
    bool is_single = false;
    // check if the node is in the community index map
    auto comm_it = communityIndexMap.find(node_idx);
    if (comm_it != communityIndexMap.end()) {
        // if the node is in the community index map, check if it is the only member of the community
        if (comm_it->second.size() == 1) {
            is_single = true;
        }
    } else {
        Rcpp::stop("Node not found in the community index map: " + std::to_string(node_idx));
    }
    return is_single;
}




