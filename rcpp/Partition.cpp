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
        //communityIndexMap.insert({community.communityIndex, community});
        communityIndexMap.emplace(community.getCommunityIndex(), community);
        // for each node in the community
        for (int node_index : community.getNodeIndices()) {
            // add the node to the node map
            nodeCommunityMap.insert({node_index, community.getCommunityIndex()});
        }
    }
    // initialize the quality
    quality = 0.0;
}

// blank constructor
Partition::Partition() {
    quality = 0.0;
    communityIndexMap = {};
    nodeCommunityMap = {};
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
    //for (int communityIndex : getCommunityIndices()) { slow
    std::vector<int> flat_set;
    std::unordered_map<int, int> new_node_community_map;
    for (auto& entry : communityIndexMap) {  // for each community in the partition
        int community_index = entry.first;
        // get the community
        auto comm = communityIndexMap.find(community_index);
        // flatten the subsets
        std::vector<int> flat_set;
        for (int nidx : comm->second.getNodeIndices()) {  // for each node in the community
            flat_set.push_back(nidx);  // add the node to the flat set
            new_node_community_map.insert({nidx, 0});  // update the node community map
        }
        // remove the community from the community index map
        // set the flat set as the community in the map
        // store as a vector of vectors!
        //comm->second.nodeIndices = std::move(flat_set);
    }
    // set the flat set partition as one community with all nodes
    setCommunityIndexMap({{0, Community(flat_set, 0)}});
    // set the node community map
    setNodeCommunityMap(new_node_community_map);
}

void Partition::addCommunity(const Community& newCommunity) {
    // Add the new community to the communities vector
    communityIndexMap.insert({newCommunity.getCommunityIndex(), newCommunity});
}

// Function to get the community of a given vertex
//int Partition::getNodeCommunity(int n_idx) const {
    // get the community of the node
//    auto node_it = nodeCommunityMap.find(n_idx);
//    if (node_it != nodeCommunityMap.end()) {
//        return node_it->second;
//    } else {
//        Rcpp::stop("Node not found in the node community map: " + std::to_string(n_idx));
//    }
//  }

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

void Partition::addNodeToCommunity(int node_index, int community_index) {
    // add a node to a community
    //Rcpp::Rcout << "P: Adding node index " << node_index << " to community " << community_index << std::endl;
    communityIndexMap.at(community_index).addNode(node_index);
    nodeCommunityMap.at(node_index) = community_index;
}

void Partition::removeNodeFromCommunity(int node_index, int community_index) {
    // remove a node from a community
    //Rcpp::Rcout << "P: Removing node index " << node_index << " from community " << community_index << std::endl;
    communityIndexMap.at(community_index).removeNode(node_index);
    nodeCommunityMap.erase(node_index);
}

void Partition::updateCommunityMembership(int node_index, int old_community_index, int new_community_index) {
    // Find the current community of the node and remove the node from the community

    auto old_comm = communityIndexMap.find(old_community_index);

    // check if the old community is in the community index map,
    if (old_comm != communityIndexMap.end()) {
        //Rcpp::Rcout << "Old community found: " << old_community_index << std::endl;
        // if the node is in the community, remove it using the removeNodeFromCommunity method
        if (old_comm->second.hasNode(node_index)) {
            Rcpp::Rcout << "Node found in the old community: " << node_index << std::endl;
            removeNodeFromCommunity(node_index, old_community_index);

        } else {
            Rcpp::stop("Node not found in the old community: " + std::to_string(node_index));
        } // I do not think that it is necessary but we can add the erase the old community to delete the old community
    } else {
        Rcpp::stop("Old community not found: " + std::to_string(old_community_index));
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
            //comm.getCommunityIndex() = old_index_to_new_index[entry.first];
            comm.setCommunityIndex(old_index_to_new_index[entry.first]);

            // add the community to the new map
            new_community_index_map.insert({comm.getCommunityIndex(), comm});
        }
        // update the community index map
        communityIndexMap = new_community_index_map;
    }
}

double Partition::calcQuality(double gamma, const Graph& G, bool print) const {
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

    for (int node_idx : G.getNodes()) {

        // get the neighbors of the node
        //if (print) {
        //    Rcpp::Rcout << "Node index: " << node_idx << std::endl;
        //}
        for (int neigh_idx : G.getNeighbors(node_idx)) {
            //Rcpp::Rcout << "----Neighbor index: " << neigh_idx << std::endl;
            // get the weight of the edge between the node and its neighbor

            //Rcpp::Rcout << "----Get weight: " << std::endl;
            double edge_weight = G.getEdgeWeights().at(node_idx).at(neigh_idx);  
            // if the neighbor is in the same community as the node
            //Rcpp::Rcout << "----Get communities: " << std::endl;
            if (print) {
                Rcpp::Rcout << "Node community: " << nodeCommunityMap.at(node_idx) << " Neighbor community: " << nodeCommunityMap.at(neigh_idx) << " has edge weight: " << edge_weight << std::endl;
            }
            if (nodeCommunityMap.at(node_idx) == nodeCommunityMap.at(neigh_idx)) {
                // add the edge weight to the quality
                quality += edge_weight;
            }
        }
    }
    if (print) {
        Rcpp::Rcout << "Get partition weights: " << std::endl;
    }
    // our graph has no self links, so we do not need to handle them for now
    std::unordered_map<int, double> p_weights = getPartitionWeights(G);

    if (print) {
        Rcpp::Rcout << "Get community indices: " << std::endl;
    }
    // for each community in the partition, get the weight of the community and remove it from the quality
    for (int c_idx : getCommunityIndices()) {

        // get the weight of the community
        double c_weight = p_weights.at(c_idx);

        // subtract the weight of the community from the quality
        quality -= gamma * c_weight * c_weight;

    }

    // print total edge weight
    //Rcpp::Rcout << "Total edge weight: " << G.totalEdgeWeight << std::endl;

    quality /= 2*G.getTotalEdgeWeight();

    // print the quality
    //Rcpp::Rcout << "Quality: " << quality << std::endl;

    return quality;   
}

void Partition::makeSingleton(const Graph& G) {
    // initialize the community index map
    std::unordered_map<int, Community> single_community_index_map;
    std::unordered_map<int, int> single_node_community_map;
    // for each node in the graph
    for (int node_idx : G.getNodes()) {
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




