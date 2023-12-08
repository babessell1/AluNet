#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"
#include "Partition.h"

/*
################################################################################
############################ PARTITION CLASS ###################################
*/

/**
 * @brief Construct a new Partition:: Partition object
 * @param communities vector of communities
 * @param G graph object
 * @note initialize the quality to 0.0
 * @note initialize the community index map and node community map
 * @note the community index map is a map of community indices to communities
 * @note the node community map is a map of node indices to community indices
**/
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

/**
 * @brief Construct a new blank Partition:: Partition object
 * @note initialize the quality to 0.0
 * @note initialize the community index map and node community map to empty maps
*/
Partition::Partition() {
    quality = 0.0;
    communityIndexMap = {};
    nodeCommunityMap = {};
}
    
/**
 * @brief get a vector of community indices in the partition
 * @return std::vector<int> : vector of community indices
 * @note might want to stop using this method in favor of new getters/setters
*/
std::vector<int> Partition::getCommunityIndices() const {
    // get community indices
    // from communityIndexMap
    std::vector<int> indices;
    for (const auto& entry : communityIndexMap) {
        indices.push_back(entry.first);
    }

    return indices;
}

/**
 * @brief Flatten the partition, i.e. set the partition to have one community with all nodes
 * @return void
 * @note last step of the optimization step in the paper
**/
void Partition::flattenPartition() {

    // for each community index in the partition
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
    }
    // set the flat set partition as one community with all nodes
    setCommunityIndexMap({{0, Community(flat_set, 0)}});
    // set the node community map
    setNodeCommunityMap(new_node_community_map);
}

/**
 * @brief Add a community to the partition
 * @param newCommunity community to add
 * @return void
 * @note add the community to the community index map
**/ 
void Partition::addCommunity(const Community& newCommunity) {

    // Add the new community to the communities vector
    communityIndexMap.insert({newCommunity.getCommunityIndex(), newCommunity});
}

/**
 * @brief Get the partition weights, i.e. the sum of weights of all nodes in each community
 * @param G graph object
 * @return std::unordered_map<int, double> : map of community indices to community weights
 * @note for each community, get the sum of weights of all nodes in the community
**/
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

/**
 * @brief Add a node to a community in the partition
 * @param node_index node index
 * @param community_index community index
 * @return void
**/ 
void Partition::addNodeToCommunity(int node_index, int community_index) {

    // add a node to a community
    communityIndexMap.at(community_index).addNode(node_index);
    nodeCommunityMap.at(node_index) = community_index;
}

/**
 * @brief Remove a node from a community in the partition
 * @param node_index node index
 * @param community_index community index
 * @return void
 * @note set the node community map to val that will throw an error if used in vector
**/
void Partition::removeNodeFromCommunity(int node_index, int community_index) {

    // remove a node from a community
    communityIndexMap.at(community_index).removeNode(node_index);

    // set the node community map to val that will throw an error if used in vector
    nodeCommunityMap.at(node_index) = -1;
}

/**
 * @brief Update the community membership of a node in the partition
 * @param node_index node index
 * @param old_community_index old community index
 * @param new_community_index new community index
 * @return void
 * @note removes the node from the old community and add it to the new community
 * @throw Rcpp::stop if the node is not found in the old community
 * @throw Rcpp::stop if the node is found in a community other than the new community
 * @throw Rcpp::stop if the old community is not found
**/
void Partition::updateCommunityMembership(int node_index, int old_community_index, int new_community_index) {
    auto old_comm = communityIndexMap.find(old_community_index);

    // check if the old community is in the community index map,
    if (old_comm != communityIndexMap.end()) {

        // if the node is in the community, remove it using the removeNodeFromCommunity method
        if (old_comm->second.hasNode(node_index)) {
            Rcpp::Rcout << "Node found in the old community: " << node_index << std::endl;
            removeNodeFromCommunity(node_index, old_community_index);

            // add the node to the new community
            addNodeToCommunity(node_index, new_community_index);

            // check if any community other than the new community has the node
            for (auto& entry : communityIndexMap) {
                int c_idx = entry.first;
                if (c_idx != new_community_index) {
                    if (entry.second.hasNode(node_index)) {
                        Rcpp::stop("Node not moved properly! " + std::to_string(c_idx));
                    }
                }
            }

        } else {
            Rcpp::stop("Node not found in the old community: " + std::to_string(node_index));
        }
    } else {
        Rcpp::stop("Old community not found: " + std::to_string(old_community_index));
    }
}

/**
 * @brief Purge empty communities from the partition, optionally renumbering the communities to be consecutive
 * @param renumber boolean to renumber the communities
 * @return void
 * @note remove empty communities from the community index map
 * @note renumber might be unstable, not tested in a long time
**/
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


/**
 * @brief Calculate the quality of the partition using the Constant Potts Model
 * @param gamma resolution parameter
 * @param G graph object
 * @param print boolean to print debug info
 * @return double : quality of the partition
 * @note CPM is: 1 / (2 * m) * sum(d(c[i], c[j]) * (a[i][j] - gamma * n[i] * n[j])),
 *       where:
 *          a[i][j] is the weight of the edge between nodes i and j
 *          d(c[i], c[j]) is 1 if c[i] == c[j] and 0 otherwise
 *          gamma is the resolution parameter
 *          n[i] is the weight of node i
 *          m is the total edge weight
 *          sum is over all pairs of nodes i and j in the graph
**/
double Partition::calcQuality(double gamma, const Graph& G, bool print) const {
    double quality = 0.0;
    for (int node_idx : G.getNodes()) {

        // get the neighbors of the node
        for (int neigh_idx : G.getNeighbors(node_idx)) {
            double edge_weight = G.getEdgeWeights().at(node_idx).at(neigh_idx);  
            if (print) {
                Rcpp::Rcout << "Node community: " << nodeCommunityMap.at(node_idx) << " Neighbor community: " << nodeCommunityMap.at(neigh_idx) << " has edge weight: " << edge_weight << std::endl;
            }
            // if the neighbor is in the same community as the node
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
    quality /= 2*G.getTotalEdgeWeight();

    return quality;   
}

/**
 * @brief Partition the graph such that each node is in its own community
 * @param G graph object
 * @return void
**/
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

/**
 * @brief Check if a node is in a singleton community (i.e. the only member of its community)
 * @param node_idx node index
 * @return bool : true if the node is in a singleton community, false otherwise
 * @throw Rcpp::stop if the node is not found in the community index map
 * @note does not check if the node exists in the graph
 * @note does not check if the node has a weight
**/
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




