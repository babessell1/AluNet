
#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/

/**
 * @brief Construct a new Community:: Community object
 * @param nodes vector of node indices
 * @param index community index
**/
Community::Community(const std::vector<int>& nodes, int index)
    : communityIndex(index), nodeIndices(std::move(nodes)) {
}


/**
 * @brief Construct a new Community:: Community object
 * @note default constructor
 * @note does not initialize any properties
**/
Community::Community() {
}

/**
 * @brief Calculate the sum of weights of all edges in the community
 * @param G graph object
 * @return double : sum of weights of all edges in the community
*/
double Community::aggregateWeights(const Graph& G) const {
    double weight_sum = 0.0;
    for (int node_index : nodeIndices) {
        // get the neighbors of the node
        std::vector<int> neighbors = G.getNeighbors(node_index);
        // for each neighbor of the node
        for (int& neighbor_index : neighbors) {
            // if the neighbor is in the community, add the weight of the edge to the sum
            weight_sum += G.getWeight(node_index, neighbor_index);
        }
    }
    if (G.getIsDirected()) {
        return weight_sum;
    }
    else {
        return weight_sum / 2.0;
    };
}

/**
 * @brief Count the number of possible edges in the community
 * @param G graph object
 * @return int : number of possible edges in the community
 * @note does not count self loops
 * @note does not count zero weight edges
 * @note for undirected graphs, divide the sum by 2
**/
int Community::countPossibleEdges(const Graph& G) const {
    // get the number of possible edges in the community, assuming undirected graph
    // do not count zero weight edges
    int n_possible_edges = 0;
    for (const int& node_index : nodeIndices) {
        // get the neighbors of the node
        std::vector<int> neighbors = G.getNeighbors(node_index);
        // count possible edges
        for (const int& neighbor_index : neighbors) {
            // if the neighbor is in the community, add the weight of the edge to the sum
            if (G.getWeight(node_index, neighbor_index) > 0.0) {
                n_possible_edges++;
            }
        }
    }
    if (G.getIsDirected()) {
        return n_possible_edges;
    }
    else {
        return n_possible_edges / 2;
    };
}

/**
 * @brief Check if the community has a node
 * @param node_index node index
 * @return bool : true if the community has the node, false otherwise
**/
bool Community::hasNode(int node_index) const {
    // check if the community has a node
    return std::find(nodeIndices.begin(), nodeIndices.end(), node_index) != nodeIndices.end();
}

/**
 * @brief Check if the community has an edge
 * @param node_index node index
 * @param neighbor_index neighbor index
 * @param G graph object
 * @return bool : true if the community has the edge, false otherwise
 * @note does not count zero weight edges
**/
bool Community::hasEdge(int node_index, int neighbor_index, const Graph& G) const {
    // check if the community has an edge
    return hasNode(node_index) && hasNode(neighbor_index) && G.getWeight(node_index, neighbor_index) > 0.0;
}

/**
 * @brief Add a node to the community
 * @param node_index node index
 * @return void
 * @note does not check if the node is already in the community
 * @note does not check if the node exists in the graph
 * @note does not check if the node has a weight
**/ 
void Community::addNode(int node_index) {
    nodeIndices.push_back(node_index);
}

/**
 * @brief Remove a node from the community
 * @param node_index node index
 * @return void
 * @note does not check if the node is already in the community
 * @note does not check if the node exists in the graph
 * @note does not check if the node has a weight
 * @note does not check if the community is empty after removing the node
**/
void Community::removeNode(int node_index) {
    nodeIndices.erase(std::remove(nodeIndices.begin(), nodeIndices.end(), node_index), nodeIndices.end());
}

/**
 * @brief Check if the community is empty
 * @return bool : true if the community is empty, false otherwise
**/ 
bool Community::isEmpty() const {
    return nodeIndices.empty();
}

/**
 * @brief Get the size of the community
 * @return size_t : size of the community
 * @note will include nodes with zero weight
**/
size_t Community::size() const {
    return nodeIndices.size();
}

/**
 * @brief Get the sum of weights of all nodes in the community
 * @param G graph object
 * @return double : sum of weights of all nodes in the community
 * @throw Rcpp::Stop if node index not found by getNodeWeights
**/
double Community::getClusterWeight(const Graph& G) const {
    // get the sum of weights of all nodes in the community
    double weight_sum = 0.0;
    std::unordered_map<int, double> node_weights = G.getNodeWeights();
    std::unordered_map<std::string, int> node_index_map = G.getNodeIndexMap();
    for (const auto& entry : node_index_map) {
        int node_index = entry.second;
        if (G.getNodeWeights().find(node_index) == G.getNodeWeights().end()) {
            Rcpp::stop("Node index not found in nodeWeights");
        }
        double weight_add = G.getNodeWeights().at(node_index);  // Using 'at' for bounds checking
        weight_sum += weight_add;
    }
    
    return weight_sum;
}
