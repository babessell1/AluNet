
#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<int>& nodes, int index)
    : communityIndex(index), nodeIndices(std::move(nodes)) {
}

// get the sum of weights of all edges in the community
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

bool Community::hasNode(int node_index) const {
    // check if the community has a node
    return std::find(nodeIndices.begin(), nodeIndices.end(), node_index) != nodeIndices.end();
}

bool Community::hasEdge(int node_index, int neighbor_index, const Graph& G) const {
    // check if the community has an edge
    return hasNode(node_index) && hasNode(neighbor_index) && G.getWeight(node_index, neighbor_index) > 0.0;
}

void Community::addNode(int nodeIndex) {
    Rcpp::Rcout << "C: Adding node index " << nodeIndex << " to community " << communityIndex << std::endl;
    nodeIndices.push_back(nodeIndex);
}

void Community::removeNode(int nodeIndex) {
    Rcpp::Rcout << "C: Removing node index " << nodeIndex << " from community " << communityIndex << std::endl;
    nodeIndices.erase(std::remove(nodeIndices.begin(), nodeIndices.end(), nodeIndex), nodeIndices.end());
}

bool Community::isEmpty() const {
    return nodeIndices.empty();
}

// get the number of nodes in the community
size_t Community::size() const {
    return nodeIndices.size();
}

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
