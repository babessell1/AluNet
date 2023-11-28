
#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<int>& nodes, int index)
    : communityIndex(index), nodeIndices(nodes) {
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
    return weight_sum;
}

int Community::countPossibleEdges(const Graph& G) const {
    // get the number of possible edges in the community, assuming undirected graph
    // do not count zero weight edges
    int n_possible_edges = 0;
    for (int node_index : nodeIndices) {
        // get the neighbors of the node
        std::vector<int> neighbors = G.getNeighbors(node_index);
        // count possible edges
        for (int& neighbor_index : neighbors) {
            // if the neighbor is in the community, add the weight of the edge to the sum
            if (G.getWeight(node_index, neighbor_index) > 0.0) {
                n_possible_edges++;
            }
        }
    }
    return n_possible_edges;
}


// get the number of nodes in the community
size_t Community::size() const {
    return nodeIndices.size();
}

double Community::getClusterWeight(const Graph& G) const {
    // get the sum of weights of all nodes in the community
    double weight_sum = 0.0;
    for (int node_index : nodeIndices) {
        double weight_add = G.nodeWeights.at(node_index);  // Using 'at' for bounds checking
        weight_sum += weight_add;
    }
    return weight_sum;
}
    