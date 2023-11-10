#include <Rcpp.h>
#include "GraphUtils.h"
#include "Leiden.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<std::vector<int>>& nodes, int index) : number_of_nodes(nodes.size()), communityIndex(index), nodeIndices(nodes) {}
/*
################################################################################
############################ PARTITION CLASS ###################################
*/
// Construct a partition based on a set of communities
Partition::Partition(const std::vector<Community>& communities) : communities(communities) {
    // Populate communityIndices based on the communities in the partition
    for (size_t i = 0; i < communities.size(); i++) {
        this->communityIndices.push_back(communities[i].communityIndex);
    }
    // Populate communityIndexMap based on the communities in the partition
    for (size_t i = 0; i < communities.size(); i++) {
        this->communityIndexMap.insert({communities[i].communityIndex, communities[i]});
    }
}

// get indices of communities in the partition (different from node indices!)
std::vector<int> Partition::getCommunityIndices() {
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
    for (int communityIndex : communityIndices) {
        // get the subsets of the community
        std::vector subsets = communityIndexMap[communityIndex].nodeIndices;
        // flatten the subsets
        std::vector<int> flat_set;
        for (std::vector<int> subset : subsets) {
            for (int nodeIndex : subset) {
                flat_set.push_back(nodeIndex);
            }
        }
        // set the flat set as the community in the map
        // store as a vector of vectors!
        communityIndexMap[communityIndex].nodeIndices = {flat_set};
    }
}

/*
size_t Partition::number_of_nodes() {
    // Implement the logic to return the number of nodes in the partition
    return graph.n; // Assuming n is the number of nodes in the graph
}

Graph Partition::get_graph() {
    // Implement the logic to return the graph object associated with the partition
    return graph;
}

void Partition::addCommunity(const Community& newCommunity) {
    // Implement logic to add a new community to the partition
    communities.push_back(newCommunity);
}

void Partition::removeCommunity(int communityIndex) {
    // Implement logic to remove a community from the partition
    // You may want to handle the case where the index is out of bounds
}

void Partition::updateCommunityMembership(int nodeIndex, int newCommunityIndex) {
    // Implement logic to update the community membership of a node
    // You may want to handle the case where the indices are out of bounds
}

*/

/*
################################################################################
############################ OPTIMIZER CLASS ###################################
*/


// Construct an optimizer based on a graph and a partition
Optimizer::Optimizer(const Graph& G, const Partition& P) : G(G), P(P) {}

/*
Partition Optimizer::moveNodesFast() {
    // Implement the Leiden move operation
    Partition updatedPartition; // Placeholder, replace with actual logic
    return updatedPartition;
}

Partition Optimizer::refinePartition() {
    // Implement creating a new partition from the existing partition
    Partition refinedPartition; // Placeholder, replace with actual logic
    return refinedPartition;
}

Partition Optimizer::mergeNodesSubset(const Community& subset) {
    // Implement refining a community
    Partition mergedPartition; // Placeholder, replace with actual logic
    return mergedPartition;
}

Graph Optimizer::aggregateGraph(const Graph& G, const Partition& P) {
    // Implement graph aggregation
    Graph aggregatedGraph; // Placeholder, replace with actual logic
    return aggregatedGraph;
}

double Optimizer::constantPotts(double gamma) const {
    double objective = 0.0;
    // Implement the constant Potts model objective function
    // Placeholder, replace with actual logic
    return objective;
}

*/
void Optimizer::optimize() {
    // Implement the Leiden algorithm iteration here
    // This involve multiple steps, such as moveNodesFast, refinePartition, mergeNodesSubset, etc.
    bool done = false;
    while (!done) {
        //moveNodesFast();
        // set done to true if partition size is equal to number of nodes
        // convert number of nodes (P.communities.size()) to int
        int num_nodes = P.communities.size();
        if (num_nodes == G.n) { 
            done = true;
        }
        //Partition P_refined = refinePartition();
        //aggregateGraph(P_refined);

        // set true no matter what for now
        done = true;

    }
    P.flattenPartition();
}


Partition initializePartition(const Graph& G) {
    std::vector<Community> communities;

    // Assign each node to its own community
    // for each node in the getNodes
    for (int nodeIndex : G.getNodes()) {
        // Construct a community with a single node
        Community community({{nodeIndex}}, nodeIndex);
        communities.push_back(community);
    }

    Partition P(communities);

    return P;
}


// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations) {

    // Create a graph from the R List
    // use listToGraph from GraphUtils.cpp
    Graph G = listToGraph(graphList);
    Partition P = initializePartition(G);

    // Create an Optimizer
    Optimizer optim(G, P);

    for (int i = 0; i < iterations; i++) {
        // Run the Leiden algorithm
        optim.optimize();
    }

    // get the communities from the partition
    std::vector<int> communities;
    for (size_t i = 0; i < P.communities.size(); i++) {
        communities.push_back(P.communities[i].communityIndex);
    }

    // Convert the vector of communities to an R List
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities
    );


    // return R integer 0 to R for now
    return result;


}
