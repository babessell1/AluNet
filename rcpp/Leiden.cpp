#include <Rcpp.h>
#include "GraphUtils.h"
#include "Leiden.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<std::vector<int>>& nodes, int index)
    : number_of_nodes(nodes.size()), communityIndex(index), nodeIndices(nodes) {}
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
        // get the community
        auto comm = communityIndexMap.find(communityIndex);
        // get the subsets of the community
        std::vector<std::vector<int>>& subsets = comm->second.nodeIndices;
        // flatten the subsets
        std::vector<int> flat_set;
        for (std::vector<int>& subset : subsets) {
            for (int nodeIndex : subset) {
                flat_set.push_back(nodeIndex);
            }
        }
        // set the flat set as the community in the map
        // store as a vector of vectors!
        comm->second.nodeIndices = {flat_set};
    }
}


void Partition::addCommunity(const Community& newCommunity) {
    // Add the new community to the communities vector
    communities.push_back(newCommunity);
    // Update the communityIndexMap with the new community's index
    communityIndexMap[newCommunity.communityIndex] = communities.size() - 1;
}

void Partition::updateCommunityMembership(int nodeIndex, int newCommunityIndex) {
    // Check if the new community index is valid
    if (communityIndexMap.find(newCommunityIndex) == communityIndexMap.end()) {
        Rcpp::stop("Invalid community index: " + std::to_string(newCommunityIndex));
    }

    // Find the current community of the node
    for (auto& community : communities) {
        auto& nodeIndices = community.nodeIndices;
        auto it = std::find(nodeIndices.begin(), nodeIndices.end(), nodeIndex);
        if (it != nodeIndices.end()) {
            // Remove the node from its current community
            nodeIndices.erase(it);
            break;
        }
    }

    // Add the node to the new community
    communities[communityIndexMap[newCommunityIndex]].nodeIndices.push_back(nodeIndex);
}

double Partition::quality(double resolution_parameter) {
    double mod = 0.0;
    for (size_t c = 0; c < this->number_of_communities(); c++) {
        double csize = this->community_size(c);
        double w = this->total_weight_in_community(c);
        double comm_possible_edges = this->graph->possible_edges(csize);
        mod += w - resolution_parameter * comm_possible_edges;
    }
    return (2.0 - this->graph->is_directed()) * mod;
}

double Partition::total_weight_in_community(size_t communityIndex) {
    double total_weight = 0.0;
    // Implement logic to calculate total weight of edges within the community
    // This will depend on how you've structured your Community and Graph classes
    return total_weight;
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
    
    return H;
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
    int communityIndex = 0;
    for (int nodeIndex : G.getNodes()) {
        // Construct a community with a single node
        Community community({{nodeIndex}}, communityIndex);
        communities.push_back(community);
        communityIndex++;
    }

    Partition P(communities);

    return P;
}


// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations) {

    // Create a graph from the R List
    // use listToGraph from GraphUtils.cpp
    Graph G = listToGraph(graphList);
    // print the number of nodes 
    Rcpp::Rcout << "Number of nodes: " << G.n << std::endl;
    Partition P = initializePartition(G);
    // print the number of communities
    Rcpp::Rcout << "Number of communities: " << P.communities.size() << std::endl;

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
