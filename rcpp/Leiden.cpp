#include <Rcpp.h>
#include "GraphUtils.h"
#include "Leiden.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<std::vector<int>>& nodes, int index)
    : communityIndex(index), nodeIndices(nodes) {}

// get the sum of weights of all edges in the community
double Community::aggregateWeights(Graph& G) {
    double weight_sum = 0.0;
    // for each subset of nodes in the community
    for (std::vector<int>& subset : nodeIndices) {
        // for each node in the subset
        for (int& nodeIndex : subset) {
            // get the neighbors of the node
            std::vector<int> neighbors = G.getNeighbors(nodeIndex);
            // for each neighbor of the node
            for (int& neighborIndex : neighbors) {
                // if the neighbor is in the community
                if (std::find(subset.begin(), subset.end(), neighborIndex) != subset.end()) {
                    // add the weight of the edge to the sum
                    weight_sum += G.getWeight(nodeIndex, neighborIndex);
                }
            }
        }
    }
    return weight_sum;
}

// get the number of nodes in the community
size_t Community::size() {
    int num_nodes = 0;
    // for each subset of nodes in the community
    for (std::vector<int>& subset : nodeIndices) {
        // add the size of the subset to the sum
        num_nodes += subset.size();
    }
    // convert num_nodes int to size
    size_t c_size = num_nodes;
    return c_size;
}
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
    for (int communityIndex : communityIndices) {
        // get the community
        auto comm = communityIndexMap.find(communityIndex);
        // get the subsets of the community
        std::vector<std::vector<int>>& subsets = comm->second.nodeIndices;
        // flatten the subsets
        std::vector<int> flat_set;
        for (std::vector<int>& subset : subsets) {
            for (int& nodeIndex : subset) {
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
    // update communityIndices
    communityIndices.push_back(newCommunity.communityIndex);
    // update communityIndexMap
    communityIndexMap.insert({newCommunity.communityIndex, newCommunity});
}

// move a node from one community to another
void Partition::updateCommunityMembership(int nodeIndex, int newCommunityIndex) {
    // Check if the new community index is valid, note that it will be a vector if not the end
    /* doesnt work
    if (communityIndexMap.find(newCommunityIndex) == communityIndexMap.end()) {
        Rcpp::stop("Invalid community index: " + std::to_string(newCommunityIndex));
    }
    */
    // Find the current community of the node, get the index of the community and remove the node from the community
    for (auto& comm : communityIndexMap) {
        auto& subsets = comm.second.nodeIndices;

        // Use iterators to avoid issues with erasing elements while iterating
        // we arent using a reference to the subset, because we are erasing elements from it
        for (auto it = subsets.begin(); it != subsets.end(); ) {
            auto& subset = *it;

            // If the node is in the subset
            auto nodeIt = std::find(subset.begin(), subset.end(), nodeIndex);
            if (nodeIt != subset.end()) {
                // Remove the node from the subset
                // this will erase, even though subset is a const reference
                // because it is a pointer to it which is not const?
                subset.erase(nodeIt);

                // If the subset is empty, remove the subset from the community
                if (subset.empty()) {
                    it = subsets.erase(it);
                } else {
                    ++it;
                }

                // Break the loop after finding and updating the node
                break;
            } else {
                // Move to the next subset
                ++it;
            }
        }
    }
    
    // Add the node to the new community
    auto comm = communityIndexMap.find(newCommunityIndex);
    if (comm != communityIndexMap.end()) {
        auto& subsets = comm->second.nodeIndices;
        subsets.push_back({nodeIndex});
    }
}

/*
double Partition::quality_messyImplementation(double resolution_parameter) {
    double mod = 0.0;
    for (size_t c = 0; c < this->number_of_communities(); c++) {
        double csize = this->community_size(c);
        double w = this->total_weight_in_community(c);
        double comm_possible_edges = this->graph->possible_edges(csize);
        mod += w - resolution_parameter * comm_possible_edges;
    }
    return (2.0 - this->graph->is_directed()) * mod;
}
*/

/*
void Partition::renumber_communities() {
    std::unordered_map<size_t, size_t> newCommunityIndices;
    size_t newIndex = 0;

    // Remove empty communities first
    communities.erase(std::remove_if(communities.begin(), communities.end(),
                                     [](const Community& community) { return community.nodeIndices.empty(); }),
                      communities.end());

    // Assign new indices to communities
    for (auto& community : communities) {
        size_t oldIndex = community.communityIndex;
        if (newCommunityIndices.find(oldIndex) == newCommunityIndices.end()) {
            newCommunityIndices[oldIndex] = newIndex++;
        }
        community.communityIndex = newCommunityIndices[oldIndex];
    }

    // Rebuild the communityIndexMap to reflect the new indexing
    communityIndexMap.clear();
    for (size_t i = 0; i < communities.size(); ++i) {
        communityIndexMap[communities[i].communityIndex] = i;
    }
}
*/
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
