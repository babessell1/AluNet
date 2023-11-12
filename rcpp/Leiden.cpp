#include <Rcpp.h>
#include "GraphUtils.h"
#include "Leiden.h"

/*
################################################################################
############################ COMMUNITY CLASS ###################################
*/
// Construct a community base on a set of node indices, and a community index
Community::Community(const std::vector<int>& nodes, int index)
    : communityIndex(index), nodeIndices(nodes) {
}

// get the sum of weights of all edges in the community
double Community::aggregateWeights(Graph& G) {
    double weight_sum = 0.0;
    for (int& nodeIndex : nodeIndices) {
        // get the neighbors of the node
        std::vector<int> neighbors = G.getNeighbors(nodeIndex);
        // for each neighbor of the node
        for (int& neighborIndex : neighbors) {
            // if the neighbor is in the community, add the weight of the edge to the sum
            weight_sum += G.getWeight(nodeIndex, neighborIndex);
        }
    }
    return weight_sum;
}

// get the number of nodes in the community
size_t Community::size() {
    return nodeIndices.size();
}
    
/*
################################################################################
############################ PARTITION CLASS ###################################
*/
// Construct a partition based on a set of communities
Partition::Partition(const std::vector<Community>& communities) : communities(communities) {
    for (const auto& community : communities) {
        // add the community index to the community vector
        communityIndices.push_back(community.communityIndex);
        // add the community to the community map
        communityIndexMap.insert({community.communityIndex, community});
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
    communities.push_back(newCommunity);
    // update communityIndices
    communityIndices.push_back(newCommunity.communityIndex);
    // update communityIndexMap
    communityIndexMap.insert({newCommunity.communityIndex, newCommunity});
}


// move a node from one community to another
void Partition::updateCommunityMembership(int nodeIndex, int newCommunityIndex) {
    // Check if the new community index is valid
    if (communityIndexMap.find(newCommunityIndex) == communityIndexMap.end()) {
        Rcpp::stop("Invalid community index: " + std::to_string(newCommunityIndex));
    }

    // Find the current community of the node and remove the node from the community
    auto it = std::find_if(communities.begin(), communities.end(),
        [nodeIndex](const Community& comm) {
            return std::find(comm.nodeIndices.begin(), comm.nodeIndices.end(), nodeIndex) != comm.nodeIndices.end();
        });

    if (it != communities.end()) {
        auto& nodeIndices = it->nodeIndices;
        auto nodeIt = std::find(nodeIndices.begin(), nodeIndices.end(), nodeIndex);
        if (nodeIt != nodeIndices.end()) {
            nodeIndices.erase(nodeIt);

            // If the community is empty, remove it from the partition
            if (nodeIndices.empty()) {
                communities.erase(it);
            }
        }
    }

    // Add the node to the new community
    auto comm = communityIndexMap.find(newCommunityIndex);
    if (comm != communityIndexMap.end()) {
        comm->second.nodeIndices.push_back(nodeIndex);
    }
} 

/*
// ######################################NOTE#####################################
// ########we should define the following two functions

// Function to get the community of a given vertex
size_t get_community_of_vertex(size_t vertex) {
  // Check if the vertex is in the map
  auto it = communityIndexMap.find(vertex);
  if (it != communityIndexMap.end()) {      
    // Return the community of the vertex
    return it->second;
  } else {
      // Handle the case where the vertex is not found
      throw std::runtime_error("Vertex not found in the partition");
    }
  }


double Partition::diff_move(size_t vertex, size_t new_community){
  // calculate the difference between moving vertex to a new community
  // and the keeping the vertex to the original community
  size_t old_community = this->get_community_of_vertex(vertex);
  
  // If the vertex is already in the new community, the difference is zero
  if (original_community == new_community) {
    return 0.0;      
  }

  // Calculate the quality of the current state
  double original_quality = this->quality();

  // Temporarily move the vertex to the new community
  this->move_vertex_to_community(vertex, new_community);
  
  // Calculate the quality after the move
  double new_quality = this->quality();

  // Move the vertex back to its original community
  this->.move_vertex_to_community(vertex, original_community);
  // The difference in quality is the new quality minus the original quality
  return new_quality - original_quality;
}

*/

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
Optimizer::Optimizer(Graph& G, Partition& P, double gamma) : G(G), P(P), gamma(gamma) {}

/*
// Leiden fast local move iteration
bool Optimizer::moveNodesFast() {
    bool update = false;  // track if the partition has changed

    std::vector<bool> stable_nodes(G.n, false); // track which nodes have not moved
    std::vector<double> cluster_weights(G.n, 0.0);
    std::vector<int> nodes_per_cluster(G.n, 0);  // track the number of nodes in each cluster (more efficient than calling size() on each cluster)
    std::vector<int> unused_clusters(G.n-1, 0);  // track which clusters are empty
    std::vector<int> n_neighboring_clusters(G.n, 0);  // track the number of neighboring clusters for each cluster


    int n_unused_clusters = 0;
    int n_unstable_nodes = G.n;

    std::vector<int> node_queue = RandomGenerator::generateRandomPermutation(G.n);

    // initialize the cluster weights and nodes per cluster
    for (const auto& entry : nodeIndexMap) {
        int nodeIndex = entry.second;
        cluster_weights[nodeIndex] = G.node_weights[nodeIndex];
        nodes_per_cluster[nodeIndex]++;
    }
    // get revered node index map second values
    std::vector<int> node_indices
    for (auto entry : nodeIndexMap) {
        node_indices.push_back(entry.second);
    }
    
    // for each node in the network (go backwards)
    for (int nidx : node_indices.reverse()) {
        // if no nodes in cluser
        if (nodes_per_cluster[nidx] == 0) {
            // add to unused clusters
            unused_clusters[n_unused_clusters] = nidx;
            n_unused_clusters++;
        }
    }

    // .... not done

}
*/

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
    // test 
    P.flattenPartition();
}

Partition initializePartition(const Graph& G) {
    std::vector<Community> communities;

    // Assign each node to its own community
    // for each node in the getNodes
    int communityIndex = 0;
    for (int nodeIndex : G.getNodes()) {
        // Construct a community with a single node
        Community community({nodeIndex}, communityIndex);
        communities.push_back(community);
        communityIndex++;
    }

    Partition P(communities);

    return P;
}

// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations) {

    // Set the resolution parameter
    double gamma = 1.0;

    // Create a graph from the R List
    // use listToGraph from GraphUtils.cpp
    Graph G = listToGraph(graphList);
    // print the number of nodes 
    Rcpp::Rcout << "Number of nodes: " << G.n << std::endl;
    Partition P = initializePartition(G);
    // print the number of communities
    Rcpp::Rcout << "Number of communities: " << P.communities.size() << std::endl;

    // Create an Optimizer
    Optimizer optim(G, P, gamma);

    // Run the Leiden algorithm for the specified number of iterations
    for (int i = 0; i < iterations; i++) {
        // Run the Leiden algorithm
        optim.optimize();
    }

    // get the communities from the partition
    std::vector<int> communities;
    for (size_t i = 0; i < optim.P.communities.size(); i++) {
        communities.push_back(optim.P.communities[i].communityIndex);
    }

    // Convert the vector of communities to an R List
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities
    );

    return result;

}
