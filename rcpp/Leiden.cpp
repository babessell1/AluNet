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
    for (int& node_index : nodeIndices) {
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

// get the number of nodes in the community
size_t Community::size() {
    return nodeIndices.size();
}
    
/*
################################################################################
############################ PARTITION CLASS ###################################
*/
// Construct a partition based on a set of communities
Partition::Partition(const std::vector<Community>& communities) {
    for (const auto& community : communities) {
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
}

void Partition::updateCommunityMembership(int node_index, int old_community_index, int new_community_index) {
    // Find the current community of the node and remove the node from the community
    auto old_comm = communityIndexMap.find(old_community_index);
    if (old_comm != communityIndexMap.end()) {
        // if the node is in the community, remove it
        auto node_it = std::find(old_comm->second.nodeIndices.begin(), old_comm->second.nodeIndices.end(), node_index);
        if (node_it != old_comm->second.nodeIndices.end()) {
            old_comm->second.nodeIndices.erase(node_it);
        } else {
            Rcpp::stop("Node not found in the old community: " + std::to_string(node_index));
        }
    } else {
        Rcpp::stop("Old community not found: " + std::to_string(old_community_index));
    }

    // Add the node to the new community
    communityIndexMap.at(new_community_index).nodeIndices.push_back(node_index);
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

void Partition::updateCommunityMembership(int node_index, int new_community_index) {
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
    for (const auto& entry : node_indexMap) {
        int node_index = entry.second;
        cluster_weights[node_index] = G.node_weights[node_index];
        nodes_per_cluster[node_index]++;
    }
    // get revered node index map second values
    std::vector<int> node_indices
    for (auto entry : node_indexMap) {
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
    int community_index = 0;
    for (int node_index : G.nodes) {
        // Construct a community with a single node
        Community community({node_index}, community_index);
        communities.push_back(community);
        community_index++;
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
    Rcpp::Rcout << "Number of communities: " << P.communityIndexMap.size() << std::endl;

    // Create an Optimizer
    Optimizer optim(G, P, gamma);

    // Run the Leiden algorithm for the specified number of iterations
    for (int i = 0; i < iterations; i++) {
        // Run the Leiden algorithm
        optim.optimize();
    }

    // move node 168 to community 126
    optim.P.updateCommunityMembership(168, 168, 126);

    // print nodes in community 126
    Rcpp::Rcout << "Nodes in community 126: " << std::endl;
    for (int node_index : optim.P.communityIndexMap.at(126).nodeIndices) {
        Rcpp::Rcout << node_index << std::endl;
    }
    // print nodes in community 168
    Rcpp::Rcout << "Nodes in community 168: " << std::endl;
    for (int node_index : optim.P.communityIndexMap.at(168).nodeIndices) {
        Rcpp::Rcout << node_index << std::endl;
    }

    // get the communities from the partition
    std::vector<int> communities;
    std::vector<int> nodes;
    for (const auto& entry : optim.P.communityIndexMap) {
        communities.push_back(entry.first);
        for (int node_index : entry.second.nodeIndices) {
            nodes.push_back(node_index);
        }
    }

    // Convert the vector of communities to an R List
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities,
        Rcpp::Named("nodes") = nodes
    );

    return result;

}
