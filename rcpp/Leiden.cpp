#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"
#include "Partition.h"
#include "Leiden.h"

#include <queue>


// [[test it]

// not finished, should define an extra function named quality()
// which counts the weights
// [[test it]]

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

double Optimizer::deltaQuality(int n_idx, int new_c_idx, double gamma) const {
  // calculate the difference between moving vertex to a new community

    int old_c_idx = P.nodeCommunityMap.at(n_idx);  // keep track of the old community

    // get quality of the partition before the move
    double old_quality = P.calcQuality(gamma, G);

    // make copy of the partition
    Partition P_new = P;

    // move the node to the new community
    P_new.updateCommunityMembership(n_idx, old_c_idx, new_c_idx);
    // purge empty communities

    // get quality of the partition after the move
    double new_quality = P_new.calcQuality(gamma, G);

    // return the difference in quality
    return new_quality - old_quality;

}

// Leiden fast local move iteration
bool Optimizer::moveNodesFast() {
    bool update = false;  // track if the partition has changed

    std::vector<bool> stable_nodes(G.n, false); // track which nodes have not moved
    std::vector<double> cluster_weights(G.n, 0.0);
    std::vector<double> edge_weights_per_cluster(G.n, 0.0);  // track the sum of edge weights per cluster
    std::vector<int> nodes_per_cluster(G.n, 0);  // track the number of nodes in each cluster (more efficient than calling size() on each cluster)
    std::vector<int> unused_clusters(G.n-1, 0);  // track which clusters are empty
    std::vector<int> n_neighboring_clusters(G.n, 0);  // track the number of neighboring clusters for each cluster

    // print starting message
    Rcpp::Rcout << "Starting move iteration" << std::endl;

    //print initial quality
    Rcpp::Rcout << "Initial quality: " << P.calcQuality(gamma, G) << std::endl;

    int n_unused_clusters = 0;
    int n_unstable_nodes = G.n;

    std::vector<int> random_nodes = RandomGenerator::generateRandomPermutation(G.n);

    std::deque<int> node_queue(G.n);
    int start = 0, end = G.n;
    for (int i = 0; i < G.n; i++) {
        node_queue[i] = random_nodes[i];
    }

    // initialize the cluster weights and nodes per cluster
    for (const auto& entry : G.nodeIndexMap) {
        int node_index = entry.second;
        cluster_weights[node_index] = G.nodeWeights[node_index];
        nodes_per_cluster[node_index]++;
    }

    // get the node indices in the graph
    std::vector<int> node_indices;
    for (auto entry : G.nodeIndexMap) {
        node_indices.push_back(entry.second);
    }
    
    // for each node in the network (go backwards)
    std::vector<int> rev_node_indices;
    for (int i = G.n - 1; i >= 0; i--) {
        rev_node_indices.push_back(node_indices[i]);
    }

    for (int nidx : rev_node_indices) {
        // print nodes per cluster
        //Rcpp::Rcout << nodes_per_cluster[nidx] << std::endl;
        // if no nodes in cluser
        if (nodes_per_cluster[nidx] == 0) {
            // add to unused clusters
            unused_clusters[n_unused_clusters] = nidx;
            n_unused_clusters++;
        }
    }

    // add an empty cluster to the end of the unused clusters if there are no empty clusters
    if (n_unused_clusters == 0) {
        int new_comm_idx = size_t(P.communityIndexMap.size());
        unused_clusters[n_unused_clusters] = new_comm_idx;
        // also add it to the community index map
        Community empty_community({}, new_comm_idx);
        P.communityIndexMap.insert({new_comm_idx, empty_community});
        n_unused_clusters++;
    }

    int counter = 0;
    int bar_width = 70; // bar width
    Rcpp::Rcout << "Progress: ";
    do {
        Rcpp::Rcout << "----------------------------------------" << std::endl;
        counter++;
        Rcpp::Rcout << "Start: " << start << std::endl;
        Rcpp::Rcout << "End: " << end << std::endl;
        Rcpp::Rcout << "n unstable nodes: " << n_unstable_nodes << std::endl;
        //Rcpp::Rcout << "n unused clusters: " << n_unused_clusters << std::endl;
        //Rcpp::Rcout << "n neighboring clusters: " << n_neighboring_clusters.size() << std::endl;
        //Rcpp::Rcout << "cluster weights: " << cluster_weights.size() << std::endl;
        //Rcpp::Rcout << "edge weights per cluster: " << edge_weights_per_cluster.size() << std::endl;
        //Rcpp::Rcout << "nodes per cluster: " << nodes_per_cluster.size() << std::endl;
        //Rcpp::Rcout << "stable nodes: " << stable_nodes.size() << std::endl;
        //Rcpp::Rcout << "node queue: " << node_queue.size() << std::endl;

        int j = node_queue[start++];  // get the next node in the queue
        // if the node is stable, skip it
        // get current community of node j
        int c_idx = P.nodeCommunityMap.at(j);

        // get the neighbors of node j
        std::vector<int> neighbors = G.getNeighbors(j);
        int n_neighboring_clusters = 1;  // track the number of neighboring clusters
        // get the neighboring clusters of node j
        std::set<int> neighboring_clusters;

        // add empty cluster to neighboring clusters
        neighboring_clusters.insert(unused_clusters[n_unused_clusters-1]);
        for (int nn_idx : neighbors) {  // neighbors of node j
            // get the community of the neighbor
            int nc_idx = P.nodeCommunityMap.at(nn_idx);  // nieghbor community index

            // if edge weight of cluster is 0
            if (edge_weights_per_cluster[nc_idx] == 0) {

                // add to neighboring clusters
                neighboring_clusters.insert(nc_idx);
                n_neighboring_clusters++;
            }
            edge_weights_per_cluster[nc_idx] += P.communityIndexMap.at(nc_idx).getClusterWeight(G);
        }

        // initialize the best cluster and the best quality
        int best_cluster = c_idx;
        double best_quality_increment = 0.0;
        // for each neighboring cluster

        for (int nc_idx : neighboring_clusters) {
            double delta_q = deltaQuality(j, nc_idx, gamma);
            // if the quality of the move is better than the best quality
            if (delta_q > best_quality_increment) {
                // update the best quality and best cluster
                best_quality_increment = delta_q;
                best_cluster = nc_idx;
            }
        }

        // mark the node as stable, remove it from the queue
        stable_nodes[j] = true;
        n_unstable_nodes--;

       if (best_cluster != c_idx) {
            // update cluster stats
            P.updateCommunityMembership(j, c_idx, best_cluster);
            // print node and new community
            Rcpp::Rcout << "Moved node: " << j << " to community: " << best_cluster << std::endl;
            cluster_weights[best_cluster] += G.nodeWeights[j];
            nodes_per_cluster[best_cluster]++;
            if (best_cluster == unused_clusters[n_unused_clusters-1]) {
                n_unused_clusters--;
            }

            // get the neighbors of node j
            std::vector<int> neighbors = G.getNeighbors(j);
            for (int nn_idx : neighbors) {
                // get community of neighbor
                int nc_idx = P.nodeCommunityMap.at(nn_idx);
                // print neighbor and neighbor community
                //Rcpp::Rcout << "Neighbor: " << nn_idx << " in community: " << nc_idx << std::endl;
                // if the neighbor is stable and neighbor's cluster is not the best cluster
                if (nc_idx != best_cluster ) { // && stable_nodes[nn_idx] == true) {
                    // print new unstable node found
                    Rcpp::Rcout << "New unstable node: " << nn_idx << std::endl;
                    stable_nodes[nn_idx] = false;
                    n_unstable_nodes++;
                    if (end == node_queue.size()) { // if the queue is full (should really set back to 0 when this happens since it is circular but having weird issues with that)
                        node_queue.resize(node_queue.size() * 1.25);
                    }
                    // add to the end of the queue
                    node_queue[end++] = nn_idx;
                }
         
            }

            update = true;
        }
        // print end of iteration
        //Rcpp::Rcout << "End of iteration: " << counter << std::endl;
        // cyclic queue
        // if the end of the queue has been reached
        if (start == end) {
            start = 0; // start at the beginning
        }

    } while (n_unstable_nodes > 0);

    if (update) {
        // purge empty clusters
        P.purgeEmptyCommunities(true);

        // print purged 
        Rcpp::Rcout << "Purged empty clusters" << std::endl;

        // print new quality
        Rcpp::Rcout << "New quality: " << P.calcQuality(gamma, G) << std::endl;

        // print largest community size
        int largest_community_size = 0;
        for (const auto& entry : P.communityIndexMap) {
            if (entry.second.nodeIndices.size() > largest_community_size) {
                largest_community_size = entry.second.nodeIndices.size();
            }
        }

    }
    return update;
}

std::vector<Community> Optimizer::getWellConnectedCommunities(const Paritition& P_ref, const Community& B) const {
    // Identify communities in the subset that are well connected, that is:
    // the edge weights of C, and the difference of the subset and C is greater than gamma * size of flat(C) * size of flat(S - C)
    std::vector<Community> well_connected_communities; // well connected communities
    // for each community in the refined partition
    for (const auto& entry : P_ref.communityIndexMap) {
        // get the community
        Community C = entry.second;
        // continue if C is not contained in B
        for (int node_index : C.nodeIndices) {
            if (B.nodeIndices.find(node_index) == B.nodeIndices.end()) {
                continue;
            }
        }
        // E(A, B-A)
        double edge_weight_B_A = 0.0;
        // calc the edge weights between B and B - A
        for (int b : B.nodeIndices) {
            if (C.nodeIndices.find(b) == C.nodeIndices.end()) { // if node is not in C
                edge_weight_B_A += G.getWeight(b, C.communityIndex);
            }
        }
        // get the sizes
        int size_A = A.nodeIndices.size();
        int size_B = B.nodeIndices.size();
        // if well connected, add to R
        if (edge_weight_B_A >= gamma * size_A * (size_B - size_A)) {
            well_connected_communities.push_back(entry.second); // add the community to R
        }
    }
    return well_connected_communities;
}

std::vector<Community> Optimizer::getWellConnectedNodes(const Community& B) const {
    // Identify nodes in the subset that are well connected, that is:
    // the edge weights of v, and the difference of the subset and v is greater than gamma * size of flat(v) * size of flat(subset - v)
    // note that v is a node so flat(v) is just 1
    std::vector<int> well_connected_nodes; // well connected nodes
    // for each node in V
    for (int v : B.nodeIndices) {
        // E(v, B-v)
        double edge_weight_B_v = 0.0;
        // calc the edge weights between B and B - v
        for (int b : B.nodeIndices) {
            if (b != v) { // if node is not v
                edge_weight_B_v += G.getWeight(b, v);
            }
        }
        // get the size of B
        int size_B = B.nodeIndices.size();
        // if well connected, add to R
        if (edge_weight_B_v >= gamma * (size_B - 1)) {
            well_connected_nodes.push_back(v); // add the node to R
        }
    }
    return well_connected_nodes;
}

Partition Optimizer::mergeNodesSubset(const Graph& G, Partition& P_ref, Community& S) {
    std::vector<int> R = getWellConnectedNodes(S); // get well connected nodes
    std::vector<Community> T = getWellConnectedCommunities(P_ref, S); // get well connected communities

    // Visit nodes in random order
    std::shuffle(R.begin(), R.end(), std::default_random_engine(std::random_device{}()));
    for (int v : R) {
        // Consider only nodes that have not yet been merged
        if (P_ref.inSingleton(v)) { // Assuming isSingleton checks if v is in a singleton community
            std::unordered_map<int, double> probabilities; // map community index to delta quality based probability
            for (const auto& C : T) { // for each well connected community contained in S
                // calculate the delta quality of moving v to C
                double delta_quality = P_ref.deltaQuality(v, C, gamma);
                if (delta_quality > 0) {  // if the quality improves
                    probabilities[c] = exp(delta_quality / temperature); // add the probability to the map
                }
            }
            // If there are possible communities to merge
            if (!delta_qualities.empty()) {
                // simmulated annealing
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(0, 1);
                int C_prime_index = -1;
                double max_delta_quality = std::numeric_limits<double>::lowest();
      
                // get sum of probabilities
                double sum_prob = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
                // Choose a random community C' with probability proportional to probabilities
                for (const auto& entry : probabilities) {
                    double prob = entry.second / sum_prob;
                    if (dis(gen) < prob) {
                        C_prime_index = entry.first;
                        break;
                    }
                }
                // If a community was chosen, move node v to community C'
                if (C_prime_index != -1) {
                    P_ref.updateCommunityMembership(v, P_ref.nodeCommunityMap.at(v), C_prime_index);
                }
            }
        }
    }
    P_ref.purgeEmptyCommunities(true); // purge empty communities
    return P_ref; // Return the modified partition
}

// this is the main idea and the functions shuold be defined.
Graph Optimizer::aggregateGraph(const Graph& G, const Partition& P) {
    size_t numCommunities = P.numberOfCommunities();
    Graph aggregatedGraph(numCommunities); 

    // Iterate through each edge in the original graph
    for (const auto& edge : G.getEdges()) { // Assuming getEdges() returns all edges in the graph
        int u = edge.first;
        int v = edge.second;
        
        int communityU = P.getCommunityOfNode(u); // Assuming getCommunityOfNode() gets the community of a node
        int communityV = P.getCommunityOfNode(v);

        // Add an edge between the communities in the aggregated graph
        // This will handle multi-edges automatically if the same edge is added multiple times
        // note that our addEdge function is string and should plug in the weight
        // so the function should be revised
        // aggregatedGraph.addEdge(communityU, communityV);
    }

    return aggregatedGraph;
}

Partition Optimizer::refinePartition() const {
    // Implement creating a new partition from the existing partition
    // copy the partition to P_refined
    Partition P_refined = P;
    P_refined.makeSingleton(G);
    // for each community in the original partition
    for (const auto& entry : P.communityIndexMap) {
        // get the community
        Community C = entry.second;
        // merge the nodes
        P_refined = mergeNodesSubset(G, P_refined, C);
    }

    return P_refined;
}

void Optimizer::optimize() {
    // Implement the Leiden algorithm iteration here
    // This involve multiple steps, such as moveNodesFast, refinePartition, mergeNodesSubset, etc.
    bool done = false;
    while (!done) {
        moveNodesFast();
        //P.updateCommunityAssignments(G);  // store assignments for plotting later
        //Partition P_refined = refinePartition();
        //aggregateGraph(P_refined);

        // set true no matter what for now
        done = true;
    }
    // test 
    P.flattenPartition();
    // print flattenedd
    Rcpp::Rcout << "Flattened partition" << std::endl;
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

    Partition P(communities, G);

    return P;
}

// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations) {

    // Set the resolution parameter
    double gamma = 1.0;
    double temperature = 0.01;

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
    
    // Run the Leiden algorithm
    bool done = false;  
    do { 
        // Run the Leiden algorithm
        optim.optimize();
        // temporary
        done = true;

    } while (!done);

    // print the number of communities
    Rcpp::Rcout << "New number of communities: " << optim.P.communityIndexMap.size() << std::endl;

    // print the nodes in each community
    //for (const auto& entry : optim.P.communityIndexMap) {
    //    Rcpp::Rcout << "Community: " << entry.first << std::endl;
    //    for (int node_index : entry.second.nodeIndices) {
    //        Rcpp::Rcout << "Node: " << node_index << std::endl;
    //    }
    //}
    
    /* OLD TESTING CODE
    //move node 168 to community 126
    optim.P.updateCommunityMembership(37, 0, 71);

    // print nodes in community 126
    Rcpp::Rcout << "Nodes in community 0: " << std::endl;
    for (int node_index : optim.P.communityIndexMap.at(0).nodeIndices) {
        Rcpp::Rcout << node_index << std::endl;
    }
    // print nodes in community 168
    Rcpp::Rcout << "Nodes in community 71: " << std::endl;
    for (int node_index : optim.P.communityIndexMap.at(71).nodeIndices) {
        Rcpp::Rcout << node_index << std::endl;
    }

    // purge empty communities
    optim.P.purgeEmptyCommunities(false);

    */

    // get the communities from the partition
    std::vector<int> communities;
    std::vector<int> nodes;
    for (const auto& entry : optim.P.communityIndexMap) {
        communities.push_back(entry.first);
        for (int node_index : entry.second.nodeIndices) {
            nodes.push_back(node_index);
        }
    }

    // creat vector of node names in the same order as the nodes vector
    std::vector<std::string> node_names;
    for (int node_index : nodes) {
        node_names.push_back(G.getNodeName(node_index));
    }

    // Convert the vector of communities to an R List
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities,
        Rcpp::Named("nodes") = nodes,
        Rcpp::Named("node_names") = node_names
    );

    return result;

}


// should define in the .h file

// Partition currentPartition;

/*
void Optimizer::refinePartition() {
  bool improved = false;
  do {
        improved = false;
        size_t numberOfNodes = currentPartition.number_of_nodes();

        for (size_t i = 0; i < numberOfNodes; ++i) {
          double bestQualityIncrease = 0.0;
            int bestCommunityIndex = -1;
            int originalCommunityIndex = currentPartition.getCommunityOfNode(i);
            auto originalQuality = calculateQuality();
            // I think that we would better define a refactor of the calculateQuality function because this would be more convenient and easy to understand hah

            // Evaluate moving node to different clusters
            size_t numberOfCommunities = currentPartition.number_of_communities();
            for (size_t j = 0; j < numberOfCommunities; ++j) {
                if (j != originalCommunityIndex) {

                    // Move node to the new community and calculate quality
                    currentPartition.updateCommunityMembership(i, j);
                    clusters[originalClusterId].removeNode(node);

                    double newQuality = calculateQuality();
                    double qualityIncrease = newQuality - originalQuality;

                    if (qualityIncrease > bestQualityIncrease) {
                        bestQualityIncrease = qualityIncrease;
                        bestCommunityIndex = j;
                    }

                    // Move node back to original community
                   currentPartition.updateCommunityMembership(i, originalCommunityIndex);
               }
           }

           // If a better community is found, move the node there
           if (bestCommunityIndex != -1) {
              currentPartition.updateCommunityMembership(i, bestCommunityIndex);
              improved = true;
          }
        }
    } while (improved);
}
*/
