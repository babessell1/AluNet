#include <Rcpp.h>
#include "GraphUtils.h"
#include "Community.h"
#include "Partition.h"
#include "Leiden.h"

#include <queue>
#include <algorithm>

/*
################################################################################
############################ OPTIMIZER CLASS ###################################
*/
// Construct an optimizer based on a graph and a partition
Optimizer::Optimizer(Graph& G, Partition& P, double gamma, double theta) : G(G), P(P), gamma(gamma), theta(theta) {
    this->communityAssignments = std::unordered_map<std::string, int>();
    // initialize the community assignments
    for (const auto& entry : G.nodeIndexMap) {
        communityAssignments.insert({entry.first, entry.second});
    }
}

void Optimizer::updateCommunityAssignments(const Partition& P, const std::unordered_map<std::string, int> original_nodeIndexMap) {
    // currently stored value in communityAssignments is the previous community assignment which is the latest node index
    // we update by matching new nodes to the previous community assignment and changing the community assignment to the new community index
    // for each node in the graph
    for (const auto& iter : original_nodeIndexMap) {
        // get the node index
        int node_idx = iter.second;
        // for each item in the community assignments
        for (auto entry : communityAssignments) {
            // if the node index matches the last community assignment
            if (entry.second == node_idx) {
                // set the community assignment to the new community index
                communityAssignments[entry.first] = P.nodeCommunityMap.at(node_idx);
            }
        }
    }
}

double Optimizer::deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const {
  // calculate the difference between moving vertex to a new community

    int old_c_idx = P.nodeCommunityMap.at(n_idx);  // keep track of the old community

    if (recalculate) {
      // get quality of the partition before the move
      double old_quality = P.calcQuality(gamma, G);
      P.quality = old_quality;
    }

    // make copy of the partition
    Partition P_new = P;

    // move the node to the new community
    P_new.updateCommunityMembership(n_idx, old_c_idx, new_c_idx);
    // purge empty communities

    // get quality of the partition after the move
    double new_quality = P_new.calcQuality(gamma, G);

    // return the difference in quality

    // print delta quality
    //Rcpp::Rcout << "Delta quality: " << new_quality - P.quality << std::endl;
    return new_quality - P.quality;

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
    size_t start = 0, end = G.n;
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
    //Rcpp::Rcout << "Progress: ";
    do {
        //Rcpp::Rcout << "----------------------------------------" << std::endl;
        counter++;
        //Rcpp::Rcout << "Start: " << start << std::endl;
        //Rcpp::Rcout << "End: " << end << std::endl;
        //Rcpp::Rcout << "n unstable nodes: " << n_unstable_nodes << std::endl;
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

        // get the quality of the partition before and after the move
        double quality = P.calcQuality(gamma, G);
        P.quality = quality;
        for (int nc_idx : neighboring_clusters) {
            double delta_q = deltaQuality(j, nc_idx, gamma, false);

            // if the quality of the move is better than the best quality, update the best quality and best cluster
            if (delta_q > best_quality_increment) {
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

                // neighbor's community is not the new best community
                if (nc_idx != best_cluster ) {
                    stable_nodes[nn_idx] = false;  // mark the neighbor as unstable
                    n_unstable_nodes++;
                    // add the neighbor to the queue
                    node_queue[end++] = nn_idx;

                    
                    //if (end == node_queue.size()) { // if the queue is full (should really set back to 0 when this happens since it is circular but having weird issues with that)
                    //    node_queue.resize(node_queue.size() * 1.25);
                    //}
                    // add to the end of the queue
                    //node_queue[end++] = nn_idx;
                }
            }
            update = true;
        }
        // cyclic queue
        // if the end of the queue has been reached
        if (start == end) {
            start = 0; // start at the beginning
        }

    } while (n_unstable_nodes > 0);

    // print new quality
    Rcpp::Rcout << "New quality: " << P.calcQuality(gamma, G) << std::endl;

    // print largest community size
    size_t largest_community_size = 0;
    for (const auto& entry : P.communityIndexMap) {
        if (entry.second.nodeIndices.size() > largest_community_size) {
            largest_community_size = entry.second.nodeIndices.size();
        }
    }

    P.purgeEmptyCommunities(true);
    return update;
}

std::vector<Community> Optimizer::getWellConnectedCommunities(const Community& B) const {
    // Identify communities in the subset that are well connected, that is:
    // the edge weights of C, and the difference of the subset and C is greater than gamma * size of flat(C) * size of flat(S - C)
    // get the sizes
    int size_B = B.nodeIndices.size();
    std::vector<Community> well_connected_communities; // well connected communities
    // for each community in the refined partition
    for (const auto& entry : P.communityIndexMap) {
        int size_C = entry.second.nodeIndices.size();
        // get the community
        Community C = entry.second;  // for each C in P
        size_t count_c_in_b = 0; // count the number of nodes in C that are in B
        double edge_weight_B_A = 0.0;
        for (int b : B.nodeIndices) {  // for each node in B (subset)
            for (int c : C.nodeIndices) {  // for each node in C
                if (b == c) {  // if node is in both B and C
                    count_c_in_b++; // increment count of nodes in C that are also in B
                } else if(G.hasEdge(b, c)) {  // if node is in B but not C
                    edge_weight_B_A += G.getWeight(b, c); // add the edge weight between B and C
                }
            }
        }
        bool is_contained = count_c_in_b == C.nodeIndices.size(); // check if C is contained by B
        if (is_contained) {
            continue;  // continue if C is not contained by B
        }
        // if well connected, add to R
        if (edge_weight_B_A >= gamma * size_C * (size_B - size_C)) {
            well_connected_communities.push_back(entry.second); // add the community to R
        }
    }
    return well_connected_communities;
}

std::vector<int> Optimizer::getWellConnectedNodes(const Community& B) const {
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
            if (b != v && G.hasEdge(b, v)) { // if node is not v
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


void Optimizer::mergeNodesSubset(Community& S) {
    std::vector<int> R = getWellConnectedNodes(S); // get well connected nodes
    std::vector<Community> T = getWellConnectedCommunities(S); // get well connected communities

    // Visit nodes in random order
    std::shuffle(R.begin(), R.end(), std::default_random_engine(std::random_device{}()));
    for (int v : R) {
        // Consider only nodes that have not yet been merged
        double quality = P.calcQuality(gamma, G);
        P.quality = quality;
        if (P.inSingleton(v)) { // Assuming isSingleton checks if v is in a singleton community
            std::unordered_map<int, double> probabilities; // map community index to delta quality based probability
            for (const auto& C : T) { // for each well connected community contained in S
                // calculate the delta quality of moving v to C
                int c_idx = C.communityIndex;
                double delta_quality = deltaQuality(v, c_idx, gamma, false);
                if (delta_quality > 0) {  // if the quality improves
                    probabilities[c_idx] = exp(delta_quality / theta); // add the probability to the map
                }
            }
            // If there are possible communities to merge
            if (!probabilities.empty()) {
                // simmulated annealing
                std::random_device rd;
                std::mt19937 gen(rd());
                std::uniform_real_distribution<> dis(0, 1);
                int C_prime_index = -1;
                //double max_delta_quality = std::numeric_limits<double>::lowest(); ??
      
                // get sum of probabilities
                double sum_prob = std::accumulate(probabilities.begin(), probabilities.end(), 0.0,
                [](double sum, const std::pair<int, double>& p) {
                        return sum + p.second;
                });
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
                    P.updateCommunityMembership(v, P.nodeCommunityMap.at(v), C_prime_index);
                }
            }
        }
    }
}

Graph Optimizer::aggregateGraph() {
    size_t num_communities = P.communityIndexMap.size();

    Graph aggregated_graph(num_communities);   // V <- P
    
    // Iterate through each edge in the original graph
    for (const auto& u : G.edgeWeights) {
        for (const auto& v : u.second) {

            int u_idx = u.first;
            int v_idx = v.first;

            double w = v.second;

            int u_comm = P.nodeCommunityMap.at(u_idx);
            int v_comm = P.nodeCommunityMap.at(v_idx);

            // TODO: HOW TO ACTUALLY HANDLE THE AGGREGATE WEIGHTS???
            // do we add them, take the max, and for self looped edges, do we add the weight to the node weight? who knows!

            // do NOT want self loops
            if (u_comm == v_comm) {
                continue;
            }

            // Add an edge between the communities if not already present, or update the weight if it is
            aggregated_graph.edgeWeights[u_comm][v_comm] += w; // Assumes edgeWeights is a suitable data structure

            // Add the node weights to the aggregated graph
            aggregated_graph.nodeWeights[u_comm] += G.nodeWeights[u_idx];
            aggregated_graph.nodeWeights[v_comm] += G.nodeWeights[v_idx];

            // Update the total edge weight
            aggregated_graph.totalEdgeWeight += w;

            // Update the number of possible edges
            if (aggregated_graph.edgeWeights[u_comm][v_comm] == 0) {
                aggregated_graph.possibleEdges++;
            }

            // Add the nodes to the aggregated graph
            aggregated_graph.nodes.push_back(u_comm);
            aggregated_graph.nodes.push_back(v_comm);

            // Update the node index map
            aggregated_graph.nodeIndexMap.insert({std::to_string(u_comm), u_comm});
            aggregated_graph.nodeIndexMap.insert({std::to_string(v_comm), v_comm});
        }
    }

    aggregated_graph.updateNodeProperties(false); // update node properties

    return aggregated_graph;
}

void Optimizer::refinePartition(const Partition& P_original) {
    // Implement creating a new partition from the existing partition
    // copy the partition to P_original
    P.makeSingleton(G);
    // for each community in the original partition
    for (const auto& entry : P_original.communityIndexMap) {
        // get the community
        Community C = entry.second;
        mergeNodesSubset(C);
    }
    // print done merging nodes
    Rcpp::Rcout << "Done merging nodes" << std::endl;
    P.purgeEmptyCommunities(true);
}

void Optimizer::optimize(int iterations) {
    // Implement the Leiden algorithm iteration here
    // This involve multiple steps, such as moveNodesFast, refinePartition, mergeNodesSubset, etc.
    bool done = false;
    int counter = 0;
    std::unordered_map<std::string, int>  og_nodeIndexMap = G.nodeIndexMap;  // store original node index map to track community assignments
    while (!done && iterations > counter) {
        bool improved = moveNodesFast();
        std::vector<int> community_indices = P.getCommunityIndices();

        // if community size is equal to number of nodes, set done to true
        counter++;
        if (static_cast<int>(community_indices.size()) == G.n || !improved) {
            done = true;

        // otherwise, refine the partition and aggregate the graph
        } else {
            Rcpp::Rcout << "debug 1 number of nodes: " << G.n << std::endl;
            Partition P_save = P;  // but maintain a copy of the partition before refinement
            Rcpp::Rcout << "Refining partition" << std::endl;
            refinePartition(P_save);
            Rcpp::Rcout << "debug 2 number of nodes: " << G.n << std::endl;

            // collapse the communities into a single node in a new graph
            Rcpp::Rcout << "Aggregating graph" << std::endl;
            this->G = aggregateGraph();
            
            // update the community assignmens
            Rcpp::Rcout << "Updating community assignments" << std::endl;
            updateCommunityAssignments(P, og_nodeIndexMap);

            // print the number of nodes in the aggregated graph
            Rcpp::Rcout << "Aggregated number of nodes: " << G.n << std::endl;

            this->P = initializePartition(G);
        }

        if (!improved) {
            // print partition could not be improved
            Rcpp::Rcout << "Partition could not be improved!" << std::endl;
        }
    }
    P.flattenPartition();
    // print flattenedd
    Rcpp::Rcout << "Flattened partition" << std::endl;
}

Partition initializePartition(Graph& G) {
    std::vector<Community> communities;
    communities.reserve(G.nodes.size());
    // Assign each node to its own community
    // for each node in the getNodes
    for (int node_index : G.nodes) {
        // Construct a community with a single node
        Community community({ node_index }, node_index); // Use the node index as the community index
        communities.push_back(community);
    }

    Partition P(communities, G);

    // Update the nodeCommunityMap in the Partition
    for (int node_index : G.nodes) {
        P.nodeCommunityMap[node_index] = node_index;
    }

    return P;
}




/*
###################################################################################
####################### R INTERFACE FUNCTIONS #####################################
*/

// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations, double gamma, double theta) {


    //double gamma = .1;
    //double theta = 0.01;

    // Create a graph from the R List
    // use listToGraph from GraphUtils.cpp
    Graph G = listToGraph(graphList);
    // print the number of nodes 
    Rcpp::Rcout << "Number of nodes: " << G.n << std::endl;
    Partition P = initializePartition(G);
    // print the number of communities
    Rcpp::Rcout << "Number of communities: " << P.communityIndexMap.size() << std::endl;

    // Create an Optimizer
    Optimizer optim(G, P, gamma, theta);
    
    // Run the Leiden algorithm
    optim.optimize(iterations);
    Rcpp::Rcout << "Final number of communities: " << optim.P.communityIndexMap.size() << std::endl;


    // print the number of communities
    Rcpp::Rcout << "Final number of communities: " << optim.P.communityIndexMap.size() << std::endl;
    
    Rcpp::Rcout << "Starting community processing..." << std::endl;
    Rcpp::List communities = Rcpp::List::create();
    int communityCount = 0;

    for (const auto& communityPair : optim.P.communityIndexMap) {
        Rcpp::Rcout << "Processing community " << communityCount << std::endl;
        Rcpp::Rcout << "Community ID: " << communityPair.first << " - Node indices: ";

        for (int nodeIndex : communityPair.second.nodeIndices) {
            Rcpp::Rcout << nodeIndex << " ";
        }
        Rcpp::Rcout << std::endl;

        Rcpp::List communityData = Rcpp::List::create(
            Rcpp::Named("nodes") = communityPair.second.nodeIndices
        );
        communities.push_back(communityData);
        communityCount++;
    }

    Rcpp::Rcout << "Processed " << communityCount << " communities." << std::endl;
    //return Rcpp::List::create(Rcpp::Named("communities") = communities);

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

    // Convert the vector of communities to an R List
    //Rcpp::List result = Rcpp::List::create(
    //    Rcpp::Named("communities") = communities,
    //    Rcpp::Named("nodes") = nodes,
    //    Rcpp::Named("node_names") = node_names
    //);

    //return result;

    //Rcpp::List output = Rcpp::List::create(
    //    Rcpp::Named("graph") = optim.G.graphToRList(),
    //    Rcpp::Named("assignments") = optim.P.communityAssignments,
    //    Rcpp::Named("quality") = optim.P.calcQuality(gamma, G)
    //);

    //return output;
    return optim.G.graphToRList(optim.communityAssignments, optim.P.calcQuality(gamma, G));

}
