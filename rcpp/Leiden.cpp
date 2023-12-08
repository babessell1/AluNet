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
    for (const auto& entry : G.getNodeIndexMap()) {
        communityAssignments.insert({entry.first, entry.second});
    }
}

void Optimizer::updateCommunityAssignments(const Partition& P, const std::unordered_map<std::string, int>& original_nodeIndexMap) {
    // currently stored value in communityAssignments is the previous community assignment which is the latest node index
    // we update by matching new nodes to the previous community assignment and changing the community assignment to the new community index
    // for each node in the graph
    for (const auto& iter : original_nodeIndexMap) {
        // get the node index
        int node_idx = iter.second;
        // for each item in the community assignments
        for (auto& entry : communityAssignments) {
            // if the node index matches the last community assignment
            if (entry.second == node_idx) {
                // set the community assignment to the new community index
                communityAssignments[entry.first] = P.getNodeCommunityIdx(node_idx);
            }
        }
    }
}

double Optimizer::deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const {
  // calculate the difference between moving vertex to a new community

    int old_c_idx = P.getNodeCommunityIdx(n_idx);  // keep track of the old community

    if (recalculate) {
      // get quality of the partition before the move
      double old_quality = P.calcQuality(gamma, G, false);
      P.setQuality(old_quality);
    }

    // make copy of the partition
    Partition P_new = P;

    // move the node to the new community
    P_new.updateCommunityMembership(n_idx, old_c_idx, new_c_idx);
    P_new.setNodeCommunity(n_idx, new_c_idx);
    // purge empty communities

    // get quality of the partition after the move
    double new_quality = P_new.calcQuality(gamma, G, false);

    // return the difference in quality

    // print delta quality
    //Rcpp::Rcout << "Delta quality: " << new_quality - P.quality << std::endl;
    return new_quality - P.getQuality();

}

// Leiden fast local move iteration
bool Optimizer::moveNodesFast() {
    bool update = false;  // track if the partition has changed

    // print intial quality
    Rcpp::Rcout << "Initial quality: " << P.calcQuality(gamma, G, false) << std::endl;

    int n_ = G.getN();
    std::vector<bool> stable_nodes(n_, false); // track which nodes have not moved
    std::vector<double> cluster_weights(n_, 0.0);
    std::vector<double> edge_weights_per_cluster(n_, 0.0);  // track the sum of edge weights per cluster
    std::vector<int> nodes_per_cluster(n_, 0);  // track the number of nodes in each cluster
    std::vector<int> unused_clusters(n_ - 1, 0);  // track which clusters are empty
    std::vector<int> n_neighboring_clusters(n_, 0);  // track the number of neighboring clusters for each cluster

    int n_unused_clusters = 0;
    int n_unstable_nodes = n_;

    std::vector<int> random_nodes = RandomGenerator::generateRandomPermutation(n_);

    std::deque<int> node_queue(n_);
    size_t start = 0, end = n_;
    for (int i = 0; i < n_; i++) {
        node_queue[i] = random_nodes[i];
    }

    // initialize the cluster weights and nodes per cluster
    for (int i = 0; i < n_; i++) {
        int comm_idx = P.getNodeCommunityIdx(i);
        cluster_weights[comm_idx] += G.getNodeWeights()[i];
        nodes_per_cluster[comm_idx]++;
    }

    // track unused clusters
    for (int i = n_ - 1; i >= 0; i--)
        if (nodes_per_cluster[i] == 0) {
            unused_clusters[n_unused_clusters] = i;
            n_unused_clusters++;
        }

    // start of the cyclical process
    do {
        int node_idx = node_queue[start++];  // get the next node in the queue
        int cluster_idx = P.getNodeCommunityIdx(node_idx);
        Rcpp::Rcout << "-------((  " << G.getNodeName(node_idx) << "  ))-------" << std::endl;
        //Rcpp::Rcout << "Next in queue (idx): " << node_idx << " - Cluster: " << cluster_idx << " iteration: " << start << std::endl;

        cluster_weights[cluster_idx] -= G.getNodeWeight(node_idx);
        nodes_per_cluster[cluster_idx]--;
        if (nodes_per_cluster[cluster_idx] == 0) {
            unused_clusters[n_unused_clusters] = cluster_idx;
            n_unused_clusters++;
        }

        //Rcpp::Rcout << "Unused clusters: " << n_unused_clusters << std::endl;

        // track neighboring clusters
        std::vector<int> neighboring_clusters(n_);
        neighboring_clusters[0] = unused_clusters[n_unused_clusters - 1];
        int n_neighboring_clusters = 1;

        // print G.getEdgeWeights().at(node_idx).size()
        //Rcpp::Rcout << "Edge weights size: " << G.getThisEdgeWeights(node_idx).size() << std::endl;
        for (const auto& neighbor_weight : G.getThisEdgeWeights(node_idx)) {
            int neighbor_cluster = P.getNodeCommunityIdx(neighbor_weight.first);
            //Rcpp::Rcout << "Neighbor cluster: " << neighbor_cluster << std::endl;
            //Rcpp::Rcout << "Edge weight, neighbor: " << neighbor_weight.second << std::endl;
            //Rcpp::Rcout << "Edge weight per cluster: " << edge_weights_per_cluster[neighbor_cluster] << std::endl;
            if (edge_weights_per_cluster[neighbor_cluster]==0) {
                neighboring_clusters[n_neighboring_clusters] = neighbor_cluster;
                n_neighboring_clusters++;
            }
            edge_weights_per_cluster[neighbor_cluster] += neighbor_weight.second;
        }

        //Rcpp::Rcout << "Neighboring clusters: " << n_neighboring_clusters << std::endl;

        // calculate the best cluster for the node to move to
        int best_cluster = cluster_idx;
        double max_delta_q = edge_weights_per_cluster[cluster_idx] - G.getNodeWeight(node_idx) * cluster_weights[cluster_idx] * gamma;
        for (int i = 0; i < n_neighboring_clusters; i++) {
            int idx = neighboring_clusters[i];
            //Rcpp::Rcout << "edge weight: " << edge_weights_per_cluster[idx] << std::endl;
            //Rcpp::Rcout << "node weight: " << G.getNodeWeight(node_idx) << std::endl;
            //Rcpp::Rcout << "cluster weight: " << cluster_weights[idx] << std::endl;
            //Rcpp::Rcout << "gamma: " << gamma << std::endl;
            // E(C’, C’) - uweight * C’weight * 𝛾
            double delta_quality = edge_weights_per_cluster[idx] - G.getNodeWeight(node_idx) * cluster_weights[idx] * gamma;
            //Rcpp::Rcout << "Delta quality: " << delta_quality << std::endl;
            if (delta_quality > max_delta_q) {
                best_cluster = idx;
                max_delta_q = delta_quality;
            }

            Rcpp::Rcout << "delta quality: " << delta_quality << std::endl;

            edge_weights_per_cluster[idx] = 0;
        }

        //Rcpp::Rcout << "Best cluster: " << best_cluster << std::endl;

        // move the node to the best cluster
        cluster_weights[best_cluster] += G.getNodeWeight(node_idx); 
        nodes_per_cluster[best_cluster]++; 

        if (best_cluster == unused_clusters[n_unused_clusters - 1]) {
            n_unused_clusters--;
        }

        //Rcpp::Rcout << "Unused clusters: " << n_unused_clusters << std::endl;

        // mark the node as stable
        stable_nodes[node_idx] = true;
        n_unstable_nodes--;

        //Rcpp::Rcout << "Unstable nodes: " << n_unstable_nodes << std::endl;

        // check if cluster is different
        if (best_cluster != cluster_idx) {
            Rcpp::Rcout << "Moving node (" << G.getNodeName(node_idx) << ")  from cluster " << cluster_idx << " to cluster: " << best_cluster << std::endl;
            Rcpp::Rcout << "with other nodes: " << std::endl;

            for ( const auto& i : P.getNodeCommunity(best_cluster).getNodeIndices() ) {
                Rcpp::Rcout << G.getNodeName(i) << std::endl;
            }
            P.updateCommunityMembership(node_idx, cluster_idx, best_cluster);
            //Rcpp::Rcout << "Node moved" << std::endl;
            update = true;
 
            // make the neighbors of the node unstable
            for (const auto& neighbor_weight : G.getThisEdgeWeights(node_idx)) {
                int neighbor_node = neighbor_weight.first; // index of neighbor          
                if (stable_nodes[neighbor_node] && P.getNodeCommunityIdx(neighbor_node) != best_cluster) {
                    stable_nodes[neighbor_node] = false;
                    n_unstable_nodes++;
                    //Rcpp::Rcout << "Adding neighbor to queue: " << neighbor_node << std::endl;
                    node_queue[(start + n_unstable_nodes < n_) ? (start + n_unstable_nodes) : (start + n_unstable_nodes - n_)] = neighbor_node;
                }
            }
        }
        start = (start < n_ - 1) ? (start + 1) : 0;

        //Rcpp::Rcout << "Start: " << start << std::endl;


    } while (n_unstable_nodes > 0);

    //Rcpp::Rcout << "Quality: " << P.calcQuality(gamma, G, false) << std::endl;

    P.purgeEmptyCommunities(false);

    return update;
}
std::vector<Community> Optimizer::getWellConnectedCommunities(const Community& B) const {
    // Identify communities in the subset that are well connected, that is:
    // the edge weights of C, and the difference of the subset and C is greater than gamma * size of flat(C) * size of flat(S - C)
    // get the sizes
    int size_B = B.getNodeIndices().size();
    std::vector<Community> well_connected_communities; // well connected communities
    // for each community in the refined partition
    for (const auto& entry : P.getCommunityIndexMap()) {
        int size_C = entry.second.getNodeIndices().size();
        // get the community
        Community C = entry.second;  // for each C in P
        size_t count_c_in_b = 0; // count the number of nodes in C that are in B
        double edge_weight_B_A = 0.0;
        for (const int& b : B.getNodeIndices()) {  // for each node in B (subset)
            for (const int& c : C.getNodeIndices()) {  // for each node in C
                if (b == c) {  // if node is in both B and C
                    count_c_in_b++; // increment count of nodes in C that are also in B
                } else if(G.hasEdge(b, c)) {  // if node is in B but not C
                    edge_weight_B_A += G.getWeight(b, c); // add the edge weight between B and C
                }
            }
        }
        bool is_contained = count_c_in_b == C.getNodeIndices().size(); // check if C is contained by B
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
    for (const int& v : B.getNodeIndices()) {
        // E(v, B-v)
        double edge_weight_B_v = 0.0;
        // calc the edge weights between B and B - v
        for (const int& b : B.getNodeIndices()) {
            if (b != v && G.hasEdge(b, v)) { // if node is not v
                edge_weight_B_v += G.getWeight(b, v);
            }
        }
        // get the size of B
        int size_B = B.getNodeIndices().size();
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
    for (const int& v : R) {
        // Consider only nodes that have not yet been merged
        double quality = P.calcQuality(gamma, G, false);
        P.setQuality(quality);
        if (P.inSingleton(v)) { // Assuming isSingleton checks if v is in a singleton community
            std::unordered_map<int, double> probabilities; // map community index to delta quality based probability
            for (const auto& C : T) { // for each well connected community contained in S
                // calculate the delta quality of moving v to C
                int c_idx = C.getCommunityIndex();
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
                // Choose a random community C' with probability proportional to probabilities  //? not good
                for (const auto& entry : probabilities) {
                    double prob = entry.second / sum_prob;
                    if (dis(gen) < prob) {
                        C_prime_index = entry.first;
                        break;
                    }
                }
                // If a community was chosen, move node v to community C'
                if (C_prime_index != -1) {
                    P.updateCommunityMembership(v, P.getNodeCommunityMap().at(v), C_prime_index);
                    P.setNodeCommunity(v, C_prime_index);
                }
            }
        }
    }
}

Graph Optimizer::aggregateGraph(Partition& P_comm) const {
    int num_communities = P_comm.getCommunityIndexMap().size();
    Graph aggregate_graph(num_communities);

    // create a map of community indices to ordered new community indices
    std::unordered_map<int, int> new_community_indices;
    int new_index = 0;
    for (const auto& entry : P.getCommunityIndexMap()) {
        int old_index = entry.first;
        new_community_indices[old_index] = new_index;
        new_index++;
    }

    // for each community in the partition
    for ( auto iter : P.getCommunityIndexMap() ) {
        int old_c_idx = iter.first;
        int c_idx = new_community_indices.at(old_c_idx);  // get the new community index

        // print community index
        //Rcpp::Rcout << "----Community index: " << c_idx << std::endl;

        std::string c_str = std::to_string(c_idx);

        //Rcpp::Rcout << "set node index" << std::endl;
        aggregate_graph.setNodeIndex(c_str, c_idx);
        //Rcpp::Rcout << "add node" << std::endl;
        aggregate_graph.addNode(c_idx);
        //Rcpp::Rcout << "set node weight" << std::endl;
        aggregate_graph.setNodeWeight(c_idx, 0);
        //Rcpp::Rcout << "init edge weight" << std::endl;
        aggregate_graph.initEdgeWeight(c_idx);

        //Rcpp::Rcout << "get nodes in community" << std::endl;

        // get the nodes in the community
        std::vector<int> nodesInCommunity = P.getNodeCommunity(old_c_idx).getNodeIndices();
        
        std::map<int, int> neighbors;  // initialize a map of neighbors
        std::map<int, double> edge_weights; // initialize a map of edge weights

        for (int node : nodesInCommunity) { // for each node in the community
            //Rcpp::Rcout << "Node: " << G.getNodeName(node) << std::endl;

            int u_comm = c_idx;   // get the community index (just a rename)

            // add the node weight to the aggregate graph
            //aggregate_graph.nodeWeights.at(c_idx) += G.getNodeWeights().at(node);
            aggregate_graph.incNodeWeight(u_comm, G.getNodeWeights().at(node));

            for (const auto& neighbour : G.getEdgeWeights()[node]) {  // for each neighbor of the node
                int v_idx = neighbour.first;  // get the neighbor index
                int v_comm = P.getNodeCommunityIdx(v_idx);  // get the neighbor community index
                std::string v_str = std::to_string(v_comm);  // convert the neighbor community to a string

                // ensure the neighbor community in the node index map
                //aggregate_graph.getNodeIndexMap()[v_str] = v_comm;   
                aggregate_graph.setNodeIndex(v_str, v_comm);

                // print neighbor index and community
                //Rcpp::Rcout << "------Neighbor index: " << v_idx << " - Neighbor community: " << v_comm << std::endl;

                double weight = neighbour.second; // get the edge weight shared between the node and the neighbor

                if (u_comm != v_comm) {  // don't want self loops
                    // check if the neighbor is already in the map
                    if (aggregate_graph.getEdgeWeights().at(u_comm) == std::unordered_map<int, double>()) {
                        // print new map
                        //Rcpp::Rcout << "------------New map" << std::endl;
                        // add ege
                        aggregate_graph.addEdge(c_str, v_str, weight);
                    } else {
                        // print existing map
                        //Rcpp::Rcout << "------------new weight" << std::endl;
                        //aggregate_graph.edgeWeights.at(c_idx)[v_comm] += weight;
                        aggregate_graph.setEdgeWeight(u_comm, v_comm, weight);
                    }
                    aggregate_graph.incTotalEdgeWeight(weight); // update the total edge weight
                } 
            }
        }
    }

    //Rcpp::stop("Done aggregating graph");

    return aggregate_graph;
}

/*
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

    aggregated_graph.updateNodeProperties(true); // update node properties

    return aggregated_graph;
}
*/

void Optimizer::refinePartition(const Partition& P_original) {
    // for each community in the original partition
    for (const auto& entry : P_original.getCommunityIndexMap()) {
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
    std::unordered_map<std::string, int>  og_nodeIndexMap = G.getNodeIndexMap();  // store original node index map to track community assignments
    while (!done && iterations > counter) {
        Rcpp::Rcout << "Number of communities before moving: " << P.getCommunityIndexMap().size() << std::endl;
        bool improved = moveNodesFast();
        std::vector<int> community_indices = P.getCommunityIndices();
        // print the number of communities in the partition
        Rcpp::Rcout << "Number of communities after moving: " << community_indices.size() << std::endl;

        // if community size is equal to number of nodes, set done to true
        counter++;
        if (static_cast<int>(community_indices.size()) == G.getN() || !improved) {
            done = true;

        // otherwise, refine the partition and aggregate the graph
        } else {
            Partition P_save = P;  // maintain a copy of the partition before refinement

            Rcpp::Rcout << "Refining partition" << std::endl;
            // set this->P to refined partition
            //this->P = initializePartition(G);
            //refinePartition(P_save);

            // collapse the communities into a single node in a new graph
            Rcpp::Rcout << "Aggregating graph" << std::endl;
            this->G = aggregateGraph(P_save);  // aggregatoin based on the refined partition
            //G.updateNodeProperties(false);  // update node properties

            // update the community assignmens
            Rcpp::Rcout << "Updating community assignments" << std::endl;
            //Rcpp::stop("Need to update community assignments");
            updateCommunityAssignments(P, og_nodeIndexMap);
            // print the number of nodes in the aggregated graph
            Rcpp::Rcout << "Aggregated number of nodes: " << G.getN() << std::endl;

            P.makeSingleton(G);  // make the partition a singleton
        }

        if (!improved) {
            // print partition could not be improved
            Rcpp::Rcout << "Partition could not be improved!" << std::endl;
        }
    }
    //P.flattenPartition();
    // print flattenedd
    //Rcpp::Rcout << "Flattened partition" << std::endl;
}

Partition initializePartition(Graph& G) {
    std::vector<Community> communities;
    communities.reserve(G.getNodes().size());
    // Assign each node to its own community
    // for each node in the getNodes
    for (int node_index : G.getNodes()) {
        // Construct a community with a single node
        Community community({ node_index }, node_index); // Use the node index as the community index
        communities.push_back(community);
    }

    Partition P(communities, G);

    // Update the nodeCommunityMap in the Partition
    for (int node_index : G.getNodes()) {
        P.setNodeCommunity(node_index, node_index);
        //P.nodeCommunityMap[node_index] = node_index;
    }

    return P;
}




/*
###################################################################################
####################### R INTERFACE FUNCTIONS #####################################
*/

// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations, double gamma, double theta) {
    Rcpp::Rcout << "Running Leiden algorithm" << std::endl;

    //double gamma = .1;
    //double theta = 0.01;

    // Create a graph from the R List
    // use listToGraph from GraphUtils.cpp
    Graph G = listToGraph(graphList);
    // print the number of nodes 
    Rcpp::Rcout << "Number of nodes: " << G.getN() << std::endl;

    Rcpp::Rcout << "Initializing partition" << std::endl;
    Partition P = initializePartition(G);
    
    // Create an Optimizer
    Rcpp::Rcout << "Initializing optimizer" << std::endl;
    Optimizer optim(G, P, gamma, theta);
    
    // Run the Leiden algorithm
    optim.optimize(iterations);
    Rcpp::Rcout << "Final number of communities: " << optim.getP().getCommunityIndexMap().size() << std::endl;
    //Rcpp::Rcout << "Final quality: " << optim.P.calcQuality(gamma, G, false) << std::endl;

    // FIX CALC QUALITY NOT WORKING AFTER AGGREGATE?

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
    //return optim.G.graphToRList(optim.communityAssignments, optim.P.calcQuality(gamma, G, false));
    std::unordered_map<std::string, int> community_assignment = optim.getCommunityAssignments();
    //final_graph = optim.getG().graphToRList(community_assignment, optim.getP().calcQuality(gamma, G, false));
    
    G = listToGraph(graphList); // reset the graph nodes to get exact ones from before (not aggregated graph)
    optim.setG(G);
    return optim.graphToRList(community_assignment, 0.0);


}
