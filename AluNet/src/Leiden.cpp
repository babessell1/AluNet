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

/**
 * @brief Construct a new Optimizer:: Optimizer object
 * @param G graph object
 * @param P partition object
 * @param gamma resolution parameter
 * @param theta temperature parameter
 * @param seed random seed (0 for random seed)
 * @note the community assignments are initialized to the node indices
**/
Optimizer::Optimizer(Graph& G, Partition& P, double gamma, double theta, int seed) : G(G), P(P), gamma(gamma), theta(theta), seed(seed) {
    this->communityAssignments = std::unordered_map<std::string, std::string>();
    this->og_communityAssignments = std::unordered_map<std::string, std::string>();
    // initialize the community assignments
    for (const auto& entry : G.getNodeIndexMap()) {
        std::string string_index = std::to_string(entry.second);
        communityAssignments.insert({entry.first, string_index});
        og_communityAssignments.insert({entry.first, string_index});
    }
    this->iteration = 0;
}

/**
 * @brief Update the community assignments based on the new partition
 * @param P partition object
 * @param original_nodeIndexMap original node index map
 * @note currently stored value in communityAssignments is the previous community assignment which is the latest node index
 * @note we update by matching new nodes to the previous community assignment and changing the community assignment to the new community index
 * @note this is critical for tracking the community assignments after the graph is aggregated
 *       which allows us to return the assigments of the original nodes after the optimization
**/
void Optimizer::updateCommunityAssignments(const Partition& P_comm, bool final_call) {
    // currently stored value in communityAssignments is the previous community assignment which is the latest node index
    // we update by matching new nodes to the previous community assignment and changing the community assignment to the new community index

    /// print nodecommunitymap of P_comm
    //Rcpp::Rcout << "Node community map of P_comm: " << std::endl;
    //for (const auto& entry : P_comm.getNodeCommunityMap()) {
        //Rcpp::Rcout << entry.first << " " << entry.second << std::endl;
    //}

    int count_replaced = 0;
    int count_nodes = 0;
    // for each node in the original community assignments
    std::unordered_map<std::string, std::string> new_community_assignments;
    for (const auto& iter : getOGCommunityAssignments()) {
        count_nodes++;
        std::string node_name = iter.first;
        // get the previous assignment
        std::string prev_assignment = iter.second;

        //Rcpp::Rcout << "Node: " << node_name << " Previous assignment: " << prev_assignment << std::endl;

        // if previous assignment is the same as the new node index and using refined partition
        if (!final_call) {
            // get the new community assignment
            for (const auto& entry : getCommunityAssignments()) {
                std::string agg_node = entry.first;
                std::string new_assignment = entry.second;

                // if the aggregated node index matches the previous assignment
                if (agg_node == prev_assignment) {
                    // replace the previous assignment with the new assignment
                    new_community_assignments[node_name] = new_assignment;
                    count_replaced++;
                    break;
                }             
            }
        } else {
            std::string new_assignment = std::to_string(P_comm.getNodeCommunityIdx(std::stoi(prev_assignment)));
            new_community_assignments[node_name] = new_assignment;
            count_replaced++;
        }
    }
    // if not all nodes were replaced, throw an error
    if (count_replaced != count_nodes) {
        Rcpp::stop("Not all nodes were replaced in the community assignments, something went wrong!");
    }
    // Replace old community assignments with the new ones
    setOGCommunityAssignments(new_community_assignments);
}

/**
 * @brief Calculate the difference in quality of the partition after moving a node to a new community
 * @param n_idx node index
 * @param new_c_idx new community index
 * @param gamma resolution parameter
 * @param recalculate whether to recalculate the quality of the partition
 * @return double : difference in quality
 * @note calculate the difference between moving vertex to a new community
 * @note recalculating is slower, but necessary if the partition quality has changed since setting the quality last
**/
double Optimizer::deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const {
    int old_c_idx = P.getNodeCommunityIdx(n_idx);  // keep track of the old community

    if (recalculate) {
      // get quality of the partition before the move
      double old_quality = P.calcQuality(gamma, G);
      P.setQuality(old_quality);
    }

    // make copy of the partition
    Partition P_new = P;

    // move the node to the new community
    P_new.updateCommunityMembership(n_idx, old_c_idx, new_c_idx);
    P_new.setNodeCommunity(n_idx, new_c_idx);

    // get quality of the partition after the move
    double new_quality = P_new.calcQuality(gamma, G);

    return new_quality - P.getQuality();
}

/**
 * @brief Move nodes around the partition to improve the quality to be as high as possible
 * @return bool : true if the partition has changed, false otherwise
 * @note quality based on Constant Potts Model
 * @note first step of the Leiden algorithm in the paper
 * @note the algorithm is based on the paper "From Louvain to Leiden: guaranteeing well-connected communities" by Traag et al.
 * @note we initialize our own data structures with map properties into vectors for faster access and iteration
 * @note the only time we change our actual data structures is when we move a node to a new community
 * @note this also allows us to use a faster calculation of delta quality instead of our method in Parition
*/
bool Optimizer::moveNodesFast() {
    bool update = false;  // track if the partition has changed
    bool print = false; // debug print flag
    int n_ = G.getN();  // number of nodes in the graph

    std::vector<bool> stable_nodes(n_, false); // track which nodes have not moved
    std::vector<double> cluster_weights(n_, 0.0);  // track the sum of node weights per cluster
    std::vector<double> edge_weights_per_cluster(n_, 0.0);  // track the sum of edge weights per cluster
    std::vector<int> nodes_per_cluster(n_, 0);  // track the number of nodes in each cluster
    std::vector<int> unused_clusters(n_ - 1, 0);  // track which clusters are empty
    std::vector<int> n_neighboring_clusters(n_, 0);  // track the number of neighboring clusters for each cluster

    int n_unused_clusters = 0;
    int n_unstable_nodes = n_;

    std::vector<int> random_nodes = G.getNodes();
    if (seed == 0) {
        std::vector<int> random_nodes = RandomGenerator::generateRandomPermutation(n_);
    } else {
        std::vector<int> random_nodes = RandomGenerator::generateRandomPermutation(n_, seed);
    }

    std::queue<int> node_queue;
    size_t start = 0, end = n_;
    for (int i = 0; i < n_; i++) {
        node_queue.push(random_nodes[i]);
    }

    // initialize the cluster weights and nodes per cluster
    for (int i = 0; i < n_; i++) {
        int comm_idx = P.getNodeCommunityIdx(i);
        cluster_weights[comm_idx] += G.getNodeWeights()[i];
        nodes_per_cluster[comm_idx]++;
    }

    // track unused clusters
    for (int i = n_ - 1; i >= 0; i--) {
        if (nodes_per_cluster[i] == 0) {
            unused_clusters[n_unused_clusters] = i;
            n_unused_clusters++;
        }
    }

    if (print) {
        Rcpp::Rcout << "Length of unused clusters: " << n_unused_clusters << std::endl;
        Rcpp::Rcout << "Length of nodes per cluster: " << nodes_per_cluster.size() << std::endl;
    }

    // start of the cyclical process
    do {
        // get the next node in the queue and its community
        int node_idx = node_queue.front();
        node_queue.pop(); 
        int cluster_idx = P.getNodeCommunityIdx(node_idx);

        if (print) {
            Rcpp::Rcout << "Queued Node " << node_idx << " is in cluster " << cluster_idx << std::endl;
        }

        // decrement the cluster weight of the node's community 
        // (we assume we will move the node to a new community)
        cluster_weights[cluster_idx] -= G.getNodeWeight(node_idx);
        nodes_per_cluster[cluster_idx]--;

        // if the cluster is empty, add it to the unused clusters
        if (nodes_per_cluster[cluster_idx] == 0) {
            unused_clusters[n_unused_clusters] = cluster_idx;
            n_unused_clusters++;
        }

        // initialize the vector of neighboring clusters
        std::vector<int> neighboring_clusters(n_);

        // add an empty cluster to the neighboring clusters
        neighboring_clusters[0] = unused_clusters[n_unused_clusters - 1]; 

        // add the neighboring clusters and the edge/node weights per cluster
        int n_neighboring_clusters = 1;  // since we added an empty cluster already
        for (const auto& neighbor_weight : G.getThisEdgeWeights(node_idx)) {

            if (print) {
                Rcpp::Rcout << "Neighbor " << neighbor_weight.first << " has edge weight " << neighbor_weight.second << std::endl;
            }

            // get the neighbor's community index
            int neighbor_cluster = P.getNodeCommunityIdx(neighbor_weight.first);

            // if the neighbor's cluster is not empty, 
            if (edge_weights_per_cluster[neighbor_cluster]==0) {

                // add the neighbor's cluster to the neighboring clusters
                neighboring_clusters[n_neighboring_clusters] = neighbor_cluster;
                n_neighboring_clusters++;  // increment the number of neighboring clusters
            }
            // add the edge weight to the neighbor's community weight
            edge_weights_per_cluster[neighbor_cluster] += neighbor_weight.second;
        }

        // initialize the best cluster as the current cluster
        int best_cluster = cluster_idx;

        // initialize the max delta quality as the current delta quality
        // note that the delta quality calculated in this fast way will 
        // not necessarily be zero!
        double max_delta_q = edge_weights_per_cluster[cluster_idx] - G.getNodeWeight(node_idx) * cluster_weights[cluster_idx] * gamma;

        if (print) {
            Rcpp::Rcout << "Max delta quality: " << max_delta_q << std::endl;
        }

        // for each neighboring cluster
        for (int i = 0; i < n_neighboring_clusters; i++) {
            int idx = neighboring_clusters[i]; // get the neighboring cluster index

            // E(C', C') - uweight * C'weight * 𝛾 
            double delta_quality = edge_weights_per_cluster[idx] - G.getNodeWeight(node_idx) * cluster_weights[idx] * gamma;
            
            if (print) {
                Rcpp::Rcout << "Neighbor cluster " << idx << " has delta quality " << delta_quality << std::endl;
            }

            // set best quality and cluster if better than current best
            if (delta_quality > max_delta_q) {
                best_cluster = idx;
                max_delta_q = delta_quality;
            }

            // reset the edge weight
            edge_weights_per_cluster[idx] = 0;
        }

        // move the node to the best cluster
        cluster_weights[best_cluster] += G.getNodeWeight(node_idx); 
        nodes_per_cluster[best_cluster]++; // increment the number of nodes in the best cluster

        bool best_is_empty = false;

        // If the best cluster is the empty cluster, decrement the number of unused clusters
        if (best_cluster == unused_clusters[n_unused_clusters - 1]) {
            n_unused_clusters--;
            best_is_empty = true;
            if (print) {
                Rcpp::Rcout << "Best cluster is empty!" << std::endl;
            }
        }

        // mark the node as stable
        stable_nodes[node_idx] = true;
        n_unstable_nodes--;

        // check if cluster is different from the current cluster
        if (best_cluster != cluster_idx) {

            if (print) {
                Rcpp::Rcout << "Moving node " << node_idx << " from cluster " << cluster_idx << " to cluster " << best_cluster << std::endl;
            }

            // if best cluster is empty and not the same as the current cluster,
            // make sure there the empty cluster is in the partition to move to
            if (best_is_empty) {
                // if the empty cluster is not in the partition, add it
                if (!P.hasCommunity(best_cluster)) {
                    Community new_empty_cluster = Community({}, best_cluster);
                    P.addCommunity(new_empty_cluster);
                }
            }

            // update the partition, flag the partition as updated
            P.updateCommunityMembership(node_idx, cluster_idx, best_cluster);
            update = true;

            if (print) {
                Rcpp::Rcout << "Node " << node_idx << " is now in cluster " << P.getNodeCommunityIdx(node_idx) << std::endl;
            }
 
            // check each neighbor of the moved node
            for (const auto& neighbor_weight : G.getThisEdgeWeights(node_idx)) {
                int neighbor_node = neighbor_weight.first; // index of neighbor  

                if (print) {
                    Rcpp::Rcout << "Checking neighbor " << neighbor_node << std::endl;
                }  

                // add neighbors to the queue if they are not in the new best cluster    
                if (stable_nodes[neighbor_node] && P.getNodeCommunityIdx(neighbor_node) != best_cluster) {
                    stable_nodes[neighbor_node] = false;
                    n_unstable_nodes++; // increment the number of unstable nodes

                    // add the neighbor to the queue
                    node_queue.push(neighbor_node);
                }
            }
        }

    // keep iterating until there are no unstable nodes
    } while (n_unstable_nodes > 0);

    if (print) {
        Rcpp::Rcout << "print dat bitch" << std::endl;
    }

    // remove all empty clusters
    P.purgeEmptyCommunities(false);

    return update;
}

/**
 * @brief Get well connected communities in the subset
 * @param B subset of the partition
 * @return std::vector<Community> : well connected communities
 * @note Identify communities in the subset that are well connected, that is:
 *       the edge weights of C, and the difference of the subset and C is greater than gamma * size of flat(C) * size of flat(S - C)
 *       get the sizes
**/
std::vector<Community> Optimizer::getWellConnectedCommunities(const Community& B) const {
    int size_B = B.getNodeIndices().size();
    std::vector<Community> well_connected_communities; // well connected communities

    // for each community in the refined partition
    for (const auto& entry : P.getCommunityIndexMap()) {

        // for each C in P
        int size_C = entry.second.getNodeIndices().size();
        Community C = entry.second;
        size_t count_c_in_b = 0; // count the number of nodes in C that are in B
        double edge_weight_B_A = 0.0;
        // for each node in B (subset)
        for (const int& b : B.getNodeIndices()) {
            // for each node in C (community of interest)
            for (const int& c : C.getNodeIndices()) { 
                if (b == c) {  // if node is in both B and C
                    count_c_in_b++; // increment count of nodes in C that are also in B
                } else if(G.hasEdge(b, c)) {  // if node is in B but not C
                    edge_weight_B_A += G.getWeight(b, c); // add the edge weight between B and C
                }
            }
        }
        // check if C is contained by B
        bool is_contained = count_c_in_b == C.getNodeIndices().size();

        //
        if (!is_contained) { //? was this wrong before or was my comment wrong?
            continue;  // continue if C is not contained by B
        }
        // if well connected, add to R
        if (edge_weight_B_A >= gamma * size_C * (size_B - size_C)) {
            well_connected_communities.push_back(entry.second); // add the community to R
        }
    }

    return well_connected_communities;
}

/**
 * @brief Get well connected nodes in the subset
 * @param B subset of the partition
 * @return std::vector<int> : well connected nodes
 * @note Identify nodes in the subset that are well connected, that is:
 *       the edge weights of v, and the difference of the subset and v is greater than gamma * size of flat(v) * size of flat(subset - v)
 *       note that v is a node so flat(v) is just 1
**/
std::vector<int> Optimizer::getWellConnectedNodes(const Community& B) const {
    std::vector<int> well_connected_nodes; // R
    int size_B = B.getNodeIndices().size();

    // Determine well-connected nodes
    for (const int& v : B.getNodeIndices()) {
        // Retrieve precomputed node weight from the graph's nodeWeights map
        double edge_weight_B_v = G.getNodeWeight(v);

        // if well connected, add to R
        if (edge_weight_B_v >= gamma * (size_B - 1)) {
            well_connected_nodes.push_back(v); // add the node to R
        }
    }

    return well_connected_nodes;
}

/**
 * @brief Merge nodes in the subset to well connected communities
 * @param S subset of the partition
 * @return void
 * @note see the paper "From Louvain to Leiden: guaranteeing well-connected communities" by Traag et al.
**/
void Optimizer::mergeNodesSubset(Community& S) {
    std::vector<int> R = getWellConnectedNodes(S); // get well connected nodes
    std::vector<Community> T = getWellConnectedCommunities(S); // get well connected communities

    // Visit nodes in random order
    if (seed != 0) {
        std::shuffle(R.begin(), R.end(), std::default_random_engine(seed));
    } else {
        std::shuffle(R.begin(), R.end(), std::default_random_engine(std::random_device()()));
    }
    for (const int& v : R) {
        // Consider only nodes that have not yet been merged
        if (P.inSingleton(v)) { // Check if v is in its own community (a singleton)

            // map community index to delta quality based probability
            std::unordered_map<int, double> log_probabilities;

            double best_delta_quality = P.calcQuality(gamma, G);
            P.setQuality(best_delta_quality);
            int best_community_index = -1;

            // for each well connected community contained in S
            for (const auto& C : T) {
                int c_idx = C.getCommunityIndex();

                // E(C', C') - uweight * C'weight * 𝛾 
                double delta_quality = C.getClusterEdgeWeights(G) - G.getNodeWeight(v) * C.getClusterWeight(G) * gamma;

                if (delta_quality > 0) {
                    log_probabilities[c_idx] = (delta_quality/theta);
                }

            }

            // If a better community was found, move node v to community C'
            if (!log_probabilities.empty()) {

                std::random_device rd;
                std::mt19937 gen(rd());  // merenne twister engine
                if (seed != 0) {
                    std::mt19937 gen(seed);
                }
                std::uniform_real_distribution<> dis(0, 1);

                int C_prime_index = -1; // initialize the new community index to nonsense value (will throw error if used)
         
                auto max_pair = *std::max_element(log_probabilities.begin(), log_probabilities.end(),
                    [](const std::pair<const int, double>& a, const std::pair<const int, double>& b) {
                        return a.second < b.second;
                    }
                );

                double max_log_prob = max_pair.second;
                double sum_exp = 0.0;

                // Choose a random community C' with probability proportional to probabilities
                for (const auto& entry : log_probabilities) {
                    double prob = exp(entry.second - max_log_prob);
                    sum_exp += prob;
                }
                double log_sum_exp = log(sum_exp) + max_log_prob; 

                for (const auto& entry : log_probabilities) {
                    double prob = exp(entry.second - log_sum_exp);
                    double roll = dis(gen);
                    if (roll < prob) {
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

/**
 * @brief Aggregate the graph based on the partition by collapsing the nodes in each 
 *        communities into a single node in a new graph. The collapsed communities
 *        are based on the refined partition, and the community assignments are based
 *        on the original partition
 * @param P_comm partition object (communities based on this partition)
 * @return Graph : aggregated graph
 * @note See the paper "From Louvain to Leiden: guaranteeing well-connected communities" by Traag et al.
 * @note P defines the aggregate graph
 * @note P_comm defines the communities in the aggregate graph
*/
void Optimizer::aggregateGraph(Partition& P_comm) {

    //reset the community assignments
    resetCommunityAssignments();

    Graph aggregate_graph;

    aggregate_graph.setIsDirected(false);

    // create a map of community indices to ordered new community indices
    // want to ensure the nodes are consecutively indexed for the 
    // next iteration of moveNodesFast so that we can
    // use the vector data structure
    std::unordered_map<int, int> new_community_indices; 
    int new_index = 0;
    for (const auto& entry : P.getCommunityIndexMap()) {
        int old_index = entry.first;
        new_community_indices[old_index] = new_index;
        aggregate_graph.addNodeName(std::to_string(new_index), new_index);
        new_index++;
    }

    // create a map of community indices to ordered new community indices (based on P_comm)
    std::unordered_map<int, int> old_community_indices; 
    int new_comm_assign_idx = 0;
    for (const auto& entry : P_comm.getCommunityIndexMap()) {
        int old_assign_index = entry.first;
        old_community_indices[old_assign_index] = new_comm_assign_idx;
        new_comm_assign_idx++;
    }

    // initalize community assignment maps
    std::unordered_map<int, Community> agg_community_index_map = P_comm.createEmptyCommunityIndexMap(); // communities based on moveNodes
    std::unordered_map<int, int> agg_node_community_map;
    std::unordered_map<std::string, int> agg_node_index_map;

    // initialize counter for the new community indices
    int new_community_assignment = -1;
    std::vector<int> used_communities;

    // for each community in the partition
    for ( auto iter : P.getCommunityIndexMap() ) {
        int old_c_idx = iter.first; // get the old community index
        int c_idx = new_community_indices.at(old_c_idx);  // get the new (consecutive) community index

        // convert the community index to a string (new node name in agg graph)
        std::string c_str = aggregate_graph.getNodeName(c_idx);

        // initialize new node in the aggregate graph
        aggregate_graph.setNodeIndex(c_str, c_idx);
        aggregate_graph.addNode(c_idx);
        aggregate_graph.setNodeWeight(c_idx, 0);
        aggregate_graph.initEdgeWeight(c_idx);

        // get the nodes in the old community
        std::vector<int> nodesInCommunity = P.getCommunity(old_c_idx).getNodeIndices();

        // initialize the community the old nodes belong to
        int old_community = -1;

        // for each node in the community
        for (int node : nodesInCommunity) {

            // get the community index (just a rename for clarity)
            int u_comm = c_idx;

            // update community assignment
            updateCommunityAssignment(std::to_string(node), c_str);
    
            // add the node weight to the aggregate graph
            aggregate_graph.incNodeWeight(u_comm, G.getNodeWeight(node));

            // for each neighbor of the node
            for (const auto& neighbour : G.getThisEdgeWeights(node)) {
                int old_v_idx = neighbour.first;  // get the neighbor index
                int v_old_comm = P.getNodeCommunityIdx(old_v_idx);  // get the neighbor community index
                int v_comm = new_community_indices.at(v_old_comm);  // get the new (consecutive) neighbor community index
                std::string v_str = aggregate_graph.getNodeName(v_comm);  // get the new (consecutive) neighbor index as a string

                // ensure the neighbor community is in the node index map
                aggregate_graph.setNodeIndex(v_str, v_comm);

                // get the edge weight shared between the node and the neighbor
                double weight = neighbour.second;

                 // don't want self loops
                if (u_comm != v_comm) {

                    // check if aggregate_graph.getEdgeWeights().at(u_comm) is empty
                    if (aggregate_graph.getEdgeWeightsFromStr(c_str).empty()) {
                        // add the neighbor to the map
                        aggregate_graph.addEdge(c_str, v_str, weight);
                    } else {
                        // update the edge weight
                        aggregate_graph.updateAggEdgeWeight(c_str, v_str, weight); // might not ever get called
                    }
                } 
            }
            // identify the community the collapsed node belong to
            int next_old_community = P_comm.getNodeCommunityIdx(node);
            if (old_community != -1) {
                if (next_old_community != old_community) {
                    Rcpp::Rcout << "New nodes best community is: " << next_old_community << std::endl;
                    Rcpp::Rcout << "but previous node in same refined community was: " << old_community << std::endl;
                    Rcpp::stop("Something went wrong in the refinement step!");
                }
            } else {
                old_community = next_old_community;
            }
        }

        // get community assigment index
        int comm_assign_idx = old_community_indices.at(old_community);

        // check if the old community is already in the used communities
        if (std::find(used_communities.begin(), used_communities.end(), old_community) == used_communities.end()) {
           
            // add the old community to the used communities
            used_communities.push_back(old_community);

            // create a new community in the aggregate graph
            Community new_community({c_idx}, comm_assign_idx);
            agg_community_index_map[comm_assign_idx] = new_community;

        } else { // if the old community is already in the used communities
            // add the node to the existing community
            agg_community_index_map.at(comm_assign_idx).addNode(c_idx);
        }
        // update the communitynew_community_assignment
        agg_node_community_map[c_idx] = comm_assign_idx;
        agg_node_index_map[c_str] = c_idx; // add the new node index to the map
    }

    // check that the aggregate graph is undirected
    if (aggregate_graph.getIsDirected()) {
        for (const auto& entry : aggregate_graph.getEdgeWeights()) {  // for each node in the aggregate graph
            for (const auto& edge : entry.second) {  // for each edge in the node
                // if there is no edge weight for the edge
                if (aggregate_graph.getEdgeWeights().at(edge.first).find(entry.first) == aggregate_graph.getEdgeWeights().at(edge.first).end()) {
                    Rcpp::Rcout << "Edge weight for edge " << edge.first << " " << entry.first << " is missing!" << std::endl;
                    Rcpp::stop("Aggregate graph is directed!");
                }
            }
        }
    }

    // set the aggregate graph properties
    P.setCommunityIndexMap(agg_community_index_map);
    P.setNodeCommunityMap(agg_node_community_map);
    aggregate_graph.setNodeIndexMap(agg_node_index_map);
    aggregate_graph.setN(); // set the number of nodes in the graph
    this->G = aggregate_graph; // set the graph to the aggregate graph
}

/**
 * @brief Refine the partition by merging nodes in each community to well connected communities
 * @param P_original original partition
 * @return void
 * @note see the paper "From Louvain to Leiden: guaranteeing well-connected communities" by Traag et al.
*/
void Optimizer::refinePartition(const Partition& P_original) {

    // for each community in the original partition
    for (const auto& entry : P_original.getCommunityIndexMap()) {
        Community C = entry.second; // get the community
        // print the community index and the number of nodes in the community
        mergeNodesSubset(C); // merge nodes in the community to well connected communities
    }
    // print the number of communities in the refined partition
    P.purgeEmptyCommunities(true); // remove empty communities
}

/**
 * @brief Main function to optimize the partition and run leiden iterations
 * @param iterations max number of iterations to run
 * @return void
 * @note see the paper "From Louvain to Leiden: guaranteeing well-connected communities" by Traag et al.
**/
void Optimizer::optimize(int iterations) {
    bool done = false; // track if the optimization is done
    int counter = 0; // track the number of iterations

    // if not converged and not at max iterations, run a Leiden iteration
    while (!done && iterations > counter) {
        counter++; // increment the iteration counter
        setIteration(counter); // set the iteration number

        //Rcpp::Rcout << "========================= ITERATION " << getIteration() << " =========================" << std::endl;
        //Rcpp::Rcout << "Number of communities before moving: " << P.getCommunityIndexMap().size() << std::endl;

        // fast local node moving step to improve the partition
        bool improved = moveNodesFast();

        std::vector<int> community_indices = P.getCommunityIndices();
        //Rcpp::Rcout << "Number of communities after moving: " << community_indices.size() << std::endl;

        // if move move node step cannot improve the partition, we are done
        if (static_cast<int>(community_indices.size()) == G.getN() || !improved) {
            done = true;

        // otherwise, refine the partition and aggregate the graph
        } else {

             // maintain a copy of the partition before refinement
            Partition P_save = P;

            // make the partition a singleton
            //Rcpp::Rcout << "Refining partition" << std::endl;
            this->P.makeSingleton(G);

            // refine the partition into a more granular and slightly random partition
            refinePartition(P_save);

            // collapse the communities into a single node in a new graph
            //Rcpp::Rcout << "Aggregating graph" << std::endl;
            aggregateGraph(P_save);
            //Rcpp::Rcout << "Aggregated number of nodes: " << G.getN() << std::endl;

            // track the community assignments for each of the original nodes
            //Rcpp::Rcout << "Updating community assignments" << std::endl;
            updateCommunityAssignments(P, false);
        }

        if (!improved) {
            // print partition could not be improved
            Rcpp::Rcout << "Partition converged after " << counter-1 << " iteration(s)" << std::endl;
        }

    }

     // update community assignments for the output
    //Rcpp::Rcout << "Updating final community assignments" << std::endl;
    updateCommunityAssignments(P, true);

}

/**
 * @brief Initialize a partition with each node in its own community
 * @param G graph object
 * @return Partition : partition object
 * @note not in the optimizer class because we want to be able to track multiple partitions 
*/
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
    }

    return P;
}

/*
###################################################################################
####################### R INTERFACE FUNCTIONS #####################################
*/

//' Takes an R List, converts it to a Graph object, and runs the Leiden algorithm
//' @param graphList R List containing the graph
//' @param iterations number of iterations to run the Leiden algorithm
//' @param gamma resolution parameter
//' @param theta temperature parameter
//' @return Rcpp::List : R List containing the graph, community assignments, and quality
//' @example runLeiden(graphList, 10, 0.1, 0.01)
//' @export
// [[Rcpp::export]]
Rcpp::List runLeiden(Rcpp::List graphList, int iterations, double gamma, double theta, int seed = 0) {
    //Rcpp::Rcout << "Running Leiden algorithm" << std::endl;

    // Create a graph from the R List
    Graph G = listToGraph(graphList);
    // print the number of nodes 
    //Rcpp::Rcout << "Number of nodes: " << G.getN() << std::endl;

    //Rcpp::Rcout << "Initializing partition" << std::endl;
    Partition P = initializePartition(G);
    
    // Create an Optimizer
    //Rcpp::Rcout << "Initializing optimizer" << std::endl;
    Optimizer optim(G, P, gamma, theta, seed);
    
    // Run the Leiden algorithm
    optim.optimize(iterations);
    //Rcpp::Rcout << "Final number of communities: " << optim.getP().getCommunityIndexMap().size() << std::endl;

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

    std::unordered_map<std::string, int> community_final_assignment = optim.getFinalCommunityAssignmentsInt();

    // reset the graph nodes to get exact ones from before (not aggregated graph)
    G = listToGraph(graphList);
    optim.setG(G);

    return optim.graphToRList(community_final_assignment, 0.0);
}
