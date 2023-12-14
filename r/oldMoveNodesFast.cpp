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
                    if (end == node_queue.size()) { // if the queue is full (should really set back to 0 when this happens since it is circular but having weird issues with that)
                        node_queue.resize(node_queue.size() * 1.25);
                    }
                    // add to the end of the queue
                    node_queue[end++] = nn_idx;
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