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
    P_new.purgeEmptyCommunities(false);

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

    int n_unused_clusters = 0;
    int n_unstable_nodes = G.n;

    std::queue<int> node_queue;
    std::vector<int> random_nodes = RandomGenerator::generateRandomPermutation(G.n);

    for (int node_index : random_nodes) {
        node_queue.push(node_index);
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

    /*
    * Iterate over the nodeOrder array in a cyclical manner. When the end
    * of the array has been reached, start again from the beginning. The
    * queue of nodes that still need to be visited is given by
    * nodeOrder[i], ..., nodeOrder[i + nUnstableNodes - 1]. Continue
    * iterating until the queue is empty.
    */
    do {
        int j = node_queue.front();
        node_queue.pop();
        // if the node is stable, skip it
        // get current community of node j
        int c_idx = P.nodeCommunityMap.at(j);

        /*
        * Identify the neighboring clusters of the currently selected
        * node, that is, the clusters with which the currently selected
        * node is connected. An empty cluster is also included in the set
        * of neighboring clusters. In this way, it is always possible that
        * the currently selected node will be moved to an empty cluster.
        */

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
        // find an empty cluster, first check if only self is in the current cluster
        // if true, then the current cluster is empty
        if (P.communityIndexMap.at(c_idx).size() == 1) {
            // add the current cluster to the neighboring clusters
            neighboring_clusters.insert(c_idx);
        } else {
            // get first empty cluster
            int empty_cluster = unused_clusters[0];
            // add the empty cluster to the neighboring clusters
            neighboring_clusters.insert(empty_cluster);
        }

        /*
        * For each neighboring cluster of the currently selected node,
        * calculate the increment of the quality function obtained by
        * moving the currently selected node to the neighboring cluster.
        * Determine the neighboring cluster for which the increment of the
        * quality function is largest. The currently selected node will be
        * moved to this optimal cluster. In order to guarantee convergence
        * of the algorithm, if the old cluster of the currently selected
        * node is optimal but there are also other optimal clusters, the
        * currently selected node will be moved back to its old cluster.
        */

        // initialize the best cluster and the best quality
        int best_cluster = c_idx;
        double best_quality_increment = 0.0;
        // for each neighboring cluster
        for (int nc_idx : neighboring_clusters) {
            // calculate the quality of the move
            double delta_q = deltaQuality(j, nc_idx, gamma);
            // if the quality of the move is better than the best quality
            if (delta_q > best_quality_increment) {
                // update the best quality and best cluster
                best_quality_increment = delta_q;
                best_cluster = nc_idx;
            }
        }
        // move selected nodet its new cluster
        // update cluster stats
        P.updateCommunityMembership(j, c_idx, best_cluster);
        cluster_weights[best_cluster] += G.nodeWeights[j];
        nodes_per_cluster[best_cluster]++;
        if (best_cluster == unused_clusters[n_unused_clusters-1]) {
            n_unused_clusters--;
        }

        // mark the node as stable, remove it from the queue
        stable_nodes[j] = true;
        n_unstable_nodes--;

        /* if the new cluster of the currently selected node is different
        from the old cluster, some further updating of statistics is 
        performed. Also, neighbors of the currently selected node that do
        not belong to the new cluster are marked as unstable and added to
        the queue. 
        */
       if (best_cluster != c_idx) {
            P.nodeCommunityMap.at(j) = best_cluster;
            // do some extra managment if best cluster is the empty cluster?
            // get the neighbors of node j
            std::vector<int> neighbors = G.getNeighbors(j);
            // get the neighboring clusters of node j
            std::set<int> neighboring_clusters;
            for (int nn_idx : neighbors) {
                // get community of neighbor
                int nc_idx = P.nodeCommunityMap.at(nn_idx);
                // if the neighbor is stable and neighbor's cluster is not the best cluster
                if (nc_idx != best_cluster && stable_nodes[nc_idx] == true) {
                    stable_nodes[nc_idx] == false;
                    n_unstable_nodes++;
                    // add to node queue with this java logic
                    // add to the end of the queue
                    node_queue.push(nn_idx);
                }
            }
            update = true;
       }
       // get next node in the queue
    } while (n_unstable_nodes > 0);

    if (update) {
        // purge empty clusters
        P.purgeEmptyCommunities(false);

    }
    return update;
}

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

*/
void Optimizer::optimize() {
    // Implement the Leiden algorithm iteration here
    // This involve multiple steps, such as moveNodesFast, refinePartition, mergeNodesSubset, etc.
    bool done = false;
    while (!done) {
        moveNodesFast();
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

    /* OLD TESTING CODE
    move node 168 to community 126
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

    // Convert the vector of communities to an R List
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities,
        Rcpp::Named("nodes") = nodes
    );

    return result;

}
