#include <Rcpp.h>
#include "GraphUtils.h"

/*
##############################################
####### GRAPH CLASS FUNCTIONS ################
*/

/**
* @brief Construct a new Graph object
* @param n Number of nodes in the graph
* @return Graph : Graph object
**/
Graph::Graph(int n) : n(n) {
    isDirected = false; 
    // initialize edge weights map
    std::unordered_map<int, std::unordered_map<int, double>> edge_weights;
    for (int i = 0; i < n; ++i) {
      edge_weights[i] = std::unordered_map<int, double>();
    }
    this->edgeWeights = edge_weights;
    // initialize node weights map
    std::unordered_map<int, double> node_weights;
    // for now, all nodes have the same weight
    for (int i = 0; i < n; ++i) {   
      node_weights[i] = 1.0;
    }
    this->nodeWeights = node_weights;
    // initialize nodes
    std::vector<int> nodes_;
    for (int i = 0; i < n; ++i) {
      nodes_.push_back(i);
    }
    this->nodes = nodes_;
    // initialize node index map
    std::unordered_map<std::string, int> nodeIndexMap;
    this->nodeIndexMap = nodeIndexMap;

    this->totalEdgeWeight = 0.0;
};

/**
 * @brief add edge to the graph
 * @param u Node u
 * @param v Node v
 * @param w Weight of the edge
 * @return void
**/
void Graph::addEdge(const std::string& u, const std::string& v, double w) {
    int u_idx = getNodeIndex(u);
    int v_idx = getNodeIndex(v);
    // add the edge to the edge weights map
    edgeWeights[u_idx][v_idx] = w;
    
    // add other way too if undirected
    if (!isDirected) {
        edgeWeights[v_idx][u_idx] = w;
    }
    // update total edge weight
    totalEdgeWeight += w;
}

/**
 * @brief get node index from node name and throw an error if not found
 * @param node Node name
 * @return int : Node index
 * @throw Rcpp::stop if node not found
**/
int Graph::getNodeIndex(const std::string& node) const {
    if (nodeIndexMap.find(node) != nodeIndexMap.end()) {
        return nodeIndexMap.at(node);
    }
    Rcpp::stop("Node not found: " + node);
}

/**
 * @brief get node name from node index and throw an error if not found
 * @param index Node index
 * @return std::string : node name
 * @throw Rcpp::stop if index not found
**/
std::string Graph::getNodeName(int index) const {
    for (const auto& entry : nodeIndexMap) {
        if (entry.second == index) {
            return entry.first;
        }
    }
    Rcpp::stop("Index not found: " + std::to_string(index));
}

/**
 * @brief reset the node property by checking nodes in the index map
 * @return void
**/
void Graph::resetNodes() {
    Rcpp::Rcout << "Resetting nodes..." << std::endl;
    // print node index map
    std::vector<int> nodes_;
    for (const auto& entry : nodeIndexMap) {
        nodes_.push_back(entry.second);
    }
    setNodes(nodes_);
}

/**
 * @brief reset the node indices to force them to be consecutive and
 *        optionally remove empty nodes from the graph
 * @param remove_empty_nodes Whether to remove empty nodes from the graph
 * @return void
**/
void Graph::updateNodeProperties(bool remove_empty_nodes = false) {
    // before updating, make sure there are no empty nodes
    if (remove_empty_nodes) {
        std::vector<int> nodes_to_remove;
        for (auto& entry : edgeWeights) {
            if (entry.second.size() == 0) {
                nodes_to_remove.push_back(entry.first);
            }
        }
        for (int node : nodes_to_remove) {
            edgeWeights.erase(node);
            nodeIndexMap.erase(getNodeName(node));
            nodeWeights.erase(node);
        }
    }

    // reassign node indices to be consecutive
    int new_index = 0;
    std::unordered_map<int, std::unordered_map<int, double>> new_edge_weights;
    std::unordered_map<int, double> new_node_weights; 
    std::unordered_map<std::string, int> new_node_index_map;
    std::unordered_map<int, int> old_node_to_new_node; // map old node index to new node inde

    // create map of new node indices from old node indices
    for (auto& entry : nodeIndexMap) {
        // map the old index to the new index (old : new)
        old_node_to_new_node[entry.second] = new_index;
        new_index++;
    }

    // use the map to update node indices in the nodeIndexMap
    for (auto& entry : old_node_to_new_node) {
        new_node_weights[entry.second] = nodeWeights.at(entry.first); // update node weights
        new_node_index_map[getNodeName(entry.first)] = entry.second; // update node index map
        for (auto& edge : edgeWeights.at(entry.first)) {
            // only update the edge weights of the connected nodes that still exist in the graph
            if (old_node_to_new_node.find(edge.first) != old_node_to_new_node.end()) {
                new_edge_weights[entry.second][old_node_to_new_node.at(edge.first)] = edge.second; // update edge weights
            }
        }
    }

    // calculate total edge weight
    double total_edge_weight = 0.0;
    for (const auto& entry : edgeWeights) {
        for (const auto& edge : entry.second) {
            total_edge_weight += edge.second;
        }
    }
    if (!isDirected) {
        total_edge_weight /= 2; // In an undirected graph, each edge is counted twice 
    }

    // update the graph properties
    totalEdgeWeight = total_edge_weight;
    edgeWeights = new_edge_weights;
    nodeWeights = new_node_weights;
    nodeIndexMap = new_node_index_map;
    resetNodes();
    setN();
}

/**
 * @brief get the neighbors of a node
 * @param n_idx Node index
 * @return std::vector<int> : neighbors of the node
**/
std::vector<int> Graph::getNeighbors(int n_idx) const {
    std::vector<int> neighbors;
    for (const auto& to : edgeWeights.at(n_idx)) {
        neighbors.push_back(to.first);
    }
    return neighbors;
}

/**
 * @brief get the weight of an edge
 * @param u Node u
 * @param v Node v
 * @return double : weight of the edge
 * @throw Rcpp::stop if edge not found
**/
double Graph::getWeight(int u, int v) const {  // find is slow, consider not checking if the edge exists
    auto u_it = edgeWeights.find(u);
    if (u_it != edgeWeights.end()) {
        auto v_it = u_it->second.find(v);
        if (v_it != u_it->second.end()) {
            return v_it->second;
        }
    }
    // otherwise, throw an error
    Rcpp::stop("Edge not found: (" + std::to_string(u) + ", " + std::to_string(v) + ")");
}

/**
 * @brief convert the graph object to an R list
 * @param community_assignments Map of original node names to community indices
 * @param quality Quality of the partition
 * @return Rcpp::List : R list representation of the graph
 * @note This should be the last function called in the R interface
 * @note should add some sanity checks here
**/
Rcpp::List Graph::graphToRList(std::unordered_map<std::string, int>& community_assignments, double quality) const {
    Rcpp::CharacterVector from;
    Rcpp::CharacterVector to;
    Rcpp::NumericVector weight;

    Rcpp::Rcout << "Writing to R List" << std::endl;

    for (int i : getNodes()) {
        // print the node
        //Rcpp::Rcout << "Node: " << getNodeName(i) << std::endl;
        std::string fromNode = getNodeName(i);
        std::unordered_map<int, double> edges = getThisEdgeWeights(i);
        // print length edges
        //Rcpp::Rcout << "Edges: " << edges.size() << std::endl;
        for (const auto& entry : edges) {
            std::string toNode = getNodeName(entry.first);
            //Rcpp::Rcout << "Edge: " << fromNode << " -> " << toNode << " : " << entry.second << std::endl;
            from.push_back(fromNode);
            to.push_back(toNode);
            weight.push_back(entry.second);
        }
    }
    Rcpp::CharacterVector node_names;
    Rcpp::NumericVector communities;
    for (const auto& entry : community_assignments) {
        node_names.push_back(entry.first);
        communities.push_back(entry.second);
    }
    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("communities") = communities,
        Rcpp::Named("node_names") = node_names,
        Rcpp::Named("quality") = quality,
        Rcpp::Named("graph") = Rcpp::List::create(
            Rcpp::Named("from") = from,
            Rcpp::Named("to") = to,
            Rcpp::Named("weight") = weight
        )
    );

    return result;
}

/**
 * @brief remove nodes with less than 'min_connections' connections
 * @param min_connections Minimum number of connections a node must have to be kept
 * @return void
**/
void Graph::removeLowConnections(int min_connections) {
    std::vector<int> nodes_to_remove;
    double removed_edge_weight = 0.0;

    // Pass through all nodes
    for ( auto& node : nodes) {
        std::vector<int> neighbors = getNeighbors(node);

        // If this node has less than 'min_connections' connections, add it to the removal list
        if (neighbors.size() < min_connections) {
            for (auto& entry : edgeWeights.at(node)) {
                removed_edge_weight += entry.second;
            }
            nodes_to_remove.push_back(node);
        }
    }

    std::sort(nodes_to_remove.begin(), nodes_to_remove.end());
    // unique will return a pointer to the end of the unique elements, so remove them
    nodes_to_remove.erase(std::unique(nodes_to_remove.begin(), nodes_to_remove.end()), nodes_to_remove.end());

    // Remove nodes from the graph
    for (auto& node : nodes_to_remove) {
        for (auto& entry : edgeWeights.at(node)) {
            // Remove the edge from the other node if it exists
            if (edgeWeights.find(entry.first) != edgeWeights.end()) {
                edgeWeights.at(entry.first).erase(node);
            }
        }
        edgeWeights.erase(node);
        nodeIndexMap.erase(getNodeName(node));
        nodeWeights.erase(node);
    }

    // Recalculate totalEdgeWeight
    totalEdgeWeight -= removed_edge_weight;
    if(!isDirected) { 
        totalEdgeWeight /= 2; // In an undirected graph, each edge is counted twice 
    }


    // Lastly, always ensure to update the node properties
    updateNodeProperties(false);

}

/**
 * @brief check if the graph has an edge between two nodes
 * @param u Node u
 * @param v Node v
 * @return bool : whether the graph has an edge between the two nodes
**/
bool Graph::hasEdge(int u, int v) const {
    auto it = edgeWeights.find(u);
    if (it != edgeWeights.end()) {
        return it->second.find(v) != it->second.end();
    }
    return false;
    
}

/**
 * @brief create custom graph object from R list
 * @param graphList R list representation of the graph
 * @return Graph : custom graph object
 * @note should add some sanity checks here
*/
Graph listToGraph(const Rcpp::List& graphList) {
    Rcpp::CharacterVector from = graphList["from"];
    Rcpp::CharacterVector to = graphList["to"];
    Rcpp::NumericVector weight = graphList["weight"];

    std::set<std::string> uniqueNodes;

    // Collect unique nodes
    Rcpp::Rcout << "Initializing graph..." << std::endl;
    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);

        uniqueNodes.insert(fromNode);
        uniqueNodes.insert(toNode);
    }

    // Initialize the graph after obtaining unique nodes
    Graph G(uniqueNodes.size());

    // Populate nodeIndexMap
    int index = 0;
    for (const std::string& node : uniqueNodes) {
        G.setNodeIndex(node, index);
        index++;
    }

    Rcpp::Rcout << "Adding edges to graph..." << std::endl;

    // Add edges to the graph & populate the edgeWeights
    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);
        double w = weight[i];

        G.addEdge(fromNode, toNode, w);
    }

    Rcpp::Rcout << "Updating node properties..." << std::endl;

    G.updateNodeProperties(false); // Updating nodes after adding all edges 
    Rcpp::Rcout << "Removing Low Connections..." << std::endl;
    //G.removeLowConnections(2); // Remove single connections

    return G;
}

/*
###################################################################################
####################### R INTERFACE FUNCTIONS #####################################
*/

/**
 * @brief Call the listToGraph function from R
 * @param graphList R list representation of the graph
 * @return Rcpp::List : R list representation of the graph and information about the graph
**/
// [[Rcpp::export]]
Rcpp::List createGraphFromList(const Rcpp::List& graphList) {
    Graph G = listToGraph(graphList); // Convert R list
    double quality = 0.0;
    std::unordered_map<std::string, int> community_assignments;
    return G.graphToRList(community_assignments, quality);
}