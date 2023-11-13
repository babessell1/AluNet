
#include "GraphUtils.h"

/*
##############################################
####### GRAPH CLASS FUNCTIONS ################
*/
Graph::Graph(int n) : n(n) {
    isDirected = false;
    possibleEdges = n * (n - 1) / 2;
    // initailize adjacency map
    std::unordered_map<int, std::vector<int>> adj;
    for (int i = 0; i < n; ++i) {
      adj[i] = std::vector<int>();
    }
    this->adj = adj;
    // initialize edge weights map
    std::unordered_map<int, std::vector<double>> edge_weights;
    for (int i = 0; i < n; ++i) {
      edge_weights[i] = std::vector<double>();
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
    std::vector<int> nodes;
    for (int i = 0; i < n; ++i) {
      nodes.push_back(i);
    }
    this->nodes = nodes;
    // initialize node index map
    std::map<std::string, int> nodeIndexMap;
    this->nodeIndexMap = nodeIndexMap;

    // calculate total edge weight
    totalEdgeWeight = 0.0;
    for (int i = 0; i < n; ++i) {
        for (double weight : edgeWeights[i]) {
            totalEdgeWeight += weight;
        }
    }
    this->totalEdgeWeight = totalEdgeWeight;
};

// add new edge to the graph, give it the two node names to connect and the weight
void Graph::addEdge(const std::string& u, const std::string& v, double w) {
    int u_idx = getNodeIndex(u);
    int v_idx = getNodeIndex(v);
    adj[u_idx].push_back(v_idx);
    edgeWeights[u_idx].push_back(w);

    if (!isDirected) { // check if graph is undirected.
        adj[v_idx].push_back(u_idx);
        edgeWeights[v_idx].push_back(w);
    }
}

// get node index from node name
int Graph::getNodeIndex(const std::string& node) const {
    if (nodeIndexMap.find(node) != nodeIndexMap.end()) {
        return nodeIndexMap.at(node);
    }
    Rcpp::stop("Node not found: " + node);
}

// get node name from node index
std::string Graph::getNodeName(int index) const {
    for (const auto& entry : nodeIndexMap) {
        if (entry.second == index) {
            return entry.first;
        }
    }
    Rcpp::stop("Index not found: " + std::to_string(index));
}

// get all node indices in the graph
std::vector<int> Graph::getNodes() const {
    std::vector<int> nodes;
    for (const auto& entry : nodeIndexMap) {
        nodes.push_back(entry.second);
    }
    return nodes;
}

void Graph::updateNodeProperties() {
    // update nodes
    nodes = getNodes();
}

// get neighbors of a node based on index
std::vector<int> Graph::getNeighbors(int n_idx) const {
    return adj.at(n_idx);
}

// get weight of an edge based on node indices
double Graph::getWeight(int u, int v) const {
    const auto& u_edges = adj.at(u);
    const auto& u_weights = edgeWeights.at(u);
    
    for (std::vector<int>::size_type i = 0; i < u_edges.size(); i++) {
        if (u_edges[i] == v) {
            return u_weights[i];
        }
    }
    Rcpp::stop("Edge not found: " + std::to_string(u) + " " + std::to_string(v));
}
// convert custom graph object to R list
Rcpp::List Graph::graphToRList() const {
    Rcpp::CharacterVector from;
    Rcpp::CharacterVector to;
    Rcpp::NumericVector weight;

    Rcpp::Rcout << "Writing to R List" << std::endl;

    for (int i : nodes) { 
        std::string fromNode = getNodeName(i);
        for (size_t j = 0; j < adj.at(i).size(); j++) {
            std::string toNode = getNodeName(adj.at(i)[j]);
            from.push_back(fromNode);
            to.push_back(toNode);
            weight.push_back(edgeWeights.at(i)[j]);
        }
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("from") = from,
        Rcpp::Named("to") = to,
        Rcpp::Named("weight") = weight
    );

    return result;
}

// trim the graph by removing nodes with only one connection
// might consider defining a minimum number of connections to keep
void Graph::removeSingleConnections() {
    // print information about nodes and edges before removal
    Rcpp::Rcout << "Nodes before removal: " << n << std::endl;

    // track total removed edge weight
    double removed_edge_weight = 0.0;

    std::vector<int> to_remove;
    for (auto it = adj.begin(); it != adj.end(); /* no increment here */) {
        // if a node has less than 2 outbound edges, mark it for removal
        if (it->second.size() < 2) {
            to_remove.push_back(it->first);
            // add the weight of the edges to be removed to the total removed edge weight
            for (double weight : edgeWeights[it->first]) {
                removed_edge_weight += weight;
            }
            it = adj.erase(it); // remove node with less than 2 edges
        } else {
            ++it;
        }
    }

    // create new maps and vectors using new indices
    std::map<std::string, int> new_node_index_map;
    std::unordered_map<int, std::vector<int>> new_adj;
    std::unordered_map<int, std::vector<double>> new_edge_weights;
    std::unordered_map<int, double> new_node_weights;
    std::vector<int> new_nodes;

    // Populate new adjacency list and other maps
    for (const auto& entry : adj) {
        int old_index = entry.first;
        int new_index = new_node_index_map.size();
        new_node_index_map[getNodeName(old_index)] = new_index;
        new_nodes.push_back(new_index);

        // Update adjacency list based on new indices
        auto& neighbors = entry.second;
        // update neighbors based on new indices
        std::vector<int> new_neighbors;
        for (int neighbor : neighbors) {
            auto it = new_node_index_map.find(getNodeName(neighbor));
            if (it != new_node_index_map.end()) {
                // update neighbors based on new indices
                new_neighbors.push_back(it->second);
            } else {
                // Handle the case where the neighbor is not found in the map
                // This may happen if the neighbor was removed
                // You can choose to ignore it or handle it based on your requirements
                // For now, I'll just pass
                //Rcpp::Rcout << "Warning: Neighbor not found in new_node_index_map: " << neighbor << std::endl;
            }

        }
        
        new_adj[new_index] = new_neighbors;
        new_node_weights[new_index] = nodeWeights[old_index];
        new_edge_weights[new_index] = edgeWeights[old_index];
    }

    // print information about new structure
    Rcpp::Rcout << "New nodes size: " << new_nodes.size() << std::endl;

    // reassign members with the new maps and vectors
    nodeIndexMap = new_node_index_map;
    nodes = new_nodes;
    adj = new_adj;
    edgeWeights = new_edge_weights;
    nodeWeights = new_node_weights;
    n = nodes.size();

    // update total edge weight
    totalEdgeWeight -= removed_edge_weight;

    // check that the set of nodes is the same as the set of keys in the adjacency list
    std::set<int> node_set(nodes.begin(), nodes.end());
    std::set<int> adj_set;
    for (const auto& entry : adj) {
        adj_set.insert(entry.first);
    }
    if (node_set != adj_set) {
        Rcpp::stop("Nodes and adjacency list do not match");
    }
    // do the same for edge weights
    std::set<int> edge_set;
    for (const auto& entry : edgeWeights) {
        edge_set.insert(entry.first);
    }
    if (node_set != edge_set) {
        Rcpp::stop("Nodes and edge weights do not match");
    }
    // do the same for node weights
    std::set<int> node_weight_set;
    for (const auto& entry : nodeWeights) {
        node_weight_set.insert(entry.first);
    }
    if (node_set != node_weight_set) {
        Rcpp::stop("Nodes and node weights do not match");
    }
}
/*
##############################################
####### R INTERFACE FUNCTIONS ################
*/
// convert R list to custom graph object
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
        G.nodeIndexMap[node] = index;
        index++;
    }

    // Add edges to the graph
    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);
        double w = weight[i];

        G.addEdge(fromNode, toNode, w);
    }

    G.updateNodeProperties();
    Rcpp::Rcout << "Removing Single Connections..." << std::endl;
    G.removeSingleConnections();

    return G;
}

// create custom graph object from R list
// [[Rcpp::export]]
Rcpp::List createGraphFromList(const Rcpp::List& graphList) {
    Graph G = listToGraph(graphList); // Convert R list
    return G.graphToRList();
}
