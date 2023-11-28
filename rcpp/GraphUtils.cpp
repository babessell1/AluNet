#include <Rcpp.h>
#include "GraphUtils.h"

/*
##############################################
####### GRAPH CLASS FUNCTIONS ################
*/
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
    std::vector<int> nodes;
    for (int i = 0; i < n; ++i) {
      nodes.push_back(i);
    }
    this->nodes = nodes;
    // initialize node index map
    std::unordered_map<std::string, int> nodeIndexMap;
    this->nodeIndexMap = nodeIndexMap;

    this->totalEdgeWeight = 0.0;
};

// add new edge to the graph, give it the two node names to connect and the weight
void Graph::addEdge(const std::string& u, const std::string& v, double w) {
    int u_idx = getNodeIndex(u);
    int v_idx = getNodeIndex(v);
    // add the edge to the edge weights map
    edgeWeights[u_idx][v_idx] = w;
    
    // add other way too if undirected
    if (!isDirected) {
        edgeWeights[v_idx][u_idx] = w;
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
    std::unordered_map<int, int> old_node_to_new_node; // map old node index to new node index
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

    
    // update Graph properties
    edgeWeights = new_edge_weights;
    nodeWeights = new_node_weights;
    nodeIndexMap = new_node_index_map;
    nodes = getNodes();
    n = nodes.size();
}


// get neighbors of a node based on keys of the edgeWeights map
std::vector<int> Graph::getNeighbors(int n_idx) const {
    std::vector<int> neighbors;

    // When the graph is directed or not, your getNeighbors can be the same.
    // You only need to check the edges for the current node
    for (const auto& to : edgeWeights.at(n_idx)) {
        neighbors.push_back(to.first);
    }
         
    return neighbors;
}


// get weight of an edge based on node indices
double Graph::getWeight(int u, int v) const {  // find is slow, consider not checking if the edge exists
    if (edgeWeights.find(u) != edgeWeights.end()) {
        if (edgeWeights.at(u).find(v) != edgeWeights.at(u).end()) {
            return edgeWeights.at(u).at(v);
        }
    }
    // otherwise, throw an error
    Rcpp::stop("Edge not found: (" + std::to_string(u) + ", " + std::to_string(v) + ")");
}

// convert custom graph object to R list
Rcpp::List Graph::graphToRList() const {
    Rcpp::CharacterVector from;
    Rcpp::CharacterVector to;
    Rcpp::NumericVector weight;

    Rcpp::Rcout << "Writing to R List" << std::endl;

    for (int i : nodes) {
        std::string fromNode = getNodeName(i);
        for (const auto& entry : edgeWeights.at(i)) {
            std::string toNode = getNodeName(entry.first);
            from.push_back(fromNode);
            to.push_back(toNode);
            weight.push_back(entry.second);
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
    updateNodeProperties(true);
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

    // Add edges to the graph & populate the edgeWeights
    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);
        double w = weight[i];

        G.addEdge(fromNode, toNode, w);
    }

    G.updateNodeProperties(true); // Updating nodes after adding all edges 
    Rcpp::Rcout << "Removing Low Connections..." << std::endl;
    G.removeLowConnections(2); // Remove single connections

    return G;
}

// create custom graph object from R list
// [[Rcpp::export]]
Rcpp::List createGraphFromList(const Rcpp::List& graphList) {
    Graph G = listToGraph(graphList); // Convert R list
    return G.graphToRList();
}
