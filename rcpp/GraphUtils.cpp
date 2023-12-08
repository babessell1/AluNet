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
    // update total edge weight
    totalEdgeWeight += w;
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
void Graph::resetNodes() {
    Rcpp::Rcout << "Resetting nodes..." << std::endl;
    // print node index map
    std::vector<int> nodes_;
    for (const auto& entry : nodeIndexMap) {
        nodes_.push_back(entry.second);
    }
    setNodes(nodes_);
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

    // calculate total edge weight
    double total_edge_weight = 0.0;
    //for (const auto& entry : new_edge_weights) {
    for (const auto& entry : edgeWeights) {
        for (const auto& edge : entry.second) {
            total_edge_weight += edge.second;
        }
    }
    if (!isDirected) {
        total_edge_weight /= 2; // In an undirected graph, each edge is counted twice 
    }

    // update Graph properties
    //edgeWeights = new_edge_weights;
    //nodeWeights = new_node_weights;
    //nodeIndexMap = new_node_index_map;
    totalEdgeWeight = total_edge_weight;
    edgeWeights = new_edge_weights;
    nodeWeights = new_node_weights;
    nodeIndexMap = new_node_index_map;
    resetNodes();
    setN();
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

// convert custom graph object to R list
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
    updateNodeProperties(false);

}

bool Graph::hasEdge(int u, int v) const {
    auto it = edgeWeights.find(u);
    if (it != edgeWeights.end()) {
        return it->second.find(v) != it->second.end();
    }
    return false;
    
}

/*
###################################################################################
####################### R INTERFACE FUNCTIONS #####################################
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

// create custom graph object from R list
// [[Rcpp::export]]
Rcpp::List createGraphFromList(const Rcpp::List& graphList) {
    Graph G = listToGraph(graphList); // Convert R list
    double quality = 0.0;
    std::unordered_map<std::string, int> community_assignments;
    return G.graphToRList(community_assignments, quality);
}