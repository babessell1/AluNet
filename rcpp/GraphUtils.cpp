#include <Rcpp.h>
using namespace Rcpp;
#include "GraphUtils.h"

Graph::Graph(int n) : n(n), adj(n), weights(n) {};

void Graph::addEdge(const std::string& u, const std::string& v, double w) {
    int uIndex = getNodeIndex(u);
    int vIndex = getNodeIndex(v);
    adj[uIndex].push_back(vIndex);
    weights[uIndex].push_back(w);
}

int Graph::getNodeIndex(const std::string& node) const {
    if (nodeIndexMap.find(node) != nodeIndexMap.end()) {
        return nodeIndexMap.at(node);
    }
    Rcpp::stop("Node not found: " + node);
}

std::string Graph::getNodeName(int index) const {
    for (const auto& entry : nodeIndexMap) {
        if (entry.second == index) {
            return entry.first;
        }
    }
    Rcpp::stop("Index not found: " + std::to_string(index));
}

std::vector<int> Graph::getNodes() const {
    std::vector<int> nodes;
    for (const auto& entry : nodeIndexMap) {
        nodes.push_back(entry.second);
    }
    return nodes;
}

Rcpp::List Graph::graphToRList() const {
    Rcpp::CharacterVector from;
    Rcpp::CharacterVector to;
    Rcpp::NumericVector weight;

    for (int i : getNodes()) {
        std::string fromNode = getNodeName(i);
        for (std::vector<int>::size_type j = 0; j < adj[i].size(); j++) {
            std::string toNode = getNodeName(adj[i][j]);
            from.push_back(fromNode);
            to.push_back(toNode);
            weight.push_back(weights[i][j]);
        }
    }

    Rcpp::List result = Rcpp::List::create(
        Rcpp::Named("from") = from,
        Rcpp::Named("to") = to,
        Rcpp::Named("weight") = weight
    );

    return result;
}

Graph listToGraph(const Rcpp::List& graphList) {
    Rcpp::CharacterVector from = graphList["from"];
    Rcpp::CharacterVector to = graphList["to"];
    Rcpp::NumericVector weight = graphList["weight"];

    std::set<std::string> uniqueNodes;

    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);

        uniqueNodes.insert(fromNode);
        uniqueNodes.insert(toNode);
    }

    Graph G(uniqueNodes.size());

    int index = 0;
    for (const std::string& node : uniqueNodes) {
        G.nodeIndexMap[node] = index;
        index++;
    }

    for (int i = 0; i < from.size(); i++) {
        std::string fromNode = Rcpp::as<std::string>(from[i]);
        std::string toNode = Rcpp::as<std::string>(to[i]);
        double w = weight[i];

        G.addEdge(fromNode, toNode, w);
    }

    return G;
}

void Graph::removeSingleConnections() {
    std::vector<int> nodesToRemove;

    for (int i : getNodes()) {
        if (adj[i].size() == 1) {
            // Node has only one connection, mark for removal
            nodesToRemove.push_back(i);
        }
    }

   for (int i : nodesToRemove) {
       // Remove node from nodeIndexMap
       nodeIndexMap.erase(getNodeName(i));

       // Remove node from adj
       adj.erase(adj.begin() + i);

       // Remove node from weights
       weights.erase(weights.begin() + i);

       // Remove node from adj and weights of other nodes
       for (int j : getNodes()) {
           for (std::vector<int>::size_type k = 0; k < adj[j].size(); k++) {
               if (adj[j][k] == i) {
                   adj[j].erase(adj[j].begin() + k);
                   weights[j].erase(weights[j].begin() + k);
               }
           }
       }
   }
}

// [[Rcpp::export]]
Rcpp::List createGraphFromMatrix(const Rcpp::NumericMatrix& edgeWeights) {
    Graph G(edgeWeights.nrow());

    for (int i = 0; i < edgeWeights.nrow(); i++) {
        for (int j = 0; j < edgeWeights.ncol(); j++) {
            if (edgeWeights(i, j) != 0) {
                G.addEdge(G.getNodeName(i), G.getNodeName(j), edgeWeights(i, j));
            }
        }
    }

    return G.graphToRList();
}

// [[Rcpp::export]]
Rcpp::List createGraphFromList(const Rcpp::List& graphList) {
    Graph G = listToGraph(graphList); // Convert R list
    G.removeSingleConnections();
    return G.graphToRList();
}

