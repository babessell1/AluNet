#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <random>
#include <unordered_map>

class Graph {
private:
    // attributes
    int n;
    std::unordered_map<std::string, int> nodeIndexMap;
    std::unordered_map<int, std::unordered_map<int, double>> edgeWeights;
    std::unordered_map<int, double> nodeWeights;
    std::vector<int> nodes;
    bool isDirected;
    int possibleEdges;
    double totalEdgeWeight;

public:
    // getters
    int getN() const{
        return n;
    }

    void addNodeName(std::string node, int index) {
        nodeIndexMap[node] = index;
    }

    std::unordered_map<std::string, int> getNodeIndexMap() const {
        return nodeIndexMap;
    }

    std::unordered_map<int, std::string> getNodeNameMap() const {
        std::unordered_map<int, std::string> node_name_map;
        for (const auto& entry : nodeIndexMap) {
            node_name_map[entry.second] = entry.first;
        }
        return node_name_map;
    }

    std::unordered_map<int, std::unordered_map<int, double>> getEdgeWeights() const {
        return edgeWeights;
    }

    std::unordered_map<int, double> getThisEdgeWeights(int u) const {
        return edgeWeights.at(u);
    }

    std::unordered_map<int, double> getNodeWeights() const {
        return nodeWeights;
    }

    double getNodeWeight(int node) const {
        return nodeWeights.at(node);
    }

    std::unordered_map<int, double> getEdgeWeightsFromStr(std::string u) {
        int u_idx = getNodeIndex(u);
        return edgeWeights.at(u_idx);
    }

    std::vector<int> getNodes() const {
        return nodes;
    }

    bool getIsDirected() const {
        return isDirected;
    }
    
    int getPossibleEdges() const {
        return possibleEdges;
    }

    double getTotalEdgeWeight() const {
        return totalEdgeWeight;
    }

    // setters
    void setN() {
        size_t n_t  = nodeIndexMap.size();
        if (n_t > INT_MAX) {
            Rcpp::stop("Number of nodes exceeds maximum integer value");
        }
        n = static_cast<int>(n_t);
    }

    void setNodeIndexMap(std::unordered_map<std::string, int> nodeIndexMap) {
        nodeIndexMap = nodeIndexMap;
    }
    
    void setNodeIndex(std::string node, int index) {
        nodeIndexMap[node] = index;
    }

    void setEdgeWeights(std::unordered_map<int, std::unordered_map<int, double>> edgeWeights) {
        this->edgeWeights = edgeWeights;
    }

    void initEdgeWeight(int u) {
        this->edgeWeights[u] = std::unordered_map<int, double>();
    }

    void setEdgeWeight(std::string u, std::string v, double w) {
        int u_idx = getNodeIndex(u);
        int v_idx = getNodeIndex(v);
        edgeWeights[u_idx][v_idx] = w;
    }

    void updateAggEdgeWeight(std::string u, std::string v, double w) {
        int u_idx = getNodeIndex(u);
        int v_idx = getNodeIndex(v);
        if (edgeWeights.at(u_idx).find(v_idx) != edgeWeights.at(u_idx).end()) {
            this->edgeWeights.at(u_idx).at(v_idx) = edgeWeights.at(u_idx).at(v_idx) + w;
        } else {
            this->edgeWeights.at(u_idx)[v_idx] = w;
        }
    }

    void setNodeWeights(std::unordered_map<int, double> nodeWeights) {
        nodeWeights = nodeWeights;
    }

    void setNodeWeight(int node, double weight) {
        nodeWeights[node] = weight;
    }

    void incNodeWeight(int node, double weight) {
        if (nodeWeights.find(node) == nodeWeights.end()) {
            this->nodeWeights[node] = 0.0;
        } else {
            this->nodeWeights[node] += weight;
        }
    }

    void setNodes(std::vector<int> nodes) {
        nodes = nodes;
    }

    void setIsDirected(bool is_directed) {
        isDirected = is_directed;
    }

    void addNode(int node) {
        nodes.push_back(node);
    }

    void incTotalEdgeWeight(double weight) {
        this->totalEdgeWeight += weight;
    }

    // methods
    Graph(int n);
    Graph();
    void addEdge(const std::string& u, const std::string& v, double w);
    int getNodeIndex(const std::string& node) const;
    std::string getNodeName(int index) const;
    void resetNodes();
    void removeLowConnections(int min_connections);
    Rcpp::List graphToRList(std::unordered_map<std::string, int>& community_assignments, double quality) const;
    std::vector<int> getNeighbors(int nodeIndex) const;
    double getWeight(int u, int v) const;
    void updateNodeProperties(bool remove_empty_nodes);
    bool hasEdge(int u, int v) const;
};

class RandomGenerator { 
public:
    static std::vector<int> generateRandomPermutation(int n) {
        std::vector<int> permutation(n);
        for (int i = 0; i < n; ++i) {
            permutation[i] = i;
        }
        // Create a random number generator
        std::random_device rd;
        std::mt19937 g(rd());

        // Use std::shuffle instead of std::random_shuffle
        std::shuffle(permutation.begin(), permutation.end(), g);
        return permutation;
    }
};

#endif
