#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <random>
#include <unordered_map>

class Graph {
public:
    int n;
    std::unordered_map<std::string, int> nodeIndexMap;
    std::unordered_map<int, std::unordered_map<int, double>> edgeWeights;
    std::unordered_map<int, double> nodeWeights;
    std::vector<int> nodes;
    bool isDirected;
    int possibleEdges;
    double totalEdgeWeight;

    Graph(int n);
    void addEdge(const std::string& u, const std::string& v, double w);
    int getNodeIndex(const std::string& node) const;
    std::string getNodeName(int index) const;
    std::vector<int> getNodes() const;
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
