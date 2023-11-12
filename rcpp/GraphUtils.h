#ifndef GRAPHUTILS_H
#define GRAPHUTILS_H

#include <Rcpp.h>
#include <vector>
#include <set>
#include <unordered_map>

class Graph {
public:
    int n;
    std::map<std::string, int> nodeIndexMap;
    std::vector<std::vector<int>> adj;
    std::vector<std::vector<double>> edge_weights;
    std::vector<int> node_weights;
    std::vector<int> nodes;
    bool isDirected;
    int possibleEdges;

    Graph(int n);
    void addEdge(const std::string& u, const std::string& v, double w);
    int getNodeIndex(const std::string& node) const;
    std::string getNodeName(int index) const;
    std::vector<int> getNodes() const;
    void removeSingleConnections();
    Rcpp::List graphToRList() const;
    std::vector<int> getNeighbors(int nodeIndex) const;
    double getWeight(int u, int v) const;
    void updateNodeProperties();
};

// listToGraph
Graph listToGraph(const Rcpp::List& graphList);

class RandomGenerator {
public:
    static std::vector<int> generateRandomPermutation(int n) {
        std::vector<int> permutation(n);
        for (int i = 0; i < n; ++i) {
            permutation[i] = i;
        }
        std::random_shuffle(permutation.begin(), permutation.end());
        return permutation;
    }
};


#endif
