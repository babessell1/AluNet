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
    std::vector<std::vector<double>> weights;

    Graph(int n);
    void addEdge(const std::string& u, const std::string& v, double w);
    int getNodeIndex(const std::string& node) const;
    std::string getNodeName(int index) const;
    std::vector<int> getNodes() const;
    void removeSingleConnections();
    Rcpp::List graphToRList() const;
    inline double Graph::possible_edges(double csize) {// Assuming an undirected graph
        return csize * (csize - 1) / 2; }

    inline bool Graph::is_directed() {// Update this based on your Graph class implementation
        return false; // or true if it's a directed graph
    }
};

// listToGraph
Graph listToGraph(const Rcpp::List& graphList);


#endif
