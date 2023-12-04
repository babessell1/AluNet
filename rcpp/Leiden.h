#ifndef LEIDEN_H
#define LEIDEN_H

#include "GraphUtils.h"

class Optimizer {
private:
    // properties
    Graph& G;
    Partition& P;
    double gamma;
    double theta;
    std::unordered_map<std::string, int> communityAssignments;

public:
    // getters
    Graph getG() const {
        return G;
    }
    Partition getP() const {
        return P;
    }
    double getGamma() const {
        return gamma;
    }
    double getTheta() const {
        return theta;
    }
    std::unordered_map<std::string, int> getCommunityAssignments() const {
        return communityAssignments;
    }

    // setters
    void setG(Graph G) {
        G = G;
    }
    void setP(Partition P) {
        P = P;
    }
    void setGamma(double gamma) {
        gamma = gamma;
    }

    // methods
    Optimizer(Graph& G, Partition& P, double gamma, double theta);
    void optimize(int iterations);
    bool moveNodesFast();
    double calcQuality(double gamma);
    double deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const;
    void refinePartition(const Partition& P_original);
    std::vector<int> getWellConnectedNodes(const Community& B) const;
    std::vector<Community> getWellConnectedCommunities(const Community& B) const;
    void mergeNodesSubset(Community& S);
    Graph aggregateGraph() const;
    void updateCommunityAssignments(const Partition& P, const std::unordered_map<std::string, int>& original_nodeIndexMap);
    Rcpp::List graphToRList(std::unordered_map<std::string, int>& community_assignments, double quality) const {
        return G.graphToRList(community_assignments, quality);
    }
};

// listToGraph
Graph listToGraph(const Rcpp::List& graphList);

Partition initializePartition(Graph& G);

#endif