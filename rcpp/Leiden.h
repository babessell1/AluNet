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
    std::unordered_map<std::string, std::string> communityAssignments;
    std::unordered_map<std::string, std::string> og_communityAssignments;
    double iteration;
    int seed;

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
    std::unordered_map<std::string, std::string> getCommunityAssignments() const {
        return communityAssignments;
    }

    std::unordered_map<std::string, std::string> getOGCommunityAssignments() const {
        return og_communityAssignments;
    }

    std::unordered_map<std::string, int> getFinalCommunityAssignmentsInt() const {
        std::unordered_map<std::string, int> community_assignments_int;
        for (const auto& entry : og_communityAssignments) {
            community_assignments_int[entry.first] = std::stoi(entry.second);
        }
        return community_assignments_int;
    }

    double getIteration() const {
        return iteration;
    }

    void resetCommunityAssignments() {
        communityAssignments.clear();
    }

    void updateCommunityAssignment(std::string node_name, std::string community_name) {
        communityAssignments[node_name] = community_name;
    }

    // setters
    void setG(Graph G) {
        this->G = G;
    }
    void setP(Partition P) {
        this->P = P;
    }
    void setGamma(double gamma) {
        this->gamma = gamma;
    }
    void setTheta(double theta) {
        theta = theta;
    }
    void setIteration(double iteration) {
        this->iteration = iteration;
    }
    void setCommunityAssignments(const std::unordered_map<std::string, std::string>& new_assignments) {
        this->communityAssignments = new_assignments;
    }
    void setOGCommunityAssignments(const std::unordered_map<std::string, std::string>& new_assignments) {
        this->og_communityAssignments = new_assignments;
    }


    // methods
    Optimizer(Graph& G, Partition& P, double gamma, double theta, int seed);
    void optimize(int iterations);
    bool moveNodesFast();
    double calcQuality(double gamma);
    double deltaQuality(int n_idx, int new_c_idx, double gamma, bool recalculate) const;
    void refinePartition(const Partition& P_original);
    std::vector<int> getWellConnectedNodes(const Community& B) const;
    std::vector<Community> getWellConnectedCommunities(const Community& B) const;
    void mergeNodesSubset(Community& S);
    void aggregateGraph(Partition& P_comm);
    void updateCommunityAssignments(const Partition& P, bool last_call);
    Rcpp::List graphToRList(std::unordered_map<std::string, int>& community_assignments, double quality) const {
        return G.graphToRList(community_assignments, quality);
    }
};

// listToGraph
Graph listToGraph(const Rcpp::List& graphList);

Partition initializePartition(Graph& G);

#endif
