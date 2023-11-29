#ifndef LEIDEN_H
#define LEIDEN_H

#include "GraphUtils.h"

class Optimizer {
public:
    // properties
    Graph& G;
    Partition& P;
    double gamma;
    double theta;

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
    Graph aggregateGraph();
    
};

Partition initializePartition(Graph& G);

#endif