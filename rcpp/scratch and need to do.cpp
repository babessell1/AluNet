// I will look at the original python file and build our own rcpp functions

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::depends(BH)]]
#include <vector>
#include <queue>

class Graph
{
private:
    int n; // number of nodes
    std::vector<std::vector<int> > adj; // adjacency list
    std::vector<std::vector<double> > weights; // weights of edges

public:
    Graph(int n) : n(n), adj(n, std::vector<int>()), weights(n, std::vector<double>()) {}

    void addEdge(int u, int v, double w)
    {
        adj[u].push_back(v);
        weights[u].push_back(w);
    }

};

/* the main leiden algorithm
while (not done){
    MoveNodesfast(G, P);
    if not done then{
        P_refined <- Refinepartition(G, P);
        G <- AggregateGraph(G, P_refined);
        P <- {{}}
    }
}
*/
// now we should define the 

