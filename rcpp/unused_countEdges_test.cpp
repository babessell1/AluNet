// import Rcpp
#include <Rcpp.h>
namespace: Rcpp;
#include <vector>
#include <set>
#include <unordered_map>

// define adjacency graph structure with weights 
class Graph {
    // adj is  size n vector of vectors of ints: adj[i] is vector of neighbors of node i 
    // weights is size n vector of vectors of doubles: weights[i] is vector of weights of edges from node i
    // n is number of nodes in graph
public:
    int n; // number of nodes
    std::vector<std::vector<int> > adj; // adjacency list
    std::vector<std::vector<double> > weights; // weights of edges
    
    Graph(int n) : n(n), adj(n), weights(n) {}
    
    void addEdge(int u, int v, double w) {
        adj[u].push_back(v); // add v to u's list of neighbors
        weights[u].push_back(w); // add weight of edge from u to v
    }
};



int countEdges(const Graph& C, const Graph& D) {
    // edge count number;
    int E = 0;
    
    
    std::unordered_map<int, std::set<int>> neighborsC;  // map from node ID to set of neighbors
    std::unordered_map<int, std::vector<double>> weightMapC;  // map from node ID to vector of edge weights
    std::unordered_map<int, std::set<int>> neighborsD;  // map from node ID to set of neighbors

    
    // Populate the maps
    for (int i = 0; i < C.n; i++) { // for each node, i in C
        for (int j = 0; j < C.adj[i].size(); j++) { // for each neighbor, j of node i
            int neighbor = C.adj[i][j];
            double weight = C.weights[i][j];
            neighborsC[i].insert(neighbor);  // add neighbor to set of neighbors
            weightMapC[i].push_back(weight);  // add weight of edge from i to neighbor
        }
    }

    for (int i = 0; i < D.n; i++) {  // for each node, i in D
        for (int j = 0; j < D.adj[i].size(); j++) {  // for each neighbor, j of node i
            int neighbor = D.adj[i][j];  // get neighbor
            neighborsD[i].insert(neighbor);  // add neighbor to set of neighbors
        }
    }

    // Count common neighbors with their weighted edges
    for (int i = 0; i < C.n; i++) {  // for each node, i in C
        for (int neighbor : neighborsC[i]) {  // for each neighbor of node i
            if (neighborsD.find(neighbor) != neighborsD.end()) {  // if neighbor is in D
                for (int k : neighborsD[neighbor]) {  // for each neighbor, k of neighbor
                    if (neighborsC[i].count(k) > 0) {  // if k is also a neighbor of i
                        int index1 = i;
                        int index2 = k;
                        int index3 = neighbor;
                        
                        for (int l = 0; l < C.adj[index1].size(); l++) {  // for each neighbor, l of node i
                            if (C.adj[index1][l] == index3) {  // if l is neighbor of i
                                E += C.weights[index1][l];  
                            }
                        }
                    }
                }
            }
        }
    }

    return E;
}


// Export the functions to be used in R
RCPP_MODULE(graph_module) {
    Rcpp::class_<Graph>("Graph")
        .constructor<int>()
        .method("addEdge", &Graph::addEdge);

    Rcpp::function("createGraph", &createGraph);
    Rcpp::function("countEdges", &countEdges);
}