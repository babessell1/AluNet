// import Rcpp
#include <Rcpp.h>

// import GraphUtils.h 
#include "GraphUtils.h"

// Export the functions to be used in R
RCPP_MODULE(graph_module) {
    Rcpp::class_<Graph>("Graph")
        .constructor<int>()
        .method("addEdge", &Graph::addEdge);

    Rcpp::function("createGraph", &createGraph);
    Rcpp::function("countEdges", &countEdges);
}