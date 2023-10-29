countEdges <- cppFunction(
    file("rcpp/countEdges.cpp", open = "r")
)