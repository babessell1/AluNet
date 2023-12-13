build <- function() {
    Rcpp.package.skeleton(
        "AluNet",
        code_files = c(
            "r/plotLeiden.R",
            "r/graphFromFitHIC.R"
        ),
        cpp_files = c(
            "rcpp/Leiden.cpp",
            "rcpp/Leiden.h",
            "rcpp/Partition.cpp",
            "rcpp/Partition.h",
            "rcpp/Community.cpp",
            "rcpp/Community.h",
            "rcpp/GraphUtils.cpp",
            "rcpp/GraphUtils.h"
        ),
        example_code = FALSE
    )
}