# AluNet

Weighted Leiden algorithm for modeling Alu-mediated enhancer-promoter interaction networks

## NOW WORKING ON

As for the implementation of the algorithm, we should work on these functions in leiden.cpp Partititon class.

```{C++}
size_t get_community_of_vertex(size_t vertex); // writen without test
double diff_move(size_t vertex, size_t new_community); // writen without test
//void removeCommunity(int communityIndex); // writen without test
//std::vector<int> getCommunityIndices(); // should define this unction
//int getCommunityIndex(int node_index); // should define this unction
//double quality(double resolution_parameter); // should define this unction
//void renumber_communities(); // writen without test
//void renumber_communities(std::vector<size_t> fixed_nodes, std::vector<size_t> fixed_membership); // writen without test
// size_t membership(size_t vertex); function that returns vertex is in which community // not defined

```

Additionally, there is more work to do with the partition class that require other classes:

```{C++}
// Graph* get_graph(); function // return the graphs of the specific partition i;
// which should define the relationship between the Partitions and Graphs function
```
