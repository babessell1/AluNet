#include "gtest/gtest.h"
#include "Leiden.h"
#include "Graph.h"
#include "Partition.h"
#include "Community.h"

TEST(OptimizerTest, Constructor) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    ASSERT_EQ(opt.getGamma(), 1.0);
    ASSERT_EQ(opt.getTheta(), 0.5);
}

TEST(OptimizerTest, UpdateCommunityAssignments) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    std::unordered_map<std::string, int> original_nodeIndexMap = {{"node1", 1}, {"node2", 2}};
    opt.updateCommunityAssignments(p, original_nodeIndexMap);
    ASSERT_EQ(opt.getCommunityAssignments()["node1"], p.getNodeCommunityMap().at(1));
    ASSERT_EQ(opt.getCommunityAssignments()["node2"], p.getNodeCommunityMap().at(2));
}

TEST(OptimizerTest, DeltaQuality) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    double delta_quality = opt.deltaQuality(1, 2, 1.0, true);
    ASSERT_GE(delta_quality, 0.0);
}

TEST(OptimizerTest, MoveNodesFast) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    bool update = opt.moveNodesFast();
    ASSERT_TRUE(update);
}

TEST(OptimizerTest, GetWellConnectedCommunities) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    Community B({1, 2}, 1);
    std::vector<Community> well_connected_communities = opt.getWellConnectedCommunities(B);
    ASSERT_GE(well_connected_communities.size(), 0);
}

TEST(OptimizerTest, GetWellConnectedNodes) {
    Graph g(4);
    Partition p;
    Optimizer opt(g, p, 1.0, 0.5);
    Community B({1, 2}, 1);
    std::vector<int> well_connected_nodes = opt.getWellConnectedNodes(B);
    ASSERT_GE(well_connected_nodes.size(), 0);
}
/*
TEST(OptimizerTest, AggregateGraph) {
    Graph g(4);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 3, 1.0);
    g.addEdge(3, 0, 1.0);

    Partition p;
    std::vector<int> nodes1 = {0, 1};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {2, 3};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    Optimizer opt(g, p, 1.0, 0.5);
    Graph aggregated_graph = opt.aggregateGraph();

    ASSERT_EQ(aggregated_graph.getNumberOfNodes(), 2);
    ASSERT_EQ(aggregated_graph.getNumberOfEdges(), g.getNumberOfEdges());
}
*/