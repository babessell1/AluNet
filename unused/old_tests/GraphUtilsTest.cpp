#include "GraphUtils.h"
#include "gtest/gtest.h"


TEST(GraphTest, Constructor) {
    Graph g(5);

    ASSERT_EQ(g.getN(), 5);
    ASSERT_EQ(g.getTotalEdgeWeight(), 0.0);
}

TEST(GraphTest, AddEdge) {
    Graph g(5);

    g.addEdge("0", "1", 1.0);
    g.addEdge("1", "2", 2.0);
    g.addEdge("2", "3", 3.0);
    g.addEdge("3", "4", 4.0);

    ASSERT_EQ(g.getTotalEdgeWeight(), 10.0);
    ASSERT_TRUE(g.hasEdge(0, 1));
    ASSERT_TRUE(g.hasEdge(1, 2));
    ASSERT_TRUE(g.hasEdge(2, 3));
    ASSERT_TRUE(g.hasEdge(3, 4));
}

TEST(GraphTest, GetNodeIndex) {
    Graph g(5);

    ASSERT_EQ(g.getNodeIndex("0"), 0);
    ASSERT_EQ(g.getNodeIndex("1"), 1);
    ASSERT_EQ(g.getNodeIndex("2"), 2);
    ASSERT_EQ(g.getNodeIndex("3"), 3);

TEST(GraphUtilsTest, RemoveLowConnections) {
    Graph g(5);

    g.addEdge("0", "1", 1.0);
    g.addEdge("1", "2", 2.0);
    g.addEdge("2", "3", 3.0);
    g.addEdge("3", "4", 4.0);

    g.removeLowConnections(2);

    ASSERT_EQ(g.getN(), 3);
    ASSERT_FALSE(g.hasNode("0"));
    ASSERT_FALSE(g.hasNode("4"));
}

TEST(GraphUtilsTest, GetNeighbors) {
    Graph g(5);

    g.addEdge("0", "1", 1.0);
    g.addEdge("0", "2", 2.0);
    g.addEdge("0", "3", 3.0);
    g.addEdge("0", "4", 4.0);

    auto neighbors = g.getNeighbors("0");

    ASSERT_EQ(neighbors.size(), 4);
    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), "1") != neighbors.end());
    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), "2") != neighbors.end());
    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), "3") != neighbors.end());
    ASSERT_TRUE(std::find(neighbors.begin(), neighbors.end(), "4") != neighbors.end());
}