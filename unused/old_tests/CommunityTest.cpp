#include "Community.h"
#include "GraphUtils.h"
#include "gtest/gtest.h"


TEST(CommunityTest, Constructor) {
    Community c;

    ASSERT_EQ(c.getCommunitySize(), 0);
}

TEST(CommunityTest, AddNode) {
    Community c;

    c.addNode("0");
    c.addNode("1");
    c.addNode("2");

    ASSERT_EQ(c.getCommunitySize(), 3);
    ASSERT_TRUE(c.hasNode("0"));
    ASSERT_TRUE(c.hasNode("1"));
    ASSERT_TRUE(c.hasNode("2"));
}

TEST(CommunityTest, RemoveNode) {
    Community c;

    c.addNode("0");
    c.addNode("1");
    c.addNode("2");

    c.removeNode("1");

    ASSERT_EQ(c.getCommunitySize(), 2);
    ASSERT_TRUE(c.hasNode("0"));
    ASSERT_FALSE(c.hasNode("1"));
    ASSERT_TRUE(c.hasNode("2"));
}

TEST(CommunityTest, AggregateWeights) {
    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 2.0);
    g.addEdge(0, 2, 3.0);

    std::vector<int> nodes = {0, 1, 2};
    int index = 1;

    Community c(nodes, index);

    ASSERT_EQ(c.aggregateWeights(g), 6.0);
}

TEST(CommunityTest, GetClusterWeight) {
    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 2.0);
    g.addEdge(0, 2, 3.0);

    std::vector<int> nodes = {0, 1, 2};
    int index = 1;

    Community c(nodes, index);

    ASSERT_EQ(c.getClusterWeight(g), 6.0);
}

TEST(CommunityTest, CountPossibleEdges) {
    Graph g(3);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 2.0);
    g.addEdge(0, 2, 3.0);

    std::vector<int> nodes = {0, 1, 2};
    int index = 1;

    Community c(nodes, index);

    ASSERT_EQ(c.countPossibleEdges(g), 3);
}

