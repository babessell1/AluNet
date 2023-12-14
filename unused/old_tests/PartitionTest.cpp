#include "Partition.h"
#include "gtest/gtest.h"

TEST(PartitionTest, Constructor) {
    Partition p;

    ASSERT_EQ(p.getPartitionSize(), 0);
}

TEST(PartitionTest, AddCommunity) {
    Partition p;

    std::vector<int> nodes1 = {0, 1, 2};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {3, 4, 5};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    ASSERT_EQ(p.getPartitionSize(), 2);
    ASSERT_TRUE(p.hasCommunity(1));
    ASSERT_TRUE(p.hasCommunity(2));
}

TEST(PartitionTest, RemoveCommunity) {
    Partition p;

    std::vector<int> nodes1 = {0, 1, 2};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {3, 4, 5};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    p.removeCommunity(1);

    ASSERT_EQ(p.getPartitionSize(), 1);
    ASSERT_FALSE(p.hasCommunity(1));
    ASSERT_TRUE(p.hasCommunity(2));
}

TEST(PartitionTest, GetCommunity) {
    Partition p;

    std::vector<int> nodes1 = {0, 1, 2};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    Community retrieved = p.getCommunity(1);

    ASSERT_EQ(retrieved.getCommunityIndex(), 1);
    ASSERT_EQ(retrieved.getNodeIndices(), nodes1);
}

TEST(PartitionTest, GetCommunities) {
    Partition p;

    std::vector<int> nodes1 = {0, 1, 2};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {3, 4, 5};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    std::vector<Community> communities = p.getCommunities();

    ASSERT_EQ(communities.size(), 2);
    ASSERT_EQ(communities[0].getCommunityIndex(), 1);
    ASSERT_EQ(communities[1].getCommunityIndex(), 2);
}

TEST(PartitionTest, CalcQuality) {
    Partition p;

    Graph g(4);
    g.addEdge(0, 1, 1.0);
    g.addEdge(1, 2, 1.0);
    g.addEdge(2, 3, 1.0);
    g.addEdge(3, 0, 1.0);

    std::vector<int> nodes1 = {0, 1};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {2, 3};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    double gamma = 1.0;
    double quality = p.calcQuality(gamma, g, false);

    ASSERT_NEAR(quality, 0.0, 0.01); // Assuming the quality is 0.0 for this partition and graph
}

TEST(PartitionTest, PurgeEmptyCommunities) {
    Partition p;

    std::vector<int> nodes1 = {0, 1};
    Community c1(nodes1, 1);
    p.addCommunity(c1);

    std::vector<int> nodes2 = {};
    Community c2(nodes2, 2);
    p.addCommunity(c2);

    p.purgeEmptyCommunities();

    ASSERT_EQ(p.getPartitionSize(), 1);
    ASSERT_TRUE(p.hasCommunity(1));
    ASSERT_FALSE(p.hasCommunity(2));
}