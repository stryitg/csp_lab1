#pragma once

#include <vector>
#include <cstdint>
#include <Graph.hpp>
#include <FordFulkerson.hpp>

class TabledGraph {
public:
    struct Node {
        size_t t;
        size_t k;
    };
    
    struct Edge {
        Node n1;
        Node n2;
    };
    
    TabledGraph(size_t T, size_t K, double vals = 0.0);    
    
    AdjacmentMatrix ToMincut() const;
    std::vector<size_t> ToLabels(const FordFulkerson::STComponents& bg_solution) const;
    
    const double& GetQ(const Node& n) const;
    double& GetQ(const Node& n);
    const double& GetG(const Edge& e) const;
    double& GetG(const Edge& e);
    
    size_t GetT() const { return m_T; }
    size_t GetK() const { return m_K; }
    
private:
    using NodesVals = std::vector<std::vector<double>>;
    using EdgesVals = std::vector<std::vector<std::vector<std::vector<double>>>>;
    
    static EdgesVals CreateEdges(size_t T, size_t K, double vals = 0.0);
    static NodesVals CreateNodes(size_t T, size_t K, double vals = 0.0);
    
    TabledGraph ToBinaryGraph() const;
    void InitNodesForBGraph(TabledGraph& bg) const;
    void InitEdgesForBGraph(TabledGraph& bg) const;
    double ComputeBinaryEdge(const Edge& e) const;
    void MakeBinaryValsPositive(TabledGraph& bg) const;
    void DecoupleLastNodes(TabledGraph& bg) const;
    void MakeUnaryValsPositive(TabledGraph& bg) const;
    AdjacmentMatrix BinaryToAdjacment(const TabledGraph& bg) const;
    void SetSTVals(AdjacmentMatrix& matrix, const TabledGraph& bg) const;
    void SetNonSTVals(AdjacmentMatrix& matrix, const TabledGraph& bg) const;
    
private:
    EdgesVals m_g;
    NodesVals m_q;
    
    const size_t m_T;
    const size_t m_K;
};