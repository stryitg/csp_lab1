#pragma once

#include <vector>
#include <unordered_map>
#include <cstdint>
#include <Graph.hpp>
#include <FordFulkerson.hpp>

class TabledGraph {
public:
    TabledGraph() = default;
    struct Node {
        size_t t;
        size_t k;
    };
    
    struct Edge {
        Node n;
        Node n_;
    };
    
    bool operator== (const TabledGraph& rhs) const {
        return rhs.m_g == m_g;
    }
    
    using Neighbors = std::unordered_map<size_t, std::vector<size_t>>;
    
    TabledGraph(size_t T, size_t K, const Neighbors& neighbors = {}, double vals = 0.0);
    
    void PropagateQs();
    
    AdjacmentMatrix ToMincut() const;
    std::vector<size_t> ToLabels(const FordFulkerson::STComponents& bg_solution) const;
    
    const double& GetQ(const Node& n) const;
    double& GetQ(const Node& n);
    const double& GetG(const Edge& e) const;
    double& GetG(const Edge& e);
    
    size_t GetT() const { return m_T; }
    size_t GetK() const { return m_K; }
    
    const std::vector<size_t>& GetNeighbors(size_t t) const { return m_neighbors.at(t); };
    
private:
    using NodesVals = std::vector<std::vector<double>>;
    using EdgesVals = std::vector<std::vector<std::vector<std::vector<double>>>>;
    
    EdgesVals CreateEdges(size_t T, size_t K, double vals = 0.0) const;
    NodesVals CreateNodes(size_t T, size_t K, double vals = 0.0) const;
    
    size_t GetMaxNeighborsSize() const;
    
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
    bool m_sparse;
    Neighbors m_neighbors;
    
    EdgesVals m_g;
    NodesVals m_q;
    
    size_t m_T;
    size_t m_K;
};