#pragma once

#include <cstdint>
#include <unordered_map>
#include <TabledGraph.hpp>

class Diffusion {
public:
    Diffusion(const TabledGraph& graph);
    
    TabledGraph Run(uint32_t iterations);
private:
    void UpdateGraph();
    std::unordered_map<size_t, double> ComputeBestEdges(size_t t, size_t k) const;
    double ComputeBestEdge(size_t t, size_t k, size_t t_) const;
    double GetAvg(const std::unordered_map<size_t, double>& edges) const;
        
    TabledGraph m_graph;
};