#pragma once

#include <random>

#include <TabledGraph.hpp>
#include <CrossingOut.hpp>

class CrossingOutPath {
public:
    CrossingOutPath(const TabledGraph& graph);
    
    std::vector<size_t> ComputePath();

private:
    bool SelectNode(size_t t);
    std::vector<size_t> GetShaffledNodes(size_t t) const;
    bool TryNodes(size_t t, const std::vector<size_t>& ks);
    bool TryNode(size_t t, size_t k);    

private:
    TabledGraph m_graph;
    CrossingOut m_crossing;
    std::vector<size_t> m_path;
    
    mutable std::random_device m_rd;
    mutable std::mt19937 m_g;
};