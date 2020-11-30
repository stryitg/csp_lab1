#include <Diffusion.hpp>
#include <algorithm>
#include <numeric>
#include <iostream>

Diffusion::Diffusion(const TabledGraph& graph) 
: m_graph(graph) {}

TabledGraph Diffusion::Run(uint32_t iterations) {
    for (uint32_t i = 0; i < iterations; ++i) {
        std::cout << "Diffustion iteration: " << i + 1 << std::endl;
        UpdateGraph();
    }
    return m_graph;
}

void Diffusion::UpdateGraph() {
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        for (size_t k = 0; k < m_graph.GetK(); ++k) {
            auto best_edges = ComputeBestEdges(t, k);
            const auto avg = GetAvg(best_edges);
            for (size_t k_ = 0; k_ < m_graph.GetK(); ++k_) {
                for (const auto t_ : m_graph.GetNeighbors(t)) {
                    if (t < t_) 
                        m_graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) += avg - best_edges[t_];
                    else
                        m_graph.GetG({{.t = t_, .k = k_}, {.t = t, .k = k}}) += avg - best_edges[t_];
                }
            }
        }
    }
}

std::unordered_map<size_t, double> Diffusion::ComputeBestEdges(size_t t, size_t k) const {
    const auto n = m_graph.GetNeighbors(t);
    std::unordered_map<size_t, double> edges(n.size());
    for (const auto t_ : n) {
        edges[t_] = ComputeBestEdge(t, k, t_); 
    }
    return edges;
}

double Diffusion::ComputeBestEdge(size_t t, size_t k, size_t t_) const {
    double max = -std::numeric_limits<double>::infinity();
    for (size_t k_ = 0; k_ < m_graph.GetK(); ++k_) {
        if (t < t_)
            max = std::max(max, m_graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}));
        else
            max = std::max(max, m_graph.GetG({{.t = t_, .k = k_}, {.t = t, .k = k}}));
    }
    return max;
}

double Diffusion::GetAvg(const std::unordered_map<size_t, double>& edges) const {
    double avg = 0;
    for (const auto [key, val] : edges) {
        avg += val;
    }
    return avg / edges.size();
}








