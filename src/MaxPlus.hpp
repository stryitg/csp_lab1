#pragma once
#include <TabledGraph.hpp>
#include <Diffusion.hpp>

class MaxPlus {
public:
    MaxPlus(const TabledGraph& graph, size_t diffusion_iter = 10, double eps_b = 1.0, double eps_e = 0.001);
    
    std::vector<size_t> ComputePath();
    
private:
    bool IsNonEmptyByCrossing(const TabledGraph& graph) const;
    std::vector<size_t> GetPath(const TabledGraph& graph) const;
    TabledGraph ToBinaryProblem(const TabledGraph& g) const;
    void SetQs(TabledGraph& g) const;
    void SetGs(TabledGraph& g) const;
    double GetMaxEdge(const TabledGraph &g, size_t t, size_t t_) const;
    
private:
    Diffusion m_diffusion;
    
    size_t m_iter;
    double m_eps_b;
    double m_eps_e;
};