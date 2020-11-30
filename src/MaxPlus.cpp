#include <MaxPlus.hpp>
#include <CrossingOutPath.hpp>
#include <Diffusion.hpp>
#include <iostream>

MaxPlus::MaxPlus(const TabledGraph& graph, size_t diffusion_iter, double eps_b, double eps_e)
: m_diffusion(graph)
, m_iter(diffusion_iter)
, m_eps_b(eps_b)
, m_eps_e(eps_e) {}

std::vector<size_t> MaxPlus::ComputePath() {
    while (m_eps_b > m_eps_e) {
        const auto graph = m_diffusion.Run(m_iter);
        while (!IsNonEmptyByCrossing(graph)) {
            m_eps_b /= 2;
            std::cout << "Current eps: " << m_eps_b << std::endl;
            if (m_eps_b < m_eps_e)
                break;
        }
        m_iter *= 1.5;
    }
    // m_eps_b *= 100;
    std::cout << "Final epsilon: " << m_eps_b << std::endl;    
    
    const auto graph = m_diffusion.Run(m_iter);
    std::vector<size_t> path;
    do {
        path = GetPath(graph);
    } while (path.empty());
    
    return path;
}

bool MaxPlus::IsNonEmptyByCrossing(const TabledGraph& graph) const {
    const auto g = ToBinaryProblem(graph);
    CrossingOut co;
    co.SetGraph(g);
    return !co.Run();
}

std::vector<size_t> MaxPlus::GetPath(const TabledGraph& graph) const {
    const auto g = ToBinaryProblem(graph);
    CrossingOutPath cop(g);
    return cop.ComputePath();
}


TabledGraph MaxPlus::ToBinaryProblem(const TabledGraph& graph) const {
    TabledGraph g = graph;
    SetQs(g);
    SetGs(g);
    return g;
}

void MaxPlus::SetQs(TabledGraph &g) const {
    for (size_t t = 0; t < g.GetT(); ++t) {
        for (size_t k = 0; k < g.GetK(); ++k) {
            g.GetQ({.t = t, .k = k}) = 1;
        }
    }
}

void MaxPlus::SetGs(TabledGraph &g) const {
    for (size_t t = 0; t < g.GetT(); ++t) {
        for (const auto t_: g.GetNeighbors(t)) {
            if (t > t_)
                continue;
            const auto max_edge = GetMaxEdge(g, t, t_);
            for (size_t k = 0; k < g.GetK(); ++k) {
                for (size_t k_ = 0; k_ < g.GetK(); ++k_) {
                    if (g.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) >= max_edge - m_eps_b)
                        g.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) = 1;
                    else
                        g.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) = 0;
                }
            }
        }
    }
}

double MaxPlus::GetMaxEdge(const TabledGraph &g, size_t t, size_t t_) const {
    auto max = -std::numeric_limits<double>::infinity();
    for (size_t k = 0; k < g.GetK(); ++k) {
        for (size_t k_ = 0; k_ < g.GetK(); ++k_) {
            max = std::max(max, g.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}));
        }
    }
    return max;
}



