#include <CrossingOutPath.hpp>
#include <CrossingOut.hpp>
#include <iostream>
#include <algorithm>

CrossingOutPath::CrossingOutPath(const TabledGraph& graph)
: m_graph(graph)
, m_path(m_graph.GetT())
, m_rd()
, m_g(m_rd()) {}

std::vector<size_t> CrossingOutPath::ComputePath() {
    std::cout << "Started path selection" << std::endl;
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        if ((t % 100) == 0)
            std::cout << "Selected for: " << t << " nodes" << std::endl;
        if (!SelectNode(t))
            return {};
    }
    std::cout << "Finished" << std::endl;
    return m_path;
}

bool CrossingOutPath::SelectNode(size_t t) {
    const auto nodes = GetShaffledNodes(t);
    return TryNodes(t, nodes);
}

std::vector<size_t> CrossingOutPath::GetShaffledNodes(size_t t) const {
    std::vector<size_t> ks;
    for (size_t k = 0; k < m_graph.GetK(); ++k) {
        if (m_graph.GetQ({.t = t, .k = k}) != 0) {
            ks.push_back(k);
        }
    }
    std::shuffle(ks.begin(), ks.end(), m_g);
    return ks;
}

bool CrossingOutPath::TryNodes(size_t t, const std::vector<size_t>& ks) {
    if (ks.size() == 1) {
        m_path[t] = ks.back();
        return true;
    }
    
    for (const auto k : ks) {
        if (TryNode(t, k)) {
            m_path[t] = k;
            return true;
        }
    }
    return false;
}

bool CrossingOutPath::TryNode(size_t t, size_t k) {
    m_graph.GetQ({.t = t, .k = k}) = 1;
    m_crossing.SetGraph(m_graph);
    if (m_crossing.Run()) {
        m_graph = m_crossing.GetGraph();
        return true;
    } else {
        m_graph.GetQ({.t = t, .k = k}) = 0;
        return false;
    }
}