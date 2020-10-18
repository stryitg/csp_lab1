#include <TabledGraph.hpp>
#include <limits>
#include <iostream>

TabledGraph::TabledGraph(size_t T, size_t K, double vals) 
: m_g(CreateEdges(T, K, vals))
, m_q(CreateNodes(T, K, vals))
, m_T(T)
, m_K(K) {}

TabledGraph::EdgesVals TabledGraph::CreateEdges(size_t T, size_t K, double vals) {
    return std::vector<std::vector<NodesVals>>(T, std::vector<NodesVals>(K, CreateNodes(T, K, vals)));
}

TabledGraph::NodesVals TabledGraph::CreateNodes(size_t T, size_t K, double vals) {
    return std::vector<std::vector<double>>(T, std::vector<double>(K, vals));
}

TabledGraph TabledGraph::ToBinaryGraph() const {
    TabledGraph bg(m_T * m_K, 2);
    InitNodesForBGraph(bg);
    InitEdgesForBGraph(bg);
    return bg;
}

void TabledGraph::InitNodesForBGraph(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K - 1; ++k) {
            const auto tb = t * m_K + k;
            bg.Get_Q({.t = tb, .k = 0}) = 0;
            bg.Get_Q({.t = tb, .k = 1}) = Get_Q({.t = t, .k = k}) - Get_Q({.t = t, .k = k + 1});
            bg.Get_G({{.t = tb, .k = 0}, {.t = tb + 1, .k = 0}}) = 0;
            bg.Get_G({{.t = tb, .k = 0}, {.t = tb + 1, .k = 1}}) = 0;
            bg.Get_G({{.t = tb, .k = 1}, {.t = tb + 1, .k = 0}}) = std::numeric_limits<double>::infinity();
            bg.Get_G({{.t = tb, .k = 1}, {.t = tb + 1, .k = 1}}) = 0;
        }
        const auto tb = t * m_K + m_K - 1;
        bg.Get_Q({.t = tb, .k = 0}) = std::numeric_limits<double>::infinity();
        bg.Get_Q({.t = tb, .k = 1}) = Get_Q({.t = t, .k = m_K - 1});
        bg.Get_G({{.t = tb - 1, .k = 0}, {.t = tb, .k = 0}}) = std::numeric_limits<double>::infinity();
        bg.Get_G({{.t = tb - 1, .k = 0}, {.t = tb, .k = 1}}) = 0;
        bg.Get_G({{.t = tb - 1, .k = 1}, {.t = tb, .k = 0}}) = std::numeric_limits<double>::infinity();
        bg.Get_G({{.t = tb - 1, .k = 1}, {.t = tb, .k = 1}}) = 0;
    }
}

void TabledGraph::InitEdgesForBGraph(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            for (size_t t_ = t + 1; t_ < m_T; ++t_) {
                for (size_t k_ = 0; k_ < m_K; ++k_) {
                    const auto tb = t * m_K + k;
                    const auto tb_ = t_ * m_K + k_;
                    bg.Get_G({{.t = tb, .k = 0}, {.t = tb_, .k = 0}}) = 0;
                    bg.Get_G({{.t = tb, .k = 0}, {.t = tb_, .k = 1}}) = 0;
                    bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}}) = 0;
                    bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 1}}) = ComputeBinaryEdge({{.t = t, .k = k}, {.t = t_, .k = k_}});
                }
            }
        }
    }
}

double TabledGraph::ComputeBinaryEdge(const Edge& e) const {
    if (e.n1.k == m_K - 1 && e.n2.k == m_K - 1) {
        return Get_G(e);
    } else if (e.n1.k == m_K - 1 && e.n2.k < m_K - 1) {
        return Get_G(e) - Get_G({.n1 = e.n1, {.t = e.n2.t, .k = e.n2.k + 1}});
    } else if (e.n1.k < m_K - 1 && e.n2.k == m_K - 1) {
        return Get_G(e) - Get_G({{.t = e.n1.t, .k = e.n1.k + 1}, .n2 = e.n2});
    } else {
        return Get_G(e) - Get_G({{.t = e.n1.t, .k = e.n1.k + 1}, .n2 = e.n2}) -
               Get_G({.n1 = e.n1, {.t = e.n2.t, .k = e.n2.k + 1}}) + Get_G({{.t = e.n1.t, .k = e.n1.k + 1}, 
                                                                            {.t = e.n2.t, .k = e.n2.k + 1}});
    }
}

AdjacmentMatrix TabledGraph::ToMinCut() const {
    auto bg = ToBinaryGraph();
    MakeBinaryValsPositive(bg);
    DecoupleLastNodes(bg);
    MakeUnaryValsPositive(bg);
    return BinaryToAdjacment(bg);
}

void TabledGraph::MakeBinaryValsPositive(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            for (size_t t_ = t + 1; t_ < m_T; ++t_) {
                for (size_t k_ = 0; k_ < m_K - 1; ++k_) {
                    const auto tb = t * m_K + k;
                    const auto tb_ = t_ * m_K + k_;
                    bg.Get_Q({.t = tb, .k = 1}) += bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 1}})
                                                  - bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}});
                    bg.Get_Q({.t = tb_, .k = 1}) += bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}})
                                                 - bg.Get_G({{.t = tb, .k = 0}, {.t = tb_, .k = 0}});
                    bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}}) = bg.Get_G({{.t = tb, .k = 0}, {.t = tb_, .k = 1}})
                                                                      + bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}})
                                                                      - bg.Get_G({{.t = tb, .k = 0}, {.t = tb_, .k = 0}})
                                                                      - bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 1}});
                }
            }
        }
    }
}

void TabledGraph::DecoupleLastNodes(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t t_ = t + 1; t_ < m_T; ++t_) {
            for (size_t k_ = 0; k_ < m_K - 1; ++k_) {
                const auto tb1 = t * m_K + m_K - 1;
                const auto tb1_ = t_ * m_K + k_;
                const auto tb2 = t * m_K + k_;
                const auto tb2_ = t_ * m_K + m_K - 1;
                auto& v1 = bg.Get_G({{.t = tb1, .k = 1}, {.t = tb1_, .k = 1}});
                bg.Get_Q({.t = tb1_, .k = 1}) += v1;
                v1 = 0;
                auto& v2 = bg.Get_G({{.t = tb2, .k = 1}, {.t = tb2_, .k = 1}});
                bg.Get_Q({.t = tb2, .k = 1}) += v2;
                v2 = 0;
            }
            const auto tb = t * m_K + m_K - 1;
            const auto tb_ = t_ * m_K + m_K - 1;
            bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 1}}) = 0;
        }
    }
}

void TabledGraph::MakeUnaryValsPositive(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            const auto tb = t * m_K + k;
            if (auto& v = bg.Get_Q({.t = tb, .k = 0}); v < 0) {
                bg.Get_Q({.t = tb, .k = 1}) -= v;
                v = 0;
            }
            if (auto& v = bg.Get_Q({.t = tb, .k = 1}); v < 0) {
                bg.Get_Q({.t = tb, .k = 0}) -= v;
                v = 0;
            }
        }
        const auto tb = t * m_K + m_K - 1;
        bg.Get_G({{.t = tb - 1, .k = 1}, {.t = tb, .k = 0}}) = 0;
    }
}

AdjacmentMatrix TabledGraph::BinaryToAdjacment(const TabledGraph& bg) const {
    const size_t size = m_T * (m_K - 1) + 2;
    AdjacmentMatrix matrix(size, std::vector<double>(size, 0));
    SetSTVals(matrix, bg);
    SetNonSTVals(matrix, bg);
    return matrix;
}

void TabledGraph::SetSTVals(AdjacmentMatrix& matrix, const TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < (m_K - 1); ++k) {
            const auto tb = t * m_K + k; 
            matrix[0][tb - t + 1] = bg.Get_Q({.t = tb, .k = 0});
            matrix[tb - t + 1][m_T * (m_K - 1) + 1] = bg.Get_Q({.t = tb, .k = 1});
        }
    }
}

void TabledGraph::SetNonSTVals(AdjacmentMatrix& matrix, const TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K - 1; ++k) {
            for (size_t t_ = t; t_ < m_T; ++t_) {
                for (size_t k_ = 0; k_ < m_K - 1; ++k_) {
                    const auto tb = t * m_K + k;
                    const auto tb_ = t_ * m_K + k_;
                    matrix[tb - t + 1][tb_ - t_ + 1] = bg.Get_G({{.t = tb, .k = 1}, {.t = tb_, .k = 0}});
                }
            }
        }
    }
}

std::vector<size_t> TabledGraph::ToPath(const FordFulkerson::STComponents& bg_solution) const {
    size_t source_c = 1;
    size_t sink_c = 0;
    const auto& source = bg_solution.s;
    const auto& sink = bg_solution.t;
    std::vector<size_t> path(m_T);
    for (size_t t = 0; t < m_T; ++t) {
        size_t count = 0;
        for (size_t k = 0; k < m_K - 1; ++k) {
            if (source_c == source.size())
                break;
            if (sink_c == sink.size() || source[source_c] < sink[sink_c]) {
                ++count;
                ++source_c;
            } else {
                ++sink_c;
            }
        }
        path[t] = m_K - 1 - count;
    }
    return path;
}

const double& TabledGraph::Get_Q(const Node& n) const {
    return m_q[n.t][n.k];
}

double& TabledGraph::Get_Q(const Node& n) {
    return const_cast<double&>(static_cast<const TabledGraph&>(*this).Get_Q(n));
}

const double& TabledGraph::Get_G(const Edge& e) const {
    return m_g[e.n1.t][e.n1.k][e.n2.t][e.n2.k];
}

double& TabledGraph::Get_G(const Edge& n) {
    return const_cast<double&>(static_cast<const TabledGraph&>(*this).Get_G(n));
}