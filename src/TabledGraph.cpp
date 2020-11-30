#include <TabledGraph.hpp>
#include <limits>
#include <iostream>
#include <cassert>
#include <algorithm>

TabledGraph::TabledGraph(size_t T, size_t K, const Neighbors& neighbors, double vals) 
: m_sparse(!neighbors.empty())
, m_neighbors(neighbors)
, m_g(CreateEdges(T, K, vals))
, m_q(CreateNodes(T, K, vals))
, m_T(T)
, m_K(K) {
}

TabledGraph::EdgesVals TabledGraph::CreateEdges(size_t T, size_t K, double vals) const {
    if (!m_sparse) {
        return std::vector<std::vector<NodesVals>>(T, std::vector<NodesVals>(K, CreateNodes(T, K, vals)));
    }
    else {
        const auto max_size = GetMaxNeighborsSize();
        return std::vector<std::vector<std::vector<std::vector<double>>>>(T, 
                std::vector<std::vector<std::vector<double>>>(max_size, 
                std::vector<std::vector<double>>(K, std::vector<double>(K, vals))));
    }
}

size_t TabledGraph::GetMaxNeighborsSize() const {
    size_t max = 0;
    for (const auto [key, val] : m_neighbors)
        max = std::max(max, val.size());
    return max;
}


TabledGraph::NodesVals TabledGraph::CreateNodes(size_t T, size_t K, double vals) const {
    return std::vector<std::vector<double>>(T, std::vector<double>(K, vals));
}

void TabledGraph::PropagateQs() {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            for (const auto t_ : m_neighbors[t]) {
                for (size_t k_ = 0; k_ < m_K; ++k_) {
                    if (t > t_) {
                        GetG({{.t = t_, .k = k_}, {.t = t, .k = k}}) += GetQ({.t = t, .k = k}) / m_neighbors[t].size();
                    } else {
                        GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) += GetQ({.t = t, .k = k}) / m_neighbors[t].size();
                    }
                }
            }
        }
    }
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
            bg.GetQ({.t = tb, .k = 0}) = 0;
            bg.GetQ({.t = tb, .k = 1}) = GetQ({.t = t, .k = k}) - GetQ({.t = t, .k = k + 1});
            bg.GetG({{.t = tb, .k = 0}, {.t = tb + 1, .k = 0}}) = 0;
            bg.GetG({{.t = tb, .k = 0}, {.t = tb + 1, .k = 1}}) = 0;
            bg.GetG({{.t = tb, .k = 1}, {.t = tb + 1, .k = 0}}) = std::numeric_limits<double>::infinity();
            bg.GetG({{.t = tb, .k = 1}, {.t = tb + 1, .k = 1}}) = 0;
        }
        const auto tb = t * m_K + m_K - 1;
        bg.GetQ({.t = tb, .k = 0}) = std::numeric_limits<double>::infinity();
        bg.GetQ({.t = tb, .k = 1}) = GetQ({.t = t, .k = m_K - 1});
        bg.GetG({{.t = tb - 1, .k = 0}, {.t = tb, .k = 0}}) = std::numeric_limits<double>::infinity();
        bg.GetG({{.t = tb - 1, .k = 0}, {.t = tb, .k = 1}}) = 0;
        bg.GetG({{.t = tb - 1, .k = 1}, {.t = tb, .k = 0}}) = std::numeric_limits<double>::infinity();
        bg.GetG({{.t = tb - 1, .k = 1}, {.t = tb, .k = 1}}) = 0;
    }
}

void TabledGraph::InitEdgesForBGraph(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            for (size_t t_ = t + 1; t_ < m_T; ++t_) {
                for (size_t k_ = 0; k_ < m_K; ++k_) {
                    const auto tb = t * m_K + k;
                    const auto tb_ = t_ * m_K + k_;
                    bg.GetG({{.t = tb, .k = 0}, {.t = tb_, .k = 0}}) = 0;
                    bg.GetG({{.t = tb, .k = 0}, {.t = tb_, .k = 1}}) = 0;
                    bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}}) = 0;
                    bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 1}}) = ComputeBinaryEdge({{.t = t, .k = k}, {.t = t_, .k = k_}});
                }
            }
        }
    }
}

double TabledGraph::ComputeBinaryEdge(const Edge& e) const {
    if (e.n.k == m_K - 1 && e.n_.k == m_K - 1) {
        return GetG(e);
    } else if (e.n.k == m_K - 1 && e.n_.k < m_K - 1) {
        return GetG(e) - GetG({.n = e.n, {.t = e.n_.t, .k = e.n_.k + 1}});
    } else if (e.n.k < m_K - 1 && e.n_.k == m_K - 1) {
        return GetG(e) - GetG({{.t = e.n.t, .k = e.n.k + 1}, .n_ = e.n_});
    } else {
        return GetG(e) - GetG({{.t = e.n.t, .k = e.n.k + 1}, .n_ = e.n_}) -
               GetG({.n = e.n, {.t = e.n_.t, .k = e.n_.k + 1}}) + GetG({{.t = e.n.t, .k = e.n.k + 1}, 
                                                                            {.t = e.n_.t, .k = e.n_.k + 1}});
    }
}

AdjacmentMatrix TabledGraph::ToMincut() const {
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
                    bg.GetQ({.t = tb, .k = 1}) += bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 1}})
                                                  - bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}});
                    bg.GetQ({.t = tb_, .k = 1}) += bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}})
                                                 - bg.GetG({{.t = tb, .k = 0}, {.t = tb_, .k = 0}});
                    bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}}) = bg.GetG({{.t = tb, .k = 0}, {.t = tb_, .k = 1}})
                                                                      + bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}})
                                                                      - bg.GetG({{.t = tb, .k = 0}, {.t = tb_, .k = 0}})
                                                                      - bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 1}});
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
                auto& v1 = bg.GetG({{.t = tb1, .k = 1}, {.t = tb1_, .k = 1}});
                bg.GetQ({.t = tb1_, .k = 1}) += v1;
                v1 = 0;
                auto& v2 = bg.GetG({{.t = tb2, .k = 1}, {.t = tb2_, .k = 1}});
                bg.GetQ({.t = tb2, .k = 1}) += v2;
                v2 = 0;
            }
            const auto tb = t * m_K + m_K - 1;
            const auto tb_ = t_ * m_K + m_K - 1;
            bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 1}}) = 0;
        }
    }
}

void TabledGraph::MakeUnaryValsPositive(TabledGraph& bg) const {
    for (size_t t = 0; t < m_T; ++t) {
        for (size_t k = 0; k < m_K; ++k) {
            const auto tb = t * m_K + k;
            if (auto& v = bg.GetQ({.t = tb, .k = 0}); v < 0) {
                bg.GetQ({.t = tb, .k = 1}) -= v;
                v = 0;
            }
            if (auto& v = bg.GetQ({.t = tb, .k = 1}); v < 0) {
                bg.GetQ({.t = tb, .k = 0}) -= v;
                v = 0;
            }
        }
        const auto tb = t * m_K + m_K - 1;
        bg.GetG({{.t = tb - 1, .k = 1}, {.t = tb, .k = 0}}) = 0;
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
            matrix[0][tb - t + 1] = bg.GetQ({.t = tb, .k = 0});
            matrix[tb - t + 1][m_T * (m_K - 1) + 1] = bg.GetQ({.t = tb, .k = 1});
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
                    matrix[tb - t + 1][tb_ - t_ + 1] = bg.GetG({{.t = tb, .k = 1}, {.t = tb_, .k = 0}});
                }
            }
        }
    }
}

std::vector<size_t> TabledGraph::ToLabels(const FordFulkerson::STComponents& bg_solution) const {
    size_t source_c = 1;
    size_t sink_c = 0;
    const auto& source = bg_solution.s;
    const auto& sink = bg_solution.t;
    std::vector<size_t> labels(m_T);
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
        labels[t] = m_K - 1 - count;
    }
    return labels;
}

const double& TabledGraph::GetQ(const Node& n) const {
    return m_q[n.t][n.k];
}

double& TabledGraph::GetQ(const Node& n) {
    return const_cast<double&>(static_cast<const TabledGraph&>(*this).GetQ(n));
}

const double& TabledGraph::GetG(const Edge& e) const {
    if (!m_sparse) {
        return m_g[e.n.t][e.n.k][e.n_.t][e.n_.k];
    } else {
        const auto iter = m_neighbors.find({e.n.t});
        assert(iter != m_neighbors.end());
        auto pos_iter = std::find(iter->second.begin(), iter->second.end(), e.n_.t);
        assert(pos_iter != iter->second.end());
        return m_g[e.n.t][pos_iter - iter->second.begin()][e.n.k][e.n_.k];
    }
    
}

double& TabledGraph::GetG(const Edge& n) {
    return const_cast<double&>(static_cast<const TabledGraph&>(*this).GetG(n));
}