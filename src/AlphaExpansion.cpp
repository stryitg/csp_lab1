#include <AlphaExpansion.hpp>
#include <algorithm>
#include <numeric>
#include <iostream>

AlphaExpansion::AlphaExpansion(const TabledGraph& graph, double eps)
: m_graph(graph)
, m_labels(m_graph.GetT())
, m_eps(eps)
, m_mincut(m_graph.GetT() + 2, std::vector<double>(m_graph.GetT() + 2))
, m_alpha_ordering(CreateRandomAlphaOrdering())
, m_alpha(m_alpha_ordering.back()) {
    m_alpha_ordering.pop_back();
    std::fill(m_labels.begin(), m_labels.end(), m_alpha);
    m_energy = ComputeEnergy();
}

std::vector<size_t> AlphaExpansion::CreateRandomAlphaOrdering() const {
    std::random_device m_rd;
    std::mt19937 m_gen(m_rd());
    std::vector<size_t> ordering(m_graph.GetK());
    std::iota(ordering.begin(), ordering.end(), 0);
    std::shuffle(ordering.begin(), ordering.end(), m_gen);
    return ordering;
}


std::vector<size_t> AlphaExpansion::ComputeLabels() {
    for (const auto alpha : m_alpha_ordering) {
        m_alpha = alpha;
        MakeAlphaExpansion();
        // if (!IsEnergyDecreased()) {
        //     break;
        // }
    }
    return m_labels;
}

void AlphaExpansion::MakeAlphaExpansion() {
    UpdateMincutProblem();
    FordFulkerson ff(m_mincut);
    const auto sink_source = ff.ComputeMincut();
    UpdateLabels(sink_source.s);
}

void AlphaExpansion::UpdateMincutProblem() {
    ResetMincut();
    SetSmoothingCosts();
    SetDataCosts();
    MakeSinkSourcePositive();
}

void AlphaExpansion::ResetMincut() {
    for (auto& row : m_mincut) {
        std::fill(row.begin(), row.end(), 0);
    }
}

void AlphaExpansion::SetSmoothingCosts() {
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        for (size_t t_ = t + 1; t_ < m_graph.GetT(); ++t_) {
            const auto source_c = ComputeSourceSmoothingCost(t, t_);
            const auto sink_c = ComputeSinkSmoothingCost(t, t_);
            m_mincut[t + 1][t_ + 1] = ComputeSmoothingCost(t, t_);
            m_mincut[0][t_ + 1] += source_c;
            m_mincut[t + 1][m_graph.GetT() + 1] += sink_c;
        }
    }
}

double AlphaExpansion::ComputeSourceSmoothingCost(size_t t, size_t t_) const {
    return m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_labels[t_]}}) -
           m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_alpha}});
}

double AlphaExpansion::ComputeSinkSmoothingCost(size_t t, size_t t_) const {
    return m_graph.GetG({{.t = t, .k = m_alpha}, {.t = t_, .k = m_alpha}}) -
           m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_alpha}});
}

double AlphaExpansion::ComputeSmoothingCost(size_t t, size_t t_) const {
    return m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_alpha}}) +
           m_graph.GetG({{.t = t, .k = m_alpha}, {.t = t_, .k = m_labels[t_]}}) -
           m_graph.GetG({{.t = t, .k = m_alpha}, {.t = t_, .k = m_alpha}}) -
           m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_labels[t_]}}); 
}

void AlphaExpansion::SetDataCosts() {
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        m_mincut[0][t + 1] += m_graph.GetQ({.t = t, .k = m_labels[t]});
        m_mincut[t + 1][m_graph.GetT() + 1] += m_graph.GetQ({.t = t, .k = m_alpha});
    }
}

void AlphaExpansion::MakeSinkSourcePositive() {
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        if (auto& v = m_mincut[0][t + 1]; v < 0) {
            m_mincut[t + 1][m_graph.GetT() + 1] -= v;
            v = 0;
        }
        if (auto& v = m_mincut[t + 1][m_graph.GetT() + 1]; v < 0) {
            m_mincut[0][t + 1] -= v;
            v = 0;
        }
    }
}

void AlphaExpansion::UpdateLabels(const std::vector<size_t>& source) {
    for (auto iter = source.begin() + 1; iter < source.end(); ++iter) {
        m_labels[*iter - 1] = m_alpha;
    }
}

bool AlphaExpansion::IsEnergyDecreased() const {
    return m_energy - ComputeEnergy() > m_eps;
}

double AlphaExpansion::ComputeEnergy() const {
    double energy = 0;
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        for (size_t t_ = t + 1; t_ < m_graph.GetT(); ++t_) {
            energy += m_graph.GetG({{.t = t, .k = m_labels[t]}, {.t = t_, .k = m_labels[t_]}});
        }
        energy += m_graph.GetQ({.t = t, .k = m_labels[t]});
    }
    return energy;
}
