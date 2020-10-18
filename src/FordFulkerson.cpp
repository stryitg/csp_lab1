#include <FordFulkerson.hpp>
#include <queue>
#include <algorithm>
#include <limits>
#include <numeric>

FordFulkerson::FordFulkerson(const AdjacmentMatrix& m)
: m_residuals(m)
, m_path(m.size(), 0)
, m_visited(m.size(), false)
, m_s(0)
, m_t(m.size() - 1) {}

FordFulkerson::STComponents FordFulkerson::ComputeMinCut() {
    while (FindAugmentingPath()) {
        UpdateResiduals();
        ResetBfs();
    }
    return GetSTComponents();
}

bool FordFulkerson::FindAugmentingPath() {
    std::queue<size_t> q;
    q.push(m_s);
    while (!q.empty()) {
        const auto node = q.front();
        q.pop();
        m_visited[node] = true;
        
        for (size_t i = 0; i < m_residuals.size(); ++i) {
            if (m_visited[i]) {
                continue;
            }
            if (m_residuals[node][i] > 0) {
                q.push(i);
                m_path[i] = node;
                if (i == m_t) {
                    return true;
                }
            }
        }
    }
    return m_visited[m_t];
}

void FordFulkerson::UpdateResiduals() {
    const auto flow = GetPathFlow();
    UpdateResiduals(flow);
}

double FordFulkerson::GetPathFlow() const {
    auto node = m_t;
    auto flow = std::numeric_limits<double>::infinity();
    while (node != m_s) {
        const auto prev = m_path[node];
        flow = std::min(flow, m_residuals[prev][node]);
        node = prev;
    }
    return flow;
}

void FordFulkerson::UpdateResiduals(double flow) {
    auto node = m_t;
    while (node != m_s) {
        const auto prev = m_path[node];
        m_residuals[prev][node] -= flow;
        m_residuals[node][prev] += flow;
        node = prev;
    }
}

void FordFulkerson::ResetBfs() {
    std::fill(m_visited.begin(), m_visited.end(), false);
    std::fill(m_path.begin(), m_path.end(), 0);
}

FordFulkerson::STComponents FordFulkerson::GetSTComponents() const {
    std::vector<size_t> vals(m_residuals.size());
    std::iota(vals.begin(), vals.end(), 0);
    std::vector<size_t> s;
    for (size_t i = 0; i < m_residuals.size(); ++i) {
        if (m_visited[i]) {
            s.push_back(i);
        }
    }
    std::vector<size_t> t(m_residuals.size() - s.size());
    std::set_difference(vals.begin(), vals.end(), s.begin(), s.end(), t.begin());
    return {.s = std::move(s), .t = std::move(t)};
}

