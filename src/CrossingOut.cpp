#include <CrossingOut.hpp>
#include <iostream>

bool CrossingOut::Run() {
    m_is_graph_changed = true;
    m_stop = false;
    
    while (m_is_graph_changed && !m_stop) {
        ResetGraphInfo();
        CrossOut();
    }
    if (m_is_not_empty)
        ++m_last_t;
    return m_is_not_empty;
}

void CrossingOut::ResetGraphInfo() {
    m_stop = false;
    m_is_graph_changed = false;
    m_is_not_empty = false;
}

void CrossingOut::CrossOut() {
    CrossOutNodes();
    CrossOutEdges();
}

void CrossingOut::CrossOutNodes() {
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        for (size_t k = 0; k < m_graph.GetK(); ++k) {
            if (m_graph.GetQ({.t = t, .k = k}) != 0 && IsToBeCrossed(t, k)) {
                if (t < m_last_t) {
                    m_stop = true;
                    m_is_not_empty = false;
                    break;
                }
                m_graph.GetQ({.t = t, .k = k}) = 0;
                m_is_graph_changed = true;
            } else if (m_graph.GetQ({.t = t, .k = k}) != 0) {
                m_is_not_empty = true;
            }
        }
        if (m_stop)
            break;
    }
    
}

bool CrossingOut::IsToBeCrossed(size_t t, size_t k) const {
    for (const auto t_ : m_graph.GetNeighbors(t)) { 
        bool dont_cross_out = false;
        for (size_t k_ = 0; k_ < m_graph.GetK(); ++k_) {
            bool edge_val = m_graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) > 0;
            if (t > t_)
                edge_val = m_graph.GetG({{.t = t_, .k = k_}, {.t = t, .k = k}}) > 0;
            dont_cross_out = dont_cross_out || ((m_graph.GetQ({.t = t_, .k = k_}) > 0) && edge_val);
                                  
        }
        if (dont_cross_out == false)
            return true;
    }
    return false;
}

void CrossingOut::CrossOutEdges() {
    if (m_stop)
        return;
        
    for (size_t t = 0; t < m_graph.GetT(); ++t) {
        for (size_t k = 0; k < m_graph.GetK(); ++k) {
            for (const auto t_ : m_graph.GetNeighbors(t)) {
                if (t > t_)
                    continue;
                for (size_t k_ = 0; k_ < m_graph.GetK(); ++k_) {
                    if (m_graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) > 0) {
                        if ((m_graph.GetQ({.t = t, .k = k}) == 0) || (m_graph.GetQ({.t = t_, .k = k_}) == 0)) {
                            m_graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) = 0;
                            m_is_graph_changed = true;
                        }
                    }
                }
            }
        }
    }
}
