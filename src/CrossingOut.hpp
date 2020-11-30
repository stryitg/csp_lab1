#pragma once

#include <TabledGraph.hpp>

class CrossingOut {
public:
    bool Run();
    
    void SetGraph(const TabledGraph& graph) { m_graph = graph; }
    const TabledGraph& GetGraph() const { return m_graph; };
    
private:
    void ResetGraphInfo();
    void CrossOut();
    void CrossOutNodes();
    bool IsToBeCrossed(size_t t, size_t k) const;
    void CrossOutEdges();
    
private:
    TabledGraph m_graph;
    
    bool m_is_graph_changed = true;
    bool m_is_not_empty;
    bool m_stop = false;
    
    size_t m_last_t = 0;
    
};