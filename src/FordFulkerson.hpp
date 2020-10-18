#pragma once

#include <Graph.hpp>

class FordFulkerson {
public:
    struct STComponents {
        std::vector<size_t> s;
        std::vector<size_t> t;
    };
    
    FordFulkerson(const AdjacmentMatrix& m);
    STComponents ComputeMinCut();
    
private:
    bool FindAugmentingPath();
    void UpdateResiduals();
    double GetPathFlow() const;
    void UpdateResiduals(double flow);
    void ResetBfs();
    STComponents GetSTComponents() const;
    
private:
    AdjacmentMatrix m_residuals;
    std::vector<size_t> m_path;
    std::vector<bool> m_visited;
    
    const size_t m_s;
    const size_t m_t;
};