#pragma once

#include <TabledGraph.hpp>
#include <FordFulkerson.hpp>
#include <random>

class AlphaExpansion {
public:
    AlphaExpansion(const TabledGraph& graph, double eps = 0.0);
    
    std::vector<size_t> ComputeLabels();
    
private:
    std::vector<size_t> CreateRandomAlphaOrdering() const;
    double ComputeEnergy() const;
    void MakeAlphaExpansion();
    void UpdateAlpha();
    void UpdateMincutProblem();
    void ResetMincut();
    void SetSmoothingCosts();
    double ComputeSourceSmoothingCost(size_t t, size_t t_) const;
    double ComputeSinkSmoothingCost(size_t t, size_t t_) const;
    double ComputeSmoothingCost(size_t t, size_t t_) const;
    void SetDataCosts();
    void MakeSinkSourcePositive();
    void UpdateLabels(const std::vector<size_t>& source);
    bool IsEnergyDecreased() const;
    
private:
    TabledGraph m_graph;
    std::vector<size_t> m_labels;
    
    double m_energy;
    double m_eps;
    
    AdjacmentMatrix m_mincut;
    
    std::vector<size_t> m_alpha_ordering;
    size_t m_alpha;
};