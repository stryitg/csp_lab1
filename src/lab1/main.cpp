#include <TabledGraph.hpp>
#include <FordFulkerson.hpp>
#include <AlphaExpansion.hpp>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <random>
#include <cassert>
#include <fstream>
#include <cmath>
#include <iomanip>

static constexpr size_t kMaxK = 5;

std::vector<std::vector<size_t>> ReadFromFile(const std::string& file_path) {
    std::ifstream fs(file_path);
    if (!fs.is_open()) {
        throw std::runtime_error("Failed to open_file " + file_path);
    }
    size_t h, w;
    fs >> h >> w;
    std::vector<std::vector<size_t>> img(h, std::vector<size_t>(w));
    for (size_t i = 0; i < h; ++i) {
        for (size_t j = 0; j < w; ++j) {
            fs >> img[i][j];
        }
    }
    return img;
}

double g(size_t k1, size_t k2) {
    const double L = 0.5;
    const double S = 1;
    const auto d = (static_cast<double>(k1) - static_cast<double>(k2)) / S;
    return L * std::log(1 + /*0.5*/ d * d);
}

double q(size_t x, size_t k) {
    const auto d = static_cast<double>(x) - static_cast<double>(k);
    return d * d;
}

void SetG(TabledGraph& graph, size_t i, size_t j, size_t h, size_t new_i, size_t new_j) {
    const auto t1 = i * h + j;
    const auto t2 = new_i * h + new_j;
    for (size_t k = 0; k < kMaxK + 1; ++k) {
        for (size_t k_ = 0; k_ < kMaxK + 1; ++k_) {
            graph.GetG({{.t = t1,.k = k}, {.t = t2, .k = k_}}) = g(k, k_);
        }
    }
}

void SetG(TabledGraph& graph, const std::vector<std::vector<size_t>>& img) {
    for (size_t i = 0; i < img.size(); ++i) {
        for (size_t j = 0; j < img[0].size(); ++j) {
            if (i >= 1) {
                SetG(graph, i, j, img.size(), i - 1, j);
            }
            if (i + 1 < img.size()) {
                SetG(graph, i, j, img.size(), i + 1, j);
            }
            if (j >= 1) {
                SetG(graph, i, j, img.size(), i, j - 1);
            }
            if (j + 1 < img[0].size()) {
                SetG(graph, i, j, img.size(), i, j + 1);
            }
        }
    }
}

void SetQ(TabledGraph& graph, const std::vector<std::vector<size_t>>& img) {
    for (size_t i = 0; i < img.size(); ++i) {
        for (size_t j = 0; j < img[0].size(); ++j) {
            for (size_t k = 0; k < kMaxK + 1; ++k) {
                const auto t = i * img.size() + j;
                graph.GetQ({.t = t, .k = k}) = q(img[i][j], k);
            }
        }
    }
}

TabledGraph CreateGraph(const std::vector<std::vector<size_t>>& img) {
    TabledGraph graph(img.size() * img[0].size(), kMaxK + 1);
    SetG(graph, img);
    SetQ(graph, img);
    return graph;
}

template <typename Distribution, typename... Args>
void AddNoise(std::vector<std::vector<size_t>>& img, Args... args) {
    std::random_device rd;
    std::mt19937 gen(rd());
    Distribution d(args...);
    for (auto& row : img) {
        for (auto& v : row) {
            v += d(gen);
            v = v % (kMaxK + 1);
        }
    }
}

struct LaplaceDistribution {
    LaplaceDistribution(double mean, double scale) 
        : m_dist(-0.5, 0.5)
        , m_mean(mean)
        , m_scale(scale) {}
    
    double sgn(double v) { return v > 0 ? 1 : -1; }    
    double operator()(std::mt19937& gen) {
        const auto v = m_dist(gen);  
        return m_mean - m_scale * sgn(v) * std::log(1 - 2 * std::abs(v)); 
    }
    std::uniform_real_distribution<double> m_dist;
    double m_mean;
    double m_scale;
};

void WriteToFile(const std::string& file_path, const std::vector<size_t>& labels, size_t h, size_t w) {
    std::ofstream os(file_path);
    if (!os.is_open()) {
        throw std::runtime_error("Failed to open_file " + file_path);
    }
    os << h << " " << w << std::endl;
    for (size_t i = 0; i < h; ++i) {
        for (size_t j = 0; j < w; ++j) {
            os << labels[i * h + j] << std::setw(5);
        }
        os << std::setw(0) << std::endl;
    }
}

std::vector<size_t> BinarizationAlgo(const TabledGraph& gr) {
    const auto m = gr.ToMincut();
    FordFulkerson ff(m);
    const auto mc = ff.ComputeMincut();
    return gr.ToLabels(mc);
}

std::vector<size_t> AlphaExpansionAlgo(const TabledGraph& gr) {
    AlphaExpansion ae(gr);
    return ae.ComputeLabels();
}

int main() {
    auto img = ReadFromFile("img.txt");
    // AddNoise<std::normal_distribution<double>>(img, static_cast<double>(kMaxK) / 2, static_cast<double>(kMaxK) / 3);
    // AddNoise<LaplaceDistribution>(img, static_cast<double>(kMaxK) / 2, static_cast<double>(kMaxK) / 3);
    // AddNoise<std::uniform_int_distribution<size_t>>(img, static_cast<double>(kMaxK) / 2, static_cast<double>(kMaxK) / 3);
    const auto gr = CreateGraph(img);
    
    const auto labels_bin = BinarizationAlgo(gr);
    WriteToFile("out_binarization.txt", labels_bin, img.size(), img[0].size());
    
    const auto labels_alpha = AlphaExpansionAlgo(gr);
    WriteToFile("out_alpha.txt", labels_alpha, img.size(), img[0].size());
    
    return 0;
}