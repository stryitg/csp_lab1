#include <iostream>
#include <array>

#include <opencv2/highgui/highgui.hpp>
#include "opencv2/imgproc/imgproc.hpp"

#include <CrossingOutPath.hpp>
#include <Diffusion.hpp>
#include <MaxPlus.hpp>

constexpr size_t kMaxShift = 5;
constexpr double kOcclusionPenalty = 100000.0;
constexpr double kAlpha = 100.0; 

struct StereoPair {
    cv::Mat l;
    cv::Mat r;
};

struct Point {
    size_t i = std::numeric_limits<size_t>::max();
    size_t j = std::numeric_limits<size_t>::max();
};

cv::Mat ReadImage(const std::string& path) {
    auto img = cv::imread(path, cv::IMREAD_GRAYSCALE);
    if (img.empty()) {
        throw std::runtime_error("Couldn't read image " + path);
    }
    return img;
}

StereoPair ReadStereoPair(const std::string& l_img_path, const std::string& r_img_path) {
    StereoPair sp{.l = ReadImage(l_img_path), .r = ReadImage(r_img_path)};
    if (sp.l.rows != sp.r.rows || sp.l.cols != sp.r.cols) {
        throw std::runtime_error("Images" + l_img_path + " and " + r_img_path + " have different dimensions");
    }
    return sp;
}

size_t PointToT(size_t i, size_t j, size_t w) {
    return i * (w - kMaxShift) + j;
}

Point TToPoint(size_t t, size_t w) {
    return {.i = t / (w - kMaxShift), .j = t % (w - kMaxShift)};
}

void SetOcclusionPenalties(TabledGraph& g) {
    for (size_t t = 0; t < g.GetT(); ++t) {
        g.GetQ({.t = t, .k = 0}) = -2.5 * kOcclusionPenalty;
    }
}

void SetDataPenalties(TabledGraph& g, const StereoPair& sp) {
    for (int i = 0; i < sp.l.rows; ++i) {
        for (int j = 0; j < (int) (sp.l.cols - kMaxShift); ++j) {
            for (size_t k = 0; k < kMaxShift + 1; ++k) {
                g.GetQ({.t = PointToT(i, j, sp.l.cols), .k = k + 1}) = 
                -(sp.r.at<uint8_t>(i, j + k) - sp.l.at<uint8_t>(i, j)) * (sp.r.at<uint8_t>(i, j + k) - sp.l.at<uint8_t>(i, j));
            }
        }
    }
}

void SetQs(TabledGraph& g, const StereoPair& sp) {
    SetOcclusionPenalties(g);
    SetDataPenalties(g, sp);
}


std::array<Point, 4> GetNeighbors(size_t i, size_t j, size_t h, size_t w) {
    std::array<Point, 4> n;
    if (i >= 1) {
        n[0].i = i - 1;
    }
    n[1].i = i;
    if (i < h - 1) {
        n[2].i = i + 1;
    }
    n[3].i = i;
    
    n[0].j = j;
    if (j < w - 1 - kMaxShift) {
        n[1].j = j + 1;
    }
    n[2].j = j;
    if (j >= 1) {
        n[3].j = j - 1;
    }
    return n;
}

void SetGs(TabledGraph&g, const StereoPair& sp) {
    for (int i = 0; i < sp.l.rows; ++i) {
        for (int j = 0; j < (int) (sp.l.cols - kMaxShift); ++j) {
            const auto n = GetNeighbors(i, j, sp.l.rows, sp.l.cols);
            for (const auto p : n) {
                if (p.i == std::numeric_limits<size_t>::max() || p.j == std::numeric_limits<size_t>::max())
                    continue;
                if (PointToT(i, j, sp.l.cols) > PointToT(p.i, p.j, sp.l.cols))
                    continue;
                for (size_t k = 0; k < kMaxShift + 1; ++k) {
                    for (size_t k_ = 0; k_ < kMaxShift + 1; ++k_) {
                        if (k == k_)
                            continue;
                        
                        auto& v = g.GetG({{.t = PointToT(i, j, sp.l.cols), .k = k + 1}, {.t = PointToT(p.i, p.j, sp.l.cols), .k = k_ + 1}});
                        if (std::max(std::abs(sp.r.at<uint8_t>(i, j + k) - sp.r.at<uint8_t>(p.i, p.j + k_)), 
                                     std::abs(sp.l.at<uint8_t>(i, j) - sp.l.at<uint8_t>(p.i, p.j))) > 8) {
                             v = -3 * kAlpha;
                        } else {
                            v = -kAlpha;
                        }
                    }
                }
            }
        }
    }
}

std::vector<size_t> ToTNeighbors(const std::array<Point, 4>& points, size_t w) {
    std::vector<size_t> n;
    for (const auto p : points) {
        if (p.i != std::numeric_limits<size_t>::max() && p.j != std::numeric_limits<size_t>::max())
            n.push_back(PointToT(p.i, p.j, w));
    }
    return n;
}

TabledGraph::Neighbors CreateNeighbors(const StereoPair& sp) {
    TabledGraph::Neighbors n;
    for (int i = 0; i < sp.l.rows; ++i) {
        for (int j = 0; j < sp.l.cols; ++j) {
            const auto points = GetNeighbors(i, j, sp.l.rows, sp.l.cols);
            n[PointToT(i, j, sp.l.cols)] = ToTNeighbors(points, sp.l.cols);
        }
    }
    return n;
}

TabledGraph CreateStereoMatchingProblem(const StereoPair& sp) {
    TabledGraph graph(sp.l.rows * (sp.l.cols - kMaxShift), kMaxShift + 2, CreateNeighbors(sp));
    std::cout << "Started setting Qs" << std::endl;
    SetQs(graph, sp);
    std::cout << "Finished Settings Qs" << std::endl;
    std::cout << "Started Settings Gs" << std::endl;
    SetGs(graph, sp);
    graph.PropagateQs();
    std::cout << "Finished Settings Gs" << std::endl;
    return graph;
}

void WriteDepthMap(const std::vector<size_t> map, const StereoPair& sp) {
    std::cout << "Writing Depth map" << std::endl;
    cv::Mat out(sp.l.rows, sp.l.cols - kMaxShift, CV_8U);
    for(size_t t = 0; t < map.size(); ++t) {
        const auto p = TToPoint(t, sp.l.cols);
        if (map[t] != 0)
            out.at<uint8_t>(p.i, p.j) = static_cast<double>(map[t] - 1) / kMaxShift * 128 + 127;
        else
            out.at<uint8_t>(p.i, p.j) = 0;
    }
    
    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    cv::imwrite("out_.png", out, compression_params);
}

int main(int argc, char** argv) {
    try {
        const auto sp = ReadStereoPair("../images/imL_small.png", "../images/imR_small.png");
        auto g = CreateStereoMatchingProblem(sp);
        MaxPlus mp(g, 500, 1, 0.00001);
        const auto path = mp.ComputePath();
        WriteDepthMap(path, sp);
    } catch (const std::exception& ex) {
        std::cout << ex.what() << std::endl;
    }
    
    return 0;
}