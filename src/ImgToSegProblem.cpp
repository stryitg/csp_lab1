#include <ImgToSegProblem.h>
#include <Eigen/LU>
#include <iostream>

ImgToSegProblem::ImgToSegProblem(const cv::Mat& img, const std::vector<Rect>& class_areas, double eps)
: m_img(img)
, m_class_areas(class_areas)
, m_eps(eps) {}

TabledGraph ImgToSegProblem::GetSegProblem() const {
    const auto classes_points = GetClassesPoints();
    const auto avg = GetClassesAvg(classes_points);
    const auto var_inv = GetClassesVarInv(classes_points, avg);
    return GetSegmentationProblem(avg, var_inv);
}

ImgToSegProblem::ClassesPoints ImgToSegProblem::GetClassesPoints() const {
    ClassesPoints classes_points(m_class_areas.size());
    for (int i = 0; i < m_img.rows; ++i) {
        for (int j = 0; j < m_img.cols; ++j) {
            if (const auto k = GetPointClass({.row = i, .col = j})) {
                classes_points[*k].push_back(CV2Eigen(m_img.at<cv::Vec3b>(i, j)));
            }
        }
    }
    return classes_points;
}

std::optional<size_t> ImgToSegProblem::GetPointClass(const Point& p) const {
    for (size_t i = 0; i <  m_class_areas.size(); ++i) {
        if (IsInRect(p, m_class_areas[i])) {
            return i;
        }
    }
    return std::nullopt;
}

bool ImgToSegProblem::IsInRect(const Point& p, const Rect& r) const {
    return (p.row >= r.up.row && p.row <= r.low.row &&
            p.col >= r.up.col && p.col <= r.low.col);
}

Eigen::Vector3d ImgToSegProblem::CV2Eigen(const cv::Vec3b& vec) const {
    return {static_cast<double>(vec[0]), static_cast<double>(vec[1]), static_cast<double>(vec[2])};
}

std::vector<Eigen::Vector3d> ImgToSegProblem::GetClassesAvg(const ClassesPoints& classes_points) const {
    std::vector<Eigen::Vector3d> avg(classes_points.size(), Eigen::Vector3d::Zero());
    for (size_t i = 0; i < classes_points.size(); ++i) {
        for (const auto& p : classes_points[i]) {
            avg[i] += p;
        }
        avg[i] /= classes_points[i].size();
    }
    return avg;
}

std::vector<Eigen::Matrix3d> ImgToSegProblem::GetClassesVarInv(const ClassesPoints& classes_points,
                                                               const std::vector<Eigen::Vector3d> avg) const {
    std::vector<Eigen::Matrix3d> var_inv(classes_points.size());
    for (size_t i = 0; i < classes_points.size(); ++i) {
        Eigen::Matrix3d var = Eigen::Matrix3d::Zero();
        for (const auto& p : classes_points[i]) {
            const auto val = p - avg[i];
            var += val * val.transpose();
        }
        var /= (classes_points[i].size() - 1);
        var_inv[i] = var.inverse();
    }
    return var_inv;
}

TabledGraph ImgToSegProblem::GetSegmentationProblem(const std::vector<Eigen::Vector3d>& avg,
                                                    const std::vector<Eigen::Matrix3d>& var_inv) const {
    const auto neigbors = GetNeighboringStructure();
    TabledGraph graph(m_img.rows * m_img.cols, m_class_areas.size(), neigbors);
    for (int i = 0; i < m_img.rows; ++i) {
        for (int j = 0; j < m_img.cols; ++j) {
            const auto t = ImgPoint2Graph({.row = i, .col = j});
            for (size_t k = 0; k < avg.size(); ++k) {
                const auto cost = GetCost(CV2Eigen(m_img.at<cv::Vec3b>(i, j)), avg[k], var_inv[k]);
                graph.GetQ({.t = t, .k = k}) = cost;
            }
            for (const auto t_ : graph.GetNeighbors(t)) {
                if (t > t_)
                    continue;
                for (size_t k = 0; k < avg.size(); ++k) {
                    for (size_t k_ = 0; k_ < avg.size(); ++k_) {
                        if (k == k_) {
                            graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) = 0;
                        } else {
                            graph.GetG({{.t = t, .k = k}, {.t = t_, .k = k_}}) = -m_eps;
                        }
                    }
                }
            }
        }
    }
    graph.PropagateQs();
    return graph;
}

TabledGraph::Neighbors ImgToSegProblem::GetNeighboringStructure() const {
    TabledGraph::Neighbors neigbors;
    for (int i = 0; i < m_img.rows; ++i) {
        for (int j = 0; j < m_img.cols; ++j) {
            const auto p = Point{.row = i, .col = j};
            neigbors[ImgPoint2Graph(p)] = GetNeighbors(p);
        }
    }
    return neigbors;
}

std::vector<size_t> ImgToSegProblem::GetNeighbors(const Point& p) const {
    std::vector<size_t> neigbors;
    if (p.row > 0) {
        neigbors.push_back(ImgPoint2Graph({.row = p.row - 1, .col = p.col}));
    }
    if (p.row + 1 < m_img.rows) {
        neigbors.push_back(ImgPoint2Graph({.row = p.row + 1, .col = p.col}));
    }
    if (p.col > 0) {
        neigbors.push_back(ImgPoint2Graph({.row = p.row, .col = p.col - 1}));
    }
    if (p.col + 1 < m_img.cols) {
        neigbors.push_back(ImgPoint2Graph({.row = p.row, .col = p.col + 1}));
    }
    return neigbors;
}

double ImgToSegProblem::GetCost(const Eigen::Vector3d& p, const Eigen::Vector3d& avg,
                                const Eigen::Matrix3d& var_inv) const {
    const auto v = p - avg;
    return -0.5 * v.transpose() * (var_inv * v);
}

cv::Mat ImgToSegProblem::SegSolutionToImg(const std::vector<size_t>& classes) const {
    cv::Mat img(m_img.rows, m_img.cols, CV_8U);
    for (size_t t = 0; t < classes.size(); ++t) {
        const auto p = GraphPoint2Img(t);
        img.at<uint8_t>(p.row, p.col) = (255 / (m_class_areas.size() - 1)) * classes[t];
    }
    return img;
}

size_t ImgToSegProblem::ImgPoint2Graph(const Point& p) const {
    return p.row * m_img.cols + p.col;
}

ImgToSegProblem::Point ImgToSegProblem::GraphPoint2Img(size_t t) const {
    return {.row = static_cast<int>(t / m_img.cols), static_cast<int>(t % m_img.cols)};
}