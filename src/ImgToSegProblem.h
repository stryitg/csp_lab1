#pragma once

#include <optional>

#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <Eigen/Core>

#include <TabledGraph.hpp>

class ImgToSegProblem {
public:
    struct Point {
        int row;
        int col;
    };

    struct Rect {
        Point up;
        Point low;
    };

    ImgToSegProblem(const cv::Mat& img, const std::vector<Rect>& class_areas, const double eps);
    TabledGraph GetSegProblem() const;
    cv::Mat SegSolutionToImg(const std::vector<size_t>& classes) const;
    
private:
    using ClassesPoints = std::vector<std::vector<Eigen::Vector3d>>;
    
    ClassesPoints GetClassesPoints() const;
    std::optional<size_t> GetPointClass(const Point& p) const;
    bool IsInRect(const Point& p, const Rect& r) const;
    Eigen::Vector3d CV2Eigen(const cv::Vec3b& vec) const;
    
    std::vector<Eigen::Vector3d> GetClassesAvg(const ClassesPoints& classes_points) const;
    std::vector<Eigen::Matrix3d> GetClassesVarInv(const ClassesPoints& classes_points,
                                                    const std::vector<Eigen::Vector3d> avg) const;
    
    TabledGraph GetSegmentationProblem(const std::vector<Eigen::Vector3d>& avg,
                                       const std::vector<Eigen::Matrix3d>& var_inv) const;
    TabledGraph::Neighbors GetNeighboringStructure() const;
    std::vector<size_t> GetNeighbors(const Point& p) const;
    double GetCost(const Eigen::Vector3d& p, const Eigen::Vector3d& avg, const Eigen::Matrix3d& var_inv) const;

    size_t ImgPoint2Graph(const Point& p) const;
    Point GraphPoint2Img(size_t t) const;
    
private:
    const cv::Mat m_img;
    const std::vector<Rect> m_class_areas;
    
    const double m_eps;
};