#include <iostream>
#include <string>
#include <ImgToSegProblem.h>
#include <MaxPlus.hpp>

cv::Mat ReadImage(const std::string& path) {
    auto img = cv::imread(path);
    if (img.empty()) {
        throw std::runtime_error("Couldn't read image " + path);
    }
    return img;
}

void WriteImage(const cv::Mat& out, const std::string& path) {
    std::vector<int> compression_params;
    compression_params.push_back(cv::IMWRITE_PNG_COMPRESSION);
    compression_params.push_back(9);
    cv::imwrite(path, out, compression_params);
}

int main() {
    try {
        const auto image = ReadImage("../images/forest_mountains.jpeg");
        // For rock.jpeg
        // ImgToSegProblem sp(image, {{.up = {.row = 0, .col = 156}, .low = {.row = 60, .col = image.cols - 1}},
        //                            {.up = {.row = 50, .col = 120}, .low = {.row = 110, .col = 150}},
        //                            {.up = {.row = 135, .col = 0}, .low = {.row = 170, .col = image.cols - 1}}}, 30);
        // For field_small.jpg
        // ImgToSegProblem sp(image, {{.up = {.row = 0, .col = 0}, .low = {.row = 50, .col = image.cols - 1}},
        //                            {.up = {.row = 50, .col = 0}, .low = {.row = image.rows - 1, .col = image.cols}}}, 30);
        // For forest_mountains.jpeg
        ImgToSegProblem sp(image, {{.up = {.row = 65, .col = 58}, .low = {.row = 86, .col = 74}},
                                   {.up = {.row = 55, .col = 4}, .low = {.row = 83, .col = 45}},
                                   {.up = {.row = 30, .col = 100}, .low = {.row = 45, .col = image.cols - 1}},
                                   {.up = {.row = 30, .col = 20}, .low = {.row = 25, .col = 88}},
                                   {.up = {.row = 0, .col = 0}, .low = {.row = 5, .col = 63}}}, 5);
        const auto graph = sp.GetSegProblem();
        MaxPlus mp(graph, 200, 1, 0.00001);
        const auto classes = mp.ComputePath();
        const auto out_img = sp.SegSolutionToImg(classes);
        WriteImage(out_img, "../images/forest_mountains_eps5.png");
    } catch (const std::exception& ex) {
        std::cout << ex.what() << std::endl;
    }
}