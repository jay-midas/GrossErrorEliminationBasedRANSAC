#pragma once

#include <data/data.h>

#include "homography.h"

#include <iostream>
#define de(x) std::cerr << #x << " is " << (x) << std::endl;

namespace Model {

class HomographyEstimator {
public:
    HomographyEstimator();
    void gaussElimination(
        Eigen::Matrix<double, 8, 9>& matrix, 
        Eigen::Matrix<double, 8, 1>& result
    ) const;
    std::optional<Homography> estimateMinimalPointModel(
        const Data::PointData& data,
        const std::vector<double>& weights
    ) const;
    std::optional<Homography> estimateNonMinimalPointModel(
        const Data::PointData& data,
        const std::vector<double>& weights
    ) const;
    std::optional<Homography> estimateFullPointModel(
        const Data::PointData& data,
        const std::vector<double>& weights
    ) const;
};

}