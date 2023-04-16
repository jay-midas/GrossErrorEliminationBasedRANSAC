#pragma once

#include <data/data.h>

#include "homography.h"

namespace Model {

class HomographyEstimator {
public:
    HomographyEstimator();
    Homography EstimatorModel(const Data::PointData& data) const ;
    void gaussElimination(
    Eigen::Matrix<double, 8, 9>& matrix_, // The matrix to which the elimination is applied
    Eigen::Matrix<double, 8, 1>& result_) const; // The resulting null-space

    bool estimateMinimalModel(
        const Data::PointData& data,
        std::vector<Homography>& models_,
        const std::vector<double>& weights) const;

    bool estimateNonMinimalModel(
        const Data::PointData& data,
        std::vector<Homography>& models_,
        const std::vector<double>& weights) const;
};

}