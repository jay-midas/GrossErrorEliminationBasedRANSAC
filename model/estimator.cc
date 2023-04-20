#include "estimator.h"


namespace Model {

bool normalizePoints(
        const Data::PointData &data, 
        Data::PointData& dstData,
        Eigen::Matrix3d &normalizing_transform_source_, 
        Eigen::Matrix3d &normalizing_transform_destination_) 
    {
        const size_t cols = data.data.size();
        double mass_point_src[2], 
            mass_point_dst[2];

        mass_point_src[0] =
            mass_point_src[1] =
            mass_point_dst[0] =
            mass_point_dst[1] =
            0.0;

        for (const auto& [lhs, rhs] : data.data)
        {
            mass_point_src[0] += lhs.x;
            mass_point_src[1] += lhs.y;
            mass_point_dst[0] += rhs.x;
            mass_point_dst[1] += rhs.y;
        }

        mass_point_src[0] /= data.data.size();
        mass_point_src[1] /= data.data.size();
        mass_point_dst[0] /= data.data.size();
        mass_point_dst[1] /= data.data.size();

        double average_distance_src = 0.0,
            average_distance_dst = 0.0;
        for (const auto& [lhs, rhs] : data.data)
        {
            const double &x1 = lhs.x;
            const double &y1 = lhs.y;
            const double &x2 = rhs.x;
            const double &y2 = rhs.y;

            const double dx1 = mass_point_src[0] - x1;
            const double dy1 = mass_point_src[1] - y1;
            const double dx2 = mass_point_dst[0] - x2;
            const double dy2 = mass_point_dst[1] - y2;

            average_distance_src += sqrt(dx1 * dx1 + dy1 * dy1);
            average_distance_dst += sqrt(dx2 * dx2 + dy2 * dy2);
        }

        average_distance_src /= data.data.size();
        average_distance_dst /= data.data.size();

        const double ratio_src = M_SQRT2 / average_distance_src;
        const double ratio_dst = M_SQRT2 / average_distance_dst;

        dstData.data.clear();
        for (const auto& [lhs, rhs] : data.data)
        {
            const double &x1 = lhs.x;
            const double &y1 = lhs.y;
            const double &x2 = rhs.x;
            const double &y2 = rhs.y;

            dstData.data.emplace_back(
                Data::Point(
                    (x1 - mass_point_src[0]) * ratio_src,
                    (y1 - mass_point_src[1]) * ratio_src
                ),
                Data::Point(
                (x2 - mass_point_dst[0]) * ratio_dst,
                (y2 - mass_point_dst[1]) * ratio_dst
                )
            );
        }

        normalizing_transform_source_ << ratio_src, 0, -ratio_src * mass_point_src[0],
            0, ratio_src, -ratio_src * mass_point_src[1],
            0, 0, 1;

        normalizing_transform_destination_ << ratio_dst, 0, -ratio_dst * mass_point_dst[0],
            0, ratio_dst, -ratio_dst * mass_point_dst[1],
            0, 0, 1;
        return true;
    }
HomographyEstimator::HomographyEstimator() = default;

void HomographyEstimator::gaussElimination(
    Eigen::Matrix<double, 8, 9>& matrix_, 
    Eigen::Matrix<double, 8, 1>& result_) const
{
    constexpr size_t _Size = 8;
    int i, j, k;
    double temp;

    for (i = 0; i < _Size; i++)                    
        for (k = i + 1; k < _Size; k++)
            if (abs(matrix_(i, i)) < abs(matrix_(k, i)))
                for (j = 0; j <= _Size; j++)
                {
                    temp = matrix_(i, j);
                    matrix_(i, j) = matrix_(k, j);
                    matrix_(k, j) = temp;
                }
    for (i = 0; i < _Size - 1; i++)            
        for (k = i + 1; k < _Size; k++)
        {
            double temp = matrix_(k, i) / matrix_(i, i);
            for (j = 0; j <= _Size; j++)
                matrix_(k, j) = matrix_(k, j) - temp * matrix_(i, j);    
        }
    for (i = _Size - 1; i >= 0; i--)                
    {                       
        result_(i) = matrix_(i, _Size);                
        for (j = i + 1; j < _Size; j++)
            if (j != i)            
                result_(i) = result_(i) - matrix_(i, j) * result_(j);
        result_(i) = result_(i) / matrix_(i, i);            
    }
}

std::optional<Homography> HomographyEstimator::estimateMinimalPointModel(
    const Data::PointData& data,
    const std::vector<double>& weights) const
{
    constexpr size_t equation_number = 2;
    assert(data.data.size() == 4);
    assert(data.data.size() == weights.size());
    Eigen::Matrix<double, 8, 9> coefficients;
    size_t row_idx = 0;
    for (size_t i = 0; i < data.data.size(); ++i) {
        const double weight = weights[i];
        const double
            & x1 = data.data[i].first.x,
            & y1 = data.data[i].first.y,
            & x2 = data.data[i].second.x,
            & y2 = data.data[i].second.y;

        const double
            minus_weight_times_x1 = -weight * x1,
            minus_weight_times_y1 = -weight * y1,
            weight_times_x2 = weight * x2,
            weight_times_y2 = weight * y2;

        coefficients(row_idx, 0) = minus_weight_times_x1;
        coefficients(row_idx, 1) = minus_weight_times_y1;
        coefficients(row_idx, 2) = -weight;
        coefficients(row_idx, 3) = 0;
        coefficients(row_idx, 4) = 0;
        coefficients(row_idx, 5) = 0;
        coefficients(row_idx, 6) = weight_times_x2 * x1;
        coefficients(row_idx, 7) = weight_times_x2 * y1;
        coefficients(row_idx, 8) = -weight_times_x2;

        ++row_idx;

        coefficients(row_idx, 0) = 0;
        coefficients(row_idx, 1) = 0;
        coefficients(row_idx, 2) = 0;
        coefficients(row_idx, 3) = minus_weight_times_x1;
        coefficients(row_idx, 4) = minus_weight_times_y1;
        coefficients(row_idx, 5) = -weight;
        coefficients(row_idx, 6) = weight_times_y2 * x1;
        coefficients(row_idx, 7) = weight_times_y2 * y1;
        coefficients(row_idx, 8) = -weight_times_y2;
        ++row_idx;
    }

    Eigen::Matrix<double, 8, 1> h;
    gaussElimination(
        coefficients,
        h);
    if (h.hasNaN())
        return std::nullopt;

    Homography model;
    model.mSymbol << h(0), h(1), h(2),
        h(3), h(4), h(5),
        h(6), h(7), 1.0;
    return model;
}

std::optional<Homography> HomographyEstimator::estimateNonMinimalPointModel(
    const Data::PointData &data, 
    const std::vector<double> &weights
) const {
        Data::PointData dstData;
        Eigen::Matrix3d normalizing_transform_source, normalizing_transform_destination;

        if (!normalizePoints(data,
            dstData,
            normalizing_transform_source, 
            normalizing_transform_destination))
            return std::nullopt;
        Model::Homography ret;
        if (auto x = estimateFullPointModel(dstData, weights); x.has_value()) {
            ret = x.value();
        } else {
            return std::nullopt;
        }
        const Eigen::Matrix3d normalizing_transform_destination_inverse = normalizing_transform_destination.inverse();
        ret.mSymbol = normalizing_transform_destination_inverse * ret.mSymbol * normalizing_transform_source;
        return ret;
}

std::optional<Homography> HomographyEstimator::estimateFullPointModel(
    const Data::PointData& data,
    const std::vector<double>& weights) const
{
    assert(weights.size() == data.data.size());
    constexpr size_t equation_number = 2;
    const size_t row_number = equation_number * data.data.size();
    Eigen::MatrixXd coefficients(row_number, 8);
    Eigen::MatrixXd inhomogeneous(row_number, 1);

    size_t row_idx = 0;

    for (size_t i = 0; i < data.data.size(); ++i)
    {
        const double weight = weights[i];
        const double
            & x1 = data.data[i].first.x,
            & y1 = data.data[i].first.y,
            & x2 = data.data[i].second.x,
            & y2 = data.data[i].second.y;

        const double
            minus_weight_times_x1 = -weight * x1,
            minus_weight_times_y1 = -weight * y1,
            weight_times_x2 = weight * x2,
            weight_times_y2 = weight * y2;

        coefficients(row_idx, 0) = minus_weight_times_x1;
        coefficients(row_idx, 1) = minus_weight_times_y1;
        coefficients(row_idx, 2) = -weight;
        coefficients(row_idx, 3) = 0;
        coefficients(row_idx, 4) = 0;
        coefficients(row_idx, 5) = 0;
        coefficients(row_idx, 6) = weight_times_x2 * x1;
        coefficients(row_idx, 7) = weight_times_x2 * y1;
        inhomogeneous(row_idx) = -weight_times_x2;
        ++row_idx;

        coefficients(row_idx, 0) = 0;
        coefficients(row_idx, 1) = 0;
        coefficients(row_idx, 2) = 0;
        coefficients(row_idx, 3) = minus_weight_times_x1;
        coefficients(row_idx, 4) = minus_weight_times_y1;
        coefficients(row_idx, 5) = -weight;
        coefficients(row_idx, 6) = weight_times_y2 * x1;
        coefficients(row_idx, 7) = weight_times_y2 * y1;
        inhomogeneous(row_idx) = -weight_times_y2;
        ++row_idx;
    }

    Eigen::Matrix<double, 8, 1> 
        h = coefficients.colPivHouseholderQr().solve(inhomogeneous);

    Homography model;
    model.mSymbol << h(0), h(1), h(2),
        h(3), h(4), h(5),
        h(6), h(7), 1.0;
    return model;
}
}