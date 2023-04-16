#include "estimator.h"

namespace Model {

HomographyEstimator::HomographyEstimator() = default;

// TODO
Homography HomographyEstimator::EstimatorModel(const Data::PointData &data) const {
    
    return Homography();
}

void HomographyEstimator::gaussElimination(
    Eigen::Matrix<double, 8, 9>& matrix_, // The matrix to which the elimination is applied
    Eigen::Matrix<double, 8, 1>& result_) const // The resulting null-space
{
    constexpr size_t _Size = 8;
    int i, j, k;
    double temp;

    //Pivotisation
    for (i = 0; i < _Size; i++)                    
        for (k = i + 1; k < _Size; k++)
            if (abs(matrix_(i, i)) < abs(matrix_(k, i)))
                for (j = 0; j <= _Size; j++)
                {
                    temp = matrix_(i, j);
                    matrix_(i, j) = matrix_(k, j);
                    matrix_(k, j) = temp;
                }

    //loop to perform the gauss elimination
    for (i = 0; i < _Size - 1; i++)            
        for (k = i + 1; k < _Size; k++)
        {
            double temp = matrix_(k, i) / matrix_(i, i);
            for (j = 0; j <= _Size; j++)
                // make the elements below the pivot elements equal to zero or elimnate the variables
                matrix_(k, j) = matrix_(k, j) - temp * matrix_(i, j);    
        }

    //back-substitution
    for (i = _Size - 1; i >= 0; i--)                
    {                       
        // result_ is an array whose values correspond to the values of x,y,z..
        result_(i) = matrix_(i, _Size);                
        //make the variable to be calculated equal to the rhs of the last equation
        for (j = i + 1; j < _Size; j++)
            if (j != i)            
                //then subtract all the lhs values except the coefficient of the variable whose value is being calculated
                result_(i) = result_(i) - matrix_(i, j) * result_(j);
        //now finally divide the rhs by the coefficient of the variable to be calculated
        result_(i) = result_(i) / matrix_(i, i);            
    }
}

bool HomographyEstimator::estimateMinimalModel(
    const Data::PointData& data,
    std::vector<Homography>& models_,
    const std::vector<double>& weights) const
{
    constexpr size_t equation_number = 2;
    const size_t sample_number_ = data.data.size();
    Eigen::Matrix<double, 8, 9> coefficients;
    size_t row_idx = 0;
    assert(4 == sample_number_);
    for (size_t i = 0; i < sample_number_; ++i)
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
        return false;

    Homography model;
    model.mSymbol << h(0), h(1), h(2),
        h(3), h(4), h(5),
        h(6), h(7), 1.0;
    models_.emplace_back(model);
    return true;
}

bool HomographyEstimator::estimateNonMinimalModel(
    const Data::PointData& data,
    std::vector<Homography>& models_,
    const std::vector<double>& weights) const
{
    constexpr size_t equation_number = 2;
    const size_t sample_number_ = data.data.size();
    const size_t row_number = equation_number * sample_number_;
    Eigen::MatrixXd coefficients(row_number, 8);
    Eigen::MatrixXd inhomogeneous(row_number, 1);

    size_t row_idx = 0;

    for (size_t i = 0; i < sample_number_; ++i)
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
    models_.emplace_back(model);
    return true;
}


}