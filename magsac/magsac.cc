#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>

#include "magsac.h"

namespace Magsac {

Magsac::Magsac() = default;

std::optional<ResultType> Magsac::magsac(
    const Data::PointData& dataSet,
    const double confidence,
    const Model::HomographyEstimator& estimator
) {

    double bestScore = -1.0;
    Model::Homography bestModel;
    maxIteration = std::numeric_limits<uint32_t>::max();
    uint32_t iteration = 0;

    std::vector<uint32_t> index(dataSet.data.size());
    std::iota(index.begin(), index.end(), 0);

    assert(dataSet.data.size() >= 4);
    

    while (iteration < maxIteration) {
        
        ++ iteration;

        std::shuffle(index.begin(), index.end(), std::default_random_engine());

        Data::PointData tmpData;
        
        for (int i = 0; i < 4; ++ i) {
            tmpData.data.emplace_back(dataSet.data[index[i]]);
        }

        auto curModel = estimator.EstimatorModel(tmpData);
        double curScore = -1.0;
        double curIteration;
        // TODO add check to remove useless model
        // TODO add sigma consensus
        getModelQuality(dataSet,curModel,curIteration, curScore);
        if (curScore < 0) {
            continue;
        }
        if (bestScore < 0 || bestScore > curScore) {
            bestScore = curScore;
            bestModel = curModel;
        }
    }

    if (bestScore < 0.0) {
        return std::nullopt;
    }
    return std::make_tuple(bestModel, bestScore);
}

void Magsac::getModelQuality(
	const Data::PointData& dataSet,
    const Model::Homography& model,
	double& iteration,
	double& score) 
{
	std::vector<double> all_residuals;

	double max_distance = 0;
	for (const auto& [lhsPoint, rhsPoint] : dataSet.data) {
        const auto descriptor_ = model.mSymbol;
        const double x1 = lhsPoint.x, y1 = lhsPoint.y, x2 = rhsPoint.x, y2 = rhsPoint.y;
        const double t1 = descriptor_(0, 0) * x1 + descriptor_(0, 1) * y1 + descriptor_(0, 2);
        const double t2 = descriptor_(1, 0) * x1 + descriptor_(1, 1) * y1 + descriptor_(1, 2);
        const double t3 = descriptor_(2, 0) * x1 + descriptor_(2, 1) * y1 + descriptor_(2, 2);
        const double d1 = x2 - (t1 / t3);
        const double d2 = y2 - (t2 / t3);

		const double residual = sqrt(d1 * d1 + d2 * d2);
		if (maximum_threshold > residual)
		{
			max_distance = std::max(max_distance, residual);
			all_residuals.emplace_back(residual);
		}
	}

    if (all_residuals.empty()) {
        score = -1.0;
        return ;
    }

	max_distance += std::numeric_limits<double>::epsilon();

	double thresholds = max_distance;
    double thresholds_squared = thresholds * thresholds;
    double thresholds_2_squared = 2.0 * thresholds_squared;

    double inliner = 0, probabilities = 1.0;
	for (const double residual : all_residuals)
	{
		double residual_squared = residual * residual;

        if (residual < thresholds)
        {
            double probability = 1.0 - residual_squared / thresholds_squared;
            inliner += 1;
            probabilities += probability;
        }
	}

	score = probabilities;
	iteration = log_confidence / log(1.0 - std::pow(inliner / all_residuals.size(), 4));
}

}

int main() {

    return 0;
}