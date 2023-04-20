#include <algorithm>
#include <cassert>
#include <numeric>
#include <random>
#include <iostream>
#include <chrono>

#define de(x) std::cerr << #x << " is " << (x) << std::endl;

#include "magsac.h"

namespace {

const double eps = 1e-14;

}

namespace Magsac {

Magsac::Magsac() : maximum_threshold(50), last_iteration_number(0) {}

std::optional<ResultType> Magsac::magsac(
    const Data::PointData& dataSet,
    const double confidence,
    const Model::HomographyEstimator& estimator
) {
	log_confidence = log(1-confidence);
    double bestScore = -1.0;
    Model::Homography bestModel;
    maxIteration = 1'000'000;
	// maxIteration = 10'000;
    uint32_t iteration = 0;

    std::vector<uint32_t> index(dataSet.data.size());
    std::iota(index.begin(), index.end(), 0);

    assert(dataSet.data.size() >= 4);

    while (iteration < maxIteration) {
        ++ iteration;
		de(iteration)
		unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();;
        std::shuffle(index.begin(), index.end(), std::default_random_engine(seed));

        Data::PointData tmpData;
        for (int i = 0; i < 4; ++ i) {
            tmpData.data.emplace_back(dataSet.data[index[i]]);
		}

        auto ret = estimator.estimateMinimalPointModel(tmpData, std::vector<double>(4, 1.0));

        Model::Homography curModel;
        if (ret.has_value()) {
            curModel = ret.value();
        } else {
            continue;
        }

        double _Score = -1.0;
		Model::Homography _Model;

		sigmaConsensus(
			dataSet,
			curModel,
			estimator,
			_Model,
			_Score
		);
		
        if (_Score >= 0 && _Score > bestScore) {
            bestScore = _Score;
            bestModel = curModel;
			maxIteration = std::min(maxIteration, static_cast<uint64_t>(last_iteration_number));
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
	score = -1.0;

	std::for_each(dataSet.data.begin(), dataSet.data.end(), [&](const auto& point) {
        if (const double residual = model.residual(point); residual < maximum_threshold) {
			max_distance = std::max(max_distance, residual);
			all_residuals.emplace_back(residual);
		}
	});

    if (all_residuals.empty()) {
        return ;
    }

	max_distance += eps;

	double thresholds = max_distance;
    // double thresholds_squared = thresholds * thresholds;
    double inliner = 0, probabilities = 1.0;
	for (const double residual : all_residuals) {
		// double residual_squared = residual * residual;
        if (residual < thresholds)
        {
            // double probability = 1.0 - residual_squared / thresholds_squared;
			double probability = 1.0 - residual / thresholds;
            inliner += 1;
            probabilities += probability;
        }
	}

	score = probabilities;
	iteration = log_confidence / log(1.0 - std::pow(inliner / dataSet.data.size(), 4));
}

void Magsac::sigmaConsensus(
    const Data::PointData& dataSet,
    const Model::Homography& model,
    const Model::HomographyEstimator& estimator,
    Model::Homography& betterModel,
    double& betterScore
) {
	constexpr double L = 1.05;
	constexpr double k = 3.64;
	constexpr double threshold_to_sigma_multiplier = 1.0 / k;
	constexpr size_t sample_size = 4;
	static auto comparator = [](std::pair<double, int> left, std::pair<double, int> right) { return left.first < right.first; };
	const int point_number = dataSet.data.size();
	double current_maximum_sigma = this->maximum_threshold;

	std::vector< std::pair<double, size_t> > all_residuals;
	all_residuals.reserve(point_number);

    for (auto point_idx = 0; point_idx < dataSet.data.size(); ++ point_idx)
    {
        auto & points = dataSet.data[point_idx];
        const double residual = model.residual(points);
        if (current_maximum_sigma > residual)
        {
            all_residuals.emplace_back(std::make_pair(residual, point_idx));
        }
    }

	std::vector<Model::Homography> sigma_models;
	std::vector<size_t> sigma_inliers;
	std::vector<double> final_weights;
	
	const size_t possible_inlier_number = all_residuals.size();
	std::sort(all_residuals.begin(), all_residuals.end(), comparator);
	current_maximum_sigma = all_residuals.back().first + eps;

	last_iteration_number = 10000;

	betterScore = -1.0;

	std::vector<double> point_weights_par(possible_inlier_number, 0);
	const double max_sigma = current_maximum_sigma;
	const auto &last_element = std::upper_bound(all_residuals.begin(), all_residuals.end(), std::make_pair(max_sigma, 0), comparator);
	const size_t sigma_inlier_number = last_element - all_residuals.begin();
	sigma_inliers.reserve(sigma_inlier_number);
	for (size_t relative_point_idx = 0; relative_point_idx < sigma_inlier_number; ++relative_point_idx) {
		sigma_inliers.emplace_back(all_residuals[relative_point_idx].second);
	}
	if (sigma_inliers.size() > sample_size)
	{
		std::vector<Model::Homography> sigma_models;
		Data::PointData tmp;
		for (auto id : sigma_inliers) {
			tmp.data.emplace_back(dataSet.data[id]);
		}
		auto ret = estimator.estimateNonMinimalPointModel(
			tmp, std::vector<double>(tmp.data.size(), 1.0)
		);
		if(ret.has_value()) {
			sigma_models.emplace_back(ret.value());
		} else {
			return ;
		}
		const double max_sigma_squared_2 = 2 * max_sigma * max_sigma;
		double residual_i_2, probability_i; 

		for (size_t relative_point_idx = 0; relative_point_idx < sigma_inliers.size(); ++relative_point_idx)
		{
			const size_t &point_idx = sigma_inliers[relative_point_idx];

			residual_i_2 = sigma_models[0].residual(dataSet.data[point_idx]);
			residual_i_2 *= residual_i_2;

			probability_i = exp(-residual_i_2 / max_sigma_squared_2);
			point_weights_par[relative_point_idx] += probability_i;
		}
	}
	sigma_inliers.clear();
	final_weights.reserve(possible_inlier_number);
	sigma_inliers.reserve(possible_inlier_number);
	for (size_t point_idx = 0; point_idx < possible_inlier_number; ++point_idx)
	{
		double weight = point_weights_par[point_idx];
		if (weight < eps)
			continue;
		sigma_inliers.emplace_back(all_residuals[point_idx].second);
		final_weights.emplace_back(weight);

	}
	if (sigma_inliers.size() < sample_size)
		return;

	Data::PointData tmp;
	for (auto id : sigma_inliers) {
		tmp.data.emplace_back(dataSet.data[id]);
	}
    auto ret = estimator.estimateNonMinimalPointModel(
		tmp, final_weights
	);
    if(!ret.has_value()) {
        return;
    }
	double best_score = -1.0, best_marginalized_iteration_number;
	auto sigma_model = ret.value();
	double marginalized_iteration_number;
	
	if(getModelQuality(dataSet, sigma_model, marginalized_iteration_number, best_score); best_score < 0) {
		return ;
	}
	if (best_score > betterScore) {
		betterScore = best_score;
		betterModel = sigma_model;
		best_marginalized_iteration_number = marginalized_iteration_number;
	}

	if (best_marginalized_iteration_number < 0 || std::isnan(best_marginalized_iteration_number)) {
		last_iteration_number = std::numeric_limits<int>::max();
	} else {
		last_iteration_number = static_cast<int>(round(best_marginalized_iteration_number));
	}
	if (sample_size <= 0) {
        betterScore = -1.0;
    }
}

}