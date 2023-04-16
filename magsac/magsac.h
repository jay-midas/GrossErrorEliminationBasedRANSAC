#include <data/data.h>
#include <model/estimator.h>
#include <optional>
#include <tuple>


namespace Magsac {

using ResultType = std::tuple<Model::Homography, double>;

class Magsac {
public:
    Magsac();
    std::optional<ResultType> magsac(
        const Data::PointData& dataSet,
        const double confidence,
        const Model::HomographyEstimator& estimator
    );
    void sigmaConsensus(
        const Data::PointData& dataSet,
        const Model::Homography& model,
        const Model::HomographyEstimator& estimator
    );
    void getModelQuality(
        const Data::PointData& dataSet,
        const Model::Homography& model,
        double &iteration,
        double &score
    );
    std::vector<double> mWeights;    
    uint64_t maxIteration;
    double maximum_threshold;
    double log_confidence;
};

}