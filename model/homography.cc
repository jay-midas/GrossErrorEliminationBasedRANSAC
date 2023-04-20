#include "homography.h"

namespace Model {

Homography::Homography() : mSymbol(Eigen::MatrixXd(3, 3)) {}

Homography::~Homography() = default;

double Homography::residual(const std::pair<Data::Point, Data::Point> &points) const {
    auto& [lhsPoint, rhsPoint] = points;
    const auto descriptor_ = mSymbol;
    const double x1 = lhsPoint.x, y1 = lhsPoint.y, x2 = rhsPoint.x, y2 = rhsPoint.y;
    const double t1 = descriptor_(0, 0) * x1 + descriptor_(0, 1) * y1 + descriptor_(0, 2);
    const double t2 = descriptor_(1, 0) * x1 + descriptor_(1, 1) * y1 + descriptor_(1, 2);
    const double t3 = descriptor_(2, 0) * x1 + descriptor_(2, 1) * y1 + descriptor_(2, 2);
    const double d1 = x2 - (t1 / t3);
    const double d2 = y2 - (t2 / t3);
    return sqrt(d1 * d1 + d2 * d2);
}
}