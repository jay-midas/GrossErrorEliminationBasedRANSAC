#include "homography.h"

namespace Model {

Homography::Homography() : mSymbol(Eigen::MatrixXd(3, 3)) {}

Homography::~Homography() = default;

}