#pragma once

#include <Eigen/Eigen>

namespace Model {

class Homography {
public:
    Homography();
    ~Homography();
    Eigen::MatrixXd mSymbol;
};

}