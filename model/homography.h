#pragma once

#include <Eigen/Eigen>

namespace Model {

class Homography {
public:
    Homography();
    ~Homography();
private:
    Eigen::MatrixXd mSymbol;
};

}