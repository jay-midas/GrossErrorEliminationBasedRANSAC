#pragma once

#include <Eigen/Eigen>
#include <data/data.h>
#include <iostream>
namespace Model {

class Homography {
public:
    Homography();
    ~Homography();
    Eigen::MatrixXd mSymbol;
    double residual(const std::pair<Data::Point, Data::Point>& points) const;
    void print() const {
        for (int i = 0; i < 3; ++ i) {
            for (int j = 0; j < 3; ++ j) {
                std::cout << mSymbol(i, j) << " ";
            }
            std::cout << std::endl;
        }
    }
};

}