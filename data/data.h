#pragma once

#include <algorithm>
#include <utility>
#include <vector>
#include <iostream>
#include <cstdio>
#include <functional>

namespace Data {

struct Point {
    double x, y, z;
    Point(double _x = 0, double _y = 0, double _z = 0);
    void print() const {
        printf("{ x=%.03f y=%.03f z=%.03f } \n", x, y, z);
    }
};

class PointData {
public:
    PointData();
    std::vector<std::pair<Point, Point>> data;
    void print() const {
        std::for_each(data.begin(), data.end(), [](const auto& x) {
            std::cout << "\n\n --- " << std::endl;
            x.first.print();
            x.second.print();
            std::cout << " --- \n\n" << std::endl;
        });
    }
};

}