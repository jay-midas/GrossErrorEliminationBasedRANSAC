#pragma once

#include <utility>
#include <vector>

namespace Data {

struct Point {
    double x, y, z;
    Point(double _x = 0, double _y = 0, double _z = 0);
};

class PointData {
public:
    PointData();
    std::vector<std::pair<Point, Point>> data;
};

}