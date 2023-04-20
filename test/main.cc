#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <chrono>

#include "data/data.h"
#include "magsac/magsac.h"
#include "model/estimator.h"

#define de(x) std::cerr << #x << " is " << (x) << std::endl;

Data::PointData GetTruePoint(const Data::PointData& data, std::vector<int>& isTrue) {
    Data::PointData mTruthPoints;
    Model::Homography minModel;
    int cnt = 0;
    while (1) {
        ++ cnt;
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();;
        std::shuffle(isTrue.begin(), isTrue.end(), std::default_random_engine(seed));
        mTruthPoints.data.clear();
        for (auto index : isTrue) {
            mTruthPoints.data.emplace_back(data.data[index]);
            if (mTruthPoints.data.size() == 4) break;
        }
        assert(mTruthPoints.data.size() == 4);
        auto ret1 = Model::HomographyEstimator().estimateMinimalPointModel(mTruthPoints, std::vector<double>(mTruthPoints.data.size(), 1.0));
        if (ret1.has_value()) {
            minModel = ret1.value();
            break;
        }
        assert(cnt <= 10'000);
    }
    mTruthPoints.data.clear();
    for (auto x : data.data) {
        if (minModel.residual(x) < 2.0) {
            mTruthPoints.data.emplace_back(x);
        }
    }
    return mTruthPoints;
}

std::vector<int> ReadData(Data::PointData& data) {
    std::fstream f;
    f.open("data.txt", std::ios_base::in);
    if (!f.is_open()) {
        exit(-1);
    }
    std::vector<int> isTrue;
    do {
        double x1, y1, z1, x2, y2, z2;
        std::string p;
        int op, i = 0;
        // while (f >> x1 >> x2 >> y1 >> y2 >> z1 >> z2 >> op) {
        while (f >> x1 >> y1 >> x2 >> y2 >> z1 >> z2 >> p >> p >> op) {
            // if (op == 0) continue;
            data.data.emplace_back(std::make_pair(Data::Point(x1, y1, z1), Data::Point(x2, y2, z2)));
            if (op) {
                isTrue.emplace_back(i);
            }
            i ++;
        }
    } while(0);
    assert(isTrue.size() >= 4);
    return isTrue;
}

int main() {
    Data::PointData data;
    std::vector<int> isTrue = ReadData(data);
    // generate truth point
    Data::PointData mTruthPoints = GetTruePoint(data, isTrue);
    // generate point
    Magsac::Magsac magsac;
    auto ret = magsac.magsac(data, 0.99, Model::HomographyEstimator());
    if (ret.has_value()) {
        auto [model, score] = ret.value();
        model.print();
        double rmse = 0; 
        for (auto point : mTruthPoints.data) {
            double err = model.residual(point);
            rmse += err * err;
        }
        rmse = sqrt(rmse / static_cast<double>(mTruthPoints.data.size()));
        std::cout  << "RMSE error = " << rmse << " px" << std::endl;
    } else {
        std::cerr << "Fail." << std::endl;
        exit(-1);
    }
    return 0;
}