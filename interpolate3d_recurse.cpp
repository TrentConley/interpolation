#include <iostream>
#include <vector>
#include <cmath>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <set>
#include <unordered_map>
#include <functional>

namespace std {
    template<>
    struct hash<vector<double>> {
        size_t operator()(const vector<double>& v) const {
            size_t seed = v.size();
            for (const auto& i : v) {
                seed ^= hash<double>()(i) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };
}

struct PointND {
    std::vector<double> coordinates;
    double value;
    PointND(const std::vector<double>& coords, double val) : coordinates(coords), value(val) {}
};

struct PointCloud {
    std::vector<PointND> points;
    std::vector<std::vector<double>> uniqueValues; // Unique values for each dimension
    size_t numDimensions; // Number of dimensions
    std::unordered_map<std::vector<double>, double> pointMap; // Point map
};

void extractUniqueValues(const std::vector<PointND>& points, std::vector<std::vector<double>>& uniqueValues) {
    size_t dimensions = points[0].coordinates.size();
    uniqueValues.resize(dimensions);

    for (size_t dim = 0; dim < dimensions; ++dim) {
        std::set<double> valueSet;
        for (const auto& point : points) {
            valueSet.insert(point.coordinates[dim]);
        }
        uniqueValues[dim].assign(valueSet.begin(), valueSet.end());
    }
}

double findLowerBound(const std::vector<double>& vec, double value) {
    auto it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end() || (it != vec.begin() && *it > value)) {
        --it;
    }
    return *it;
}

double getValueAtPoint(const PointCloud& points, const std::vector<double>& queryCoords) {
    auto it = points.pointMap.find(queryCoords);
    if (it != points.pointMap.end()) {
        return it->second;
    }
    throw std::runtime_error("Value not found for given point");
}

double interpolateRecursive(const std::vector<double>& queryPoint, const PointCloud& points, size_t dim) {
    if (dim == 0) {
        return getValueAtPoint(points, queryPoint);
    }

    double lower = findLowerBound(points.uniqueValues[dim - 1], queryPoint[dim - 1]);
    double upper = *std::upper_bound(points.uniqueValues[dim - 1].begin(), points.uniqueValues[dim - 1].end(), lower);

    std::vector<double> lowerCoords = queryPoint;
    std::vector<double> upperCoords = queryPoint;
    lowerCoords[dim - 1] = lower;
    upperCoords[dim - 1] = upper;

    double lowerValue = interpolateRecursive(lowerCoords, points, dim - 1);
    double upperValue = interpolateRecursive(upperCoords, points, dim - 1);

    return (upper - queryPoint[dim - 1]) / (upper - lower) * lowerValue + (queryPoint[dim - 1] - lower) / (upper - lower) * upperValue;
}

double interpolate(const std::vector<double>& queryPoint, const PointCloud& points) {
    return interpolateRecursive(queryPoint, points, points.numDimensions);
}

std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open the file!\n";
        return data;
    }

    std::string line, cell;
    std::getline(file, line); // Skip the header
    while (std::getline(file, line)) {
        std::vector<double> row;
        std::stringstream lineStream(line);

        while (std::getline(lineStream, cell, ',')) {
            row.push_back(std::stod(cell));
        }
        data.push_back(row);
    }

    file.close();
    std::cout << "CSV file read successfully.\n";
    return data;
}

void populatePointCloud(PointCloud& pointCloud, const std::vector<std::vector<double>>& data) {
    pointCloud.points.clear();
    pointCloud.pointMap.clear();
    pointCloud.numDimensions = data[0].size() - 1; // Last column is the dependent variable

    for (const auto& row : data) {
        std::vector<double> coords(row.begin(), row.end() - 1);
        double value = row.back();
        pointCloud.points.emplace_back(coords, value);
        pointCloud.pointMap[coords] = value;
    }

    extractUniqueValues(pointCloud.points, pointCloud.uniqueValues);
    std::cout << "PointCloud populated successfully with " << data.size() << " points.\n";
}

int main() {
    std::string filename = "testing_coeffs_2.csv"; 
    std::vector<std::vector<double>> data = readCSV(filename);

    std::cout << "Populating PointCloud...\n";
    PointCloud pointCloud;
    populatePointCloud(pointCloud, data);

    std::vector<std::vector<double>> queryPoints = {
        {0, 0, 0, 0}, {0.001, 0, 0, 0}, {30.001, 0, 0, 0}, {2, 15, 0, 0},
        {2, 13.43, 13, 12}, {5.67, 8.91, 14.32, 3.45}, {10.23, 7.89, 2.34, 19.56},
        {4.56, 12.34, 6.78, 9.01}, {11.11, 13.13, 14.14, 15.15}, {16.16, 17.17, 18.18, 19.19},
        {1.23, 2.34, 3.45, 4.56}, {7.89, 8.90, 9.01, 10.12}, {13.14, 14.15, 15.16, 16.17},
        {18.19, 19.20, 0.21, 1.22}
    };

    for (const auto& queryPoint : queryPoints) {
        std::cout << "\nQuery point: ";
        for (const auto& val : queryPoint) {
            std::cout << val << " ";
        }
        std::cout << "\n";

        auto start = std::chrono::high_resolution_clock::now();
        double result = interpolate(queryPoint, pointCloud);
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed = end - start;

        std::cout << "Interpolated value: " << result << "\n";
        std::cout << "Time taken: " << elapsed.count() << " seconds.\n";
    }

    std::cout << "All interpolations completed.\n";
    return 0;
}