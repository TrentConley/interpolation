#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <chrono>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>
#include <set>
#include <unordered_map>

#include <functional> // Add this include

// Add this custom hash function before the PointCloud struct
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

// Define a custom point cloud data structure
struct PointCloud {
    std::vector<PointND> pts;
    std::vector<std::vector<double>> uniqueVals; // Unique values for each dimension
    size_t nDims; // Number of dimensions
    std::unordered_map<std::vector<double>, double> point_map; // New point map
};

// Helper function to extract unique values from points
void get_unique_values(const std::vector<PointND>& points, std::vector<std::vector<double>>& unique_values) {
    size_t dimensions = points[0].coordinates.size();
    unique_values.resize(dimensions);

    for (size_t dim = 0; dim < dimensions; ++dim) {
        std::set<double> value_set;
        for (const auto& point : points) {
            value_set.insert(point.coordinates[dim]);
        }
        unique_values[dim].assign(value_set.begin(), value_set.end());
    }
}

// Helper function to find the largest value less than or equal to the query point
double find_lower_bound(const std::vector<double>& vec, double value) {
    auto it = std::lower_bound(vec.begin(), vec.end(), value);
    if (it == vec.end() || (it != vec.begin() && *it > value)) {
        --it;
    }
    return *it;
}

// Helper function to find the value at a given point
double find_value_at_point(const PointCloud& points, const std::vector<double>& query_coords) {
    auto it = points.point_map.find(query_coords);
    if (it != points.point_map.end()) {
        return it->second;
    }
    throw std::runtime_error("Value not found for given point");
}

// Recursive interpolation function
double interpolate_recursive(const std::vector<double>& query_point, const PointCloud& points, size_t dim) {
    
    if (dim == 0) {
        return find_value_at_point(points, query_point);
    }

    double lower = find_lower_bound(points.uniqueVals[dim - 1], query_point[dim - 1]);
    double upper = *std::upper_bound(points.uniqueVals[dim - 1].begin(), points.uniqueVals[dim - 1].end(), lower);

    std::vector<double> lower_coords = query_point;
    std::vector<double> upper_coords = query_point;
    lower_coords[dim - 1] = lower;
    upper_coords[dim - 1] = upper;

    double lower_value = interpolate_recursive(lower_coords, points, dim - 1);
    double upper_value = interpolate_recursive(upper_coords, points, dim - 1);
    // std::cout << "Lower bound for dimension " << dim-1 << ": " << lower << std::endl;
    // std::cout << "Upper bound for dimension " << dim-1 << ": " << upper << std::endl;

    return (upper - query_point[dim - 1]) / (upper - lower) * lower_value + (query_point[dim - 1] - lower) / (upper - lower) * upper_value;
}

double interpolate(const std::vector<double>& query_point, const PointCloud& points) {

    return interpolate_recursive(query_point, points, points.nDims);
}

// Function to read a CSV file
std::vector<std::vector<double>> readCSV(const std::string& filename) {
    std::vector<std::vector<double>> data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Could not open the file!\n";
        return data;
    }

    std::string line, cell;
    // Skip the header
    std::getline(file, line);
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

// Function to populate the PointCloud structure from data
void populatePointCloud(PointCloud& pointCloud, const std::vector<std::vector<double>>& data) {
    pointCloud.pts.clear();
    pointCloud.point_map.clear(); // Clear the point map
    pointCloud.nDims = data[0].size() - 1; // Last column is the dependent variable

    for (const auto& row : data) {
        std::vector<double> coords(row.begin(), row.end() - 1);
        double value = row.back();
        pointCloud.pts.emplace_back(coords, value);
        pointCloud.point_map[coords] = value; // Add to point map
    }

    // Determine the unique values for each dimension using the helper function
    get_unique_values(pointCloud.pts, pointCloud.uniqueVals);

    std::cout << "PointCloud populated successfully with " << data.size() << " points.\n";
}

int main() {
    // Use dummy data instead of reading from CSV file
    std::string filename = "testing_coeffs_2.csv"; // Replace with your CSV file path

    std::vector<std::vector<double>> data = readCSV(filename);

    // Generate a point cloud based on the dummy data
    std::cout << "Populating PointCloud...\n";
    PointCloud pointCloud;
    populatePointCloud(pointCloud, data);

    // Print the PointCloud
    // std::cout << "PointCloud contents:\n";
    // for (const auto& point : pointCloud.pts) {
    //     std::cout << "Coordinates: [";
    //     for (size_t i = 0; i < point.coordinates.size(); ++i) {
    //         std::cout << point.coordinates[i];
    //         if (i < point.coordinates.size() - 1) {
    //             std::cout << ", ";
    //         }
    //     }
    //     std::cout << "], Value: " << point.value << "\n";
    // }

    std::vector<double> values;
    for (const auto& point : pointCloud.pts) {
        values.push_back(point.value);
    }

    // Define multiple query points (number of dimensions is one less than the total number of columns in pointCloud)
    std::vector<std::vector<double>> queryPoints = {
        {0, 0, 0, 0},  //    -0.0185284  
        {0.001, 0, 0, 0},  //    -0.0185284  
        {30.001, 0, 0, 0},  //    -0.0185284  
        {2, 15, 0, 0},      //  -0.0145791
        {2, 13.43, 13, 12},        //-0.0234039
        {5.67, 8.91, 14.32, 3.45},// -0.0155981
        {10.23, 7.89, 2.34, 19.56}, //-0.00652823
        {4.56, 12.34, 6.78, 9.01}, // -0.0132203
        {11.11, 13.13, 14.14, 15.15}, // -0.00353603
        {16.16, 17.17, 18.18, 19.19}, // -0.00494464
        {1.23, 2.34, 3.45, 4.56}, // -0.0194466
        {7.89, 8.90, 9.01, 10.12}, // -0.00528553
        {13.14, 14.15, 15.16, 16.17}, // -0.00132123
        {18.19, 19.20, 0.21, 1.22} // 0.00182494
    };
    // Perform interpolation for each query point
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

        // Output the result
        std::cout << "Interpolated value: " << result << "\n";
        std::cout << "Time taken: " << elapsed.count() << " seconds.\n";
    }

    std::cout << "All interpolations completed.\n";

    return 0;
}
