#include <iostream>
#include <vector>
#include <algorithm>

// Function to perform 1D linear interpolation
double interpolate(double query_point, const std::vector<double>& points, const std::vector<double>& values) {
    if (points.size() != values.size() || points.empty()) {
        throw std::invalid_argument("Points and values vectors must have the same non-zero length.");
    }

    // Find the lower bound index for the query point
    auto low_idx = std::lower_bound(points.begin(), points.end(), query_point);

    // If the query point is out of bounds, return the nearest value
    if (low_idx == points.begin()) {
        return values.front();
    }
    if (low_idx == points.end()) {
        return values.back();
    }

    // Indices for interpolation
    size_t idx_upper = std::distance(points.begin(), low_idx);
    std::cout << "idx upper:  " << idx_upper << std::endl;
    size_t idx_lower = idx_upper - 1;
    std::cout << "idx lower:  " << idx_lower << std::endl;

    // Get the surrounding points and their corresponding values
    double x0 = points[idx_lower];
    double x1 = points[idx_upper];
    double y0 = values[idx_lower];
    double y1 = values[idx_upper];

    // Perform linear interpolation
    double interpolated_value = y0 + (y1 - y0) * (query_point - x0) / (x1 - x0);

    return interpolated_value;
}

int main() {
    // Define points and their corresponding values
    std::vector<double> points = {0.0, 1.0, 2.0, 3.0, 4.0};
    std::vector<double> values = {0.0, 1.0, 4.0, 9.0, 16.0}; // Example values: y = x^2

    // Define query points
    std::vector<double> query_points = {-1, 0.2, 1.9, 2.5, 3.5, 4.2};

    // Perform interpolation for each query point
    for (const auto& query_point : query_points) {
        try {
            double value = interpolate(query_point, points, values);
            std::cout << "Interpolated value at " << query_point << " is " << value << std::endl;
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    return 0;
}
