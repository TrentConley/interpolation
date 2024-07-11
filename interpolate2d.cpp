#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <set>

struct Point2D {
    double x, y;
    Point2D(double x, double y) : x(x), y(y) {}
};

// Helper function to extract unique values from points
void get_unique_values(const std::vector<Point2D>& points, std::vector<double>& unique_x, std::vector<double>& unique_y) {
    std::set<double> x_set;
    std::set<double> y_set;
    
    for (const auto& point : points) {
        x_set.insert(point.x);
        y_set.insert(point.y);
    }

    unique_x.assign(x_set.begin(), x_set.end());
    unique_y.assign(y_set.begin(), y_set.end());
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
double find_value_at_point(const std::vector<Point2D>& points, const std::vector<double>& values, double x_val, double y_val) {
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].x == x_val && points[i].y == y_val) {
            return values[i];
        }
    }
    throw std::runtime_error("Value not found for given point");
}

// Helper function to print the unique coordinates
void print_unique_coordinates(const std::vector<double>& unique_x, const std::vector<double>& unique_y) {
    std::cout << "Unique x coordinates: ";
    for (const auto& x : unique_x) {
        std::cout << x << " ";
    }
    std::cout << std::endl;

    std::cout << "Unique y coordinates: ";
    for (const auto& y : unique_y) {
        std::cout << y << " ";
    }
    std::cout << std::endl;
}

// Helper function to print the interpolation points
void print_interpolation_points(const std::vector<double>& unique_x, const std::vector<double>& unique_y, size_t x1_idx, size_t y1_idx, size_t x2_idx, size_t y2_idx) {
    std::cout << "Points used for interpolation:" << std::endl;
    std::cout << "Lower left: (" << unique_x[x1_idx] << ", " << unique_y[y1_idx] << ")" << std::endl;
    std::cout << "Lower right: (" << unique_x[x2_idx] << ", " << unique_y[y1_idx] << ")" << std::endl;
    std::cout << "Upper left: (" << unique_x[x1_idx] << ", " << unique_y[y2_idx] << ")" << std::endl;
    std::cout << "Upper right: (" << unique_x[x2_idx] << ", " << unique_y[y2_idx] << ")" << std::endl;
}

// Function to perform 2D linear interpolation
double interpolate(const Point2D& query_point, const std::vector<Point2D>& points, const std::vector<double>& values) {
    std::vector<double> unique_x, unique_y;
    get_unique_values(points, unique_x, unique_y);

    // Optional: print unique x and y coordinates
    // print_unique_coordinates(unique_x, unique_y);

    double x_lower = find_lower_bound(unique_x, query_point.x);
    double y_lower = find_lower_bound(unique_y, query_point.y);

    // Print x_lower and y_lower values
    std::cout << "x_lower: " << x_lower << std::endl;
    std::cout << "y_lower: " << y_lower << std::endl;

    if (x_lower == unique_x.back() || y_lower == unique_y.back()) {
        throw std::out_of_range("Query point is out of bounds");
    }

    auto x_lower_it = std::find(unique_x.begin(), unique_x.end(), x_lower);
    auto y_lower_it = std::find(unique_y.begin(), unique_y.end(), y_lower);

    size_t x1_idx = std::distance(unique_x.begin(), x_lower_it);
    size_t y1_idx = std::distance(unique_y.begin(), y_lower_it);
    size_t x2_idx = x1_idx + 1;
    size_t y2_idx = y1_idx + 1;

    // Print the 2D points used for the vertices of the hypercube
    print_interpolation_points(unique_x, unique_y, x1_idx, y1_idx, x2_idx, y2_idx);

    // Set x and y for query point
    double x = query_point.x;
    double y = query_point.y;

    // Get the surrounding points and their corresponding values
    double x1 = unique_x[x1_idx];
    double x2 = unique_x[x2_idx];
    double y1 = unique_y[y1_idx];
    double y2 = unique_y[y2_idx];

    double q11 = find_value_at_point(points, values, x1, y1);
    double q21 = find_value_at_point(points, values, x2, y1);
    double q12 = find_value_at_point(points, values, x1, y2);
    double q22 = find_value_at_point(points, values, x2, y2);

    std::cout << "Lower left val: " << q11 << std::endl;
    std::cout << "Lower right val: " << q21 << std::endl;
    std::cout << "Upper left val: " << q12 << std::endl;
    std::cout << "Upper right val: " << q22 << std::endl;

    // x direction first
    double f_x_y1 = (x2 - x) * q11 / (x2 - x1) + (x - x1) * q21 / (x2 - x1);
    double f_x_y2 = (x2 - x) * q12 / (x2 - x1) + (x - x1) * q22 / (x2 - x1);

    // y direction next
    double interpolated_value = (y2 - y) * f_x_y1 / (y2 - y1) + (y - y1) * f_x_y2 / (y2 - y1);

    return interpolated_value;
}

int main() {
    // Define points on a 2D grid and their corresponding values
    std::vector<Point2D> points = { 
        Point2D(0.0, 0.0), // 0
        Point2D(1.0, 0.0), // 1
        Point2D(2.0, 0.0), // 4
        Point2D(0.0, 1.0), // 1
        Point2D(1.0, 1.0), // 2
        Point2D(2.0, 1.0), // 5
        Point2D(0.0, 2.0), // 4
        Point2D(1.0, 2.0), // 5
        Point2D(2.0, 2.0)  // 8
    };
    std::vector<double> values = { 
        0.0, 1.0, 4.0, 
        1.0, 2.0, 5.0, 
        4.0, 5.0, 8.0
    }; // Example values on a 3x3 grid

    // Define query points
    std::vector<Point2D> query_points = { 
        Point2D(0.5, 0.5), Point2D(1.5, 0.5), 
        Point2D(1.5, 1.5), Point2D(0.7, 1.85)
    };

    // Perform interpolation for each query point
    for (const auto& query_point : query_points) {
        try {
            double value = interpolate(query_point, points, values);
            std::cout << "Interpolated value at (" << query_point.x << ", " << query_point.y << ") is " << value << "\n" << std::endl;
        } catch (const std::invalid_argument& e) {
            std::cerr << e.what() << std::endl;
        } catch (const std::out_of_range& e) {
            std::cerr << e.what() << std::endl;
        } catch (const std::runtime_error& e) {
            std::cerr << e.what() << std::endl;
        }
    }

    return 0;
}
