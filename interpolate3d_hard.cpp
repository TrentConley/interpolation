#include <iostream>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <set>

struct Point3D {
    double x, y, z;
    Point3D(double x, double y, double z) : x(x), y(y), z(z) {}
};

// Helper function to extract unique values from points
void get_unique_values(const std::vector<Point3D>& points, std::vector<double>& unique_x, std::vector<double>& unique_y, std::vector<double>& unique_z) {
    std::set<double> x_set;
    std::set<double> y_set;
    std::set<double> z_set;

    
    for (const auto& point : points) {
        x_set.insert(point.x);
        y_set.insert(point.y);
        z_set.insert(point.z);
    }

    unique_x.assign(x_set.begin(), x_set.end());
    unique_y.assign(y_set.begin(), y_set.end());
    unique_z.assign(z_set.begin(), z_set.end());
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
double find_value_at_point(const std::vector<Point3D>& points, const std::vector<double>& values, double x_val, double y_val, double z_val) {
    for (size_t i = 0; i < points.size(); ++i) {
        if (points[i].x == x_val && points[i].y == y_val && points[i].z == z_val) {
            return values[i];
        }
    }
    throw std::runtime_error("Value not found for given point");
}

// Helper function to print the unique coordinates
void print_unique_coordinates(const std::vector<double>& unique_x, const std::vector<double>& unique_y, const std::vector<double>& unique_z) {
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

    std::cout << "Unique z coordinates: ";
    for (const auto& z : unique_z) {
        std::cout << z << " ";
    }
    std::cout << std::endl;
}

// Helper function to print the interpolation points
void print_interpolation_points(const std::vector<double>& unique_x, const std::vector<double>& unique_y, const std::vector<double>& unique_z, size_t x1_low_idx, size_t y1_low_idx, size_t z1_low_idx) {
    
    std::cout << "Points used for interpolation:" << std::endl;
    std::cout << "111: (" << unique_x[x1_low_idx] << ", " << unique_y[y1_low_idx] << ", " << unique_y[z1_low_idx] << ")" << std::endl;
    std::cout << "112: (" << unique_x[x1_low_idx] << ", " << unique_y[y1_low_idx] << ", " << unique_y[z1_low_idx + 1] << ")" << std::endl;
    std::cout << "121: (" << unique_x[x1_low_idx] << ", " << unique_y[y1_low_idx + 1] << ", " << unique_y[z1_low_idx] << ")" << std::endl;
    std::cout << "122: (" << unique_x[x1_low_idx] << ", " << unique_y[y1_low_idx + 1] << ", " << unique_y[z1_low_idx + 1] << ")" << std::endl;
    std::cout << "211: (" << unique_x[x1_low_idx+1] << ", " << unique_y[y1_low_idx] << ", " << unique_y[z1_low_idx] << ")" << std::endl;
    std::cout << "212: (" << unique_x[x1_low_idx+1] << ", " << unique_y[y1_low_idx] << ", " << unique_y[z1_low_idx+1] << ")" << std::endl;
    std::cout << "221: (" << unique_x[x1_low_idx+1] << ", " << unique_y[y1_low_idx + 1] << ", " << unique_y[z1_low_idx] << ")" << std::endl;
    std::cout << "222: (" << unique_x[x1_low_idx+1] << ", " << unique_y[y1_low_idx + 1] << ", " << unique_y[z1_low_idx+1] << ")" << std::endl;
}

// Function to perform 3D linear interpolation
double interpolate(const Point3D& query_point, const std::vector<Point3D>& points, const std::vector<double>& values) {
    std::vector<double> unique_x, unique_y, unique_z;
    get_unique_values(points, unique_x, unique_y, unique_z);

    // Optional: print unique x, y, and z coordinates
    print_unique_coordinates(unique_x, unique_y, unique_z);

    double x_lower = find_lower_bound(unique_x, query_point.x);
    double y_lower = find_lower_bound(unique_y, query_point.y);
    double z_lower = find_lower_bound(unique_z, query_point.z);

    // Print x_lower, y_lower, and z_lower values
    std::cout << "x_lower: " << x_lower << std::endl;
    std::cout << "y_lower: " << y_lower << std::endl;
    std::cout << "z_lower: " << z_lower << std::endl;

    if (x_lower == unique_x.back() || y_lower == unique_y.back() || z_lower == unique_z.back()) {
        throw std::out_of_range("Query point is out of bounds");
    }

    auto x_lower_it = std::find(unique_x.begin(), unique_x.end(), x_lower);
    auto y_lower_it = std::find(unique_y.begin(), unique_y.end(), y_lower);
    auto z_lower_it = std::find(unique_z.begin(), unique_z.end(), z_lower);

    size_t x1_idx = std::distance(unique_x.begin(), x_lower_it);
    size_t y1_idx = std::distance(unique_y.begin(), y_lower_it);
    size_t z1_idx = std::distance(unique_z.begin(), z_lower_it);
    size_t x2_idx = x1_idx + 1;
    size_t y2_idx = y1_idx + 1;
    size_t z2_idx = z1_idx + 1;

    // Print the 3D points used for the vertices of the hypercube
    print_interpolation_points(unique_x, unique_y, unique_z, x1_idx, y1_idx, z1_idx);

    // Set x, y, and z for query point
    double x = query_point.x;
    double y = query_point.y;
    double z = query_point.z;

    // Get the surrounding points and their corresponding values
    double x1 = unique_x[x1_idx];
    double x2 = unique_x[x2_idx];
    double y1 = unique_y[y1_idx];
    double y2 = unique_y[y2_idx];
    double z1 = unique_y[z1_idx];
    double z2 = unique_y[z2_idx];

    double q111 = find_value_at_point(points, values, x1, y1, z1);
    double q112 = find_value_at_point(points, values, x1, y1, z2);
    double q121 = find_value_at_point(points, values, x1, y2, z1);
    double q122 = find_value_at_point(points, values, x1, y2, z2);
    double q211 = find_value_at_point(points, values, x2, y1, z1);
    double q212 = find_value_at_point(points, values, x2, y1, z2);
    double q221 = find_value_at_point(points, values, x2, y2, z1);
    double q222 = find_value_at_point(points, values, x2, y2, z2);

    std::cout << "111 val: " << q111 << std::endl;
    std::cout << "112 val: " << q112 << std::endl;
    std::cout << "121 val: " << q121 << std::endl;
    std::cout << "122 val: " << q122 << std::endl;
    std::cout << "211 val: " << q211 << std::endl;
    std::cout << "212 val: " << q212 << std::endl;
    std::cout << "221 val: " << q221 << std::endl;
    std::cout << "222 val: " << q222 << std::endl;


    // Interpolation in x-direction
    double f_y1_z1 = (x2 - x) / (x2 - x1) * q111 + (x - x1) / (x2 - x1) * q211;
    double f_y1_z2 = (x2 - x) / (x2 - x1) * q112 + (x - x1) / (x2 - x1) * q212;
    double f_y2_z1 = (x2 - x) / (x2 - x1) * q121 + (x - x1) / (x2 - x1) * q221;
    double f_y2_z2 = (x2 - x) / (x2 - x1) * q122 + (x - x1) / (x2 - x1) * q222;

    // Interpolation in y-direction
    double f_z1 = (y2 - y) / (y2 - y1) * f_y1_z1 + (y - y1) / (y2 - y1) * f_y2_z1;
    double f_z2 = (y2 - y) / (y2 - y1) * f_y1_z2 + (y - y1) / (y2 - y1) * f_y2_z2;

    // Interpolation in z-direction
    double interpolated_value = (z2 - z) / (z2 - z1) * f_z1 + (z - z1) / (z2 - z1) * f_z2;

    return interpolated_value;
}

int main() {
    // Define points on a 3D grid and their corresponding values
    std::vector<Point3D> points = { 
        Point3D(0.0, 0.0, 0.0), Point3D(1.0, 0.0, 0.0), Point3D(2.0, 0.0, 0.0),
        Point3D(0.0, 1.0, 0.0), Point3D(1.0, 1.0, 0.0), Point3D(2.0, 1.0, 0.0),
        Point3D(0.0, 2.0, 0.0), Point3D(1.0, 2.0, 0.0), Point3D(2.0, 2.0, 0.0),
        Point3D(0.0, 0.0, 1.0), Point3D(1.0, 0.0, 1.0), Point3D(2.0, 0.0, 1.0),
        Point3D(0.0, 1.0, 1.0), Point3D(1.0, 1.0, 1.0), Point3D(2.0, 1.0, 1.0),
        Point3D(0.0, 2.0, 1.0), Point3D(1.0, 2.0, 1.0), Point3D(2.0, 2.0, 1.0),
        Point3D(0.0, 0.0, 2.0), Point3D(1.0, 0.0, 2.0), Point3D(2.0, 0.0, 2.0),
        Point3D(0.0, 1.0, 2.0), Point3D(1.0, 1.0, 2.0), Point3D(2.0, 1.0, 2.0),
        Point3D(0.0, 2.0, 2.0), Point3D(1.0, 2.0, 2.0), Point3D(2.0, 2.0, 2.0)
    };
    std::vector<double> values = { 
        0.0, 1.0, 4.0, 
        1.0, 2.0, 5.0, 
        4.0, 5.0, 8.0,
        1.0, 2.0, 5.0, 
        2.0, 3.0, 6.0, 
        5.0, 6.0, 9.0,
        4.0, 5.0, 8.0, 
        5.0, 6.0, 9.0, 
        8.0, 9.0, 12.0
    }; // Example values on a 3x3x3 grid

    // Define query points
    std::vector<Point3D> query_points = { 
        Point3D(0.5, 0.5, 0.5), Point3D(1.5, 0.5, 0.5), 
        Point3D(1.5, 1.5, 1.5), Point3D(0.7, 1.85, 1.2),
        Point3D(1.0, 1.0, 1.0)
    };

    // Perform interpolation for each query point
    for (const auto& query_point : query_points) {
        try {
            double value = interpolate(query_point, points, values);
            std::cout << "Interpolated value at (" << query_point.x << ", " << query_point.y << ", " << query_point.z << ") is " << value << "\n" << std::endl;
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