//
// Created by Charles Du on 1/16/21.
//

#include <CLI/CLI.hpp>
#include <string>
#include <limits>
#include <cmath>

//#include <fstream>
#include <iostream>
#include <CGAL/property_map.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Shape_detection/Efficient_RANSAC.h>
// Type declarations.
typedef CGAL::Exact_predicates_inexact_constructions_kernel  Kernel;
typedef Kernel::FT                                           FT;
typedef std::pair<Kernel::Point_3, Kernel::Vector_3>         Point_with_normal;
typedef std::vector<Point_with_normal>                       Pwn_vector;
typedef CGAL::First_of_pair_property_map<Point_with_normal>  Point_map;
typedef CGAL::Second_of_pair_property_map<Point_with_normal> Normal_map;
typedef CGAL::Shape_detection::Efficient_RANSAC_traits
        <Kernel, Pwn_vector, Point_map, Normal_map>             Traits;
typedef CGAL::Shape_detection::Efficient_RANSAC<Traits> Efficient_ransac;
typedef CGAL::Shape_detection::Cone<Traits>             Cone;
typedef CGAL::Shape_detection::Cylinder<Traits>         Cylinder;
typedef CGAL::Shape_detection::Plane<Traits>            Plane;
typedef CGAL::Shape_detection::Sphere<Traits>           Sphere;
typedef CGAL::Shape_detection::Torus<Traits>            Torus;


int main(int argc, char** argv) {
    struct {
        std::string point_file;
        std::string primitive_type = "all";
        // Set probability to miss the largest primitive at each iteration.
        double probability = -1;
        // Detect shapes with at least xxx points.
        size_t min_points = std::numeric_limits<std::size_t>::max();
        // Set maximum Euclidean distance between a point and a shape.
        double epsilon = -1;
        // Set maximum Euclidean distance between points to be clustered.
        double cluster_epsilon = -1;
        // Set maximum normal deviation.
        // normal_threshold < dot(surface_normal, point_normal);
        double normal_threshold = -1;
        // repeat efficient ransac to obtain more accurate fitting
        size_t num_repeat = 1;
    } args;

    CLI::App app{"Efficient RANSAC demo"};
    app.add_option("-P,--point-file", args.point_file, "input point/normal file")
            ->required();
    app.add_option("-T,--primitive-type", args.primitive_type, "Primitive type")
            ->required();
    app.add_option("-p,--probability", args.probability,
                   "probability to miss the largest primitive at each iteration");
    app.add_option("-m,--min_points", args.min_points,
                   "Detect shapes with at least xxx points");
    app.add_option("-e,--epsilon", args.epsilon,
                   "maximum Euclidean distance between a point and a shape");
    app.add_option("-c,--cluster_epsilon", args.cluster_epsilon,
                   "maximum Euclidean distance between points to be clustered");
    app.add_option("-n,--normal_threshold", args.normal_threshold,
                   "maximum normal deviation");
    app.add_option("-r,--num_repeat", args.num_repeat, "repeat efficient ransac");
    CLI11_PARSE(app, argc, argv);

    // Points with normals.
    Pwn_vector points;

    std::ifstream stream(args.point_file);
    if (!stream ||
        !CGAL::read_xyz_points(
                stream,
                std::back_inserter(points),
                CGAL::parameters::point_map(Point_map()).
                        normal_map(Normal_map()))) {
        std::cerr << "Error: cannot read file " << args.point_file << "!" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << points.size() << " points" << std::endl;

    // Instantiate shape detection engine.
    Efficient_ransac ransac;
    // Provide input data.
    ransac.set_input(points);

    // Register shapes for detection.
    if (args.primitive_type == "plane") {
        ransac.add_shape_factory<Plane>();
    } else if (args.primitive_type == "sphere") {
        ransac.add_shape_factory<Sphere>();
    } else if (args.primitive_type == "cylinder") {
        ransac.add_shape_factory<Cylinder>();
    } else if (args.primitive_type == "cone") {
        ransac.add_shape_factory<Cone>();
    } else if (args.primitive_type == "torus") {
        ransac.add_shape_factory<Torus>();
    } else {
        ransac.add_shape_factory<Plane>();
        ransac.add_shape_factory<Sphere>();
        ransac.add_shape_factory<Cylinder>();
        ransac.add_shape_factory<Cone>();
        ransac.add_shape_factory<Torus>();
    }


    // Set parameters for shape detection.
    Efficient_ransac::Parameters parameters;
    // Set probability to miss the largest primitive at each iteration.
    if (args.probability != -1) {
        parameters.probability = args.probability;
    }
    std::cout << "probability: " << parameters.probability << std::endl;
    // Detect shapes with at least min_points points.
    if (args.min_points != std::numeric_limits<std::size_t>::max()) {
        parameters.min_points = args.min_points;
    }
    std::cout << "min_points: " << parameters.min_points << std::endl;
    // Set maximum Euclidean distance between a point and a shape.
    if (args.epsilon != -1) {
        parameters.epsilon = args.epsilon;
    }
    std::cout << "epsilon: " << parameters.epsilon << std::endl;
    // Set maximum Euclidean distance between points to be clustered.
    if (args.cluster_epsilon != -1) {
        parameters.cluster_epsilon = args.cluster_epsilon;
    }
    std::cout << "cluster_epsilon: " << parameters.cluster_epsilon << std::endl;
    // Set maximum normal deviation.
    // normal_threshold < dot(surface_normal, point_normal);
    if (args.normal_threshold != -1) {
        parameters.normal_threshold = args.normal_threshold;
    }
    std::cout << "normal_threshold: " << parameters.normal_threshold << std::endl;

    std::cout << "num_repeat: " << args.num_repeat << std::endl;

    double min_dist = std::numeric_limits<double>::max();
    for (int iter = 0; iter <args.num_repeat; ++iter) {
        std::cout << "~~~~~~~~~~~~ iter " << iter << " ~~~~~~~~~~~~" << std::endl;

        ransac.detect(parameters);

        // Print number of detected shapes and unassigned points.
        std::cout << ransac.shapes().end() - ransac.shapes().begin()
                  << " detected shapes, "
                  << ransac.number_of_unassigned_points()
                  << " unassigned points." << std::endl;

        Efficient_ransac::Shape_range shapes = ransac.shapes();
        Efficient_ransac::Shape_range::iterator it = shapes.begin();

        while (it != shapes.end()) {
            // Get specific parameters depending on the detected shape.
            // Print the parameters of the detected shape.
//            std::cout << "----------------" << std::endl;
//            std::cout << (*it)->info();

            // iterate through all input points
            FT max_distance = 0;
            for (const auto &p : points) {
                FT dist = CGAL::sqrt((*it)->squared_distance(p.first));
                if (dist > max_distance) {
                    max_distance = dist;
                }
            }
//            std::cout << " max distance: " << max_distance << std::endl;

            //
            if (max_distance < min_dist) {
                min_dist = max_distance;
                std::cout << "----------------" << std::endl;
                std::cout << (*it)->info();
                std::cout << " maxDist: " << max_distance << std::endl;
            }

            // Proceed with the next detected shape.
            it++;
        }
    }

    return EXIT_SUCCESS;



}
