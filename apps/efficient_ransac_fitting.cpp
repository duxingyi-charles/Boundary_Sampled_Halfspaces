//
// Created by Charles Du on 1/17/21.
//
#include <CLI/CLI.hpp>
#include <string>
#include <limits>

#include <Eigen/Dense>

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

bool ransac_fit(Pwn_vector &points, double error_bound, const std::string &primitive_type,
                double probability, size_t min_points, double epsilon, double cluster_epsilon,
                double normal_threshold, size_t num_repeat, bool verbose,
        //output
                std::string &detected_shape_str) {
    // Instantiate shape detection engine.
    Efficient_ransac ransac;
    // Provide input data.
    ransac.set_input(points);

    // Register shapes for detection.
    if (primitive_type == "plane") {
        ransac.add_shape_factory<Plane>();
    } else if (primitive_type == "sphere") {
        ransac.add_shape_factory<Sphere>();
    } else if (primitive_type == "cylinder") {
        ransac.add_shape_factory<Cylinder>();
    } else if (primitive_type == "cone") {
        ransac.add_shape_factory<Cone>();
    } else if (primitive_type == "torus") {
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
    if (probability != -1) {
        parameters.probability = probability;
    }
    // Detect shapes with at least min_points points.
    if (min_points != std::numeric_limits<std::size_t>::max()) {
        parameters.min_points = min_points;
    }
    // Set maximum Euclidean distance between a point and a shape.
    if (epsilon != -1) {
        parameters.epsilon = epsilon;
    }
    // Set maximum Euclidean distance between points to be clustered.
    if (cluster_epsilon != -1) {
        parameters.cluster_epsilon = cluster_epsilon;
    }
    // Set maximum normal deviation.
    // normal_threshold < dot(surface_normal, point_normal);
    if (normal_threshold != -1) {
        parameters.normal_threshold = normal_threshold;
    }
    if (verbose) {
        std::cout << "probability: " << parameters.probability << std::endl;
        std::cout << "min_points: " << parameters.min_points << std::endl;
        std::cout << "epsilon: " << parameters.epsilon << std::endl;
        std::cout << "cluster_epsilon: " << parameters.cluster_epsilon << std::endl;
        std::cout << "normal_threshold: " << parameters.normal_threshold << std::endl;
        std::cout << "num_repeat: " << num_repeat << std::endl;
    }


    double min_dist = std::numeric_limits<double>::max();
    bool detect_succeed = false;
    for (size_t iter = 0; iter <num_repeat; ++iter) {
        if (verbose) {
            std::cout << "~~~~~~~~~~~~ iter " << iter << " ~~~~~~~~~~~~" << std::endl;
        }

        ransac.detect(parameters);

        // Print number of detected shapes and unassigned points.
        if (verbose) {
            std::cout << ransac.shapes().end() - ransac.shapes().begin()
                      << " detected shapes, "
                      << ransac.number_of_unassigned_points()
                      << " unassigned points." << std::endl;
        }

        Efficient_ransac::Shape_range shapes = ransac.shapes();
        Efficient_ransac::Shape_range::iterator it = shapes.begin();

        while (it != shapes.end()) {
            // compute max distance between input points and the current shape
            FT max_distance = 0;
            for (const auto &p : points) {
                FT dist = CGAL::sqrt((*it)->squared_distance(p.first));
                if (dist > max_distance) {
                    max_distance = dist;
                }
            }

            // update best fitting result
            if (max_distance < error_bound && max_distance < min_dist) {
                bool is_shape_valid = true;
                // special treatment for torus: we don't allow minor_radius >= major_radius
                std::string info = (*it)->info();
                if (info.find("torus") != std::string::npos) {
                    std::size_t pos1 = info.find("major radius");
                    std::size_t pos2 = info.find("minor radius");
                    std::size_t pos3 = info.find("#Pts");
//                    std::cout << info << std::endl;
                    std::string str = info.substr(pos1+15, pos2 - (pos1+15));
                    double major_radius = std::stod(str);
//                    std::cout << str << std::endl;
//                    std::cout << major_radius << std::endl;
                    str = info.substr(pos2+15, pos3 - (pos2+15));
                    double minor_radius = std::stod(str);
//                    std::cout << str << std::endl;
//                    std::cout << minor_radius << std::endl;
                    if (minor_radius >= major_radius) {
                        is_shape_valid = false;
                    }
                }
                //
                if (is_shape_valid) {
                    detect_succeed = true;
                    min_dist = max_distance;
                    detected_shape_str = (*it)->info();
                    if (verbose) {
                        std::cout << "----------------" << std::endl;
                        std::cout << (*it)->info();
                        std::cout << " maxDist: " << max_distance << std::endl;
                    }
                }
            }

            // Proceed with the next detected shape.
            it++;
        }
    }

    return detect_succeed;
}

bool export_primitive(const std::string &filename, const std::string &primitive_str) {
    std::ofstream fout(filename, std::ofstream::out);
    if (!fout.good()) {
        std::cout << "Can not create output file " << filename << std::endl;
        return false;
    }

    fout << primitive_str << std::endl;

    fout.close();
    std::cout << "export_primitive finish: " << filename << std::endl;
    return true;
}


int main(int argc, char** argv) {
    // parse command line
    struct {
        std::string point_normal_file;
        std::string type;
        double error_bound;
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
        // verbose printout
        bool verbose = true;
        // ransac output
        std::string output_primitive_file;
    } args;

    CLI::App app{"Mixed Fitting demo"};
    app.add_option("-P,--point-normal-file", args.point_normal_file,
                   "input point/normal file")->required();
    app.add_option("-T,--primitive-type", args.type, "primtive type")->required();
    app.add_option("-E,--error-bound", args.error_bound,
                   "Error bound for fitting")->required();
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
    app.add_option("-r,--num_repeat", args.num_repeat,
                   "repeat efficient ransac multiple times");
    app.add_option("-v,--verbose", args.verbose,
                   "verbose printout");
    app.add_option("output_primitive_file", args.output_primitive_file,
                   "Output primtive file")
            ->required();
    CLI11_PARSE(app, argc, argv);



    // import points with normals
    Pwn_vector points;

    std::ifstream stream(args.point_normal_file);
    if (!stream ||
        !CGAL::read_xyz_points(
                stream,
                std::back_inserter(points),
                CGAL::parameters::point_map(Point_map()).
                        normal_map(Normal_map()))) {
        std::cerr << "Error: cannot read file " << args.point_normal_file << "!" << std::endl;
        return EXIT_FAILURE;
    }
    std::cout << points.size() << " points" << std::endl;

    // efficient ransac fitting
    std::cout << "ransac " << args.type << " fitting..." << std::endl;
    std::string detected_shape_str;
    if (ransac_fit(points,args.error_bound,args.type,args.probability,args.min_points,
                   args.epsilon,args.cluster_epsilon,args.normal_threshold,
                   args.num_repeat,args.verbose,detected_shape_str)){
        std::cout << detected_shape_str << std::endl;
        export_primitive(args.output_primitive_file, detected_shape_str);
        return EXIT_SUCCESS;
    } else {
        std::cout << "ransac fitting failed: no shape detected." << std::endl;
        return EXIT_FAILURE;
    }


}
