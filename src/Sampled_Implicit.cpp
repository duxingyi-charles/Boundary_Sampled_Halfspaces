//
// Created by Charles Du on 10/9/20.
//

#include "Sampled_Implicit.h"
#include <iostream>
#include <fstream>

bool Sampled_Implicit::import_xyz(const std::string &filename, std::vector<Point> &pts) {

    std::ifstream reader(filename.data(), std::ofstream::in);

    if (!reader.good()) {
        std::cout << "Can not open the file " << filename << std::endl;
        return false;
    }else {
        std::cout << "Reading: "<<filename<< std::endl;
    }

    // first number: dimension (2 or 3)
    int dim;
    reader >> dim;
    // read point coordinates
    if (dim != 3) {
        std::cout << "Can't handle non-3D points." << std::endl;
        reader.close();
        return false;
    }
    pts.clear();
    double x,y,z;
    while ((reader >> x >> y >> z)) {
        pts.emplace_back(Point(x,y,z));
    }

    reader.close();
    return true;
}

bool Sampled_Implicit::export_xyz(const std::string &filename, const std::vector<Point> &pts) {
//    if (pts.empty()) {
//        std::cout << "Vector of points is empty." << std::endl;
//        return false;
//    }
    std::ofstream fout(filename, std::ofstream::out);
    if (!fout.good()) {
        std::cout << "Can not create output file " << filename << std::endl;
        return false;
    }

    //precision of output
    fout.precision(std::numeric_limits<double>::max_digits10);

    int dim = pts[0].size();
    fout << dim << std::endl;

    for (const auto &p : pts) {
        for (int i = 0; i < dim; ++i) {
            fout << p(i) << " ";
        }
        fout << std::endl;
    }

    fout.close();
    std::cout << "export_xyz finish: " << filename << std::endl;
    return true;
}