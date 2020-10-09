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
    while(!reader.eof()){
        reader >> x >> y >> z;
        pts.emplace_back(Point(x,y,z));
    }

    reader.close();
    return true;
}