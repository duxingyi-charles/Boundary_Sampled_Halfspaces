//
// Created by Charles Du on 8/14/20.
//

#ifndef PSI_GRID_H
#define PSI_GRID_H

#include "Sampled_Implicit.h"
#include <vector>
#include <string>

typedef std::pair<int,int> Edge;

class Grid {
public:

    Grid() = default;
    // construct a single cuboid grid with lower corner p_min and upper corner p_max
    Grid(const Point& p_min, const Point& p_max);

    ~Grid() = default;

    void compute_arrangement(const Sampled_Implicit &);

    bool export_grid(const std::string &filename);

private:
    // vertices
    std::vector<Point> V;
    // edges
    std::vector<Edge>  E;
    // faces
    std::vector<std::vector<int>> F;
    // cells
    std::vector<std::vector<int>> C;

    //

};

#endif //PSI_GRID_H
