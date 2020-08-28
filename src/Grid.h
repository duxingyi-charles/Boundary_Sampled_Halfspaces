//
// Created by Charles Du on 8/14/20.
//

#ifndef PSI_GRID_H
#define PSI_GRID_H

#include "Sampled_Implicit.h"
#include <vector>

typedef std::pair<int,int> Edge;

class Grid {
public:

    Grid() = default;
    ~Grid() = default;

    void compute_arrangement(const Sampled_Implicit &);


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
