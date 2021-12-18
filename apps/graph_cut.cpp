#include <BSH.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <limits>
#include <string>

constexpr size_t INVALID = std::numeric_limits<size_t>::max();

struct ArrangementGraph
{
    struct Patch
    {
        size_t adj_cell_0 = INVALID;
        size_t adj_cell_1 = INVALID;

        double weight = -1;
        size_t implicit_id = INVALID;
        size_t num_samples = 0;
        std::vector<size_t> adj_patches;
    };

    size_t num_cells;
    std::vector<Patch> patches;
};

void print_graph(const ArrangementGraph& g)
{
    std::cout << "num_cells: " << g.num_cells << std::endl;
    const size_t num_patches = g.patches.size();
    for (size_t i = 0; i < num_patches; i++) {
        const auto& p = g.patches[i];
        std::cout << "(" << p.adj_cell_0 << ", " << p.adj_cell_1 << "): " << p.weight
                  << " func_id: " << p.implicit_id << " #samples: " << p.num_samples
                  << " adj_patches: ";
        for (auto adj_p : p.adj_patches) {
            std::cout << " " << adj_p;
        }
        std::cout << std::endl;
    }
}

ArrangementGraph load_graph(const std::string& filename)
{
    std::ifstream fin(filename.c_str());
    nlohmann::json info;
    fin >> info;

    ArrangementGraph g;
    g.num_cells = info["num_cells"];

    const auto& info_patches = info["patches"];
    size_t num_patches = info_patches.size();
    g.patches.reserve(num_patches);
    for (size_t i = 0; i < num_patches; i++) {
        const auto& patch = info_patches[i];
        ArrangementGraph::Patch p;
        p.adj_cell_0 = patch["adj_cell_0"].get<size_t>();
        p.adj_cell_1 = patch["adj_cell_1"].get<size_t>();
        p.weight = patch["weight"].get<size_t>();
        p.implicit_id = patch["implicit_id"].get<size_t>();
        p.num_samples = patch["num_samples"].get<size_t>();

        size_t num_adj_patches = patch["adj_patches"].size();
        p.adj_patches.reserve(num_adj_patches);
        for (auto adj_p : patch["adj_patches"]) {
            p.adj_patches.push_back(adj_p);
        }

        g.patches.push_back(std::move(p));
    }

    return g;
}

std::vector<bool> graph_cut(const ArrangementGraph& g)
{
    const size_t num_cells = g.num_cells;
    const size_t num_patches = g.patches.size();

    std::vector<double> P_dist;
    std::vector<std::vector<int>> P_samples;
    std::vector<std::vector<int>> P_block;
    std::vector<std::vector<int>> P_sign;
    std::vector<std::vector<int>> B_patch;
    std::vector<std::vector<int>> P_Adj_same;
    std::vector<std::vector<int>> P_Adj_diff;
    int topK = 1;
    bool consider_adj_diff = true;
    int max_search_count = std::numeric_limits<int>::max();

    P_dist.reserve(num_patches);
    P_samples.reserve(num_patches);
    P_block.reserve(num_patches);
    P_sign.reserve(num_patches);
    B_patch.resize(num_cells);

    size_t sample_count = 0;
    for (size_t i = 0; i < num_patches; i++) {
        const auto& p = g.patches[i];
        P_dist.push_back(p.weight);
        P_block.push_back({static_cast<int>(p.adj_cell_0), static_cast<int>(p.adj_cell_1)});
        P_sign.push_back({1, -1});

        P_samples.emplace_back();
        auto& samples = P_samples.back();
        for (size_t j=0; j<p.num_samples; j++) {
            samples.push_back(sample_count);
            sample_count++;
        }

        B_patch[p.adj_cell_0].push_back(i);
        B_patch[p.adj_cell_1].push_back(i);
    }

    P_Adj_same.resize(num_patches);
    P_Adj_diff.resize(num_patches);
    for (size_t i = 0; i < num_patches; i++) {
        const auto& p = g.patches[i];
        const auto& adj_patches = p.adj_patches;
        for (auto j : adj_patches) {
            if (g.patches[i].implicit_id == g.patches[j].implicit_id) {
                P_Adj_same[i].push_back(j);
            } else {
                P_Adj_diff[i].push_back(j);
            }
        }
    }

    std::vector<bool> B_label;
    std::vector<bool> P_label;
    std::vector<int>  P_prohibited;
    BSH::connected_graph_cut(P_dist,
        P_samples,
        P_block,
        P_sign,
        B_patch,
        P_Adj_same,
        P_Adj_diff,
        topK,
        consider_adj_diff,
        max_search_count,
        B_label,
        P_label,
        P_prohibited);

    return B_label;
}

void save_labels(const std::string& filename, const std::vector<bool>& labels) {
    std::ofstream fout(filename.c_str());
    nlohmann::json info;
    info["labels"] = labels;
    fout << info;
}

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input_graph.json output_label.json" << std::endl;
        return 1;
    }

    auto g = load_graph(argv[1]);
    print_graph(g);

    auto labels = graph_cut(g);
    save_labels(argv[2], labels);

    return 0;
}
