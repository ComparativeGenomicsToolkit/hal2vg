// Small positive and deliberately damaged PackedGraphs for the verifier tests.
#include "bdsg/packed_graph.hpp"
#include <fstream>
#include <string>
#include <vector>
using namespace std;
using namespace handlegraph;

int main(int argc, char** argv) {
    if (argc != 3) return 2;
    const string mode = argv[1];
    bdsg::PackedGraph graph;
    vector<handle_t> forward;
    auto add = [&](const string& sequence, nid_t id) { return graph.create_handle(sequence, id); };
    if (mode == "split" || mode == "circular_split") {
        forward = {add("AC", 1), add("gT", 11), add("NN", 2), add("A", 12), add("t", 3), add("Gc", 13)};
    } else if (mode == "renumber_flip") {
        forward = {graph.flip(add("AcGT", 999999999)), add("NNA", 5), graph.flip(add("gCa", 999999))};
    } else if (mode == "loop" || mode == "loop_split") {
        if (mode == "loop") forward = {add("AT", 1)};
        else forward = {add("A", 30), add("T", 10)};
    } else if (mode == "orientation" || mode == "orientation_bad") {
        forward = {add("NN", 1)};
    } else {
        forward = {add(mode == "mutation" ? "TCgT" : "ACgT", 1), add("NNA", 2), add("tGc", 3)};
    }
    const bool loop = mode == "loop" || mode == "loop_split";
    const bool circular = mode == "circular" || mode == "circular_split";
    auto append = [&](const string& name, const vector<handle_t>& steps, bool is_circular) {
        auto path = graph.create_path_handle(name, is_circular);
        for (auto handle : steps) graph.append_step(path, handle);
        for (size_t j = 1; j < steps.size(); ++j) {
            const auto a = graph.get_id(steps[j - 1]), b = graph.get_id(steps[j]);
            if (!(mode == "missing_edge" && ((a == 2 && b == 3) || (a == 3 && b == 2))))
                graph.create_edge(steps[j - 1], steps[j]);
        }
        if (is_circular && !steps.empty()) graph.create_edge(steps.back(), steps.front());
    };
    vector<handle_t> query;
    if (mode == "sharing_loss") {
        for (size_t j = forward.size(); j; --j)
            query.push_back(graph.flip(add(graph.get_sequence(forward[j - 1]), 100 + j)));
    } else {
        for (auto j = forward.rbegin(); j != forward.rend(); ++j) query.push_back(graph.flip(*j));
    }
    if (mode == "orientation_bad") query = forward;
    if (loop) {
        const auto one = forward, reverse_one = query;
        forward.insert(forward.end(), one.begin(), one.end());
        query.insert(query.end(), reverse_one.begin(), reverse_one.end());
    }
    append("a_ref", forward, circular);
    append("z_hap", query, circular);
    append("empty", {}, false);
    if (mode == "orphan_node") add("A", 200);
    if (mode == "unused_edge") graph.create_edge(forward.front(), forward.back());
    if (mode == "hairpin") {
        append("hairpin", {forward[0], graph.flip(forward[0])}, false);
    }
    ofstream output(argv[2], ios::binary);
    graph.serialize(output);
}
