#include "low_memory_plan.hpp"
#include <cassert>
#include <iostream>
#include <set>
#include <random>

using namespace hal2vg_low_memory;

static std::set<size_t> group_members(const std::vector<uint8_t>& key, size_t member) {
    std::set<size_t> result;
    for (size_t i = 0; i < key.size(); ++i) {
        if (key[i] == key[member]) result.insert(i);
    }
    return result;
}

int main() {
    // One unchanged group and one singleton must continue while other members split.
    Builder example(5);
    example.append_run({0, 0, 1, 1, 2}, 2);
    example.append_run({0, 0, 1, 2, 3}, 5);
    example.append_run({0, 0, 1, 1, 2}, 8);
    Plan p = example.finish();
    assert(p.node_for(0, 0) == p.node_for(0, 2));
    assert(p.node_for(4, 0) == p.node_for(4, 2));
    assert(p.node_for(2, 0) != p.node_for(2, 2));
    assert(p.nodes[p.node_for(4, 0)].begin == 0 && p.nodes[p.node_for(4, 0)].end == 8);

    // Independent oracle: at each position, shared nodes must mean exactly the
    // input equivalence classes; node continuation means equal adjacent member sets.
    std::mt19937 random(20261025);
    for (size_t members : {size_t(1), size_t(3), size_t(5), size_t(9), size_t(17),
                           size_t(33), size_t(65), size_t(129), size_t(256)}) {
        std::vector<std::vector<uint8_t>> keys;
        Builder builder(members);
        for (size_t run = 0; run < 41; ++run) {
            std::vector<uint8_t> key(members);
            int canonical[256]; std::fill(canonical, canonical + 256, -1);
            unsigned groups = 0;
            for (size_t i = 0; i < members; ++i) {
                unsigned base = run == 0 ? i : random() % 9;
                if (canonical[base] < 0) canonical[base] = groups++;
                key[i] = canonical[base];
            }
            keys.push_back(key);
            builder.append_run(key, run + 1);
        }
        Plan plan = builder.finish();
        for (size_t run = 0; run < keys.size(); ++run) {
            for (size_t i = 0; i < members; ++i) {
                const auto node = plan.node_for(i, run);
                assert(plan.nodes[node].begin <= run && plan.nodes[node].end > run);
                for (size_t j = 0; j < members; ++j) {
                    assert((node == plan.node_for(j, run)) == (keys[run][i] == keys[run][j]));
                }
                if (run != 0) {
                    assert((node == plan.node_for(i, run - 1)) ==
                        (group_members(keys[run], i) == group_members(keys[run - 1], i)));
                }
            }
        }
    }
    std::cout << "low-memory plans: grouping, singleton continuation and bit packing PASS\n";
}
