#ifndef HAL2VG_LOW_MEMORY_PLAN_HPP
#define HAL2VG_LOW_MEMORY_PLAN_HPP

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

namespace hal2vg_low_memory {
using std::size_t;

// A node is shared by exactly the same original block members throughout its
// interval. packed is populated lazily by the normal handle-graph conversion.
struct Node {
    uint64_t begin;
    uint64_t end;
    uint64_t packed = 0;
};

struct Plan {
    size_t members = 0;
    size_t runs = 0;
    unsigned bits = 1;
    std::vector<Node> nodes;
    std::vector<uint64_t> run_offsets{0};
    std::vector<uint64_t> run_nodes;
    std::vector<uint8_t> labels; // member-major, bits bits per label

    uint64_t node_for(size_t member, size_t run) const {
        if (member >= members || run >= runs) {
            throw std::out_of_range("invalid low-memory split-plan coordinate");
        }
        const size_t bit = (member * runs + run) * bits;
        const size_t byte = bit / 8;
        const unsigned shift = bit % 8;
        unsigned group = labels[byte] >> shift;
        if (shift + bits > 8) {
            group |= static_cast<unsigned>(labels[byte + 1]) << (8 - shift);
        }
        group &= (1u << bits) - 1;
        return run_nodes[run_offsets[run] + group];
    }

    uint64_t storage_bytes() const {
        return nodes.capacity() * sizeof(Node) + run_offsets.capacity() * sizeof(uint64_t) +
            run_nodes.capacity() * sizeof(uint64_t) + labels.capacity();
    }
};

class Builder {
public:
    explicit Builder(size_t members) {
        if (members == 0) {
            throw std::invalid_argument("a split plan needs at least one member");
        }
        plan.members = members;
    }

    // key numbers the equivalence classes at every position in [previous_end,end).
    // Labels are canonical dense numbers, as in hal2vg's grouping_at().
    void append_run(const std::vector<uint8_t>& key, uint64_t end) {
        if (key.size() != plan.members || end <= previous_end) {
            throw std::invalid_argument("invalid low-memory grouping run");
        }
        const size_t groups = static_cast<size_t>(*std::max_element(key.begin(), key.end())) + 1;
        maximum_groups = std::max(maximum_groups, groups);
        const size_t none = std::numeric_limits<size_t>::max();
        std::vector<size_t> counts(groups, 0), parent(groups, none);
        std::vector<bool> consistent(groups, true);
        for (size_t i = 0; i < key.size(); ++i) {
            const size_t g = key[i];
            ++counts[g];
            if (!previous_key.empty()) {
                if (parent[g] == none) {
                    parent[g] = previous_key[i];
                } else if (parent[g] != previous_key[i]) {
                    consistent[g] = false;
                }
            }
        }
        std::vector<uint64_t> current_nodes(groups);
        for (size_t g = 0; g < groups; ++g) {
            if (counts[g] == 0) {
                throw std::invalid_argument("non-dense low-memory group labels");
            }
            // Equal member sets continue a node even if OTHER groups changed.
            // This also keeps singleton intervals continuous, without pinching.
            if (!previous_key.empty() && consistent[g] &&
                counts[g] == previous_counts[parent[g]]) {
                current_nodes[g] = previous_nodes[parent[g]];
                plan.nodes[current_nodes[g]].end = end;
            } else {
                current_nodes[g] = plan.nodes.size();
                plan.nodes.push_back({previous_end, end, 0});
            }
        }
        row_labels.insert(row_labels.end(), key.begin(), key.end());
        plan.run_nodes.insert(plan.run_nodes.end(), current_nodes.begin(), current_nodes.end());
        plan.run_offsets.push_back(plan.run_nodes.size());
        ++plan.runs;
        previous_key = key;
        previous_counts = std::move(counts);
        previous_nodes = std::move(current_nodes);
        previous_end = end;
    }

    bool has_disagreement() const { return maximum_groups > 1; }

    Plan finish() {
        if (plan.runs == 0) {
            throw std::invalid_argument("empty low-memory split plan");
        }
        plan.bits = 1;
        while ((1u << plan.bits) < maximum_groups) {
            ++plan.bits;
        }
        if (plan.members > std::numeric_limits<size_t>::max() / plan.runs ||
            plan.members * plan.runs > (std::numeric_limits<size_t>::max() - 7) / plan.bits) {
            throw std::length_error("low-memory split plan exceeds address space");
        }
        const size_t cells = plan.members * plan.runs;
        plan.labels.assign((cells * plan.bits + 7) / 8, 0);
        // Transpose once: subsequent path walks read each member's labels in order.
        for (size_t run = 0; run < plan.runs; ++run) {
            for (size_t member = 0; member < plan.members; ++member) {
                const unsigned value = row_labels[run * plan.members + member];
                const size_t bit = (member * plan.runs + run) * plan.bits;
                const size_t byte = bit / 8;
                const unsigned shift = bit % 8;
                plan.labels[byte] |= value << shift;
                if (shift + plan.bits > 8) {
                    plan.labels[byte + 1] |= value >> (8 - shift);
                }
            }
        }
        return std::move(plan);
    }

private:
    Plan plan;
    std::vector<uint8_t> row_labels;
    std::vector<uint8_t> previous_key;
    std::vector<size_t> previous_counts;
    std::vector<uint64_t> previous_nodes;
    uint64_t previous_end = 0;
    size_t maximum_groups = 0;
};

} // namespace hal2vg_low_memory
#endif
