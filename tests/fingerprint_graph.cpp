// Full PackedGraph comparison by path spelling and canonical homology maps.
// No per-base graph objects or full sequence exports are materialized.
#include "bdsg/packed_graph.hpp"
#include <handlegraph/util.hpp>
#include <openssl/evp.h>
#include <jansson.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

using namespace std;
using namespace handlegraph;

class Digest {
public:
    Digest() : context(EVP_MD_CTX_create()) {
        if (!context || EVP_DigestInit_ex(context, EVP_sha256(), nullptr) != 1)
            throw runtime_error("cannot initialize SHA256");
        pending.reserve(65536);
    }
    ~Digest() { EVP_MD_CTX_destroy(context); }
    Digest(const Digest&) = delete;
    Digest& operator=(const Digest&) = delete;
    void byte(unsigned char value) {
        pending.push_back(static_cast<char>(value));
        if (pending.size() == 65536) flush();
    }
    void integer(uint64_t value) {
        for (unsigned i = 0; i < 8; ++i) { byte(value & 255); value >>= 8; }
    }
    string finish() {
        flush();
        unsigned char bytes[EVP_MAX_MD_SIZE]; unsigned size = 0;
        if (EVP_DigestFinal_ex(context, bytes, &size) != 1 || size != 32)
            throw runtime_error("SHA256 finalization failed");
        const char* hex = "0123456789abcdef";
        string result;
        for (unsigned i = 0; i < size; ++i) {
            result.push_back(hex[bytes[i] >> 4]); result.push_back(hex[bytes[i] & 15]);
        }
        return result;
    }
private:
    EVP_MD_CTX* context;
    string pending;
    void flush() {
        if (!pending.empty() && EVP_DigestUpdate(context, pending.data(), pending.size()) != 1)
            throw runtime_error("SHA256 update failed");
        pending.clear();
    }
};

static char upper(char base) {
    return base >= 'a' && base <= 'z' ? base - ('a' - 'A') : base;
}
static char complement(char base) {
    char result;
    switch (upper(base)) {
        case 'A': result = 'T'; break; case 'C': result = 'G'; break;
        case 'G': result = 'C'; break; case 'T': result = 'A'; break;
        case 'N': result = 'N'; break; case 'R': result = 'Y'; break;
        case 'Y': result = 'R'; break; case 'K': result = 'M'; break;
        case 'M': result = 'K'; break; case 'S': result = 'S'; break;
        case 'W': result = 'W'; break; case 'B': result = 'V'; break;
        case 'D': result = 'H'; break; case 'H': result = 'D'; break;
        case 'V': result = 'B'; break;
        default: throw runtime_error("unsupported DNA character in graph");
    }
    return base >= 'a' && base <= 'z' ? result + ('a' - 'A') : result;
}

struct NodeInfo {
    uint64_t length = 0;
    uint64_t anchor_path = numeric_limits<uint64_t>::max();
    uint64_t anchor_start = 0;
    bool anchor_reverse = false;
    bool cached = false;
    string sequence;
};

class Nodes {
public:
    explicit Nodes(const bdsg::PackedGraph& graph) {
        const uint64_t count = graph.get_node_count();
        minimum = count ? graph.min_node_id() : 0;
        const uint64_t range = count ? graph.max_node_id() - minimum + 1 : 0;
        dense = range <= count * 2 + 1024;
        if (dense) array.resize(range); else sparse.reserve(count);
        graph.for_each_handle([&](const handle_t& handle) {
            const nid_t id = graph.get_id(handle);
            NodeInfo& info = dense ? array[id - minimum] : sparse[id];
            info.length = graph.get_length(handle);
            if (info.length == 0) throw runtime_error("zero-length node is not supported");
        });
    }
    NodeInfo& at(nid_t id) {
        if (dense) {
            if (id < minimum || static_cast<uint64_t>(id - minimum) >= array.size() ||
                array[id - minimum].length == 0) throw runtime_error("path visits a missing node");
            return array[id - minimum];
        }
        return sparse.at(id);
    }
private:
    nid_t minimum;
    bool dense;
    vector<NodeInfo> array;
    unordered_map<nid_t, NodeInfo> sparse;
};

struct Edge {
    uint64_t a, b;
    bool operator==(const Edge& other) const { return a == other.a && b == other.b; }
};
struct EdgeHash {
    size_t operator()(const Edge& edge) const {
        uint64_t value = edge.a + 0x9e3779b97f4a7c15ULL;
        value ^= edge.b + (value << 6) + (value >> 2);
        value ^= value >> 30; value *= 0xbf58476d1ce4e5b9ULL;
        value ^= value >> 27; value *= 0x94d049bb133111ebULL;
        return value ^ (value >> 31);
    }
};
static Edge edge_key(const bdsg::PackedGraph& graph, handle_t a, handle_t b) {
    pair<uint64_t, uint64_t> forward{as_integer(a), as_integer(b)};
    pair<uint64_t, uint64_t> reverse{as_integer(graph.flip(b)), as_integer(graph.flip(a))};
    const auto key = min(forward, reverse);
    return {key.first, key.second};
}

static void text(json_t* object, const char* key, const string& value) {
    if (json_object_set_new(object, key, json_string(value.c_str())) != 0)
        throw runtime_error("cannot encode manifest field");
}
static void number(json_t* object, const char* key, uint64_t value) {
    // Decimal strings also preserve unsigned sentinel values above INT64_MAX.
    text(object, key, to_string(value));
}
static void output(json_t* object) {
    if (json_dumpf(object, stdout, JSON_COMPACT | JSON_SORT_KEYS) != 0 || fputc('\n', stdout) == EOF)
        throw runtime_error("cannot write manifest");
    json_decref(object);
}

struct Interval {
    uint64_t query = 0, length = 0, anchor_path = 0, anchor = 0;
    bool reverse = false;
};
static bool continues(const Interval& old, const Interval& next) {
    if (old.length == 0 || old.anchor_path != next.anchor_path || old.reverse != next.reverse ||
        old.query + old.length != next.query) return false;
    return old.reverse ? old.anchor >= old.length && old.anchor - old.length == next.anchor
                       : old.anchor + old.length == next.anchor;
}
static void record(Digest& digest, const Interval& interval) {
    digest.integer(interval.query); digest.integer(interval.length);
    digest.integer(interval.anchor_path); digest.integer(interval.anchor);
    digest.byte(interval.reverse ? 1 : 0);
}

int main(int argc, char** argv) {
    try {
        if (argc != 2) throw runtime_error("usage: hal2vg-fingerprint input.vg > manifest.jsonl");
        ifstream input(argv[1], ios::binary);
        if (!input) throw runtime_error("cannot open input graph");
        bdsg::PackedGraph graph;
        graph.deserialize(input);
        if (input.bad()) throw runtime_error("error reading input graph");
        vector<pair<string, path_handle_t>> paths;
        graph.for_each_path_handle([&](const path_handle_t& path) {
            paths.emplace_back(graph.get_path_name(path), path);
        });
        sort(paths.begin(), paths.end(), [](const auto& a, const auto& b) { return a.first < b.first; });
        for (size_t i = 1; i < paths.size(); ++i)
            if (paths[i - 1].first == paths[i].first) throw runtime_error("duplicate path name");
        Nodes nodes(graph);
        unordered_map<Edge, bool, EdgeHash> covered;
        graph.for_each_edge([&](const edge_t& edge) { covered.emplace(edge_key(graph, edge.first, edge.second), false); });
        const uint64_t edges = covered.size();
        uint64_t total = 0;
        cerr << "[verify] loaded nodes=" << graph.get_node_count() << " edges=" << edges
             << " paths=" << paths.size() << "; anchor/coverage pass" << endl;
        for (size_t p = 0; p < paths.size(); ++p) {
            uint64_t offset = 0; bool have_previous = false;
            handle_t first, previous;
            auto visit_edge = [&](handle_t a, handle_t b) {
                auto found = covered.find(edge_key(graph, a, b));
                if (found == covered.end()) throw runtime_error("path transition is missing its graph edge");
                found->second = true;
            };
            graph.for_each_step_in_path(paths[p].second, [&](const step_handle_t& step) {
                const handle_t handle = graph.get_handle_of_step(step);
                NodeInfo& info = nodes.at(graph.get_id(handle));
                if (info.anchor_path == numeric_limits<uint64_t>::max()) {
                    info.anchor_path = p; info.anchor_start = offset;
                    info.anchor_reverse = graph.get_is_reverse(handle);
                }
                if (have_previous) visit_edge(previous, handle); else first = handle;
                previous = handle; have_previous = true; offset += info.length;
            });
            if (have_previous && graph.get_is_circular(paths[p].second)) visit_edge(previous, first);
            total += offset;
            if (p % 100 == 0) cerr << "[verify] anchors " << p << '/' << paths.size() << endl;
        }
        for (const auto& entry : covered)
            if (!entry.second) throw runtime_error("graph has an edge not used by any path");
        graph.for_each_handle([&](const handle_t& handle) {
            if (nodes.at(graph.get_id(handle)).anchor_path == numeric_limits<uint64_t>::max())
                throw runtime_error("graph has a node not used by any path");
        });
        covered.clear(); covered.rehash(0);
        cerr << "[verify] spelling/homology pass" << endl;
        for (size_t p = 0; p < paths.size(); ++p) {
            const auto path = paths[p].second;
            Digest raw, uppercase, homology;
            Interval pending;
            uint64_t offset = 0, intervals = 0;
            graph.for_each_step_in_path(path, [&](const step_handle_t& step) {
                const handle_t handle = graph.get_handle_of_step(step);
                const nid_t id = graph.get_id(handle);
                NodeInfo& info = nodes.at(id);
                const bool reverse = graph.get_is_reverse(handle);
                string temporary;
                if (!info.cached && info.length <= 64) {
                    info.sequence = graph.get_sequence(graph.get_handle(id)); info.cached = true;
                }
                if (!info.cached) temporary = graph.get_sequence(graph.get_handle(id));
                const string& sequence = info.cached ? info.sequence : temporary;
                for (uint64_t j = 0; j < info.length; ++j) {
                    const char base = reverse ? complement(sequence[info.length - 1 - j]) : sequence[j];
                    raw.byte(base); uppercase.byte(upper(base));
                }
                const bool opposite = reverse != info.anchor_reverse;
                Interval current{offset, info.length, info.anchor_path,
                    info.anchor_start + (opposite ? info.length - 1 : 0), opposite};
                if (continues(pending, current)) pending.length += current.length;
                else {
                    if (pending.length) { record(homology, pending); ++intervals; }
                    pending = current;
                }
                offset += info.length;
            });
            if (pending.length) { record(homology, pending); ++intervals; }
            json_t* row = json_object();
            text(row, "type", "path"); text(row, "name", paths[p].first);
            number(row, "length", offset); number(row, "homology_intervals", intervals);
            text(row, "sequence_sha256", raw.finish()); text(row, "uppercase_sha256", uppercase.finish());
            text(row, "homology_sha256", homology.finish());
            json_t* metadata = json_object();
            number(metadata, "sense", static_cast<uint64_t>(graph.get_sense(path)));
            text(metadata, "sample", graph.get_sample_name(path)); text(metadata, "locus", graph.get_locus_name(path));
            number(metadata, "haplotype", graph.get_haplotype(path)); number(metadata, "phase_block", graph.get_phase_block(path));
            const auto range = graph.get_subrange(path);
            number(metadata, "subrange_start", range.first); number(metadata, "subrange_end", range.second);
            json_object_set_new(metadata, "circular", json_boolean(graph.get_is_circular(path)));
            json_object_set_new(row, "metadata", metadata); output(row);
            if (p % 100 == 0) cerr << "[verify] fingerprints " << p << '/' << paths.size() << endl;
        }
        json_t* summary = json_object();
        text(summary, "type", "summary"); number(summary, "format_version", 1);
        number(summary, "node_count", graph.get_node_count()); number(summary, "edge_count", edges);
        number(summary, "path_count", paths.size()); number(summary, "total_bases", total);
        number(summary, "unused_nodes", 0); number(summary, "unused_edges", 0);
        json_object_set_new(summary, "complete", json_true()); output(summary);
        if (fflush(stdout) != 0) throw runtime_error("cannot flush manifest");
        return 0;
    } catch (const exception& error) {
        cerr << "hal2vg-fingerprint: " << error.what() << endl;
        return 1;
    }
}
