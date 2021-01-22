// Chop regions (from BED File) out of paths in vg graphs, creating subpath names and cutting out nodes or parts of nodes
// Assumes that:
// - regions don't overlap (error otherwise)
// - regions are only ever part of at most one path each (assert false otherwise)

//#define debug

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>
#include <unistd.h>
#include <getopt.h>

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"

using namespace std;
using namespace handlegraph;
using namespace bdsg;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> <bed>" << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl;
}    

static unordered_map<string, vector<pair<int64_t, int64_t>>> load_bed(istream& bed_stream);
static unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream);
static vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems);
static void chop_path_intervals(MutablePathMutableHandleGraph* graph,
                                const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals);
static vector<handle_t> chop_path(MutablePathMutableHandleGraph* graph,
                                  path_handle_t path_handle,
                                  const vector<pair<int64_t, int64_t>>& intervals);
// Create a subpath name (todo: make same function in vg consistent (it only includes start))
static inline string make_subpath_name(const string& path_name, size_t offset, size_t length) {
    return path_name + "[" + std::to_string(offset) + ":" + std::to_string(offset + length) + "]";
}

int main(int argc, char** argv) {
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "h",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'h':
        case '?':
            /* getopt_long already printed an error message. */
            help(argv);
            exit(1);
            break;
        default:
            abort ();
        }
    }

    if (argc <= 1) {
        help(argv);
        return 1;
    }

    // Parse the positional argument
    if (optind >= argc) {
        cerr << "[chop-vg-paths] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    if (optind != argc - 2) {
        cerr << "[chop-vg-paths] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    string graph_path = argv[optind++];
    string bed_path = argv[optind++];

    ifstream bed_stream(bed_path);
    if (!bed_stream) {
        cerr << "[chop-vg-paths] error: Unable to open input BED file " << bed_path << endl;
        return 1;
    }
    unordered_map<string, vector<pair<int64_t, int64_t>>> bed_intervals = load_bed(bed_stream);
    bed_stream.close();
    
    ifstream graph_stream(graph_path);
    if (!graph_stream) {
        cerr << "[chop-vg-paths] error: Unable to open input graph " << graph_path << endl;
        return 1;
    }    
    unique_ptr<MutablePathMutableHandleGraph> graph = load_graph(graph_stream);
    graph_stream.close();

    chop_path_intervals(graph.get(), bed_intervals);

    dynamic_cast<SerializableHandleGraph*>(graph.get())->serialize(cout);

    return 0;
}

unordered_map<string, vector<pair<int64_t, int64_t>>> load_bed(istream& bed_stream) {
    // load bed
    unordered_map<string, vector<pair<int64_t, int64_t>>> intervals;
    string buffer;
    while (getline(bed_stream, buffer)) {
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);
        if (toks.size() >= 3) {
            string& name = toks[0];
            int64_t start = stol(toks[1]);
            int64_t end = stol(toks[2]);            
            intervals[name].push_back(make_pair(start, end));
        }
    }
    // verify bed
    for (auto& seq_intervals : intervals) {
        sort(seq_intervals.second.begin(), seq_intervals.second.end(),
             [](const pair<int64_t, int64_t>& b1, const pair<int64_t, int64_t>& b2) {
                 return b1.first < b2.first || (b1.first == b2.first && b1.second < b2.second);
             });
        for (size_t i = 1; i < seq_intervals.second.size(); ++i) {
            if (seq_intervals.second[i].first < seq_intervals.second[i-1].second) {
                cerr << "Overlapping bed intervals found:\n"
                     << " " << seq_intervals.first << "\t"
                     << seq_intervals.second[i-1].first << "\t"
                     << seq_intervals.second[i-1].second << endl
                     << " " << seq_intervals.first << "\t"
                     << seq_intervals.second[i].first << "\t"
                     << seq_intervals.second[i].second << endl
                     << "These are not supported.  Please clean up (ex with bedools merge) first" << endl;
                exit(1);
            }
        }
    }
    return intervals;
}

unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream) {

    char magic_bytes[4];
    graph_stream.read(magic_bytes, 4);
    uint32_t magic_number = ntohl(*((uint32_t*) magic_bytes));
    graph_stream.clear();
    graph_stream.seekg(0, ios::beg);

    if (magic_number != ODGI().get_magic_number()) {
      cerr << "Only ODGI supported until this bug is fixed: https://github.com/vgteam/libbdsg/issues/94" << endl;
      exit(1);
    }
    
    MutablePathMutableHandleGraph* graph;
    if (magic_number == PackedGraph().get_magic_number()) {
        graph = new PackedGraph();
    } else if (magic_number == HashGraph().get_magic_number()) {
        graph = new HashGraph();
    } else if (magic_number == ODGI().get_magic_number()) {
        graph = new ODGI();
    } else {
        cerr << "Unable to parse input graph with magic number " << magic_number << endl;
        exit(1);
    }
    dynamic_cast<SerializableHandleGraph*>(graph)->deserialize(graph_stream);

    return unique_ptr<MutablePathMutableHandleGraph>(graph);
}

vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems) {
    size_t start = string::npos;
    for (size_t i = 0; i < s.size(); ++i) {
        if (delims.find(s[i]) != string::npos) {
            if (start != string::npos && i > start) {
                elems.push_back(s.substr(start, i - start));
            }
            start = string::npos;
        } else if (start == string::npos) {
            start = i;
        }
    }
    if (start != string::npos && start < s.size()) {
        elems.push_back(s.substr(start, s.size() - start));
    }
    return elems;
}

void chop_path_intervals(MutablePathMutableHandleGraph* graph,
                         const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals) {

    // keep some stats to print
    size_t chopped_paths = 0;
    size_t chopped_nodes = 0;
    size_t chopped_bases = 0;
    
    // careful not to iterate and chop, as we could hit new subpaths made
    vector<path_handle_t> path_handles;    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });
    
    for (auto path_handle : path_handles) {
        string path_name = graph->get_path_name(path_handle);
        auto it = bed_intervals.find(path_name);
        bool was_chopped = false;
        if (it != bed_intervals.end()) {
            vector<handle_t> chopped_handles = chop_path(graph, path_handle, it->second);
            if (!chopped_handles.empty()) {
                graph->destroy_path(path_handle);
                for (handle_t handle : chopped_handles) {
                    assert(graph->steps_of_handle(handle).empty());
                    chopped_bases += graph->get_length(handle);
                    was_chopped = true;
                    ++chopped_nodes;                    
                    dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(handle);
                }
            }
        }
        if (was_chopped) {
            ++chopped_paths;
        }
    }
    cerr << "[clip-vg] : Clipped "
         << chopped_bases << " bases from "
         << chopped_nodes << " nodes in "
         << chopped_paths << " paths" << endl;
}

vector<handle_t> chop_path(MutablePathMutableHandleGraph* graph,
                           path_handle_t path_handle,
                           const vector<pair<int64_t, int64_t>>& intervals) {

    // get the breakpoints
    set<int64_t> breakpoints;
    for (const pair<int64_t, int64_t>& interval : intervals) {
        breakpoints.insert(interval.first);
        breakpoints.insert(interval.second); // we're cutting before offset, so the open coordinate is what we want
    }

    // to be safe, don't cut and iterate at the same time, so load up steps here
    vector<handle_t> steps;
    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
            steps.push_back(graph->get_handle_of_step(step_handle));
        });
    
    // cut the nodes to ensure breakpoints at node boundaries
    int64_t offset = 0;
    for (auto handle : steps) {
        int64_t len = graph->get_length(handle);
        // find breakpoints in node
        vector<size_t> cut_points;
        for (auto i = breakpoints.lower_bound(offset); i != breakpoints.end() && *i - offset < len; ++i) {
            int64_t cut_point = *i - offset;
            // libbdsg is buggy and can't accept cutpoints on ends on reverse strand
            if (cut_point > 0 && cut_point < len) {
                cut_points.push_back(cut_point);
            }
        }
        // chop the node
#ifdef debug
        if (!cut_points.empty()) {
            cerr << "dividing node_id=" << graph->get_id(handle) << ":" << graph->get_is_reverse(handle) << " seq=" << graph->get_sequence(handle)
                 << " for path " << graph->get_path_name(path_handle) << " at cut points:";
            for (auto cp : cut_points) {
                cerr << " " << cp;
            }
            cerr << endl;
        }
#endif
#ifdef debug
        vector<handle_t> pieces = graph->divide_handle(handle, cut_points) ;
        for (auto piece : pieces) {
            cerr << " piece " << graph->get_id(piece) << ":" << graph->get_is_reverse(piece) << " " << graph->get_sequence(piece) << endl;
        }
#else
        graph->divide_handle(handle, cut_points);
#endif
        offset += len;
    }
    
    steps.clear();
    int64_t original_path_length = offset;
    vector<handle_t> chopped_handles;
    offset = 0;
    step_handle_t current_step = graph->path_begin(path_handle);
#ifdef debug
    cerr << "init step to " << graph->get_id(graph->get_handle_of_step(current_step)) << ":" << graph->get_is_reverse(graph->get_handle_of_step(current_step))
         << " seq=" <<graph->get_sequence(graph->get_handle_of_step(current_step)) << endl;
#endif
    vector<path_handle_t> subpaths;
    
    // cut out a subpath and make a new path out of it
    function<void(int64_t)> cut_to = [&](int64_t end_offset) {
#ifdef debug
        cerr << "\ncut_to " << end_offset << " where current offset is " << offset << endl;
#endif
        vector<handle_t> steps;
        int64_t start_offset = offset;
        int64_t path_length = 0;
        while (offset < end_offset && current_step != graph->path_end(path_handle)) {
            handle_t handle = graph->get_handle_of_step(current_step);
            steps.push_back(handle);
#ifdef debug
            cerr << " pushing subpath step " << graph->get_id(handle) << " len=" << graph->get_length(handle) << endl;
#endif
            offset += graph->get_length(handle);
            current_step = graph->get_next_step(current_step);
            path_length += graph->get_length(handle);
        }
#ifdef debug
        cerr << "start offset=" << start_offset << " path length=" << path_length << " end offset=" << end_offset << endl;
#endif
        assert(start_offset + path_length == end_offset);

        if (path_length > 0) {
            path_handle_t subpath_handle = graph->create_path_handle(make_subpath_name(graph->get_path_name(path_handle), start_offset, path_length));
            for (auto step : steps) {
                graph->append_step(subpath_handle, step);
            }
            subpaths.push_back(subpath_handle);
        }
    };

    for (size_t i = 0; i < intervals.size(); ++i) {
        if (intervals[i].first > offset) {
            // cut everythign left of the interval
            cut_to(intervals[i].first);
        }
        // scan past the interval
        while (offset < intervals[i].second && current_step != graph->path_end(path_handle)) {
            handle_t handle = graph->get_handle_of_step(current_step);
            offset += graph->get_length(handle);
            current_step = graph->get_next_step(current_step);
#ifdef debug
            cerr << "deleting " << graph->get_id(handle) << endl;
#endif
            chopped_handles.push_back(handle);
        }
    }

    // cut the last bit
    if (offset < original_path_length) {
        cut_to(original_path_length);
    }
    
    return chopped_handles;    
}
