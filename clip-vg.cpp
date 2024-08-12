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

using namespace std;
using namespace handlegraph;
using namespace bdsg;

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph>" << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "options: " << endl
       << "    -b, --bed FILE            Intervals to clip in BED format" << endl
       << "    -m, --min-length N        Only clip paths of length < N" << endl
       << "    -u, --max-unaligned N     Clip out unaligned regions of length > N" << endl
       << "    -a, --anchor PREFIX       If set, consider regions not aligned to a path with PREFIX unaligned (with -u)" << endl
       << "    -e, --ref-prefix STR      Forwardize (but don't clip) paths whose name begins with STR" << endl
       << "    -c, --allow-cycle         Do not fail with error when reference cycle detected" << endl
       << "    -f, --force-clip          Don't abort with error if clipped node overlapped by multiple paths" << endl
       << "    -r, --name-replace S1>S2  Replace (first occurrence of) S1 with S2 in all path names" << endl
       << "    -n, --no-orphan-filter    Don't filter out new subpaths that don't align to anything" << endl
       << "    -d, --drop-path PREFIX    Remove all paths with given PREFIX, and all nodes that are on no other paths (done after other filters)" << endl
       << "    -L, --leave-aligned       When used in conjunction with -d, paths are prserved if they align to a non-dropped path" << endl 
       << "    -o, --out-bed FILE        Save all clipped intervals here" << endl
       << "    -p, --progress            Print progress" << endl
       << endl;
}    

static unordered_map<string, vector<pair<int64_t, int64_t>>> load_bed(istream& bed_stream, const string& ref_prefix);
static unordered_map<string, vector<pair<int64_t, int64_t>>> find_unaligned(const PathHandleGraph* graph, int64_t max_unaligned,
                                                                            const string& ref_prefix, const string& anchor_prefix);
static unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream);
static vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems);
static void chop_path_intervals(MutablePathMutableHandleGraph* graph,
                                const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals,
                                bool force_clip, bool orphan_filter,
                                const string& ref_prefix,
                                bool progress);
static pair<unordered_set<handle_t>, vector<path_handle_t>> chop_path(MutablePathMutableHandleGraph* graph,
                                                                      path_handle_t path_handle,
                                                                      const vector<pair<int64_t, int64_t>>& intervals);
static void replace_path_name_substrings(MutablePathMutableHandleGraph* graph, const vector<string>& to_replace,
                                         bool progress);
static void forwardize_paths(MutablePathMutableHandleGraph* graph, const string& ref_prefix, bool allow_ref_cycles, bool progress);
static vector<unordered_set<nid_t>> weakly_connected_components(const HandleGraph* graph);
static void drop_paths(MutablePathMutableHandleGraph* graph, const string& drop_prefix, bool leave_aligned, bool progress);

static unordered_map<string, vector<pair<int64_t, int64_t>>> get_path_intervals(const PathHandleGraph* graph);

static unordered_map<string, vector<pair<int64_t, int64_t>>> get_clipped_intervals(
    const unordered_map<string, vector<pair<int64_t, int64_t>>>& input_intervals,
    const unordered_map<string, vector<pair<int64_t, int64_t>>>& output_intervals);

// Create a subpath name (todo: make same function in vg consistent (it only includes start))
static inline string make_subpath_name(const string& path_name, size_t offset, size_t length) {
    PathSense sense;
    string sample;
    string locus;
    size_t haplotype;
    size_t phase_block;
    subrange_t subrange;
    PathMetadata::parse_path_name(path_name, sense, sample, locus, haplotype, phase_block, subrange);
    subrange.first = subrange != PathMetadata::NO_SUBRANGE ? subrange.first : 0;
    subrange.first += offset;
    subrange.second = subrange.first + length;
    return PathMetadata::create_path_name(sense, sample, locus, haplotype, phase_block, subrange);
}

int main(int argc, char** argv) {

    string bed_path;
    int64_t min_length = 0;
    int64_t max_unaligned = 0;
    string anchor_prefix;
    string ref_prefix;
    bool allow_ref_cycles = false;
    size_t input_count = 0;
    bool force_clip = false;
    bool orphan_filter = true;
    bool progress = false;
    vector<string> replace_list;
    string drop_prefix;
    bool leave_aligned_drop_paths = false;
    string out_bed_path;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"bed", required_argument, 0, 'b'},
            {"min-length", required_argument, 0, 'm'},
            {"max-unaligned", required_argument, 0, 'u'},
            {"anchor", required_argument, 0, 'a'},
            {"ref-prefix", required_argument, 0, 'e'},
            {"allow-cycle", no_argument, 0, 'c'},
            {"force-clip", no_argument, 0, 'f'},
            {"name-replace", required_argument, 0, 'r'},
            {"no-orphan_filter", no_argument, 0, 'n'},
            {"drop-prefix", required_argument, 0, 'd'},
            {"leave-aligned", no_argument, 0, 'L'},
            {"out-bed", required_argument, 0, 'o'},
            {"progress", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpb:m:u:a:e:cfnr:d:Lo:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'b':
            bed_path = optarg;
            ++input_count;
            break;
        case 'm':
            min_length = stol(optarg);
            ++input_count;
            break;
        case 'u':
            max_unaligned = stol(optarg);
            ++input_count;
            break;
        case 'a':
            anchor_prefix = optarg;
            break;
        case 'e':
            ref_prefix = optarg;
            break;
        case 'c':
            allow_ref_cycles = true;
            break;
        case 'f':
            force_clip = true;
            break;
        case 'n':
            orphan_filter = false;
            break;
        case 'r':
            replace_list.push_back(optarg);
            break;
        case 'd':
            drop_prefix = optarg;
            break;
        case 'L':
            leave_aligned_drop_paths = true;
            break;
        case 'o':
            out_bed_path = optarg;
            break;
        case 'p':
            progress = true;
            break;
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
        cerr << "[clip-vg] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    if (optind != argc - 1) {
        cerr << "[clip-vg] error: too many arguments" << endl;
        help(argv);
        return 1;
    }

    if (input_count > 1) {
        cerr << "[clip-vg] error: at most one of -b, -m or -u can be used at a time" << endl;
        return 1;
    }

    if (input_count == 0 && replace_list.empty() && ref_prefix.empty()) {
        cerr <<  "[clip-vg] error: at east one of -b, -m, -u, -e or -r must be specified" << endl;
        return 1;
    }
    if (!anchor_prefix.empty() && max_unaligned <= 0) {
        cerr << "[clip-vg] error: -a cannot be used without -u" << endl;
        return 1;
    }
    
    if (leave_aligned_drop_paths && drop_prefix.empty()) {
        cerr << "[clip-vg] error: -L can only be used with -d" << endl;
        return 1;
    }

    string graph_path = argv[optind++];
    ifstream graph_stream(graph_path);
    if (!graph_stream) {
        cerr << "[clip-vg] error: Unable to open input graph " << graph_path << endl;
        return 1;
    }    
    unique_ptr<MutablePathMutableHandleGraph> graph = load_graph(graph_stream);
    graph_stream.close();
    if (progress) {
        cerr << "[clip-vg]: Loaded graph" << endl;
    }

    unordered_map<string, vector<pair<int64_t, int64_t>>> input_graph_intervals;
    if (!out_bed_path.empty()) {
        input_graph_intervals = get_path_intervals(graph.get());
        if (progress) {
            cerr << "[clip-vg]: Graph has " << input_graph_intervals.size() << " paths." << endl;
        }
    }

    unordered_map<string, vector<pair<int64_t, int64_t>>> bed_intervals;

    if (!bed_path.empty()) {
        ifstream bed_stream(bed_path);
        if (!bed_stream) {
            cerr << "[clip-vg] error: Unable to open input BED file " << bed_path << endl;
            return 1;
        }
        bed_intervals = load_bed(bed_stream, ref_prefix);
    } else if (min_length != 0) {
        // apply min length to all paths to get intervals
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                if (ref_prefix.empty() || path_name.substr(0, ref_prefix.length()) != ref_prefix) {
                    int64_t path_length = 0;
                    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                            path_length += graph->get_length(graph->get_handle_of_step(step_handle));
                            return path_length < min_length;
                        });
                    if (path_length < min_length) {
                        bed_intervals[path_name].push_back(make_pair(0, path_length));
                    }
                }
            });
    } else if (max_unaligned != 0) {
        // apply max unaligned length to all paths
        if (progress) {
            cerr << "[clip-vg]: Finding unaligned intervals >= " << max_unaligned
                 << " using anchor prefix " << anchor_prefix << " and ref prefix " << ref_prefix << endl;
        }
        bed_intervals = find_unaligned(graph.get(), max_unaligned, ref_prefix, anchor_prefix);
    }
    
    if (progress) {
        size_t num_intervals = 0;
        for (auto& bi : bed_intervals) {
            num_intervals += bi.second.size();
        }
        cerr << "[clip-vg]: Loaded " << num_intervals << " BED intervals over " << bed_intervals.size() << " sequences" << endl;
    }

    if (!bed_intervals.empty()) {
        chop_path_intervals(graph.get(), bed_intervals, force_clip, orphan_filter, ref_prefix, progress);
    }

    if (!ref_prefix.empty()) {
        forwardize_paths(graph.get(), ref_prefix, allow_ref_cycles, progress);  
    }
    
    if (!replace_list.empty()) {
        replace_path_name_substrings(graph.get(), replace_list, progress);
    }

    if (!drop_prefix.empty()) {
        drop_paths(graph.get(), drop_prefix, leave_aligned_drop_paths, progress);
    }

    if (!out_bed_path.empty()) {
        unordered_map<string, vector<pair<int64_t, int64_t>>> output_graph_intervals = get_path_intervals(graph.get());
#ifdef debug
        for (const auto& xx : output_graph_intervals) {
            cerr << " got output intervals " << xx.first << " count = " << xx.second.size() << endl;
        }
#endif
        unordered_map<string, vector<pair<int64_t, int64_t>>> clipped_graph_intervals = get_clipped_intervals(input_graph_intervals, output_graph_intervals);
        ofstream out_bed_file(out_bed_path);
        size_t icount = 0;
        for (const auto& pi : clipped_graph_intervals) {
            for (const auto& i : pi.second) {
                out_bed_file << pi.first << "\t" << i.first << "\t" << i.second << "\n";
                ++icount;
            }
        }
        out_bed_file.flush();
        if (progress) {
            cerr << "[clip-vg]: Outputted " << icount << " clipped intervals to " << out_bed_path << endl;
        }
    }

    dynamic_cast<SerializableHandleGraph*>(graph.get())->serialize(cout);

    return 0;
}

unordered_map<string, vector<pair<int64_t, int64_t>>> load_bed(istream& bed_stream, const string& ref_prefix) {
    // load bed
    unordered_map<string, vector<pair<int64_t, int64_t>>> intervals;
    string buffer;
    while (getline(bed_stream, buffer)) {
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);
        if (toks.size() >= 3) {
            string& name = toks[0];
            if (ref_prefix.empty() || name.substr(0, ref_prefix.length()) != ref_prefix) {
                int64_t start = stol(toks[1]);
                int64_t end = stol(toks[2]);
                intervals[name].push_back(make_pair(start, end));
            }
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

unordered_map<string, vector<pair<int64_t, int64_t>>> find_unaligned(const PathHandleGraph* graph, int64_t max_unaligned,
                                                                     const string& ref_prefix, const string& anchor_prefix) {
    unordered_map<string, vector<pair<int64_t, int64_t>>> intervals;

    // anchor-prefix means we consider a node unaligned if it doesn't align to a path with that prefix
    // to do this check, we need a table of nodes on these paths:
    unordered_set<nid_t> minigraph_nodes;
    if (!anchor_prefix.empty()) {
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                string path_name = graph->get_path_name(path_handle);
                if (path_name.compare(0, anchor_prefix.length(), anchor_prefix) == 0) {
                    graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        minigraph_nodes.insert(graph->get_id(graph->get_handle_of_step(step_handle)));
                    });
                }
            });
    }
    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (ref_prefix.empty() || path_name.substr(0, ref_prefix.length()) != ref_prefix) {
                int64_t offset = 0;
                int64_t start = -1;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        int64_t len = (int64_t)graph->get_length(handle);
                        bool aligned = minigraph_nodes.count(graph->get_id(handle));
                        if (!aligned && anchor_prefix.empty()) {
                            graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle_2) {
                                if (graph->get_path_handle_of_step(step_handle_2) != path_handle) {
                                    aligned = true;
                                }
                                return !aligned;
                            });
                        }
                        // start an unaligned interval
                        if (start < 0 && aligned == false) {
                            start = offset;
                        }
                        // end an unaligned interval
                        if (aligned == true) {
                            if (start >= 0 && offset - start > max_unaligned) {
                                intervals[path_name].push_back(make_pair(start, offset));
                            }
                            start = -1;
                        }
                        offset += len;
                    });
                if (start >= 0 && offset - start > max_unaligned) {
                    intervals[path_name].push_back(make_pair(start, offset));
                }
            }
        });
    return intervals;
}


unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream) {

    char magic_bytes[4];
    graph_stream.read(magic_bytes, 4);
    uint32_t magic_number = ntohl(*((uint32_t*) magic_bytes));
    graph_stream.clear();
    graph_stream.seekg(0, ios::beg);

    MutablePathMutableHandleGraph* graph;
    if (magic_number == PackedGraph().get_magic_number()) {
        graph = new PackedGraph();
    } else if (magic_number == HashGraph().get_magic_number()) {
        graph = new HashGraph();
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
                         const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals,
                         bool force_clip, bool orphan_filter, const string& ref_prefix,
                         bool progress) {

    // keep some stats to print
    size_t chopped_paths = 0;
    size_t chopped_nodes = 0;
    size_t chopped_bases = 0;
    
    // careful not to iterate and chop, as we could hit new subpaths made
    vector<path_handle_t> path_handles;    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });

    // when force_clip is true, store handles here to given them second chance at destruction
    // after all paths are deleted
    unordered_set<nid_t> to_destroy;
    // newly created subpaths
    vector<path_handle_t> subpaths;
    // paths to destroy (faster to do in single api call)
    vector<path_handle_t> paths_to_destroy;
    
    for (auto path_handle : path_handles) {
        string path_name = graph->get_path_name(path_handle);
        auto it = bed_intervals.find(path_name);
        bool was_chopped = false;
        if (it != bed_intervals.end()) {
            if (progress) {
                cerr << "[clip-vg]: Clipping " << it->second.size() << " intervals from path " << path_name << endl;
            }
            auto chopped_handles_subpaths = chop_path(graph, path_handle, it->second);
            auto& chopped_handles = chopped_handles_subpaths.first;
            subpaths.insert(subpaths.end(), chopped_handles_subpaths.second.begin(), chopped_handles_subpaths.second.end());
            if (!chopped_handles.empty()) {
#ifdef debug
                cerr << "adding path to destroy list" << graph->get_path_name(path_handle) << endl;
#endif
                paths_to_destroy.push_back(path_handle);
                for (handle_t handle : chopped_handles) {
                    if (force_clip) {
                        to_destroy.insert(graph->get_id(handle));
                    } else {
                        vector<step_handle_t> steps = graph->steps_of_handle(handle);
                        bool aligned = false;
                        for (size_t i = 0; i < steps.size() && !aligned; ++i) {
                            string other_path_name = graph->get_path_name(graph->get_path_handle_of_step(steps[i]));
                            if (path_name.substr(0, path_name.rfind("[")) !=
                                other_path_name.substr(0, other_path_name.rfind("["))) {
                                aligned = true;
                            }
                        }
                        if (!aligned) {
                            chopped_bases += graph->get_length(handle);
                            was_chopped = true;
                            ++chopped_nodes;
                            to_destroy.insert(graph->get_id(handle));
                        } else {
                            cerr << "[clip-vg]: Unable to clip node " << graph->get_id(handle) << ":" << graph->get_is_reverse(handle)
                                 << " in path " << path_name << " because it is found in the following other paths:\n";
                            for (step_handle_t step : graph->steps_of_handle(handle)) {
                                cerr <<"\t" << graph->get_path_name(graph->get_path_handle_of_step(step)) << endl;
                            }
                            cerr << " Use the -f option to not abort in this case" << endl;
                            exit(1);
                        }
                    }
                }
            }
        }
        if (was_chopped) {
            ++chopped_paths;
        }
    }

    // delete all the paths
#ifdef debug
    cerr << "destroying " << paths_to_destroy.size() << " paths" << endl;
#endif
    graph->destroy_paths(paths_to_destroy);
    paths_to_destroy.clear();

    // trim out fragments between clipped regions that would otherwise be left disconnected from the graph
    size_t removed_subpath_count = 0;
    size_t removed_subpath_base_count = 0;
    size_t removed_component_count = 0;
    size_t removed_component_base_count = 0;
    if (orphan_filter) {
        for (path_handle_t subpath_handle : subpaths) {
            bool connected = false;
            graph->for_each_step_in_path(subpath_handle, [&](step_handle_t step_handle) {
                    connected = graph->steps_of_handle(graph->get_handle_of_step(step_handle)).size() > 1;
                    return !connected;
                });
            if (!connected) {
                graph->for_each_step_in_path(subpath_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        to_destroy.insert(graph->get_id(handle));
                        removed_subpath_base_count += graph->get_length(handle);
                    });
                paths_to_destroy.push_back(subpath_handle);
                if (progress) {
                    cerr << "[clip-vg]: Removing orphaned subpath " << graph->get_path_name(subpath_handle) << endl;
                }
                ++removed_subpath_count;
            }
        }
        graph->destroy_paths(paths_to_destroy);

        // use the reference path prefix (if given) to clip out components that aren't anchored to it
        // (this would take care of above filter, but we leave that one as it's not dependent on path name)
        if (!ref_prefix.empty()) {
            vector<unordered_set<nid_t>> components = weakly_connected_components(graph);
            for (auto& component : components) {
                bool ref_anchored = false;
                for (auto ni = component.begin(); !ref_anchored && ni != component.end(); ++ni) {
                    vector<step_handle_t> steps = graph->steps_of_handle(graph->get_handle(*ni));
                    for (size_t si = 0; !ref_anchored && si < steps.size(); ++si) {
                        string step_path_name = graph->get_path_name(graph->get_path_handle_of_step(steps[si]));
                        if (step_path_name.substr(0, ref_prefix.length()) == ref_prefix) {
                            ref_anchored = true;
                        }
                    }
                }
                if (!ref_anchored) {
                    ++removed_component_count;
                    for (auto node_id : component) {
                        handle_t node_handle = graph->get_handle(node_id);
                        removed_component_base_count += graph->get_length(node_handle);
                        // destroy here instead of adding to to_destroy, becuase we don't care
                        // if there are paths or not (so don't require -f)
                        dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(node_handle);
                        if (to_destroy.count(node_id)) {
                            to_destroy.erase(node_id);
                        }
                    }
                }
            }
        }
    }
    
    for (nid_t nid : to_destroy) {
        assert(graph->has_node(nid));
        handle_t handle = graph->get_handle(nid);
        if (graph->steps_of_handle(handle).empty()) {
            chopped_bases += graph->get_length(handle);
            ++chopped_nodes;
            dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(handle);
#ifdef debug
            cerr << "force destroying handle " << graph->get_id(handle) << ":" << graph->get_is_reverse(handle) << endl;
#endif
        }
    }
    
    if (progress) {
        cerr << "[clip-vg]: Clipped "
             << chopped_bases << " bases from "
             << chopped_nodes << " nodes";
        if (!force_clip) {
            cerr << " in " << chopped_paths << " paths";
        }
        cerr << endl;        
        if (removed_subpath_count > 0) {
            cerr << "[clip-vg]: " << removed_subpath_count << " orphaned subpaths were removed with total "
                 << removed_subpath_base_count << " bases" << endl;
        }
        if (removed_component_count > 0) {
            cerr << "[clip-vg]: " << removed_component_count << " orphaned connected components were removed with total "
                 << removed_component_base_count << " bases" << endl;
        }
    }
}

pair<unordered_set<handle_t>, vector<path_handle_t>> chop_path(MutablePathMutableHandleGraph* graph,
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
        if (!cut_points.empty()) {
#ifdef debug
            cerr << "dividing node_id=" << graph->get_id(handle) << ":" << graph->get_is_reverse(handle) << " seq=" << graph->get_sequence(handle)
                 << " for path " << graph->get_path_name(path_handle) << " at cut points:";
            for (auto cp : cut_points) {
                cerr << " " << cp;
            }
            cerr << endl;
#endif        
            size_t total_pieces_length = 0;
            vector<handle_t> pieces = graph->divide_handle(handle, cut_points) ;
            for (size_t i = 0; i < pieces.size(); ++i) {
                handle_t& piece = pieces[i];
                size_t piece_length = graph->get_length(piece);
                if (i == 0) {
                    assert(piece_length == cut_points[0]);
                } else if (i < pieces.size() - 1) {
                    assert(piece_length == cut_points[i] - cut_points[i-1]);
                }
                total_pieces_length += piece_length;
#ifdef debug
                cerr << " piece " << graph->get_id(piece) << ":" << graph->get_is_reverse(piece) << " " << graph->get_sequence(piece)
                     << " tlen=" << total_pieces_length << "/" << len << endl;
#endif
            }
            // bugs in divide-handle turning out to be a real issue.  add this sanity check to catch them early.
            assert(total_pieces_length == (size_t)len);
        }
        offset += len;
    }
    
    steps.clear();
    int64_t original_path_length = offset;
    unordered_set<handle_t> chopped_handles;
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
#ifdef debug
                cerr << " pushing subpath step " << graph->get_id(step) << ":" << graph->get_is_reverse(step)
                     << " len=" << graph->get_length(step) <<  " to " << graph->get_path_name(subpath_handle) << endl;
#endif
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
            cerr << "adding to delete set: " << graph->get_id(handle) << endl;
#endif
            chopped_handles.insert(handle);
        }
    }

    // cut the last bit
    if (offset < original_path_length) {
        cut_to(original_path_length);
    }
    
    return make_pair(chopped_handles, subpaths);    
}

void replace_path_name_substrings(MutablePathMutableHandleGraph* graph, const vector<string>& to_replace,
                                  bool progress) {
    // parse the strings
    vector<pair<string, string>> replace;
    for (const string& repstring : to_replace) {
        size_t sep = repstring.find('>');
        if (sep == string::npos || sep == 0 || sep == repstring.length() - 1) {
            cerr << "[clip-vg]: Unable to find separator '>' in " << repstring << ". Replacement must be"
                 << " specified with \"s1>s2\"" << endl;
            exit(1);
        }
        replace.push_back(make_pair(repstring.substr(0, sep), repstring.substr(sep + 1)));
        if (replace.back().first == replace.back().second) {
            replace.pop_back();
        }
    }

    size_t replacement_count = 0;
    size_t path_count = 0;
    // take care to not modify path handles while iterating path handles, just in case
    vector<string> path_names;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_names.push_back(graph->get_path_name(path_handle));
        });
    vector<path_handle_t> paths_to_destroy;
    for (string& path_name : path_names) {
        path_handle_t path_handle = graph->get_path_handle(path_name);
        bool changed = false;
        for (auto& rep : replace) {
            size_t p = path_name.find(rep.first);
            if (p != string::npos) {
                path_name.replace(p, rep.first.length(), rep.second);
                ++replacement_count;
                changed = true;
            }
        }
        if (changed) {
            ++path_count;
            if (graph->has_path(path_name)) {
                cerr << "[clip-vg] error: cannot change name of path from " << graph->get_path_name(path_handle) << " to "
                     << path_name << " because the latter already exists in the graph" << endl;
                exit(1);
            }
            path_handle_t new_path_handle = graph->create_path_handle(path_name, graph->get_is_circular(path_handle));
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    graph->append_step(new_path_handle, graph->get_handle_of_step(step_handle));
                });
            paths_to_destroy.push_back(path_handle);
        }
    }
    graph->destroy_paths(paths_to_destroy);
    
    if (progress) {
        cerr << "[clip-vg]: Replaced " << replacement_count << " substrings in " << path_count << " path names" << endl;
    }
}

void forwardize_paths(MutablePathMutableHandleGraph* graph, const string& ref_prefix, bool allow_ref_cycles, bool progress) {
    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (path_name.substr(0, ref_prefix.length()) == ref_prefix) {
                size_t fw_count = 0;
                size_t total_steps = 0;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        ++total_steps;
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        if (graph->get_is_reverse(handle)) {
                            vector<step_handle_t> steps = graph->steps_of_handle(handle);
                            size_t ref_count = 0;
                            for (step_handle_t step : steps) {
                                if (graph->get_path_handle_of_step(step) == path_handle) {
                                    ++ref_count;
                                }
                                if (ref_count > 1) {
                                    break;
                                }
                            }
                            if (ref_count > 1) {
                                if (allow_ref_cycles) {
                                    // todo: should be able to forwardize ref cycle if all steps are reverse
                                    return;
                                } else {
                                    cerr << "[clip-vg] error: Cycle detected in reference path " << path_name << " at node " << graph->get_id(handle) << endl;
                                exit(1);
                                }
                            }
                            handle_t flipped_handle = graph->create_handle(graph->get_sequence(handle));
			    graph->follow_edges(handle, true, [&](handle_t prev_handle) {
                                    if (graph->get_id(prev_handle) != graph->get_id(handle)) {
                                        graph->create_edge(prev_handle, flipped_handle);
                                    }
			      });
			    graph->follow_edges(handle, false, [&](handle_t next_handle) {
                                    if (graph->get_id(handle) != graph->get_id(next_handle)) {
                                        graph->create_edge(flipped_handle, next_handle);
                                    }
			      });
                            // self-loop cases we punted on above:
                            if (graph->has_edge(handle, handle)) {
                                graph->create_edge(flipped_handle, flipped_handle);
                            }
                            if (graph->has_edge(handle, graph->flip(handle))) {
                                graph->create_edge(flipped_handle, graph->flip(flipped_handle));                                
                            }
                            if (graph->has_edge(graph->flip(handle), handle)) {
                                graph->create_edge(graph->flip(flipped_handle), flipped_handle);
                            }
                            for (step_handle_t step : steps) {
                                step_handle_t next_step = graph->get_next_step(step);
                                handle_t new_handle = graph->get_is_reverse(graph->get_handle_of_step(step)) ? flipped_handle :
                                    graph->flip(flipped_handle);
                                graph->rewrite_segment(step, next_step, {new_handle});
                            }
                            ++fw_count;
                            assert(graph->steps_of_handle(handle).empty());
                            dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(handle);
                        }
                    });
                if (fw_count > 0 && progress) {
                    cerr << "[clip-vg]: Forwardized " << fw_count << " / " << total_steps << " steps in reference path " << path_name << endl;
                }
            }
        });

    if (!allow_ref_cycles) {
        // do a check just to be sure
        graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (path_name.substr(0, ref_prefix.length()) == ref_prefix) {
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    handle_t handle = graph->get_handle_of_step(step_handle);
                    if (graph->get_is_reverse(handle)) {
                        cerr << "[clip-vg] error: Failed to fowardize node " << graph->get_id(handle) << " in path " << path_name << endl;
                        exit(1);
                    }
                });
            }
        });
    }
}

// this is pasted from libhandlegraph
// todo: update libhandlegraph to version that contains algorithms!!!
vector<unordered_set<nid_t>> weakly_connected_components(const HandleGraph* graph) {
    vector<unordered_set<nid_t>> to_return;
    
    // This only holds locally forward handles
    unordered_set<handle_t> traversed;
    
    graph->for_each_handle([&](const handle_t& handle) {
        
        // Only think about it in the forward orientation
        auto forward = graph->forward(handle);
        
        if (traversed.count(forward)) {
            // Already have this node, so don't start a search from it.
            return;
        }
        
        // The stack only holds locally forward handles
        vector<handle_t> stack{forward};
        to_return.emplace_back();
        while (!stack.empty()) {
            handle_t here = stack.back();
            stack.pop_back();
            
            traversed.insert(here);
            to_return.back().insert(graph->get_id(here));
            
            // We have a function to handle all connected handles
            auto handle_other = [&](const handle_t& other) {
                // Again, make it forward
                auto other_forward = graph->forward(other);
                
                if (!traversed.count(other_forward)) {
                    stack.push_back(other_forward);
                }
            };
            
            // Look at edges in both directions
            graph->follow_edges(here, false, handle_other);
            graph->follow_edges(here, true, handle_other);
            
        }
    });
    return to_return;
}

// this was written to filter out minigraph-only nodes from the graph (and all the minigraph paths)
// this used to get done by hal2vg, but it's useful to keep them around so that -u option will work
// better.  ie this way, a path that's private to a sample but still in the minigraph will be kept
// because it aligns to two paths whereas if minigraph paths weren't in, it'd be deleted
void drop_paths(MutablePathMutableHandleGraph* graph, const string& drop_prefix, bool leave_aligned, bool progress) {

    unordered_set<nid_t> to_destroy;
    size_t removed_path_count = 0;
    size_t removed_base_count = 0;

    // careful not to iterate and chop, as we could hit new subpaths made
    vector<path_handle_t> path_handles;    
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            path_handles.push_back(path_handle);
        });

    for (path_handle_t& path_handle : path_handles) {
        string path_name = graph->get_path_name(path_handle);
        if (path_name.compare(0, drop_prefix.length(), drop_prefix) == 0) {
            // we've found a path with the given prefix: now destroy all handles that don't touch
            // any path *without* the prefix
            size_t offset = 0;
            vector<pair<int64_t, int64_t>> intervals;
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    handle_t handle = graph->get_handle_of_step(step_handle);
                    bool has_other_path = false;
                    size_t len = graph->get_length(handle);
                    graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle_2) {
                            path_handle_t other_path_handle = graph->get_path_handle_of_step(step_handle_2);
                            if (other_path_handle != path_handle &&
                                graph->get_path_name(other_path_handle).compare(0, drop_prefix.length(), drop_prefix) != 0) {
                                has_other_path = true;
                                return false;
                            }
                            return true;
                        });
                    if (!has_other_path) {
                        to_destroy.insert(graph->get_id(handle));
                        if (!intervals.empty() && offset == intervals.back().second) {
                            intervals.back().second += len;
                        } else {
                            intervals.push_back(make_pair(offset, offset + len));
                        }
                    }
                    offset += len;
                });
            if (leave_aligned && !intervals.empty()) {
                // chop the path (to keep fragments for any non-destroyed nodes)
                chop_path(graph, path_handle, intervals);
            }
            if (!leave_aligned || !intervals.empty()) {
                // detroy the path            
                graph->destroy_path(path_handle);
                removed_path_count++;
            }
        }
        }

    // destroy the nodes
    size_t removed_node_count = to_destroy.size();
    for (nid_t node_id : to_destroy) {
        handle_t node_handle = graph->get_handle(node_id);
        removed_base_count += graph->get_length(node_handle);
        dynamic_cast<DeletableHandleGraph*>(graph)->destroy_handle(node_handle);        
    }

    if (progress) {
        cerr << "[clip-vg]: Drop prefix removed " << removed_base_count << " bases from " << removed_node_count << " nodes in " << removed_path_count << " paths" << endl;        
    }
}
unordered_map<string, vector<pair<int64_t, int64_t>>> get_path_intervals(const PathHandleGraph* graph) {
    unordered_map<string, vector<pair<int64_t, int64_t>>> path_intervals;
    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            size_t path_len = 0;
            graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                    path_len += graph->get_length(graph->get_handle_of_step(step_handle));
                });
            subrange_t path_range = graph->get_subrange(path_handle);
            int64_t start = path_range == PathMetadata::NO_SUBRANGE ? 0 : path_range.second;
            vector<pair<int64_t, int64_t>>& intervals = path_intervals[path_name];        
            intervals.push_back(make_pair(start, start + path_len));
        });

    for (auto& pi : path_intervals) {
        std::sort(pi.second.begin(), pi.second.end(), [](const pair<int64_t, int64_t>& i1, const pair<int64_t, int64_t>& i2) {
                return i1.first < i2.first || (i1.first == i2.first && i1.second < i2.second);
            });
    }

    return path_intervals;
}

static unordered_map<string, vector<pair<int64_t, int64_t>>> get_clipped_intervals(
    const unordered_map<string, vector<pair<int64_t, int64_t>>>& input_intervals,
    const unordered_map<string, vector<pair<int64_t, int64_t>>>& output_intervals) {

    unordered_map<string, vector<pair<int64_t, int64_t>>> clipped_intervals;

    for (const auto& input_pi : input_intervals) {
        const string& path_name = input_pi.first;
        const auto& in_intervals = input_pi.second;
        if (!output_intervals.count(path_name)) {
            // path doesn't appear in output -> everything was clipped
            clipped_intervals[path_name].insert(clipped_intervals[path_name].end(), in_intervals.begin(), in_intervals.end());
#ifdef debug
            cerr << "clippin everything for " << path_name << endl;
#endif
        } else {
#ifdef debug
            cerr << "doin frag clip for " << path_name << endl;
#endif
            const auto& out_intervals = output_intervals.at(path_name);
            // note: clipping here is fairly simple because output intervals are a subset
            // and nothing overlaps
            int64_t j = 0;
            int64_t k = 0;
            for (int64_t i = 0; i < in_intervals.size(); ++i) {
                // j: first out_interval that's not completely left of in_interval;
                for (; j < out_intervals.size() && out_intervals[j].second < in_intervals[i].first; ++j);
                // k: first out_interval that's completely right of in_interval
                for (k = j; k < out_intervals.size() && out_intervals[k].first < in_intervals[i].second; ++k);
                // [j, k) all out_intervals that are sub intervals of in_interval[i]
                vector<pair<int64_t, int64_t>>&  gap_intervals = clipped_intervals[path_name];
                if (k == j) {
                    // nothing overlaps -- everything clipped
                    gap_intervals.push_back(in_intervals[i]);
                } else {
                    // left of first interval                    
                    if (out_intervals[j].first > in_intervals[i].first) {
                        gap_intervals.push_back(make_pair(in_intervals[i].first, out_intervals[j].first));
                    }
                    // gaps between out_intervals
                    for (int64_t l = 0; l < (int64_t)out_intervals.size() - 2; ++l) {
                        if (out_intervals[l+1].first > out_intervals[l].second) {
                            gap_intervals.push_back(make_pair(out_intervals[l].second, out_intervals[l+1].first));
                        }
                    }
                    // right of last interval
                    if (in_intervals[i].second > out_intervals.back().second) {
                        gap_intervals.push_back(make_pair(out_intervals.back().second, in_intervals[i].second));
                    }
                }
            }
        }
    }

    return clipped_intervals;
}
