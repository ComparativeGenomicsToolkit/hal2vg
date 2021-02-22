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
  cerr << "usage: " << argv[0] << " [options] <graph>" << endl
       << "Chop out path intervals from a vg graph" << endl
       << endl
       << "options: " << endl
       << "    -b, --bed FILE            Intervals to clip in BED format" << endl
       << "    -m, --min-length N        Only clip paths of length < N" << endl
       << "    -u, --max-unaligned N     Clip out unaligned regions of length > N" << endl
       << "    -e, --ref-prefix STR      Ignore paths whose name begins with STR" << endl
       << "    -f, --force-clip          Don't abort with error if clipped node overlapped by multiple paths" << endl
       << "    -r, --name-replace S1>S2  Replace (first occurrence of) S1 with S2 in all path names" << endl
       << "    -n, --no-orphan-filter    Don't filter out new subpaths that don't align to anything" << endl
       << "    -p, --progress            Print progress" << endl
       << endl;
}    

static unordered_map<string, vector<pair<int64_t, int64_t>>> load_bed(istream& bed_stream, const string& ref_prefix);
static unordered_map<string, vector<pair<int64_t, int64_t>>> find_unaligned(const PathHandleGraph* graph, int64_t max_unaligned,
                                                                            const string& ref_prefix);
static unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream);
static vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems);
static void chop_path_intervals(MutablePathMutableHandleGraph* graph,
                                const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals,
                                bool force_clip, bool orphan_filter,
                                bool progress);
static pair<unordered_set<handle_t>, vector<path_handle_t>> chop_path(MutablePathMutableHandleGraph* graph,
                                                                      path_handle_t path_handle,
                                                                      const vector<pair<int64_t, int64_t>>& intervals);
static void replace_path_name_substrings(MutablePathMutableHandleGraph* graph, const vector<string>& to_replace,
                                         bool progress);
// Create a subpath name (todo: make same function in vg consistent (it only includes start))
static inline string make_subpath_name(const string& path_name, size_t offset, size_t length) {
    return path_name + "[" + std::to_string(offset) + "-" + std::to_string(offset + length) + "]";
}

int main(int argc, char** argv) {

    string bed_path;
    int64_t min_length = 0;
    int64_t max_unaligned = 0;
    string ref_prefix;
    size_t input_count = 0;
    bool force_clip = false;
    bool orphan_filter = true;
    bool progress = false;
    vector<string> replace_list;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"bed", required_argument, 0, 'b'},
            {"min-length", required_argument, 0, 'm'},
            {"max-unaligned", required_argument, 0, 'u'},
            {"ref-prefix", required_argument, 0, 'e'},
            {"force-clip", no_argument, 0, 'f'},
            {"name-replace", required_argument, 0, 'r'},
            {"no-orphan_filter", no_argument, 0, 'n'},
            {"progress", no_argument, 0, 'p'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpb:m:u:e:fnr:",
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
        case 'e':
            ref_prefix = optarg;
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

    if (input_count == 0 && replace_list.empty()) {
        cerr <<  "[clip-vg] error: at east one of -b, -m, -u or -r must be specified" << endl;
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
        bed_intervals = find_unaligned(graph.get(), max_unaligned, ref_prefix);
    }
    
    if (progress) {
        size_t num_intervals = 0;
        for (auto& bi : bed_intervals) {
            num_intervals += bi.second.size();
        }
        cerr << "[clip-vg]: Loaded " << num_intervals << " BED intervals over " << bed_intervals.size() << " sequences" << endl;
    }

    if (!bed_intervals.empty()) {
        chop_path_intervals(graph.get(), bed_intervals, force_clip, orphan_filter, progress);
    }

    if (!replace_list.empty()) {
        replace_path_name_substrings(graph.get(), replace_list, progress);
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
                                                                     const string& ref_prefix) {
    unordered_map<string, vector<pair<int64_t, int64_t>>> intervals;

    graph->for_each_path_handle([&](path_handle_t path_handle) {
            string path_name = graph->get_path_name(path_handle);
            if (ref_prefix.empty() || path_name.substr(0, ref_prefix.length()) != ref_prefix) {
                int64_t offset = 0;
                int64_t start = -1;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        int64_t len = (int64_t)graph->get_length(handle);
                        bool aligned = false;
                        graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle_2) {
                                if (graph->get_path_handle_of_step(step_handle_2) != path_handle) {
                                    aligned = true;
                                }
                                return !aligned;
                            });
                        // start an unaligned interval
                        if (start < 0 && aligned == false) {
                            start = offset;
                        }
                        // end an unaligned interval
                        if (aligned == true) {
                            if (start >= 0 && offset + len - start > max_unaligned) {
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
                         const unordered_map<string, vector<pair<int64_t, int64_t>>>& bed_intervals,
                         bool force_clip, bool orphan_filter,
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
                cerr << "destroying path " << graph->get_path_name(path_handle) << endl;
#endif
                graph->destroy_path(path_handle);
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

    // trim out fragments between clipped regions that would otherwise be left disconnected from the graph
    size_t removed_subpath_count = 0;
    size_t removed_subpath_base_count = 0;
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
                graph->destroy_path(subpath_handle);
                if (progress) {
                    cerr << "[clip-vg]: Removing orphaned subpath " << graph->get_path_name(subpath_handle) << endl;
                }
                ++removed_subpath_count;
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
            graph->destroy_path(path_handle);
        }
    }
    
    if (progress) {
        cerr << "[clip-vg]: Replaced " << replacement_count << " substrings in " << path_count << " path names" << endl;
    }
}
