// Filter big deletions from PAF.  These can sometimes arise from minigraph split alignments.  They are rare, but can really
// mess up topology of graph if even one gets into cactus.

// 1) Estimate anchors along reference path for every node in input graph
// 2) Scan every query in PAF in order, and look at target blocks
// 3) Use table from 1) in order to estimate the distances between target blocks
// 4) If two conesecutive target blocks span too big of a distance, delete the smaller block

//#define debug

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>
#include <unistd.h>
#include <getopt.h>
#include <omp.h>

#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"

#include "IntervalTree.h"
#include "paf.hpp"

using namespace std;
using namespace handlegraph;
using namespace bdsg;

struct Anchor {
    path_handle_t path_handle;
    int64_t max_offset;
    int64_t min_offset;
};

struct PafDelta {
    int64_t delta;
    int64_t ref_delta;
    int64_t query_delta;
    int64_t ref_overlap_size;
    int64_t prev_ref_start;
    int64_t prev_ref_end;
    int64_t cur_ref_start;
    int64_t cur_ref_end;
};
    
static unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream);
static pair<unordered_map<string, nid_t>, unordered_map<nid_t, string>> load_trans(const string& trans_path);
static unordered_map<nid_t, Anchor> index_graph(const PathHandleGraph* graph,
                                                const string& ref_prefix);
static unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>> index_deletions(const PathHandleGraph* graph, const unordered_map<nid_t, Anchor>& index);
static pair<vector<PafLine>, unordered_map<int64_t, int64_t>> load_paf(ifstream& paf_file);
static int64_t for_each_query_block(const vector<PafLine>& paf_lines, const vector<bool>& masking,
                                    function<void(int64_t, int64_t)> visit_block);
static int64_t check_delta(int64_t max_deletion_threshold, int64_t max_insertion_threshold, const PafDelta& paf_delta, double overlap_threshold);
static PafDelta get_delta(path_handle_t ref_path, const PafLine& prev_paf, const PafLine& cur_paf,
                          const unordered_map<string, nid_t>& mg_to_vg, const unordered_map<nid_t, Anchor>& ref_index,
                          const unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>>& ref_deletions);
                                                             
void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <graph.vg> <trans> <aln.paf> -d <del-threshold>\n" << endl
         << "Use distances from graph to filter out implied deletions from PAF (cigars not considered, only blocks)" << endl
         << "  <graph.vg> : minigraph as obtained from vg convert -g graph.gfa" << endl
         << "  <trans> : node translation from vg convert -g -T" << endl
         << "  <aln.paf> : paf alignment from cactus-graphmap" << endl
         << endl
         << "options: " << endl
         << "    -d  --del-threshold F     only remove deletions greater than this. if < 1, then interpreted as fraction of reference path size" << endl
         << "    -i, --ins-threshold F     Like <del-threshold>, but applied to insertions instead of deletions [-1]" << endl
         << "    -m, --max-filter F        If F*<threshold> matches need to be pulled apart to resolve a single deletion, just leave it alone [1]" << endl
         << "    -r, --ref-prefix STR      Only consider paths whose names start with STR" << endl
         << "    -p, --progress            Print progress" << endl
         << "    -o, --filter-off-ref      Filter mappings that aren't in dominant ref" << endl        
         << "    -v, --verbose             Print deletions" << endl
         << "    -t, --threads N           number of threads to use (used only for indexing graph) [default: all available]" << endl
       << endl;
}    

int main(int argc, char** argv) {

    string ref_prefix;
    bool progress = false;
    bool verbose = false;
    bool keep_off_ref = true;
    // only filter deletions that don't overlap an existing deletion by at least this much
    // (doesn't seem to a factor -- most big deletions not in minigraph)
    double overlap_threshold = 0.5;
    double filter_threshold = 1.0;
    double max_insertion = -1.0;
    double max_deletion = -1.0;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"del-threshold", required_argument, 0, 'd'},            
            {"ins-threshold", required_argument, 0, 'i'},
            {"max-filter", required_argument, 0, 'm'},            
            {"ref-prefix", required_argument, 0, 'r'},
            {"filter-off-ref", no_argument, 0, 'o'},
            {"help", no_argument, 0, 'h'},            
            {"progress", no_argument, 0, 'p'},
            {"verbose", no_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "d:i:m:r:khpvt:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'd':
            max_deletion = stof(optarg);
            break;                                                
        case 'i':
            max_insertion = stof(optarg);
            break;                                    
        case 'm':
            filter_threshold = stof(optarg);
            break;
        case 'r':
            ref_prefix = optarg;
            break;
        case 'o':
            keep_off_ref = false;
            break;
        case 'v':
            verbose = true;
            break;
        case 'p':
            progress = true;
            break;
        case 't':
        {
            int num_threads = stoi(optarg);
            if (num_threads <= 0) {
                cerr << "[filter-paf-deletions] error: Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                exit(1);
            }
            omp_set_num_threads(num_threads);
            break;
        }                        
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

    if (argc <= 3) {
        cerr << "[filter-paf-deletions] error: too few arguments\n" << endl;
        help(argv);
        return 1;
    }

    // Parse the positional argument
    if (optind >= argc) {
        cerr << "[filter-paf-deletions] error: too few arguments\n" << endl;
        help(argv);
        return 1;
    }

    if (optind != argc - 3) {
        cerr << "[filter-paf-deletions] error: too many arguments\n" << endl;
        help(argv);
        return 1;
    }

    if (max_deletion <= 0 && max_insertion <= 0) {
        cerr << "[filter-paf-deletions] error: at least one of -d or -i must be set to positive value" << endl;
        return 1;
    }
    
    string graph_path = argv[optind++];
    string trans_path = argv[optind++];
    string paf_path = argv[optind++];

    // load the graph
    ifstream graph_stream(graph_path);
    if (!graph_stream) {
        cerr << "[filter-paf-deletions] error: Unable to open input graph " << graph_path << endl;
        return 1;
    }    
    unique_ptr<PathHandleGraph> graph = load_graph(graph_stream);
    graph_stream.close();
    if (progress) {
        cerr << "[filter-paf-deletions]: Loaded graph" << endl;
    }

    // load the minigraph <-> vg id translation table (because our PAF is expressed in terms of the minigraph
    // ids but we lose them when converting to vg.)
    unordered_map<string, nid_t> mg_to_vg;
    unordered_map<nid_t, string> vg_to_mg;
    std::tie(mg_to_vg, vg_to_mg) = load_trans(trans_path);

    if (progress) {
        cerr << "[filter-paf-deletions]: Loaded " << mg_to_vg.size() << " translations." << endl;
    }

    // open the paf
    ifstream paf_file(paf_path);
    if (!paf_file) {
        cerr << "[filter-paf-deletions] error: Unable to open PAF" << endl;
        return 1;
    }

    // load the paf into memory
    vector<PafLine> paf_lines;
    unordered_map<int64_t, int64_t> orig_to_sorted;
    std::tie(paf_lines, orig_to_sorted) = load_paf(paf_file);
    if (progress) {
        cerr << "[filter-paf-deletions]: Loaded " << paf_lines.size() << " paf lines" << endl;
    }

    // index the minigraph
    // this maps each node in the graph to a (maximal) reference interval
    unordered_map<nid_t, Anchor> ref_index = index_graph(graph.get(), ref_prefix);
    if (progress) {
        cerr << "[filter-paf-deletions]: Created reference path index" << endl;
    }

    unordered_map<path_handle_t, int64_t> ref_path_to_length;
    if ((max_deletion > 0 && max_deletion <= 1.) || (max_insertion > 0. && max_insertion <= 1.)) {
        graph->for_each_path_handle([&](path_handle_t path_handle) {
                int64_t len = 0;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        len += graph->get_length(graph->get_handle_of_step(step_handle));
                    });
                ref_path_to_length[path_handle] = len;
            });
        if (progress) {
            cerr << "[filter-paf-deletions]: Computed lengths for " << ref_path_to_length.size() << " reference paths" << endl;
        }
    }

    unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>> ref_deletions = index_deletions(graph.get(), ref_index);
    if (progress) {
        cerr << "[filter-paf-deletions]: Created reference deletion index" << endl;
    }

#ifdef debug
    for (auto fam : ref_index) {
        cerr << fam.first << " -> " << graph->get_path_name(fam.second.path_handle) << " " << fam.second.min_offset << " " << fam.second.max_offset << endl;
    }
#endif

    // we have everything needed to filter the paf
    vector<bool> filtered_lines(paf_lines.size(), false);
    int64_t filtered_line_total = 0;
    int64_t filtered_line_it = 0;
    int64_t filtered_match_total = 0;
    vector<bool> off_ref_filtered_lines(paf_lines.size(), false);
    int64_t off_ref_filtered_line_it = 0;
    int64_t off_ref_filtered_match_total = 0;
    int64_t iteration = 0;

    do {
        filtered_line_it = 0;
        for_each_query_block(paf_lines, filtered_lines, [&](int64_t block_start, int64_t block_end) {
                assert(!filtered_lines[block_start] && !filtered_lines[block_end]);
                // get some stats about the block
                unordered_map<path_handle_t, int64_t> ref_path_sizes;
                int64_t total_matches = 0;
                for (int64_t i = block_start; i <= block_end; ++i) {
                    if (!filtered_lines[i]) {
                        const PafLine& paf = paf_lines[i];
                        nid_t target_id = mg_to_vg.at(paf.target_name);
                        const Anchor& anchor = ref_index.at(target_id);
                        ref_path_sizes[anchor.path_handle] += paf.num_matching;
                        total_matches += paf.num_matching;
                    } else {
                        assert(iteration > 0);
                    }
                }
                if (total_matches == 0) {
                    // whole block was filtered, nothing to be done
                    return;
                }
                // find the number one reference path by match coverage
                // todo: what about tie?
                path_handle_t ref_path;
                int64_t ref_path_size = -1;
                for (const auto& rps : ref_path_sizes) {
                    if (rps.second > ref_path_size) {
                        ref_path_size = rps.second;
                        ref_path = rps.first;
                    }
                }

                // support fracitonal thresholds which apply to path length
                int64_t max_deletion_threshold = max_deletion;
                int64_t max_insertion_threshold = max_insertion;
                if (max_deletion > 0. && max_deletion <= 1.) {
                    max_deletion_threshold = max_deletion * ref_path_to_length.at(ref_path);
                }
                if (max_insertion > 0. && max_insertion <= 1.) {
                    max_insertion_threshold = max_insertion * ref_path_to_length.at(ref_path);
                }

                // mask out everything off this path
                // get rid of all off-reference path mappings right away
                // note: these are flagged to be ignored but not actually filtered
                // unless keep_off_ref is set to false, then they are removed
                int64_t off_ref_total = 0;
                int64_t off_ref_match_total = 0;
                for (int64_t i = block_start; i <= block_end; ++i) {
                    if (!filtered_lines[i]) {
                        nid_t cur_target_id = mg_to_vg.at(paf_lines[i].target_name);
                        const Anchor& cur_anchor = ref_index.at(cur_target_id);
                        if (cur_anchor.path_handle != ref_path) {
                            off_ref_filtered_lines[i] = true;
                            ++off_ref_total;
                            off_ref_match_total += paf_lines[i].num_matching;
                        }
                    }
                }
                off_ref_filtered_line_it += off_ref_total;
                off_ref_filtered_match_total += off_ref_match_total;
                if (!keep_off_ref && verbose && off_ref_total > 0) {
                    cerr << "[filter-paf-deletions]: filtered " << off_ref_total << " lines with " << off_ref_match_total << " bases "
                         << " because they did not map to reference sequence " << graph->get_path_name(ref_path) << " all in block "
                         << "\n  I=" << block_start <<": " << paf_lines[block_start]
                         << "\n  J=" << block_end << ": " << paf_lines[block_end] << endl << endl;
                }

                // try to find a gap that exceeds the length
                int64_t prev_idx = -1;
                // these are the boundaries of the deletions in the block
                // the deletion is between cut_point[i] and cut_point[i] - 1
                vector<int64_t> cut_points;
                vector<int64_t> cut_sizes;
                for (int64_t i = block_start; i <= block_end; ++i) {
                    if (filtered_lines[i] || off_ref_filtered_lines[i]) {
                        continue;
                    }
                    if (prev_idx == -1) {
                        prev_idx = i;
                        continue;
                    }
                    // if we got this far that means we're on the path and we have a prev on the path too
                    // do a rough delta check
                    assert(prev_idx < i);
                    const PafLine& cur_paf = paf_lines[i];
                    const PafLine& prev_paf = paf_lines[prev_idx];
                    PafDelta paf_delta = get_delta(ref_path, prev_paf, cur_paf, mg_to_vg, ref_index, ref_deletions);
                    
                    int64_t checked_delta = check_delta(max_deletion_threshold, max_insertion_threshold, paf_delta, overlap_threshold);

                    if (checked_delta != 0) {
                        if (verbose) {
                            cerr << "[filter-paf-deletions]: detected " << (checked_delta > 0 ? "deletion" : "insertion")
                                 << " of size " << (int64_t)abs(paf_delta.delta) << " with overlap " << paf_delta.ref_overlap_size
                                 << " on ref path " << graph->get_path_name(ref_path) << " with cur anchor ("
                                 << paf_delta.cur_ref_start << ", " << paf_delta.cur_ref_end << ") and prev anchor (" << paf_delta.prev_ref_start << ", "
                                 << paf_delta.prev_ref_end << ") and threshold " << max_deletion_threshold
                                 << " on following paf line:\n  I=" << (prev_idx) <<": " << prev_paf
                                 << "\n  J=" << i << ": " << cur_paf << endl << endl;
                        }
                        cut_points.push_back(i);
                        cut_sizes.push_back(paf_delta.delta);
                    } 
                
                    prev_idx = i;
                }

                // greedy heuristic: for every deletion cut point, we try scanning forward and backward to find the minimum
                // resolving block of lines whose removal solves the deletion
                for (int64_t j = 0; j < cut_points.size(); ++j) {
                    if (filtered_lines[cut_points[j]]) {
                        continue;
                    }
                    // go backward                        
                    int64_t backward_matches = 0;
                    int64_t backward_candidate = -1; // last *unfiltered* line scanning left
                    int64_t prev_idx = -1;
                    for (int64_t k = cut_points[j] - 1; k >= block_start && backward_candidate == -1; --k) {
                        if (!filtered_lines[k] && !off_ref_filtered_lines[k]) {
                            if (prev_idx == -1) {
                                prev_idx = k;
                            }
                            const PafLine& prev_paf = paf_lines[k];
                            const PafLine& cur_paf = paf_lines[cut_points[j]];
                            PafDelta paf_delta = get_delta(ref_path, prev_paf, cur_paf, mg_to_vg, ref_index, ref_deletions);
                            int64_t checked_delta = check_delta(max_deletion_threshold, max_insertion_threshold, paf_delta, overlap_threshold);
                            if (checked_delta == 0) {
                                backward_candidate = k;
                            } else {
                                backward_matches += paf_lines[k].num_matching;
                            }
                        }
                    }
                    if (backward_candidate == -1) {
                        // need to delete block_start too
                        backward_candidate = block_start - 1;        
                    }

                    // go forward
                    int64_t forward_matches = 0;
                    int64_t forward_candidate = -1; // last *unfiltered* line scanning right
                    for (int64_t k = cut_points[j]; k <= block_end && forward_candidate == -1 && prev_idx != -1; ++k) {
                        if (!filtered_lines[k] && !off_ref_filtered_lines[k]) {
                            const PafLine& prev_paf = paf_lines[prev_idx];
                            const PafLine& cur_paf = paf_lines[k];
                            int64_t max_deletion_threshold = max_deletion;
                            int64_t max_insertion_threshold = max_insertion;
                            PafDelta paf_delta = get_delta(ref_path, prev_paf, cur_paf, mg_to_vg, ref_index, ref_deletions);
                            int64_t checked_delta = check_delta(max_deletion_threshold, max_insertion_threshold, paf_delta, overlap_threshold);
                            if (checked_delta == 0) {
                                forward_candidate = k;
                            } else {
                                forward_matches += paf_lines[k].num_matching;
                            }
                        }
                    }
                    if (forward_candidate == -1) {
                        // need to delete block_end too
                        forward_candidate = block_end + 1;
                    }

                    assert(backward_candidate >= block_start - 1 && forward_candidate <= block_end + 1);

                    int64_t min_segment_start = -1;
                    int64_t min_segment_end = -1;
                    int64_t min_segment_matches = -1;
                    if (backward_matches < forward_matches) {
                        min_segment_start = backward_candidate + 1; // want to filter right of candidate (not include)
                        min_segment_end = cut_points[j] - 1;
                        min_segment_matches = backward_matches;
                    } else {
                        min_segment_start = cut_points[j];
                        min_segment_end = forward_candidate - 1; // want to filter left of candidate (not include)
                        min_segment_matches = forward_matches;                        
                    }

                    int64_t max_matches_deleted = filter_threshold * cut_sizes[j];
                    if (j == 0) {
                        if (min_segment_end < min_segment_start) {
                            cerr << "ms me " << min_segment_start << " " << min_segment_end << " bm fm " << backward_matches << " " << forward_matches
                                 << " bc fc " << backward_candidate << " " << forward_candidate << endl;
                        }
                        assert(min_segment_start <= min_segment_end);
                    }
                    if (min_segment_matches > 0) {
                  
                        if (min_segment_matches <= max_matches_deleted) {
                            int64_t lines_in_segment = 0;
                            for (int64_t k = min_segment_start; k <= min_segment_end; ++k) {
                                if (!filtered_lines[k]) {
                                    filtered_lines[k] = true;
                                    ++filtered_line_it;
                                    filtered_match_total += paf_lines[k].num_matching;
                                    ++lines_in_segment;
                                }
                            }

                            if (verbose) {
                                cerr << "[filter-paf-deletions]: filtering " << lines_in_segment << " PAF lines between (inclusively)\n  I="
                                     << min_segment_start << ": " << paf_lines[min_segment_start]
                                     << "\n  J=" << min_segment_end << ":  " << paf_lines[min_segment_end]
                                     << "\nfor a total of " << min_segment_matches << " matches" << endl << endl;
                            }
                        } else {
                            if (verbose) {
                                cerr << "[filter-paf-deletions]: leaving in PAF lines between (inclusively)\n  I="
                                     << min_segment_start << ": " << paf_lines[min_segment_start]
                                     << "\n  J=" << min_segment_end << ":  " << paf_lines[min_segment_end]
                                     << "\nfor a total of " << min_segment_matches << " matches, which exceeds deletion threshold of " << max_matches_deleted
                                     << endl << endl;
                            }

                        }
                    }
                }
            });
        if (!keep_off_ref) {
            filtered_line_it += off_ref_filtered_line_it;
        }
        
        if (progress) {
            cerr << "[filter-paf-deletions]: Iteration " << iteration << ": Found " << filtered_line_it << " lines to filter" << endl;
        }
        ++iteration;
        filtered_line_total += filtered_line_it;
        if (!keep_off_ref) {
            filtered_line_total += off_ref_filtered_line_it;
        }
    } while (filtered_line_it > 0);

    if (!keep_off_ref) {
        filtered_match_total += off_ref_filtered_match_total;
    }
    if (progress) {
        cerr << "[filter-paf-deletions]: Filtering out " << filtered_line_total << " paf lines totaling " << filtered_match_total << " matches" << endl;
    }

    // output the unfiltered lines
    paf_file.clear();
    paf_file.seekg(0, ios::beg) ;
    string buffer;
    for (int64_t line_no = 0; line_no < filtered_lines.size(); ++line_no) {
        int64_t sorted_line_no = orig_to_sorted.at(line_no);
        const auto& ret = getline(paf_file, buffer);
        assert(ret);

        // sanity check:
        PafLine paf_line = parse_paf_line(buffer);
        assert(paf_line.query_name == paf_lines[sorted_line_no].query_name);
        assert(paf_line.query_start == paf_lines[sorted_line_no].query_start);
        assert(paf_line.query_end == paf_lines[sorted_line_no].query_end);
        assert(paf_line.target_name == paf_lines[sorted_line_no].target_name);
        assert(paf_line.target_start == paf_lines[sorted_line_no].target_start);
        assert(paf_line.target_end == paf_lines[sorted_line_no].target_end);
                
        if (filtered_lines[sorted_line_no] == false && (keep_off_ref || off_ref_filtered_lines[sorted_line_no] == false)) {
            cout << buffer << "\n";
        }
    }
    cout << flush;

    return 0;
}

static string strip_prefix(const string& name) {
    if (name.compare(0, 3, "id=") == 0) {
        size_t p = name.find('|', 3);
        assert(p != string::npos);
        return name.substr(p + 1);
    }
    return name;
}

unordered_map<nid_t, Anchor> index_graph(const PathHandleGraph* graph,
                                         const string& ref_prefix) {

    // start by making a path position index
    // minigraph assumption: no more than one path per handle!
    unordered_map<handle_t, int64_t> position_index;
    graph->for_each_path_handle([&](path_handle_t path_handle) {            
            if (graph->get_path_name(path_handle).compare(0, ref_prefix.length(), ref_prefix) == 0) {
                size_t offset = 0;
                graph->for_each_step_in_path(path_handle, [&](step_handle_t step_handle) {
                        handle_t handle = graph->get_handle_of_step(step_handle);
                        size_t len = graph->get_length(handle);
                        assert(len > 0);
                        assert(!position_index.count(handle));
                        position_index[handle] = offset;
                        position_index[graph->flip(handle)] = offset + len - 1;
                        offset += len;
                    });
            }
        });

    if (position_index.empty()) {
        cerr << "[filter-paf-deletions] error: no reference path found" << endl;
        exit(0);
    }
    
    vector<unordered_map<nid_t, Anchor>> thread_results(get_thread_count());
    
    // really slow brute-force relies on minigraph not having too many nodes
    graph->for_each_handle([&](handle_t handle) {
            unordered_map<nid_t, Anchor>& ref_index = thread_results[omp_get_thread_num()];

            // find all reference nodes that are connected via BFS
            unordered_set<handle_t> context;
            unordered_set<handle_t> ref_handles;
            vector<handle_t> cur_handles = {handle};
            while (!cur_handles.empty()) {
                vector<handle_t> next_handles;
                for (auto& h : cur_handles) {
                    if (!context.count(h)) {
                        context.insert(h);
                        if (position_index.count(h)) {
                            // dead-end on reference
                            ref_handles.insert(h);
                        } else {
                            graph->follow_edges(h, false, [&](handle_t n) {
                                    next_handles.push_back(n);
                                });
                            graph->follow_edges(h, true, [&](handle_t p) {
                                    next_handles.push_back(p);
                                });
                        }
                    }
                }
                cur_handles = std::move(next_handles);
            }

            // update the index with reference offsets
            unordered_set<path_handle_t> ref_path_set;
            int64_t min_ref_offset = numeric_limits<int64_t>::max();
            int64_t max_ref_offset = -1;
            for (handle_t ref_handle : ref_handles) {
                vector<step_handle_t> steps = graph->steps_of_handle(ref_handle);
                assert(steps.size() == 1);
                path_handle_t ref_path_handle = graph->get_path_handle_of_step(steps.back());
                ref_path_set.insert(ref_path_handle);
                // assumption: only one reference path in component
                // (fair for minigraph, but may need to do better than prefix for path selection)
                assert(ref_path_set.size() == 1);
                int64_t ref_offset = position_index.at(ref_handle);
                int64_t ref_offset_rev = position_index.at(graph->flip(ref_handle));
                min_ref_offset = std::min(min_ref_offset, min(ref_offset, ref_offset_rev));
                max_ref_offset = std::max(max_ref_offset, max(ref_offset, ref_offset_rev));
                assert(max_ref_offset >= min_ref_offset);
            }
            if (!ref_path_set.empty()) {
                assert(ref_path_set.size() == 1);
                Anchor& anchor = ref_index[graph->get_id(handle)];
                anchor.path_handle = *ref_path_set.begin();
                anchor.min_offset = min_ref_offset;
                anchor.max_offset = max_ref_offset;
                assert(anchor.max_offset >= anchor.min_offset);
            }
        }, true);

    // merge up the indexes
    for (size_t i = 1; i < thread_results.size(); ++i) {
        for (const auto& id_anchor : thread_results[i]) {
            thread_results[0][id_anchor.first] = id_anchor.second;
        }
        thread_results[i].clear();
    }

    return thread_results[0];
}

pair<vector<PafLine>, unordered_map<int64_t, int64_t>> load_paf(ifstream& paf_file) {

    vector<pair<int64_t, PafLine>> numbered_lines;
    string buffer;
    int64_t line_no = 0;
    while (getline(paf_file, buffer)) {
        PafLine paf_line = parse_paf_line(buffer);
        // dont use this
        paf_line.cigar = "";
        numbered_lines.push_back(make_pair(line_no++, paf_line));
    }
    std::sort(numbered_lines.begin(), numbered_lines.end(), [&](const pair<int64_t, PafLine>& p1, const pair<int64_t, PafLine>& p2) {
            return p1.second.query_name < p2.second.query_name ||
                                          (p1.second.query_name == p2.second.query_name && p1.second.query_start < p2.second.query_start);
        });

    vector<PafLine> paf_lines;
    unordered_map<int64_t, int64_t> orig_to_sorted;
    for (int64_t i = 0; i < numbered_lines.size(); ++i) {
        paf_lines.push_back(numbered_lines[i].second);
        orig_to_sorted[numbered_lines[i].first] = i;
    }
    
    return make_pair(paf_lines, orig_to_sorted);
}

int64_t for_each_query_block(const vector<PafLine>& paf_lines, const vector<bool>& filtered_lines,
                             function<void(int64_t, int64_t)> visit_block) {
    if (paf_lines.empty()) {
        assert(false);
    }
    int64_t block_start = -1;
    int64_t block_end = -1;
    string prev_query;
    int64_t num_visits = 0;
    for (int64_t i = 0; i < paf_lines.size(); ++i) {
        if (filtered_lines[i]) {
            continue;
        }
        const PafLine& paf = paf_lines[i];
        if (block_start == -1) {
            block_start = i;
        } else if (paf.query_name != prev_query) {
            assert(!prev_query.empty());
            if (block_start > -1) {
                // visit the previous block
                visit_block(block_start, block_end);
            }
            ++num_visits;
            //start a new block
            block_start = i;
        }
        // update end of current block
        block_end = i;
        prev_query = paf.query_name;
    }

    if (block_end != -1) {
        // visit last block if present
        visit_block(block_start, block_end);
        ++num_visits;
    }
    return num_visits;
}

unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>> index_deletions(const PathHandleGraph* graph, const unordered_map<nid_t, Anchor>& index) {

    vector<unordered_map<path_handle_t, vector<Interval<int64_t, int64_t>>>> thread_deletions(get_thread_count());

    // get approximate deletion intervals using the index
    graph->for_each_edge([&](edge_t edge) {
            const Anchor& a1 = index.at(graph->get_id(edge.first));
            const Anchor& a2 = index.at(graph->get_id(edge.second));
            if (a1.path_handle == a2.path_handle) {
                Interval<int64_t, int64_t> interval(0, 0, 0);
                if (a1.min_offset < a2.min_offset) {
                    interval.start = a1.max_offset;
                    interval.stop = a2.min_offset;
                } else {
                    interval.start = a2.max_offset;
                    interval.stop = a1.min_offset;                
                }
                interval.value = interval.stop - interval.start;
                if (interval.value > 1) {
                    thread_deletions[omp_get_thread_num()][a1.path_handle].push_back(interval);
                }
            }
        }, true);

    for (size_t i = 1; i < thread_deletions.size(); ++i) {
        for (const auto& pi : thread_deletions[i]) {
            for (const auto& interval : pi.second) {
                thread_deletions[0][pi.first].push_back(interval);
            }
        }
        thread_deletions[i].clear();
    }

    unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>> path_to_tree;
    for (const auto& pi : thread_deletions[0]) {
        path_to_tree[pi.first] = IntervalTree<int64_t, int64_t>(pi.second);
    }
    return path_to_tree;
}


pair<unordered_map<string, nid_t>, unordered_map<nid_t, string>> load_trans(const string& trans_path) {
    ifstream trans_file(trans_path);
    if (!trans_file) {
        cerr << "[filter-paf-deletions] error: Unable to load trans file" << endl;
        exit(1);
    }

    unordered_map<string, nid_t> mg_to_vg;
    unordered_map<nid_t, string> vg_to_mg;

    string buffer;
    while (getline(trans_file, buffer)) {
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);
        assert(toks.size() == 3 && toks[0] == "T");
        string& mg_name = toks[1];
        bool has_prefix = mg_name.compare(0, 3, "id=") == 0;
        mg_to_vg[mg_name] = stol(toks[2]);
        // hack to support prefixed or not minigraph
        // just by keeping both versions in the map no matter what
        // todo: parameterize prefix name
        if (has_prefix) {
            mg_to_vg[strip_prefix(mg_name)] = stol(toks[2]);
        } else {
            mg_to_vg["id=_MINIGRAPH_|" + mg_name] = stol(toks[2]);
        }
        vg_to_mg[stol(toks[2])] = mg_name;
    }

    return make_pair(mg_to_vg, vg_to_mg);    
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

// return -Delta for insertion, +Delta for deletion and 0 if it doesn't pass thresholds
int64_t check_delta(int64_t max_deletion_threshold, int64_t max_insertion_threshold, const PafDelta& paf_delta, double overlap_threshold) {
    int64_t ret = 0;
    if (paf_delta.delta != 0) {
        // note: paf_delta.delta is ref-query, so deletions are positive and insertions are negative
        if (paf_delta.delta < 0 && max_insertion_threshold > 0 && -paf_delta.delta > max_insertion_threshold) {
            ret = paf_delta.delta;
        } else if (paf_delta.delta > 0 && max_deletion_threshold > 0 && paf_delta.delta > max_deletion_threshold  &&
                  abs((double)paf_delta.ref_overlap_size / paf_delta.delta) < overlap_threshold) {
            ret = paf_delta.delta;
        }
    }
    return ret;
}

PafDelta get_delta(path_handle_t ref_path, const PafLine& prev_paf, const PafLine& cur_paf,
                   const unordered_map<string, nid_t>& mg_to_vg, const unordered_map<nid_t, Anchor>& ref_index,
                   const unordered_map<path_handle_t, IntervalTree<int64_t, int64_t>>& ref_deletions) {
    
    PafDelta paf_delta;
        
    paf_delta.query_delta = cur_paf.query_start - prev_paf.query_end; // not abs because sorted
    
    nid_t prev_target_id = mg_to_vg.at(prev_paf.target_name);
    const Anchor& prev_anchor = ref_index.at(prev_target_id);

    nid_t cur_target_id = mg_to_vg.at(cur_paf.target_name);
    const Anchor& cur_anchor = ref_index.at(cur_target_id);
                    
    // todo : verify 
    paf_delta.cur_ref_start = cur_anchor.min_offset + cur_paf.target_start;
    paf_delta.cur_ref_end = cur_anchor.max_offset - (cur_paf.target_len - cur_paf.target_end);
    paf_delta.prev_ref_start = prev_anchor.min_offset + prev_paf.target_start;
    paf_delta.prev_ref_end = prev_anchor.max_offset - (prev_paf.target_len - prev_paf.target_end);

    int64_t cur_ref_start = paf_delta.cur_ref_start;
    int64_t cur_ref_end = paf_delta.cur_ref_end;
    int64_t prev_ref_start = paf_delta.prev_ref_start;
    int64_t prev_ref_end = paf_delta.prev_ref_end;
        
    // sort the ref intervals
    if (cur_ref_start < prev_ref_start) {
        swap(cur_ref_start, prev_ref_start);
        swap(cur_ref_end, prev_ref_end);
    }
    paf_delta.ref_delta = cur_ref_start - prev_ref_end;
                                    
    paf_delta.delta = paf_delta.ref_delta > 0 ? paf_delta.ref_delta - paf_delta.query_delta : 0;

    paf_delta.ref_overlap_size = 0;
    if (paf_delta.delta > 0) {                     
        if (ref_deletions.count(ref_path)) {
            vector<Interval<int64_t, int64_t>> overlaps = ref_deletions.at(ref_path).findOverlapping(prev_ref_end, cur_ref_start);
            for (const auto& overlap : overlaps) {
                int64_t intersection_start = max(prev_ref_end, overlap.start);
                int64_t intersection_stop = min(cur_ref_start, overlap.stop);
                paf_delta.ref_overlap_size = max(paf_delta.ref_overlap_size, intersection_stop - intersection_start);                                
            }
        }
    }

    return paf_delta;
}
