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

using namespace std;
using namespace handlegraph;
using namespace bdsg;

struct Anchor {
    path_handle_t path_handle;
    int64_t max_offset;
    int64_t min_offset;
};

static unique_ptr<MutablePathMutableHandleGraph> load_graph(istream& graph_stream);
static pair<unordered_map<string, nid_t>, unordered_map<nid_t, string>> load_trans(const string& trans_path);
static vector<string> &split_delims(const string &s, const string& delims, vector<string> &elems);
static unordered_map<nid_t, Anchor> index_graph(const PathHandleGraph* graph,
                                                const string& ref_prefix);
static IntervalTree<int64_t, int64_t> index_deletions(const PathHandleGraph* graph, const unordered_map<nid_t, Anchor>& index);
static int64_t filter_paf(const PathHandleGraph* graph, ifstream& paf_file, const unordered_map<nid_t, Anchor>& ref_index,
                          const IntervalTree<int64_t, int64_t>& ref_deletions,
                          const unordered_map<string, nid_t>& mg_to_vg, int64_t max_deletion,
                          double overlap_threshold, vector<bool>& filtered_lines, bool progress, bool verbose);

void help(char** argv) {
    cerr << "usage: " << argv[0] << " [options] <graph.vg> <trans> <aln.paf> <reference-prefix> <threshold>\n" << endl
         << "Use distances from graph to filter out implied deletions from PAF (cigars not considered, only blocks)" << endl
         << "  <graph.vg> : minigraph as obtained from vg convert -g graph.gfa" << endl
         << "  <trans> : node translation from vg convert -g -T" << endl
         << "  <aln.paf> : paf alignment from cactus-graphmap" << endl
         << "  <reference-prefix> : prefix of reference path(s) in graph" 
         << endl
         << "options: " << endl
         << "    -p, --progress            Print progress" << endl
         << "    -v, --verbose             Print deletions" << endl
         << "    -t, --threads N           number of threads to use (used only for computing snarls) [default: all available]" << endl
       << endl;
}    

int main(int argc, char** argv) {

    bool progress = false;
    bool verbose = false;
    // only filter deletions that don't overlap an existing deletion by at least this much
    // (doesn't seem to a factor -- most big deletions not in minigraph)
    double overlap_threshold = 0.5;
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"progress", no_argument, 0, 'p'},
            {"verbose", no_argument, 0, 'v'},
            {"threads", required_argument, 0, 't'},
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hpvt:",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
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

    if (argc <= 5) {
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

    if (optind != argc - 5) {
        cerr << "[filter-paf-deletions] error: too many arguments\n" << endl;
        help(argv);
        return 1;
    }
    
    string graph_path = argv[optind++];
    string trans_path = argv[optind++];
    string paf_path = argv[optind++];
    string ref_prefix = argv[optind++];
    int64_t max_deletion = stol(argv[optind++]);

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

    // index the minigraph
    // this maps each node in the graph to a (maximal) reference interval
    unordered_map<nid_t, Anchor> ref_index = index_graph(graph.get(), ref_prefix);
    if (progress) {
        cerr << "[filter-paf-deletions]: Created reference path index" << endl;
    }

    IntervalTree<int64_t, int64_t> ref_deletions = index_deletions(graph.get(), ref_index);
    if (progress) {
        cerr << "[filter-paf-deletions]: Created reference deletion index" << endl;
    }

#ifdef debug
    for (auto fam : ref_index) {
        cerr << fam.first << " -> " << graph->get_path_name(fam.second.path_handle) << " " << fam.second.min_offset << " " << fam.second.max_offset << endl;
    }
#endif

    // we have everything needed to filter the paf
    vector<bool> filtered_lines;
    int64_t filtered_line_total = 0;
    int64_t filtered_line_it = 0;
    int64_t iteration = 0;
    do {
        filtered_line_it = filter_paf(graph.get(), paf_file, ref_index, ref_deletions, mg_to_vg, max_deletion,
                                      overlap_threshold, filtered_lines, progress, verbose);
        filtered_line_total += filtered_line_it;
        if (progress) {
            cerr << "[filter-paf-deletions]: Iteration " << iteration << ": Found " << filtered_line_it << " lines to filter" << endl;
        }
        ++iteration;
    } while (filtered_line_it > 0);

    if (progress) {
        cerr << "[filter-paf-deletions]: Filtering out " << filtered_line_total << " paf lines" << endl;
    }

    // output the unfiltered lines
    paf_file.clear();
    paf_file.seekg(0, ios::beg) ;
    string buffer;
    for (int64_t line_no = 0; line_no < filtered_lines.size(); ++line_no) {
        const auto& ret = getline(paf_file, buffer);
        assert(ret);
        if (filtered_lines[line_no] == false) {
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

int64_t filter_paf(const PathHandleGraph* graph, ifstream& paf_file, const unordered_map<nid_t, Anchor>& ref_index,
                   const IntervalTree<int64_t, int64_t>& ref_deletions,
                   const unordered_map<string, nid_t>& mg_to_vg, int64_t max_deletion,
                   double overlap_threshold, vector<bool>& filtered_lines, bool progress, bool verbose) {

    paf_file.clear();
    paf_file.seekg(0, ios::beg) ;
    string buffer;

    unordered_set<string> query_set;
    string prev_query_name;
    int64_t prev_query_start;
    int64_t prev_query_end;
    nid_t prev_target_id = 0;
    int64_t prev_ref_start;
    int64_t prev_ref_end;

    // store matches in current block (run of mappings from same query that don't span a max_deletion)
    vector<int64_t> line_to_block; // map paf line number to a block (this vector contains offsets to below vector)
    vector<int64_t> blocks;
    unordered_set<int64_t> cut_blocks; // offsets in above that mark block before max_deletion

    size_t line_no = 0;
    size_t unfiltered_line_no = 0;
    while (getline(paf_file, buffer)) {
        if (unfiltered_line_no >= filtered_lines.size()) {
            // first pass
            filtered_lines.push_back(false);
        } else if (filtered_lines[unfiltered_line_no]) {
            ++unfiltered_line_no;
            continue;
        }
        vector<string> toks;
        split_delims(buffer, "\t\n", toks);

        // pretty silly not to have a paf parser in this repo
        const string& query_name = toks[0];
        int64_t query_start = stol(toks[2]);
        int64_t query_end = stol(toks[3]);
        string target_name = strip_prefix(toks[5]);
        int64_t target_length = stol(toks[6]);
        int64_t target_start = stol(toks[7]);
        int64_t target_end = stol(toks[8]);
        int64_t matches = stol(toks[9]);

        // map the target interval onto the reference path
        if (!mg_to_vg.count(target_name)) {
            cerr << "[filter-paf-deletions] error: target name from PAF not found in translation: " << target_name << endl;
            exit(1);
        }
        nid_t target_id = mg_to_vg.at(target_name);
        const Anchor& anchor = ref_index.at(target_id);
        int64_t ref_start = anchor.min_offset + target_start;
        int64_t ref_end = anchor.max_offset - (target_length - target_end);

#ifdef debug
        cerr << "anchor " << target_id << " = " << anchor.min_offset << " , " << anchor.max_offset << endl;
        cerr << "ref start = " << anchor.min_offset  << " + " << target_start << "   ref_end = " << anchor.max_offset << " - ("
             << target_length << " - " << target_end << ");" << endl;
#endif
        bool continues_block = query_name == prev_query_name;
        if (continues_block) {

            if (prev_query_start > query_start) {
                cerr << "Input paf not in order.  Sort with \"sort -k 1,1 -k 3,3n\" first.  Offending line:\n" << buffer << endl;
                exit(1);
            }
            int64_t query_delta = query_start - prev_query_end;
            int64_t ref_delta = ref_start - prev_ref_end;
            int64_t delta = ref_delta - query_delta;

            int64_t ref_overlap_size = 0;
            if (delta > 1) {
                vector<Interval<int64_t, int64_t>> overlaps = ref_deletions.findOverlapping(prev_query_end, query_start);
                for (const auto& overlap : overlaps) {
                    int64_t intersection_start = max(prev_query_end, overlap.start);
                    int64_t intersection_stop = min(query_start, overlap.stop);
                    ref_overlap_size = max(ref_overlap_size, intersection_stop - intersection_start);
                }
            }

            if (delta > max_deletion && ((double)ref_overlap_size / delta < overlap_threshold)) {
                if (verbose) {
                    cerr << "[filter-paf-deletions]: detected deletion of size " << delta << " with overlap " << ref_overlap_size <<
                        " on following paf line:"
                         << " \n  " << buffer << endl;
                }
                continues_block = false;
                cut_blocks.insert(blocks.size() - 1);
            }
        }

        // update the current block may adding the matches (or assigng a new block if necessary)
        if (!continues_block) {
            blocks.push_back(0);
        }
        assert(!blocks.empty());
        blocks.back() += matches;
        line_to_block.push_back(blocks.size() - 1);

        prev_query_name = query_name;
        prev_query_start = query_start;
        prev_query_end = query_end;
        prev_target_id = target_id;
        prev_ref_start = ref_start;
        prev_ref_end = ref_end;
        
        ++line_no;
        ++unfiltered_line_no;
    }

    // second pass to do filtering now that we have block sizes
    paf_file.clear();
    paf_file.seekg(0, ios::beg) ;

    size_t filtered_line_count = 0;
        
    line_no = 0;
    unfiltered_line_no = 0;
    while (getline(paf_file, buffer)) {
        if (filtered_lines[unfiltered_line_no]) {
            ++unfiltered_line_no;
            continue;
        }

        bool filter_out = false;        
        int64_t cur_block = line_to_block[line_no];
        if (cur_block > 0) {
            int64_t prev_block = cur_block - 1;
            if (cut_blocks.count(prev_block) && blocks[prev_block] >= blocks[cur_block]) {
                filter_out = true;
            }
        }
        if (!filter_out && cur_block < blocks.size() - 1) {
            int64_t next_block = cur_block + 1;
            if (cut_blocks.count(cur_block) && blocks[next_block] > blocks[cur_block]) {
                filter_out =true;
            }
        }

        if (filter_out) {
            filtered_lines[unfiltered_line_no] = true;
            ++filtered_line_count;
        }
        ++unfiltered_line_no;
        ++line_no;
    }

    return filtered_line_count;
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
    
IntervalTree<int64_t, int64_t> index_deletions(const PathHandleGraph* graph, const unordered_map<nid_t, Anchor>& index) {

    vector<vector<Interval<int64_t, int64_t>>> thread_deletions(get_thread_count());

    // get approximate deletion intervals using the index
    graph->for_each_edge([&](edge_t edge) {
            const Anchor& a1 = index.at(graph->get_id(edge.first));
            const Anchor& a2 = index.at(graph->get_id(edge.second));
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
                thread_deletions[omp_get_thread_num()].push_back(interval);
            }            
        }, true);

    for (size_t i = 1; i < thread_deletions.size(); ++i) {
        for (const auto& interval : thread_deletions[i]) {
            thread_deletions[0].push_back(interval);
        }
        thread_deletions[i].clear();
    }
    IntervalTree<int64_t, int64_t> tree(thread_deletions[0]);
    return tree;
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
        mg_to_vg[toks[1]] = stol(toks[2]);
        vg_to_mg[stol(toks[2])] = toks[1];
    }

    return make_pair(mg_to_vg, vg_to_mg);    
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

