// Count the number of bases that aren't in a given reference sample.
// Print the table of results stratisfied by number of covering samples
// Assume's current cactus convertion of Sample.Haplotype.Contig

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

using namespace std;
using namespace handlegraph;
using namespace bdsg;

static unique_ptr<PathHandleGraph> load_graph(istream& graph_stream) {

    char magic_bytes[4];
    graph_stream.read(magic_bytes, 4);
    uint32_t magic_number = ntohl(*((uint32_t*) magic_bytes));
    graph_stream.clear();
    graph_stream.seekg(0, ios::beg);

    PathHandleGraph* graph;
    if (magic_number == PackedGraph().get_magic_number()) {
        graph = new PackedGraph();
    } else if (magic_number == HashGraph().get_magic_number()) {
        graph = new HashGraph();
    } else {
        cerr << "Unable to parse input graph with magic number " << magic_number << endl;
        exit(1);
    }
    dynamic_cast<SerializableHandleGraph*>(graph)->deserialize(graph_stream);

    return unique_ptr<PathHandleGraph>(graph);
}

void help(char** argv) {
  cerr << "usage: " << argv[0] << " [options] <graph> [graph] [graph] [...]" << endl
       << "Count nodes and bp in graph covered by different sample counts\n" 
       << "Assumes SAMPLE.HAPLOTYPE.CONTIG path name format" << endl
       << endl
       << "options: " << endl
       << "    -r, --reference           Include counts of nodes that are not present in the given reference sample prefix" << endl
       << "    -i, --ignore              Completely ignore all paths with given prefix [default: _MINIGRAPH_]" << endl
       << "    -t, --threads             Number of threads [default: all]" << endl
       << "    -s, --separator           Use this separator for tokenizing path name. Haplotype key will be first 2 tokens (or all tokens if fewer than 2) [default=.]" << endl
       << "    -p, --progress            Print progress" << endl 
       << endl;
}    

// returns SAMPLE.HAPLOTYPE
// todo: vg/bdsg in progress of adpoting conventions / api
// to manage stuff like this -- should switch to using that
const string& get_sample_name(const PathHandleGraph* graph, path_handle_t path_handle,
                              unordered_map<path_handle_t, string>& name_map,
                              char separator) {
    if (!name_map.count(path_handle)) {
        string path_name = graph->get_path_name(path_handle);
        string sample;
        int dots = 0;
        for (int64_t i = 0; i < path_name.length(); ++i) {
            if (path_name[i] == separator) {
                ++dots;
            }
            if (dots == 2) {
                break;
            }
            sample.push_back(path_name[i]);
        }
        name_map[path_handle] = sample;        
    }
    return name_map.at(path_handle);
}

int main(int argc, char** argv) {

    string ref_sample;
    string ignore_sample = "_MINIGRAPH_";
    char separator = '.';
    bool progress = false;
    
    int c;
    optind = 1; 
    while (true) {

        static const struct option long_options[] = {
            {"help", no_argument, 0, 'h'},
            {"ref-sample", required_argument, 0, 'r'},
            {"ignore", required_argument, 0, 'i'},
            {"separator", required_argument, 0, 's'},
            {"threads", required_argument, 0, 't'},
            {"progress", no_argument, 0, 'p'},            
            {0, 0, 0, 0}
        };

        int option_index = 0;

        c = getopt_long (argc, argv, "hr:s:i:t:p",
                         long_options, &option_index);

        // Detect the end of the options.
        if (c == -1)
            break;

        switch (c)
        {
        case 'r':
            ref_sample = optarg;
            break;
        case 'i':
            ignore_sample = optarg;
            break;
        case 's':
            assert(strlen(optarg) == 1);
            separator = optarg[0];
            break;
        case 't':
            {
                int num_threads = stoi(optarg);
                if (num_threads <= 0) {
                    cerr << "error:[count-vg-hap-depth] Thread count (-t) set to " << num_threads << ", must set to a positive integer." << endl;
                    exit(1);
                }
                omp_set_num_threads(num_threads);
                break;                
            }
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
        cerr << "[count-vg-hap-depth] error: too few arguments" << endl;
        help(argv);
        return 1;
    }

    // depth stats (one per thread)
    vector<vector<int64_t>> depth_base_counts(get_thread_count());
    vector<vector<int64_t>> depth_nfree_base_counts(get_thread_count());
    vector<vector<int64_t>> depth_node_counts(get_thread_count());
    vector<vector<int64_t>> depth_base_counts_nonref(get_thread_count());
    vector<vector<int64_t>> depth_nfree_base_counts_nonref(get_thread_count());    
    vector<vector<int64_t>> depth_node_counts_nonref(get_thread_count());    

    // do counts for each graph arg
    while(optind < argc) {

        string graph_path = argv[optind++];
        ifstream graph_stream(graph_path);
        if (!graph_stream) {
            cerr << "[count-vg-hap-depth] error: Unable to open input graph " << graph_path << endl;
            return 1;
        }    
        unique_ptr<PathHandleGraph> graph = load_graph(graph_stream);
        graph_stream.close();
        if (progress) {
            cerr << "[count-vg-hap-depth]: Loaded graph" << endl;
        }

        // path handle to sample key (one per thread)
        vector<unordered_map<path_handle_t, string>> name_maps(get_thread_count());
    
        if (progress) {
            cerr << "[count-vg-hap-depth]: Calculating coverage with " << depth_base_counts.size() << " threads" << endl;
        }

        graph->for_each_handle([&](handle_t handle) {
                int64_t t = omp_get_thread_num();
                // collect all the samples that step on the node
                set<string> sample_set;
                bool ref = false;
                graph->for_each_step_on_handle(handle, [&](step_handle_t step_handle) {
                        const string& sample_name = get_sample_name(graph.get(), graph->get_path_handle_of_step(step_handle), name_maps[t], separator);
                        if (ignore_sample.empty() || sample_name.compare(0, ignore_sample.length(), ignore_sample) != 0) {
                            if (!ref && sample_name.compare(0, ref_sample.length(), ref_sample) == 0) {
                                ref = true;
                            }
                            sample_set.insert(sample_name);
                        }
                    });
                // update the total coverage
                int64_t coverage = sample_set.size();
                if (depth_base_counts[t].size() <= coverage) {
                    depth_base_counts[t].resize(coverage + 1, 0);
                    depth_node_counts[t].resize(coverage + 1, 0);
                    depth_nfree_base_counts[t].resize(coverage + 1, 0);
                }
                int64_t node_len = graph->get_length(handle);
                int64_t num_ns = 0;
                string node_seq = graph->get_sequence(handle);
                for (auto c : node_seq) {
                    if (c == 'N' || c == 'n') {
                        ++num_ns;
                    }
                }
                depth_base_counts[t][coverage] += node_len;
                depth_nfree_base_counts[t][coverage] += node_len - num_ns;
                depth_node_counts[t][coverage] += 1;
                            
                if (!ref && !ref_sample.empty()) {
                    // update the nonref coverage 
                    int64_t coverage = sample_set.size();
                    if (depth_base_counts_nonref[t].size() <= coverage) {
                        depth_base_counts_nonref[t].resize(coverage + 1, 0);
                        depth_node_counts_nonref[t].resize(coverage + 1, 0);
                        depth_nfree_base_counts_nonref[t].resize(coverage + 1, 0);
                    }
                    depth_base_counts_nonref[t][coverage] += node_len;
                    depth_nfree_base_counts_nonref[t][coverage] += node_len - num_ns;                    
                    depth_node_counts_nonref[t][coverage] += 1;
                }
            },
            true);
    }

    // make sure all tables have same size
    size_t max_size = 0;
    for (int64_t t = 0; t < get_thread_count(); ++t) {
        max_size = std::max(max_size, depth_base_counts[t].size());
        max_size = std::max(max_size, depth_base_counts_nonref[t].size());
    }
    for (int64_t t = 0; t < get_thread_count(); ++t) {
        if (depth_base_counts[t].size() < max_size) {
            depth_base_counts[t].resize(max_size, 0);
            depth_nfree_base_counts[t].resize(max_size, 0);
            depth_node_counts[t].resize(max_size, 0);
        }
        if (depth_base_counts_nonref[t].size() < max_size) {
            depth_base_counts_nonref[t].resize(max_size, 0);
            depth_nfree_base_counts_nonref[t].resize(max_size, 0);
            depth_node_counts_nonref[t].resize(max_size, 0);
        }
        assert(depth_base_counts[t].size() == max_size);
        assert(depth_nfree_base_counts[t].size() == max_size);
        assert(depth_node_counts[t].size() == max_size);
        assert(depth_base_counts_nonref[t].size() == max_size);
        assert(depth_nfree_base_counts_nonref[t].size() == max_size);
        assert(depth_node_counts_nonref[t].size() == max_size);
    }
    
    if (progress) {
        cerr << "[count-vg-hap-depth]: Merging data from different threads" << endl;
    }
    
    // merge up the threads
    for (int64_t t = 1; t < get_thread_count(); ++t) {
        for (int64_t coverage = 0; coverage < depth_base_counts[t].size(); ++coverage) {
            assert(depth_base_counts[0].size() > coverage);
            depth_base_counts[0][coverage] += depth_base_counts[t][coverage];
            depth_nfree_base_counts[0][coverage] += depth_nfree_base_counts[t][coverage];
            depth_node_counts[0][coverage] += depth_node_counts[t][coverage];

            if (!ref_sample.empty()) {
                assert(depth_base_counts_nonref[0].size() > coverage);
                depth_base_counts_nonref[0][coverage] += depth_base_counts_nonref[t][coverage];
                depth_nfree_base_counts_nonref[0][coverage] += depth_nfree_base_counts_nonref[t][coverage];
                depth_node_counts_nonref[0][coverage] += depth_node_counts_nonref[t][coverage];                
            }
        }
    }

    // there's almost certainly an stl one-line for this.. oh well
    function<vector<int64_t>(const vector<int64_t>&)> get_cumul = [](const vector<int64_t>& v) {
        int64_t tot = 0;
        vector<int64_t> cumul(v.size(), 0);
        for (int64_t i = 0; i < v.size(); ++i) {
            tot += v[i];
            cumul[i] = tot;
        }
        return cumul;
    };
    function<vector<int64_t>(const vector<int64_t>&)> get_lumuc = [](const vector<int64_t>& v) {
        int64_t tot = 0;
        vector<int64_t> cumul(v.size(), 0);
        for (int64_t i = v.size() - 1; i >= 0; --i) {
            tot += v[i];
            cumul[i] = tot;
        }
        return cumul;
    };
    
    // keep cumulative counts while we're at it
    // cumulate from 0
    vector<int64_t> node_counts_cumul = get_cumul(depth_node_counts[0]);
    vector<int64_t> base_counts_cumul = get_cumul(depth_base_counts[0]);
    vector<int64_t> nfree_base_counts_cumul = get_cumul(depth_nfree_base_counts[0]);    
    vector<int64_t> node_counts_nonref_cumul = get_cumul(depth_node_counts_nonref[0]);
    vector<int64_t> base_counts_nonref_cumul = get_cumul(depth_base_counts_nonref[0]);
    vector<int64_t> nfree_base_counts_nonref_cumul = get_cumul(depth_nfree_base_counts_nonref[0]);

    //cumulate from end
    vector<int64_t> node_counts_lumuc = get_lumuc(depth_node_counts[0]);
    vector<int64_t> base_counts_lumuc = get_lumuc(depth_base_counts[0]);
    vector<int64_t> nfree_base_counts_lumuc = get_lumuc(depth_nfree_base_counts[0]);
    vector<int64_t> node_counts_nonref_lumuc = get_lumuc(depth_node_counts_nonref[0]);
    vector<int64_t> base_counts_nonref_lumuc = get_lumuc(depth_base_counts_nonref[0]);
    vector<int64_t> nfree_base_counts_nonref_lumuc = get_lumuc(depth_nfree_base_counts_nonref[0]);

    // print the results
    cout << "hap-depth"
         << "\t" << "nodes" << "\t" << "bases" << "\t" << "non-n-bases"
         << "\t" << "nodes-cumul" << "\t" <<"bases-cumul" << "\t" << "non-n-bases-cumul"
         << "\t" << "nodes-cumul-rev" << "\t" << "bases-cumul-rev" << "\t" << "non-n-bases-cumul-rev";
    if (!ref_sample.empty()) {
        cout << "\t" << "nodes-nonref" << "\t" << "bases-nonref" << "\t" << "non-n-bases-nonref"
             << "\t" << "nodes-cumul-nonref" << "\t" << "bases-cumul-nonref" << "\t" << "non-n-bases-cumul-nonref"
             << "\t" << "nodes-cumul-rev-nonref" << "\t" << "bases-cumul-rev-nonref" << "\t" << "non-n-bases-cumul-rev-nonref";
    }
    cout << endl;

    for (int64_t coverage = 0; coverage < depth_base_counts[0].size(); ++coverage) {        
        cout << coverage
             << "\t" << depth_node_counts[0][coverage] << "\t" << depth_base_counts[0][coverage] << "\t" << depth_nfree_base_counts[0][coverage]
             << "\t" << node_counts_cumul[coverage] << "\t" << base_counts_cumul[coverage] << "\t" << nfree_base_counts_cumul[coverage]
             << "\t" << node_counts_lumuc[coverage] << "\t" << base_counts_lumuc[coverage] << "\t" << nfree_base_counts_lumuc[coverage];
        if (!ref_sample.empty()) {
            cout << "\t" << depth_node_counts_nonref[0][coverage] << "\t" << depth_base_counts_nonref[0][coverage] << "\t" << depth_nfree_base_counts_nonref[0][coverage]
                 << "\t" << node_counts_nonref_cumul[coverage] << "\t" << base_counts_nonref_cumul[coverage] << "\t" << nfree_base_counts_nonref_cumul[coverage]
                 << "\t" << node_counts_nonref_lumuc[coverage] << "\t" << base_counts_nonref_lumuc[coverage] << "\t" << nfree_base_counts_nonref_lumuc[coverage];
        }
        cout << "\n";
    }
    
    return 0;
}
