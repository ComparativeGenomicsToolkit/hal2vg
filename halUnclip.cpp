/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// Convert clipped sequences (like chr1_sub_110000_22220000) back to their original states

//#define debug

#include <cstdlib>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>
#include <unordered_set>

#include "hal.h"
#include "commonC.h"
#include "bioioC.h"
#include "subpaths.h"

using namespace std;
using namespace hal;

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("inFile", "input HAL file");
    optionsParser->addArgument("seqFile", "cactus-style seqfile. 1st col=genome name, 2nd col=(original) fasta file.  only local paths supported");
    optionsParser->addArgument("outFile", "output HAL file");
    optionsParser->addOptionFlag("progress",
                                 "show progress",
                                 false);
    optionsParser->setDescription("Fill back clipped sequence (removed by cactus-preprocess) using the original fasta files"
        ". Star trees only");
}

static vector<string> split_delims(const string &s, const string& delims) {
    vector<string> elems;
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

// c++ wrapper for sonlib -- load fasta file into memory
static unordered_map<string, pair<size_t, string>> read_fasta(const string& fa_path) {
    FILE* fa_file = fopen(fa_path.c_str(), "r");
    if (!fa_file) {
        cerr << "Unable to open fastat file: " << fa_path << endl;
        exit(1);
    }

    List* seqs = constructEmptyList(0, free);
    List* seq_lens = constructEmptyList(0, free);
    List* seq_names = constructEmptyList(0, free);

    fastaRead(fa_file, seqs, seq_lens, seq_names);

    unordered_map<string, pair<size_t, string>> fa_info; 
    for (int64_t i = 0; i < seqs->length; ++i) {
        string name = (char*)seq_names->list[i];
        size_t len = (size_t)listGetInt(seq_lens, i);
        string seq = (char*)seqs->list[i];
        fa_info[name] = make_pair(len, seq);
    }

    destructList(seqs);
    destructList(seq_lens);
    destructList(seq_names);

    return fa_info;
}

// do a pass over the seqfile to get the total lengths of every sequence
static unordered_map<string, size_t> get_dimensions_from_seqfile(const string& seqfile_path) {
    unordered_map<string, size_t> seq_map;
    
    ifstream seqfile(seqfile_path);
    if (!seqfile) {
        cerr << "[halUnclip]: Unable to open seqfile: " << seqfile_path << endl;
        exit(1);
    }

    string buffer;
    while (getline(seqfile, buffer)) {
        vector<string> toks = split_delims(buffer, " \t");
        if (toks.size() == 2) {
            string name = toks[0];
            string fa_path = toks[1];
            unordered_map<string, pair<size_t, string>> fa_info = read_fasta(fa_path);
            for (auto& fi : fa_info) {
                seq_map[name + "." + fi.first] = fi.second.first;
            }
        }
    }
    
    return seq_map;
}

static unordered_map<string, vector<Sequence::Info>> get_filled_dimensions(AlignmentConstPtr alignment, const unordered_map<string, size_t>& seq_d, bool progress) {

    unordered_map<string, vector<Sequence::Info>> dim_map;

    
    vector<string> names = alignment->getChildNames(alignment->getRootName());
    names.push_back(alignment->getRootName());

    for (const string& name : names) {
        const Genome* genome = alignment->openGenome(name);
        vector<Sequence::Info>&  dimensions = dim_map[name];
        if (progress) {
            cerr << "[halUnclip]: Scanning dimensions of genome " << genome->getName() << endl;
        }

        // map base name to sequence fragments
        unordered_map<string, vector<const Sequence*>> frag_map;

        // pass 1, map all hal sequences back to their base name and check that they correspond to a fasta sequence
        for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            string sequence_name = sequence->getName();
            string parsed_name = parse_subpath_name(sequence_name);
            string full_name = genome->getName() + "." + parsed_name;
            size_t fa_len = sequence->getSequenceLength();
            if (name != alignment->getRootName()) {
                if (!seq_d.count(full_name)) {
                    cerr << "[halUnclip]: Unable to find sequence (from HAL) " << full_name << " in dimension map from input fasta" << endl;
                    exit(1);
                }
                fa_len = seq_d.at(full_name);
            }
            if (parsed_name == sequence_name && sequence->getSequenceLength() != fa_len) {
                cerr << "[halUnclip]: Sequence " << full_name << " has len=" << fa_len << " in fasta but len=" << sequence->getSequenceLength() << " in hal" << endl;
                exit(1);
            }
            if (parsed_name != sequence_name && sequence->getSequenceLength() > fa_len) {
                cerr << "[halUnclip]: Sequence " << sequence->getFullName() << " has len=" << fa_len << " in fasta but len=" << sequence->getSequenceLength() << " in hal" << endl;
                exit(1);
            }

            frag_map[parsed_name].push_back(sequence);
        }

        // pass 2: compute the dimensions for each base sequence
        for (auto& nf : frag_map) {
            const string& base_name = nf.first;
            string full_name = genome->getName() + "." + base_name;
            vector<const Sequence*>& frags = nf.second;
            size_t fa_len;
            if (name == alignment->getRootName()) {
                assert(frags.size() == 1);
                fa_len = frags[0]->getSequenceLength();
            } else {
                fa_len = seq_d.at(full_name);
            }
            // sort the fragments by start position
            map<size_t, const Sequence*> start_to_frag;
            for (const Sequence* frag : frags) {
                int64_t start;
                string parsed_name = parse_subpath_name(frag->getName(), &start);
                if (start == -1) {
                    start = 0;
                    assert(frags.size() == 1);
                }
                start_to_frag[start] = frag;
            }

            // count the top segments
            size_t top = 0;
            // count the gaps (separate counter just for debugging)
            size_t gaps = 0;
            if (start_to_frag.begin()->first > 0) {
                // gap in front
                ++gaps;
            }
            for (auto i = start_to_frag.begin(); i != start_to_frag.end(); ++i) {
                auto next = i;
                ++next;
                if (next != start_to_frag.end()) {
                    if (i->first + i->second->getSequenceLength() < next->first) {
                        // gap in middle
                        ++gaps;
                    }
                }
                top += i->second->getNumTopSegments();                
            }
            if (start_to_frag.rbegin()->first + start_to_frag.rbegin()->second->getSequenceLength() < fa_len) {
                // gap in back
                ++gaps;
            }
            if (name == alignment->getRootName()) {
                dimensions.push_back(Sequence::Info(base_name, fa_len, 0, frags[0]->getNumBottomSegments()));
            } else {
                dimensions.push_back(Sequence::Info(base_name, fa_len, top + gaps, 0));
            }
        }
        
        alignment->closeGenome(genome);
    }

    return dim_map;
}

static void copy_and_fill(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment, const unordered_map<string, size_t>& seq_dims) {

    // root gets the same dimensions
    // 
    
}

    
int main(int argc, char** argv) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(&optionsParser);
    string in_hal_path;
    string out_hal_path;
    string seqfile_path;
    bool progress;
    try {
        optionsParser.parseOptions(argc, argv);
        in_hal_path = optionsParser.getArgument<string>("inFile");
        seqfile_path = optionsParser.getArgument<string>("seqFile");
        out_hal_path = optionsParser.getArgument<string>("outFile");
        progress = optionsParser.getFlag("progress");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    // load the input genome
    if (progress) {
        cerr << "[halUnclip]: Opening input alignment" << endl;
    }
    AlignmentConstPtr in_alignment(openHalAlignment(in_hal_path, &optionsParser, READ_ACCESS));

    // and the output genome
    if (progress) {
        cerr << "[halUnclip]: Creating output alignment object" << endl;
    }    
    AlignmentPtr out_alignment(openHalAlignment(out_hal_path, &optionsParser, READ_ACCESS | WRITE_ACCESS | CREATE_ACCESS));

    // and load the fasta sequence sizes from the seqfile
    if (progress) {
        cerr << "[halUnclip]: Reading fasta dimensions from seqfile" << endl;
    }
    unordered_map<string, size_t> seq_dims = get_dimensions_from_seqfile(seqfile_path);

    if (progress) {
        cerr << "[halUnclip]: Computing new hal dimensions" << endl;
    }
    unordered_map<string, vector<Sequence::Info>> dimensions = get_filled_dimensions(in_alignment, seq_dims, progress);

    // set up the size of each genome, staring with the root
    string root_name = in_alignment->getRootName();
    Genome* root_genome = out_alignment->addRootGenome(root_name);
    root_genome->setDimensions(dimensions.at(root_name));
    for (auto& kv : dimensions) {
        if (kv.first != root_name) {
            Genome* leaf_genome = out_alignment->addLeafGenome(kv.first, root_name, 1);
            leaf_genome->setDimensions(kv.second);
        }
    }
    
    // copy over the filled graph
    if (progress) {
        cerr << "[halUnclip]: Copying and filling the graph" << endl;
    }
    copy_and_fill(in_alignment, out_alignment, seq_dims);

    // add back the fasta sequences
    if (progress) {
        cerr << "[halUnclip]: Adding fasta sequences" << endl;
    }    
    //add_fasta_sequences(out_alignment, seqfile_path);

    if (progress) {
        cerr << "[halUnclip]: Writing output alignment" << endl;
    }
     
    return 0;
}

