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
    optionsParser->addOptionFlag("validate",
                                 "run a (non-exhaustive) check on the output",
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

    fclose(fa_file);

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

    // copy root exactly as is
    vector<Sequence::Info>& root_dims = dim_map[alignment->getRootName()];
    const Genome* root_genome = alignment->openGenome(alignment->getRootName());
    for (SequenceIteratorPtr seqIt = root_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        root_dims.push_back(Sequence::Info(sequence->getName(), sequence->getSequenceLength(), sequence->getNumTopSegments(), sequence->getNumBottomSegments()));
    }
    
    vector<string> names = alignment->getChildNames(alignment->getRootName());

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
            if (!seq_d.count(full_name)) {
                if (parsed_name != sequence_name) {
                    cerr << "[halUnclip]: Unable to find sequence (from HAL) " << full_name << " in dimension map from input fasta" << endl;
                    exit(1);
                }
            } else {
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
            fa_len = seq_d.at(full_name);
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
            dimensions.push_back(Sequence::Info(base_name, fa_len, top + gaps, 0));
        }
        
        alignment->closeGenome(genome);
    }

    return dim_map;
}

static void copy_and_fill(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment, const unordered_map<string, size_t>& seq_dims) {

    // just copy the root
    const Genome* in_root_genome = in_alignment->openGenome(in_alignment->getRootName());
    Genome* out_root_genome = out_alignment->openGenome(in_alignment->getRootName());
    BottomSegmentIteratorPtr in_botit = in_root_genome->getBottomSegmentIterator();
    BottomSegmentIteratorPtr out_botit = out_root_genome->getBottomSegmentIterator();
    assert(in_root_genome->getNumBottomSegments() == out_root_genome->getNumBottomSegments());
    assert(in_root_genome->getNumChildren() == out_root_genome->getNumChildren());
    for (size_t i = 0; i < in_root_genome->getNumBottomSegments(); ++i) {
        out_botit->bseg()->setCoordinates(in_botit->bseg()->getStartPosition(), in_botit->bseg()->getLength());
        // don't set child indexes, they will get done in both directions by the leaves.
        for (size_t j = 0; j < in_root_genome->getNumChildren(); ++j) {
            out_botit->bseg()->setChildIndex(j, NULL_INDEX);
            out_botit->bseg()->setChildReversed(j, false);
        }
        out_botit->bseg()->setTopParseIndex(NULL_INDEX);
        in_botit->toRight();
        out_botit->toRight();
    }
    
    vector<string> names = in_alignment->getChildNames(in_alignment->getRootName());

    for (const string& name : names) {
        const Genome* in_genome = in_alignment->openGenome(name);
        Genome* out_genome = out_alignment->openGenome(name);
        hal_index_t out_child_no = out_root_genome->getChildIndex(out_genome);

        // map base name to sequence fragments
        // todo: same thing done in above funciton -- generalize?
        unordered_map<string, vector<const Sequence*>> frag_map;

        // pass 1, map all hal sequences back to their base name and check that they correspond to a fasta sequence
        for (SequenceIteratorPtr seqIt = in_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            string sequence_name = sequence->getName();
            string parsed_name = parse_subpath_name(sequence_name);
            frag_map[parsed_name].push_back(sequence);
        }

        // pass 2, copy each sequence fragment by fragment
        for (auto& nf : frag_map) {
            const string& base_name = nf.first;
            // the one output sequence that corresponds to the list of fragments in the input
            Sequence* out_sequence = out_genome->getSequence(base_name);
            assert(out_sequence != nullptr);
            TopSegmentIteratorPtr out_top = out_sequence->getTopSegmentIterator();
            
            int64_t cur_pos = 0;
            int64_t frag_start = -1;
            vector<const Sequence*>& frags = nf.second;
            for (size_t i = 0; i < frags.size(); ++i) {
                const Sequence* in_sequence_frag = frags[i];
                parse_subpath_name(in_sequence_frag->getName(), &frag_start);
                if (frag_start == -1) {
                    frag_start = 0;
                }
                if (frag_start > cur_pos + 1) {
                    // need to add a gap *before* this fragment                    
                    out_top->tseg()->setCoordinates(cur_pos, out_sequence->getStartPosition() + frag_start - cur_pos);
                    out_top->tseg()->setParentIndex(NULL_INDEX);
                    out_top->tseg()->setNextParalogyIndex(NULL_INDEX);
                    out_top->tseg()->setBottomParseIndex(NULL_INDEX);
                    cur_pos += out_top->tseg()->getLength();
                    out_top->toRight();
                }
                // copy the fragment.  note that the ancestor coordinates haven't changed
                // any, so those coordinates can go directly
                TopSegmentIteratorPtr frag_top = in_sequence_frag->getTopSegmentIterator();
                int64_t out_top_offset = out_top->tseg()->getArrayIndex();
#ifdef debug
                cerr << "frag " << in_sequence_frag->getFullName() << " has " << in_sequence_frag->getNumTopSegments() << " topsegs which will map to range "
                     << out_sequence->getTopSegmentIterator()->tseg()->getArrayIndex() << " - "
                     << (out_sequence->getTopSegmentIterator()->tseg()->getArrayIndex() + in_sequence_frag->getNumTopSegments()) << endl;
#endif
                for (size_t i = 0; i < in_sequence_frag->getNumTopSegments(); ++i) {
                    out_top->tseg()->setCoordinates(out_sequence->getStartPosition() + cur_pos, frag_top->tseg()->getLength());
                    out_top->tseg()->setParentIndex(frag_top->tseg()->getParentIndex());
                    out_top->tseg()->setParentReversed(frag_top->tseg()->getParentReversed());
                    if (frag_top->tseg()->hasNextParalogy()) {
                        out_top->tseg()->setNextParalogyIndex(frag_top->tseg()->getNextParalogyIndex() + out_top_offset);
                    } else {
                        out_top->tseg()->setNextParalogyIndex(NULL_INDEX);
                    }
                    if (frag_top->tseg()->hasParent()) {
                        out_botit->toParent(out_top);
                        out_botit->bseg()->setChildIndex(out_child_no, out_top->tseg()->getArrayIndex());
                        out_botit->bseg()->setChildReversed(out_child_no, out_top->tseg()->getParentReversed());
                    }
                    out_top->tseg()->setBottomParseIndex(NULL_INDEX);
                    
                    cur_pos += out_top->tseg()->getLength();         
                    frag_top->toRight();
                    out_top->toRight();
                }            
            }
            if (cur_pos < (int64_t)out_sequence->getSequenceLength()) {
                // needto add a gap *after* the last fragment
                out_top->tseg()->setCoordinates(out_sequence->getStartPosition() + cur_pos, (int64_t)out_sequence->getSequenceLength() - cur_pos);
                out_top->tseg()->setParentIndex(NULL_INDEX);
                out_top->tseg()->setNextParalogyIndex(NULL_INDEX);
                out_top->tseg()->setBottomParseIndex(NULL_INDEX);
                cur_pos += out_top->tseg()->getLength();                
                out_top->toRight();
            }
            assert(cur_pos == (int64_t)out_sequence->getSequenceLength());

        }
        in_alignment->closeGenome(in_genome);
        out_alignment->closeGenome(out_genome);
    }
}

// go in and rewerite the sequences from the fasta
void add_fasta_sequences(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment, const string& seqfile_path) {
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
            Genome* genome = out_alignment->openGenome(name);
            assert(genome != nullptr);
            for (auto& fi : fa_info) {
                Sequence* sequence = genome->getSequence(fi.first);
                if (sequence != nullptr) {
                    assert(sequence->getSequenceLength() == fi.second.first);
                    sequence->setString(fi.second.second);
                }
            }
        }
    }

    const Genome* in_root_genome = in_alignment->openGenome(in_alignment->getRootName());
    Genome* out_root_genome = out_alignment->openGenome(in_alignment->getRootName());
    for (SequenceIteratorPtr seqIt = in_root_genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence* in_sequence = seqIt->getSequence();
        Sequence* out_sequence = out_root_genome->getSequence(in_sequence->getName());
        in_sequence->getString(buffer);
        out_sequence->setString(buffer);
    }    
}


// root->leaf alignments are consistent
static void validate_alignments(AlignmentConstPtr in_alignment, AlignmentPtr out_alignment) {

    validateAlignment(out_alignment.get());

    const Genome* in_root_genome = in_alignment->openGenome(in_alignment->getRootName());
    Genome* out_root_genome = out_alignment->openGenome(in_alignment->getRootName());
    assert(in_root_genome->getNumBottomSegments() == out_root_genome->getNumBottomSegments());
    assert(in_root_genome->getNumChildren() == out_root_genome->getNumChildren());
    // we go by genome (instead of segment) to hopefully be cache-friendlier
    for (size_t j = 0; j < in_root_genome->getNumChildren(); ++j) {
        const Genome* in_genome = in_root_genome->getChild(j);
        Genome* out_genome = out_root_genome->getChild(j);
        BottomSegmentIteratorPtr in_botit = in_root_genome->getBottomSegmentIterator();
        BottomSegmentIteratorPtr out_botit = out_root_genome->getBottomSegmentIterator();
        TopSegmentIteratorPtr in_topit = in_genome->getTopSegmentIterator();
        TopSegmentIteratorPtr out_topit = out_genome->getTopSegmentIterator();
        for (size_t i = 0; i < in_genome->getNumBottomSegments(); ++i) {
            in_topit->toChild(in_botit, j);
            out_topit->toChild(out_botit, j);
            
            string s1, s2;
            if (j == 0) {
                in_botit->getString(s1);
                out_botit->getString(s2);
                assert(s1 == s2);
            }
            in_topit->getString(s1);
            out_topit->getString(s2);
            assert(s1 == s2);

            string in_seq_name = in_topit->tseg()->getSequence()->getName();
            string out_seq_name = out_topit->tseg()->getSequence()->getName();
            int64_t start;
            string in_base_name = parse_subpath_name(in_seq_name, &start);
            assert(in_base_name == out_seq_name);
            assert(in_topit->getReversed() == out_topit->getReversed());
            if (!in_topit->getReversed()) {
                // punt on reverse check for now
                assert(in_topit->getStartPosition() + start == out_topit->getStartPosition());
            }

            in_botit->toRight();
            out_botit->toRight();
        }
        in_alignment->closeGenome(in_genome);
        out_alignment->closeGenome(out_genome);        
    }
}

int main(int argc, char** argv) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(&optionsParser);
    string in_hal_path;
    string out_hal_path;
    string seqfile_path;
    bool progress;
    bool validate;
    try {
        optionsParser.parseOptions(argc, argv);
        in_hal_path = optionsParser.getArgument<string>("inFile");
        seqfile_path = optionsParser.getArgument<string>("seqFile");
        out_hal_path = optionsParser.getArgument<string>("outFile");
        progress = optionsParser.getFlag("progress");
        validate = optionsParser.getFlag("validate");
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
    for (auto& kv : dimensions) {
        if (kv.first != root_name) {
            Genome* leaf_genome = out_alignment->addLeafGenome(kv.first, root_name, 1);
            leaf_genome->setDimensions(kv.second);
            if (progress) {
                cerr << "[halUnclip]: Adding leaf genome " << kv.first << " with length " << leaf_genome->getSequenceLength() << " and " << leaf_genome->getNumTopSegments() << " top segments" << endl;
            }
        }
    }

    // important to set root dimensions after adding leaves so bottom segments have right number of slots
    root_genome->setDimensions(dimensions.at(root_name));
    if (progress) {
        cerr << "[halUnclip]: Adding root genome " << root_name << " with length " << root_genome->getSequenceLength() << " and " << root_genome->getNumBottomSegments() << " bottom segments" << endl;
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
    add_fasta_sequences(in_alignment, out_alignment, seqfile_path);

    if (validate) {
        if (progress) {
            cerr << "[halUnclip]: Validating alignment" << endl;
        }
        validate_alignments(in_alignment, out_alignment);
    }    
            
    if (progress) {
        cerr << "[halUnclip]: Writing output alignment" << endl;
    }
     
    return 0;
}

