/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// Merge chromosome HAL files into one big one.  Only star trees with same root name supported (ie what comes out of cactus-align-batch). 

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

using namespace std;
using namespace hal;

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("inFiles", "comma-separated (only way in HAL parser!) list of input HAL files to merge");
    optionsParser->addArgument("outFile", "output HAL file");
    optionsParser->addOptionFlag("progress",
                                 "show progress",
                                 false);
    optionsParser->setDescription("Merge chromosome HALs into combined file.  Ancestral sequences are renamed as needed to avoid conflicts"
        ". Star trees only");
}

// we expect to see the same ancestor sequence names in multiple input files. we uniqify them by adding
// .i to them where i is the file's position in the input.
static string anc_seq_name(const string& seq_name, size_t idx) {
    return seq_name + ".hmc" + to_string(idx);
}
// undo the above
static string orig_seq_name(const string& seq_name) {
    return seq_name.substr(0, seq_name.rfind(".hmc"));
}

// get the dimensions from all genomes in all input files
static pair<string, unordered_map<string, vector<Sequence::Info>>> get_hal_dimensions(CLParser* optionsParser,
                                                                                      const vector<string>& hal_paths) {
    // genome -> dimensions (covering all input)
    unordered_map<string, vector<Sequence::Info>> dimensions;
    // to check uniqueness
    unordered_set<string> sequence_names;

    string root_name;
    for (size_t i = 0; i < hal_paths.size(); ++i) {
        const string& hal_path = hal_paths[i];
        
        // open the hal file
        AlignmentConstPtr alignment(openHalAlignment(hal_path, optionsParser, READ_ACCESS));
    
        // for every genome
        vector<string> genome_names = alignment->getChildNames(alignment->getRootName());
        genome_names.push_back(alignment->getRootName());
        if (root_name.empty()) {
            root_name = alignment->getRootName();
        } else if (alignment->getRootName() != root_name) {
            throw hal_exception("Root mismatch: " + root_name + " vs " + alignment->getRootName());        
        }
        for (const string& genome_name : genome_names) {
            const Genome* genome = alignment->openGenome(genome_name);
            vector<Sequence::Info>& genome_dimensions = dimensions[genome_name];
            // for every sequence
            for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
                const Sequence *sequence = seqIt->getSequence();
                // add a little suffix to make ancestral sequences unique
                string seq_name = genome->getParent() ? sequence->getName() : anc_seq_name(sequence->getName(), i);
                genome_dimensions.emplace_back(seq_name,
                                               sequence->getSequenceLength(),
                                               sequence->getNumTopSegments(),
                                               sequence->getNumBottomSegments());
                string full_name = genome_name + "." + seq_name;
                if (sequence_names.count(full_name)) {
                    throw hal_exception("Duplicate sequence name found: " + full_name);
                } else {
                    sequence_names.insert(full_name);
                }
            }
            alignment->closeGenome(genome);        
        }
    }
    return make_pair(root_name, dimensions);
}

// append each input hal to the out_alignment, one after another.  all arrays are copied over,
// but need to be adjsuted to reflect their new offsets.
static void merge_hals(CLParser* optionsParser, AlignmentPtr out_alignment, const vector<string>& in_paths, bool progress) {

    // keep track of where we are in the output    
    vector<size_t> top_offsets(out_alignment->getChildNames(out_alignment->getRootName()).size(), 0);
    size_t bot_offset = 0;
    
    for (size_t i = 0; i < in_paths.size(); ++i) {
        AlignmentConstPtr in_alignment(openHalAlignment(in_paths[i], optionsParser, READ_ACCESS));
        const Genome* in_root = in_alignment->openGenome(in_alignment->getRootName());
        Genome* out_root = out_alignment->openGenome(in_alignment->getRootName());
        assert(in_root->getName() == out_root->getName());
        size_t in_root_degree = in_root->getNumChildren();
        size_t out_root_degree = out_root->getNumChildren();
        vector<const Genome*> in_genomes = {in_root};
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
            in_genomes.push_back(in_alignment->openGenome(in_child_name));
        }

        // copy the dna sequence by sequence
        for (const Genome* in_genome : in_genomes) {
            if (progress) {
                cerr << "[halMergeChroms]: copying dna for " << in_genome->getName() << " from " << in_paths[i] << endl;
            }
            Genome* out_genome = out_alignment->openGenome(in_genome->getName());
            for (SequenceIteratorPtr in_si = in_genome->getSequenceIterator(); !in_si->atEnd(); in_si->toNext()) {
                const Sequence* in_sequence = in_si->getSequence();
                string out_seq_name = in_genome->getParent() ? in_sequence->getName() : anc_seq_name(in_sequence->getName(), i);
                Sequence* out_sequence = out_genome->getSequence(out_seq_name);
                DnaIteratorPtr in_di = in_sequence->getDnaIterator(0);
                DnaIteratorPtr out_di = out_sequence->getDnaIterator(0);
                assert(in_sequence->getSequenceLength() == out_sequence->getSequenceLength());
                string dna;
                in_sequence->getString(dna);
                out_sequence->setString(dna);
            }
        }

        // make a child index map (in -> out) for the root genome
        // assume : all genomes in in_genome present in out_genome
        vector<size_t> in_ci_to_out_ci(in_root->getNumChildren());
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
            in_ci_to_out_ci.at(in_root->getChildIndex(in_alignment->openGenome(in_child_name))) =
                out_root->getChildIndex(out_alignment->openGenome(in_child_name));
        }

        // copy over the bottom segments of the root
        if (progress) {
            cerr << "[halMergeChroms]: copying bottom segments for " << in_root->getName() << " from " << in_paths[i]
                 << " with bseg offset " << bot_offset << endl;
        }
        BottomSegmentIteratorPtr in_bi = in_root->getBottomSegmentIterator(0);
        BottomSegmentIteratorPtr out_bi = out_root->getBottomSegmentIterator(bot_offset);
        for (;!in_bi->atEnd(); in_bi->toRight(), out_bi->toRight()) {
            // set the segment in the root genome
            assert(out_bi->bseg()->getArrayIndex() == in_bi->bseg()->getArrayIndex() + bot_offset);
            assert(out_bi->bseg()->getNumChildren() >= in_bi->bseg()->getNumChildren());
            out_bi->bseg()->setTopParseIndex(NULL_INDEX);
            // determine the sequence-relative coordinate in the input
            const Sequence* in_sequence = in_bi->bseg()->getSequence();
            int64_t in_start_coord = in_bi->bseg()->getStartPosition() - in_sequence->getStartPosition();
            assert(in_start_coord >= 0 && in_start_coord < in_sequence->getSequenceLength());
            // set the sequence relative coordinate in the output
            const Sequence* out_sequence = out_root->getSequence(anc_seq_name(in_sequence->getName(), i));
            int64_t out_start_coord = out_sequence->getStartPosition() + in_start_coord;
            assert(out_start_coord >= 0 && out_start_coord < out_root->getSequenceLength());
            out_bi->bseg()->setCoordinates(out_start_coord, in_bi->bseg()->getLength());
            // set the segment in the child genomes
            for (size_t out_ci = 0; out_ci < out_root_degree; ++out_ci) {
                out_bi->bseg()->setChildIndex(out_ci,  NULL_INDEX);
            }
            for (size_t in_ci = 0; in_ci < in_root_degree; ++in_ci) {
                size_t out_ci = in_ci_to_out_ci.at(in_ci);
                assert(out_ci < out_bi->bseg()->getNumChildren());
                if (in_bi->bseg()->hasChild(in_ci)) {
                    out_bi->bseg()->setChildIndex(out_ci,  in_bi->bseg()->getChildIndex(in_ci) + top_offsets[out_ci]);
                    out_bi->bseg()->setChildReversed(out_ci, in_bi->bseg()->getChildReversed(in_ci));
                }
            }
        }

        // for every child genome, copy over the top segments
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
            if (progress) { 
                cerr << "[halMergeChroms]: copying top segments for " << in_child_name << " from " << in_paths[i] << endl;
            }
            const Genome* in_child = in_alignment->openGenome(in_child_name);
            Genome* out_child = out_alignment->openGenome(in_child_name);
            
            size_t in_ci = in_root->getChildIndex(in_child);
            size_t out_ci = in_ci_to_out_ci[in_ci];
            size_t top_offset = top_offsets[out_ci];
            TopSegmentIteratorPtr in_ti = in_child->getTopSegmentIterator(0);
            TopSegmentIteratorPtr out_ti = out_child->getTopSegmentIterator(top_offsets[out_ci]);

            for (;!in_ti->atEnd(); in_ti->toRight(), out_ti->toRight()) {
                // set the segment in the child genome
                assert(out_ti->tseg()->getArrayIndex() == in_ti->tseg()->getArrayIndex() + top_offset);
                if (in_ti->tseg()->hasParent()) {
                    out_ti->tseg()->setParentIndex(in_ti->tseg()->getParentIndex() + bot_offset);
                    out_ti->tseg()->setParentReversed(in_ti->tseg()->getParentReversed());
                } else {
                    out_ti->tseg()->setParentIndex(NULL_INDEX);
                }
                out_ti->tseg()->setBottomParseIndex(NULL_INDEX);
                // determine the sequence-relative coordinate in the input
                const Sequence* in_sequence = in_ti->tseg()->getSequence();
                int64_t in_start_coord = in_ti->tseg()->getStartPosition() - in_sequence->getStartPosition();
                // set the sequence relative coordinate in the output
                const Sequence* out_sequence = out_child->getSequence(in_sequence->getName());
                int64_t out_start_coord = out_sequence->getStartPosition() + in_start_coord;
                out_ti->tseg()->setCoordinates(out_start_coord, in_ti->tseg()->getLength());
                // set the paralogy edge
                if (in_ti->tseg()->hasNextParalogy()) {
                    out_ti->tseg()->setNextParalogyIndex(in_ti->tseg()->getNextParalogyIndex() + top_offset);
                } else {
                    out_ti->tseg()->setNextParalogyIndex(NULL_INDEX);
                }
            }
        }

        // update the offsets to move past the current alignment in all genomes
        bot_offset += in_root->getNumBottomSegments();
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
            const Genome* in_child = in_alignment->openGenome(in_child_name);
            size_t in_ci = in_root->getChildIndex(in_child);
            size_t out_ci = in_ci_to_out_ci[in_ci];
            top_offsets[out_ci] += in_child->getNumTopSegments();
        }
    }
}
    
int main(int argc, char** argv) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(&optionsParser);
    string in_hal_paths;
    string out_hal_path;
    bool progress;
    try {
        optionsParser.parseOptions(argc, argv);
        in_hal_paths = optionsParser.getArgument<string>("inFiles");
        out_hal_path = optionsParser.getArgument<string>("outFile");
        progress = optionsParser.getFlag("progress");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    vector<string> in_paths = chopString(in_hal_paths, ",");
    
    // map genome -> dimensions for each input alignment
    if (progress) {
        cerr << "[halMergeChroms]: Scanning dimensions of " << in_paths.size() << " input files." << endl;
    }
    pair<string, unordered_map<string, vector<Sequence::Info>>> rd = get_hal_dimensions(&optionsParser, in_paths);
    string& root_name = rd.first;
    unordered_map<string, vector<Sequence::Info>>& dimensions = rd.second;

    // create the new file
    if (progress) {
        cerr << "[halMergeChroms]: Creating empty alignment: " << out_hal_path << endl;
    }
    AlignmentPtr alignment(openHalAlignment(out_hal_path, &optionsParser, READ_ACCESS | WRITE_ACCESS | CREATE_ACCESS));

    // set up the size of each genome, staring with the root
    Genome* root_genome = alignment->addRootGenome(root_name);
    for (auto& kv : dimensions) {
        if (kv.first != root_name) {
            Genome* leaf_genome = alignment->addLeafGenome(kv.first, root_name, 1);
            leaf_genome->setDimensions(kv.second);
        }
    }
    // important to set root dimensions after adding leaves so bottom segments have right number of slots
    root_genome->setDimensions(dimensions.at(root_name));

    // copy over over everything
    merge_hals(&optionsParser, alignment, in_paths, progress);

    if (progress) {
      cerr << "[halMergeChroms]: Writing merged alignment" << endl;
    }
     
    return 0;
}

