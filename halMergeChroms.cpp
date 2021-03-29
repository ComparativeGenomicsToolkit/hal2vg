/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// Merge chromosome HAL files into one big one.  Only star trees supported, with same root name supported. 

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
    optionsParser->setDescription("Merge chromosome HALs into combined file.  Ancestral sequences are renamed as needed to avoid conflicts"
        ". Star trees only");
}

// get the dimensions from a file
static pair<string, unordered_map<string, vector<Sequence::Info>>> get_hal_dimensions(CLParser* optionsParser, const string& hal_path) {

    // open the hal file
    AlignmentConstPtr alignment(openHalAlignment(hal_path, optionsParser, READ_ACCESS));

    // our output map
    unordered_map<string, vector<Sequence::Info>> dimensions;
    
    // for every genome
    vector<string> genome_names = alignment->getChildNames(alignment->getRootName());
    genome_names.push_back(alignment->getRootName());
    for (const string& genome_name : genome_names) {
        const Genome* genome = alignment->openGenome(genome_name);
        vector<Sequence::Info>& genome_dimensions = dimensions[genome_name];
        // for every sequence
        for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
            const Sequence *sequence = seqIt->getSequence();
            genome_dimensions.emplace_back(sequence->getName(),
                                           sequence->getSequenceLength(),
                                           sequence->getNumTopSegments(),
                                           sequence->getNumBottomSegments());
        }
        alignment->closeGenome(genome);        
    }

    return make_pair(alignment->getRootName(), dimensions);
}

// append each input hal to the out_alignment, one after another.  all arrays are copied over,
// but need to be adjsuted to reflect their new offsets.
static void merge_hals(CLParser* optionsParser, AlignmentPtr out_alignment, const vector<string>& in_paths) {

    // keep track of where we are in the output    
    vector<size_t> top_offsets(out_alignment->getChildNames(out_alignment->getRootName()).size(), 0);
    size_t bot_offset = 0;
    
    for (size_t i = 0; i < in_paths.size(); ++i) {
        AlignmentConstPtr in_alignment(openHalAlignment(in_paths[i], optionsParser, READ_ACCESS));
        const Genome* in_root = in_alignment->openGenome(in_alignment->getRootName());
        Genome* out_root = out_alignment->openGenome(in_alignment->getRootName());
        assert(in_root->getName() == out_root->getName());
        size_t in_root_degree = in_root->getNumChildren();
        vector<const Genome*> in_genomes = {in_root};

        // make a child index map (in -> out) for the root genome
        // assume : all genomes in in_genome present in out_genome
        vector<size_t> in_ci_to_out_ci(in_root->getNumChildren());
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
            in_ci_to_out_ci.at(in_root->getChildIndex(in_alignment->openGenome(in_child_name))) =
                out_root->getChildIndex(out_alignment->openGenome(in_child_name));
        }

        // copy over the bottom segments of the root
        BottomSegmentIteratorPtr in_bi = in_root->getBottomSegmentIterator(0);
        BottomSegmentIteratorPtr out_bi = out_root->getBottomSegmentIterator(bot_offset);
        for (;!in_bi->atEnd(); in_bi->toRight(), out_bi->toRight()) {
            // set the segment in the root genome
            assert(out_bi->bseg()->getArrayIndex() == in_bi->bseg()->getArrayIndex() + bot_offset);
            out_bi->bseg()->setTopParseIndex(NULL_INDEX);
            // determine the sequence-relative coordinate in the input
            const Sequence* in_sequence = in_bi->bseg()->getSequence();
            int64_t in_start_coord = in_bi->bseg()->getStartPosition() - in_sequence->getStartPosition();
            // set the sequence relative coordinate in the output
            const Sequence* out_sequence = out_root->getSequence(in_sequence->getName());
            int64_t out_start_coord = out_sequence->getStartPosition() + in_start_coord;
            out_bi->bseg()->setCoordinates(out_start_coord, in_bi->bseg()->getLength());

            // set the segment in the child genomes
            for (size_t in_ci = 0; in_ci < in_root_degree; ++in_ci) {
                size_t out_ci = in_ci_to_out_ci[in_ci];
                out_bi->bseg()->setChildIndex(out_ci, in_bi->bseg()->getChildIndex(in_ci) + top_offsets[out_ci]);
                out_bi->bseg()->setChildReversed(out_ci, in_bi->bseg()->getChildReversed(in_ci));
            }
        }

        // for every child genome, copy over the top segments
        for (const string& in_child_name : in_alignment->getChildNames(in_root->getName())) {
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
                out_ti->tseg()->setParentIndex(in_ti->tseg()->getParentIndex() + bot_offset);
                out_ti->tseg()->setBottomParseIndex(NULL_INDEX);
                // determine the sequence-relative coordinate in the input
                const Sequence* in_sequence = in_ti->tseg()->getSequence();
                int64_t in_start_coord = in_ti->tseg()->getStartPosition() - in_sequence->getStartPosition();
                // set the sequence relative coordinate in the output
                const Sequence* out_sequence = out_root->getSequence(in_sequence->getName());
                int64_t out_start_coord = out_sequence->getStartPosition() + in_start_coord;
                out_ti->tseg()->setCoordinates(out_start_coord, in_ti->tseg()->getLength());
                // set the paralogy edge
                if (in_ti->tseg()->hasNextParalogy()) {
                    out_ti->tseg()->setNextParalogyIndex(in_ti->tseg()->getNextParalogyIndex() + top_offset);
                } else {
                    out_ti->tseg()->setNextParalogyIndex(NULL_INDEX);
                }
            }
            in_genomes.push_back(in_child);
        }

        // copy the dna sequence by sequence
        for (const Genome* in_genome : in_genomes) {
            Genome* out_genome = out_alignment->openGenome(in_genome->getName());
            for (SequenceIteratorPtr in_si = in_genome->getSequenceIterator(); !in_si->atEnd(); in_si->toNext()) {
                const Sequence* in_sequence = in_si->getSequence();
                Sequence* out_sequence = out_genome->getSequence(in_sequence->getName());
                string dna;
                in_sequence->getString(dna);
                out_sequence->setString(dna);
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
    try {
        optionsParser.parseOptions(argc, argv);
        in_hal_paths = optionsParser.getArgument<string>("inFiles");
        out_hal_path = optionsParser.getArgument<string>("outFile");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }

    vector<string> in_paths = chopString(in_hal_paths, ",");
    
    // map genome -> dimensions for each input alignment
    vector<unordered_map<string, vector<Sequence::Info>>> dimensions_list;
    string root_name;
    for (const string& in_path : in_paths) {
        auto rd = get_hal_dimensions(&optionsParser, in_path);
        dimensions_list.push_back(rd.second);
        if (root_name.empty()) {
            root_name = rd.first;
            assert(in_path == in_paths[0]);
        } else if (rd.first != root_name) {
            throw hal_exception("Root mismatch: " + rd.first + " vs " + root_name);
        }
    }

    // check uniqueness and sum up
    // todo: resolve ancestor name collisions    
    unordered_set<string> sequence_names;
    unordered_map<string, vector<Sequence::Info>> total_dimensions;    
    for (const unordered_map<string, vector<Sequence::Info>>& dimensions : dimensions_list) {
        for (const auto& kv : dimensions) {
            for (const Sequence::Info& sequence_info : kv.second) {
                string name = kv.first + "." + sequence_info._name;
                if (sequence_names.count(name)) {
                    throw hal_exception("Conflict: sequence name found in more than one file: " + name);
                }
                sequence_names.insert(name);
            }
            total_dimensions[kv.first].insert(total_dimensions[kv.first].begin(), kv.second.begin(), kv.second.end());
        }
    }

    // create the new file
    AlignmentPtr alignment(openHalAlignment(out_hal_path, &optionsParser, READ_ACCESS | WRITE_ACCESS | CREATE_ACCESS));

    // set up the size of each genome
    Genome* root_genome = alignment->addRootGenome(root_name);
    root_genome->setDimensions(total_dimensions.at(root_name));    
    for (const auto& gd : total_dimensions) {
        if (gd.first != root_name) {
            Genome* leaf_genome = alignment->addLeafGenome(gd.first, root_name, 1);
            leaf_genome->setDimensions(gd.second);
        }
    }

    // copy over over everything
    merge_hals(&optionsParser, alignment, in_paths);

     
    return 0;
}

