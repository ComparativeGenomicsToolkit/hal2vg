/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// This file was created by merging hal2sg.cpp and sg2vg.cpp with
// a small amount of glue for the interface. 

#define debug

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>

#include "stPinchGraphs.h"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"
#include "hal.h"

using namespace std;
using namespace hal;
using namespace handlegraph;
using namespace bdsg;
using namespace handlegraph;

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("halFile", "input hal file");
    optionsParser->addOption("rootGenome", 
                             "process only genomes in clade with specified root"
                             " (HAL root if empty)", 
                             "\"\"");
    optionsParser->addOptionFlag("onlySequenceNames",
                                 "use only sequence names for output names.  By "
                                 "default, the UCSC convention of "
                                 "Genome.Sequence is used",
                                 false);
    optionsParser->addOptionFlag("keepCase",
                                 "don't convert all nucleotides to upper case",
                                 false);
    optionsParser->addOption("outputFormat",
                             "output graph format in {pg, hg, odgi} [default=pg]",
                             "pg");

    optionsParser->setDescription("Convert HAL alignment to handle graph");

}

static void add_genome_threads(const Genome* genome,
                               stPinchThreadSet* threads,
                               vector<string>& IDToName,
                               unordered_map<const Sequence*, int64_t>& seqToID,
                               bool fullNames);

static void pinch_genome(const Genome* genome,
                         stPinchThreadSet* threads,
                         unordered_map<const Sequence*, int64_t>& seqToID);

static void pinch_to_handle(const Genome* genome,
                            stPinchThreadSet* threadSet,
                            const vector<string>& IDToName,
                            int64_t& seqID,
                            unordered_map<stPinchBlock*, nid_t>& blockToNode,
                            MutablePathMutableHandleGraph& graph,
                            bool fullNames);

int main(int argc, char** argv) {
    CLParser optionsParser;
    initParser(&optionsParser);
    string halPath;
    string rootGenomeName;
    bool fullNames;
    bool keepCase;
    string outputFormat;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        rootGenomeName = optionsParser.getOption<string>("rootGenome");
        fullNames = !optionsParser.getFlag("onlySequenceNames");
        keepCase = optionsParser.getFlag("keepCase");
        outputFormat = optionsParser.getOption<string>("outputFormat");
        if (outputFormat != "pg" && outputFormat != "hg" && outputFormat != "odgi") {
            throw hal_exception("--outputFormat must be one of {pg, hg, odgi}");
        }
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
    try {
        AlignmentConstPtr alignment(openHalAlignment(halPath, &optionsParser));
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("input hal alignmenet is empty");
        }
        
        // root is specified either by the parameter or as the alignment root
        // by default
        const Genome* rootGenome = NULL;
        if (rootGenomeName != "\"\"") {
            rootGenome = alignment->openGenome(rootGenomeName);
        } else {
            rootGenome = alignment->openGenome(alignment->getRootName());
        }
        if (rootGenome == NULL) {
            throw hal_exception(string("Root genome, ") + rootGenomeName + 
                                ", not found in alignment");
        }

        // map Sequence pointers to integers (assumes sequence pointers stable within hal)
        vector<string> IDToName;
        unordered_map<const Sequence*, int64_t> seqToID;
        
        // start up our pinch graph
        stPinchThreadSet* threadSet = stPinchThreadSet_construct();
        
        const Genome* parentGenome = rootGenome;
        string parentName = rootGenome->getName();
        add_genome_threads(parentGenome, threadSet, IDToName, seqToID, fullNames);

        vector<string> childs = alignment->getChildNames(rootGenome->getName());
        deque<string> queue(childs.begin(), childs.end());

        while (!queue.empty()) {
            string childName = queue.front();
            queue.pop_front();
            const Genome* childGenome = alignment->openGenome(childName);
            string parentName = alignment->getParentName(childName);
            if (parentName != parentGenome->getName()) {
                alignment->closeGenome(parentGenome);
                parentGenome = childGenome->getParent();
                add_genome_threads(parentGenome, threadSet, IDToName, seqToID, fullNames);
            }
            
            add_genome_threads(childGenome, threadSet, IDToName, seqToID, fullNames);
            pinch_genome(childGenome, threadSet, seqToID);
                
            childs = alignment->getChildNames(childName);
            for (int i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }

            alignment->closeGenome(childGenome);
        }

        alignment->closeGenome(parentGenome);
        seqToID.clear();

        // make a handle graph
        unique_ptr<MutablePathMutableHandleGraph> graph;
        if (outputFormat == "pg") {
            graph = unique_ptr<MutablePathMutableHandleGraph>(new PackedGraph());
        } else if (outputFormat == "hg") {
            graph = unique_ptr<MutablePathMutableHandleGraph>(new HashGraph());
        } else if (outputFormat == "odgi") {
            graph = unique_ptr<MutablePathMutableHandleGraph>(new ODGI());
        } else {
            assert(false);
        }

        // keep track of where blocks fit into the handle graph
        unordered_map<stPinchBlock*, nid_t> blockToNode;

        // start iterating over the genomes again in order to export to handle graph
        //
        // IMPORTANT:  This code relies on the the genomes being visited in the
        // exact same order as above
        int64_t seqID = 0;
            
        queue.push_back(parentName);
        while (!queue.empty()) {
            string genomeName = queue.front();
            queue.pop_front();
            const Genome* genome = alignment->openGenome(genomeName);
                        
            pinch_to_handle(genome, threadSet, IDToName, seqID, blockToNode, *graph, fullNames);

            childs = alignment->getChildNames(genomeName);
            for (int i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }
            alignment->closeGenome(genome);
        }

        // write out the graph
        dynamic_cast<SerializableHandleGraph*>(graph.get())->serialize(cout);
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        exit(1);
    }
     
    return 0;
}

// Add every sequence from the genome into the pinch graph
void add_genome_threads(const Genome* genome,
                       stPinchThreadSet* threads,
                       vector<string>& IDToName,
                       unordered_map<const Sequence*, int64_t>& seqToID,
                       bool fullNames) {
    
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        hal_size_t seqLen = sequence->getSequenceLength();
        string name = fullNames ? sequence->getFullName() : sequence->getName();
        // update lookups to map hal sequence to numeric id
        int64_t seqID = IDToName.size(); 
        seqToID[sequence] = seqID;
        IDToName.push_back(name);
        // add to thread set
#ifdef debug
        cerr << "Adding sequence " << name << " as thread " << seqID << " with length " << seqLen << endl;
#endif
        stPinchThreadSet_addThread(threads, seqID, 0, seqLen);
    }
}

// Use exact pairwise alginments from genome to its parent to make the pinch graph
void pinch_genome(const Genome* genome, stPinchThreadSet* threads,
                  unordered_map<const Sequence*, int64_t>& seqToID) {

    TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIt = genome->getParent()->getBottomSegmentIterator();

    // avoid thread set lookups
    const Sequence* topSeq = nullptr;
    const Sequence* botSeq = nullptr;
    stPinchThread* topThread = nullptr;
    stPinchThread* botThread = nullptr;
    string topString;
    string botString;
    
    for (; not topIt->atEnd(); topIt->toRight()) {
        if (topIt->tseg()->hasParent()) {
            botIt->toParent(topIt);

            int64_t topID = seqToID[topIt->tseg()->getSequence()];
            int64_t botID = seqToID[botIt->bseg()->getSequence()];

            if (topIt->tseg()->getSequence() != topSeq) {
                topSeq = topIt->tseg()->getSequence();
                topThread = stPinchThreadSet_getThread(threads, topID);
            }
            if (botIt->bseg()->getSequence() != botSeq) {
                botSeq = botIt->bseg()->getSequence();
                botThread = stPinchThreadSet_getThread(threads, botID);
            }

            topIt->getString(topString);
            botIt->getString(botString);

#ifdef debug
            cerr << "pinching " << endl
                 << "   " << *topIt << endl
                 << "  -" << topString << endl
                 << "   " << *botIt << endl
                 << "  -" << botString << endl;
#endif

            int64_t first_match = -1;
            int64_t last_match = -1;
            for (int64_t i = 0; i < topString.length(); ++i) {
                if (std::toupper(topString[i]) == std::toupper(botString[i])) {
                    if (first_match == -1) {
                        first_match = i;
                    }
                    last_match = i;
                }
                if (std::toupper(topString[i]) != std::toupper(botString[i]) || i == topString.length() - 1) {
                    if (last_match >= first_match && first_match >= 0) {
                        hal_index_t length = last_match - first_match + 1;
                        hal_index_t start1 = topIt->tseg()->getStartPosition() + first_match - topSeq->getStartPosition();
                        hal_index_t start2;
                        if (!botIt->getReversed()) {
                            start2 = botIt->bseg()->getStartPosition() + first_match - botSeq->getStartPosition();
                        } else {
                            start2 = botIt->bseg()->getEndPosition() - first_match - length + 1 - botSeq->getStartPosition();
                        }
#ifdef debug
                        cerr << " inserting (fm=" << first_match <<",lm=" << last_match << ", s1=" << start1 << ",s2=" << start2 << ",l=" << length
                             << ", hl1=" << topSeq->getSequenceLength() << ",hl2=" << botSeq->getSequenceLength() << ",pl1=" << stPinchThread_getLength(topThread)
                             << ", pl2=" << stPinchThread_getLength(botThread) << ", rev=" << botIt->getReversed() << endl
                             << "   " << topString.substr(first_match, length) << endl;
#endif
                        stPinchThread_pinch(topThread,
                                            botThread,
                                            start1,
                                            start2,
                                            length,
                                            !botIt->getReversed());
                    }
                    first_match = -1;
                    last_match = -1;
                }
            }            
        }
    }
}

// create nodes and edges for a genome using the pinch graph
void pinch_to_handle(const Genome* genome,
                     stPinchThreadSet* threadSet,
                     const vector<string>& IDToName,
                     int64_t& seqID,
                     unordered_map<stPinchBlock*, nid_t>& blockToNode,
                     MutablePathMutableHandleGraph& graph,
                     bool fullNames) {

    // iterate over the sequences of the genome
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext(), ++seqID) {
        const Sequence *sequence = seqIt->getSequence();
        string seqName = fullNames ? sequence->getFullName() : sequence->getName();        
        assert(IDToName[seqID] == seqName);
        stPinchThread* thread = stPinchThreadSet_getThread(threadSet, seqID);

        // create the path
        path_handle_t pathHandle = graph.create_path_handle(seqName);
        string pathString;
        
        // iterate over the segments of the sequence
        stPinchSegment* prevSeg = NULL;
        stPinchBlock* prevBlock = NULL;
        bool prevRev = false;
        handle_t prevHandle;
        stPinchSegment* lastSeg = stPinchThread_getLast(thread);
        hal_index_t segStart = 0;
        string seqString;
        for (stPinchSegment* seg = stPinchThread_getFirst(thread); ;
             seg = stPinchSegment_get3Prime(seg)) {

            // get the segment's block.  note that if it's not aligned to anything, it will have no block
            stPinchBlock* block = stPinchSegment_getBlock(seg);
            bool reversed = block != NULL && stPinchSegment_getBlockOrientation(seg) == 0;
            handle_t handle;

            // have we already converted this block?
            auto bi = blockToNode.find(block);
            if (bi == blockToNode.end()) {
                // no: it is a new block
                sequence->getSubString(seqString, segStart, stPinchSegment_getLength(seg));
                if (reversed) {
                    // we always work in block-relative orientation
                    reverseComplement(seqString);
                }
#ifdef debug
                handle = graph.create_handle(seqString);
                if (block != NULL) {
                    blockToNode[block] = graph.get_id(handle);
                } 
                cerr << "created node " << graph.get_id(handle) << " for block " << block << " from " << sequence->getFullName() << " at " << segStart << endl;
#endif
            } else {
                // yes: we can find it in the table
                handle = graph.get_handle(bi->second);
#ifdef debug
                cerr << "found node " << graph.get_id(handle) << " for block " << block << " from " << sequence->getFullName() << " at " << segStart << endl;
#endif                
            }
            assert(!graph.get_is_reverse(handle));
            if (reversed) {
                graph.flip(handle);
            }
                   
            // wire up the edge to previous
            if (prevSeg != NULL) {
#ifdef debug
                cerr << "creating edge from " << graph.get_id(prevHandle) << ":" << graph.get_is_reverse(prevHandle) << " -> "
                     << graph.get_id(handle) << ":" << graph.get_is_reverse(handle) << endl;
#endif
                graph.create_edge(prevHandle, handle);
            }

            // add the node to the path
            graph.append_step(pathHandle, handle);
            pathString += graph.get_sequence(handle);

            prevRev = reversed;
            prevBlock = block;
            prevSeg = seg;
            prevHandle = handle;
            
            segStart += stPinchSegment_getLength(seg);
            
            if (seg == lastSeg) {
                break;
            }
        }

        // make sure the path we added is the same as the hal
        string halPathString;
        sequence->getString(halPathString);
        if (pathString.length() != halPathString.length()) {
            throw runtime_error("Incorrect length in coverted path for " + seqName + ": " + std::to_string(pathString.length()) +
                                ". Should be: " + std::to_string(halPathString.length()));
        }
        vector<size_t> mismatches;
        for (size_t i = 0; i < halPathString.size(); ++i) {
            if (toupper(pathString[i]) != toupper(halPathString[i])) {
                mismatches.push_back(i);
            }
        }
        if (!mismatches.empty()) {
            stringstream msg;
            msg << mismatches.size() << " mismatches found in converted path for " << seqName << ":\n";
            for (size_t i = 0; i < mismatches.size() && i < 10; ++i) {
                msg << " path[" << mismatches[i] << "]=" << pathString[mismatches[i]] << ". should be " << halPathString[mismatches[i]] << "\n";
            }
            throw runtime_error(msg.str());
        }                    
    }
}
