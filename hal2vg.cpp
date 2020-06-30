/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// This file was created by merging hal2sg.cpp and sg2vg.cpp with
// a small amount of glue for the interface. 

//#define debug

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
    optionsParser->addOption("outputFormat",
                             "output graph format in {pg, hg, odgi} [default=pg]",
                             "pg");

    optionsParser->setDescription("Convert HAL alignment to handle graph");

}

static void add_genome_threads(const Genome* genome,
                               stPinchThreadSet* threads,
                               vector<string>& IDToName,
                               unordered_map<string, int64_t>& nameToID,
                               bool fullNames);

static void pinch_genome(const Genome* genome,
                         stPinchThreadSet* threads,
                         unordered_map<string, int64_t>& nameToID,
                         bool fullNames);

static void pinch_to_handle(const Genome* genome,
                            stPinchThreadSet* threadSet,
                            const vector<string>& IDToName,
                            const unordered_map<string, int64_t>& nameToID,
                            unordered_map<stPinchBlock*, nid_t>& blockToNode,
                            MutablePathMutableHandleGraph& graph,
                            bool fullNames);

int main(int argc, char** argv) {
    CLParser optionsParser;
    initParser(&optionsParser);
    string halPath;
    string rootGenomeName;
    bool fullNames;
    string outputFormat;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        rootGenomeName = optionsParser.getOption<string>("rootGenome");
        fullNames = !optionsParser.getFlag("onlySequenceNames");
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
        if (rootGenomeName == "\"\"") {
            rootGenomeName = alignment->getRootName();
        } 

        // map Sequence pointers to integers (assumes sequence pointers stable within hal)
        vector<string> IDToName;
        unordered_map<string, int64_t> nameToID;
        
        // start up our pinch graph
        stPinchThreadSet* threadSet = stPinchThreadSet_construct();
        
        const Genome* parentGenome = nullptr;
        string parentName;

        deque<string> queue = {rootGenomeName};

        while (!queue.empty()) {
            string genomeName = queue.front();
            queue.pop_front();
            const Genome* genome = alignment->openGenome(genomeName);
            // add the genome sequences as threads
            add_genome_threads(genome, threadSet, IDToName, nameToID, fullNames);

            string curParent = alignment->getParentName(genomeName);
            if (!curParent.empty()) {
                // load up the parent genome if it's not already open, taking care
                // to only ever have one parent open at a time
                if (curParent != parentName) {
                    if (parentGenome != nullptr) {
                        alignment->closeGenome(parentGenome);
                    }
                    parentName = curParent;
                    parentGenome = alignment->openGenome(parentName);
                }

                // pinch the child with its parent
                cerr << "pinching " << genome->getName() << endl;
                pinch_genome(genome, threadSet, nameToID, fullNames);
            }

            // recurse on children                
            vector<string> childs = alignment->getChildNames(genomeName);
            for (size_t i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }

            // todo: this logic not very efficient for normal (ie non-star trees)
            alignment->closeGenome(genome);

        }

        if (parentGenome != nullptr) {
            alignment->closeGenome(parentGenome);
        }

        // clean up the pinch graph
        stPinchThreadSet_joinTrivialBoundaries(threadSet);

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
        queue = {rootGenomeName};
        while (!queue.empty()) {
            string genomeName = queue.front();
            queue.pop_front();
            const Genome* genome = alignment->openGenome(genomeName);

            pinch_to_handle(genome, threadSet, IDToName, nameToID, blockToNode, *graph, fullNames);

            vector<string> childs = alignment->getChildNames(genomeName);
            for (size_t i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }
            alignment->closeGenome(genome);
        }

        stPinchThreadSet_destruct(threadSet);

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
                       unordered_map<string, int64_t>& nameToID,
                       bool fullNames) {
    
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        hal_size_t seqLen = sequence->getSequenceLength();
        string name = fullNames ? sequence->getFullName() : sequence->getName();
        // update lookups to map hal sequence to numeric id
        int64_t seqID = IDToName.size(); 
        nameToID[name] = seqID;
        IDToName.push_back(name);
        // add to thread set
#ifdef debug
        cerr << "Adding sequence " << name << " as thread " << seqID << " with length " << seqLen << endl;
#endif
        stPinchThreadSet_addThread(threads, seqID, 0, seqLen);
    }
}

// Use exact pairwise alginments from genome to its parent to make the pinch graph
void pinch_genome(const Genome* genome,
                  stPinchThreadSet* threads,
                  unordered_map<string, int64_t>& nameToID,
                  bool fullNames) {

    TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIt = genome->getParent()->getBottomSegmentIterator();

    // avoid thread set lookups
    const Sequence* topSeq = nullptr;
    const Sequence* botSeq = nullptr;
    stPinchThread* topThread = nullptr;
    stPinchThread* botThread = nullptr;
    string topString;
    string botString;

    // merge up consecutive segments for fewer pinches
    stPinchThread* prevTopThread = nullptr;
    stPinchThread* prevBotThread = nullptr;
    hal_index_t prevStart1 = -1;
    hal_index_t prevStart2 = -1;
    hal_index_t prevLength = -1;
    bool prevReversed = false;
    
    for (; not topIt->atEnd(); topIt->toRight()) {
        if (topIt->tseg()->hasParent()) {
            botIt->toParent(topIt);

            // todo: lots of string lookups
            int64_t topID = nameToID[fullNames ? topIt->tseg()->getSequence()->getFullName() : topIt->tseg()->getSequence()->getName()];
            int64_t botID = nameToID[fullNames ? botIt->bseg()->getSequence()->getFullName() : botIt->bseg()->getSequence()->getName()];

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
                 << "  " << topString << endl
                 << "   " << *botIt << endl
                 << "  " << botString << endl;
#endif

            int64_t first_match = -1;
            int64_t last_match = -1;
            for (int64_t i = 0; i < (int64_t)topString.length(); ++i) {
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
                             << ", pl2=" << stPinchThread_getLength(botThread) << ", rev=" << botIt->getReversed()
                             << " sp1g=" << (start1 + topSeq->getStartPosition()) << " sp2g=" << (start2 + botSeq->getStartPosition()) << endl
                             << "   " << topString.substr(first_match, length) << endl;
#endif
                        // are we dealing with two consectuive segments? 
                        bool canMerge = topThread == prevTopThread &&
                           botThread == prevBotThread &&
                           start1 == prevStart1 + prevLength &&
                           botIt->getReversed() == prevReversed &&
                           ((!prevReversed && start2 == prevStart2 + prevLength) ||
                            (prevReversed && start2 + length == prevStart2));

                        if (canMerge) {
                            // if consecutive, just merge
                            prevLength += length;
                            if (botIt->getReversed()) {
                                prevStart2 = start2;
                            }
                        } else {
                            // otherwise
                            if (prevTopThread != nullptr) {
                                // pinch the last segment
                                stPinchThread_pinch(prevTopThread,
                                                    prevBotThread,
                                                    prevStart1,
                                                    prevStart2,
                                                    prevLength,
                                                    !prevReversed);
                            }
                            // and update our previous
                            prevTopThread = topThread;
                            prevBotThread = botThread;
                            prevStart1 = start1;
                            prevStart2 = start2;
                            prevLength = length;
                            prevReversed = botIt->getReversed();
                        }
                                                
                    }
                    first_match = -1;
                    last_match = -1;
                }
            }            
        }
    }
    // do that last pinch
    if (prevTopThread != nullptr) {
        stPinchThread_pinch(prevTopThread,
                            prevBotThread,
                            prevStart1,
                            prevStart2,
                            prevLength,
                            !prevReversed);
    }
}

// create nodes and edges for a genome using the pinch graph
void pinch_to_handle(const Genome* genome,
                     stPinchThreadSet* threadSet,
                     const vector<string>& IDToName,
                     const unordered_map<string, int64_t>& nameToID,
                     unordered_map<stPinchBlock*, nid_t>& blockToNode,
                     MutablePathMutableHandleGraph& graph,
                     bool fullNames) {

    // iterate over the sequences of the genome
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        string seqName = fullNames ? sequence->getFullName() : sequence->getName();
        int64_t seqID = nameToID.find(seqName)->second;
        stPinchThread* thread = stPinchThreadSet_getThread(threadSet, seqID);

        // create the path
        path_handle_t pathHandle = graph.create_path_handle(seqName);
        string pathString;
        
        // iterate over the segments of the sequence
        stPinchSegment* prevSeg = nullptr;
        handle_t prevHandle;
        stPinchSegment* lastSeg = stPinchThread_getLast(thread);
        hal_index_t segStart = 0;
        string seqString;
        for (stPinchSegment* seg = stPinchThread_getFirst(thread); ;
             seg = stPinchSegment_get3Prime(seg)) {

            // get the segment's block.  note that if it's not aligned to anything, it will have no block
            stPinchBlock* block = stPinchSegment_getBlock(seg);
            bool reversed = block != nullptr && stPinchSegment_getBlockOrientation(seg) == 0;
            handle_t handle;

            // get the segment's dna sequence from the hal
            sequence->getSubString(seqString, segStart, stPinchSegment_getLength(seg));
            if (reversed) {
                // we always work in block-relative orientation
                reverseComplement(seqString);
            }
            
            // have we already converted this block?
            auto bi = blockToNode.find(block);
            if (bi == blockToNode.end()) {
                // no: it is a new block
                handle = graph.create_handle(seqString);
                if (block != nullptr) {
                    blockToNode[block] = graph.get_id(handle);
                }
#ifdef debug                
                cerr << "created node " << graph.get_id(handle) << " for block " << block << " from " << sequence->getFullName() << " at " << segStart
                     << " rev=" << reversed << " len=" << seqString.length()
                     << endl;
                cerr << "node seq " << graph.get_sequence(handle) << endl;
#endif
            } else {
                // yes: we can find it in the table
                handle = graph.get_handle(bi->second);
#ifdef debug
                cerr << "found node " << graph.get_id(handle) << " for block " << block << " from " << sequence->getFullName() << " at " << segStart
                     << " rev=" << reversed << " len=" << seqString.length()
                     << endl;
                cerr << "node seq " << graph.get_sequence(handle) << endl;
                cerr << "my substring " << seqString << endl;
#endif
            }
            assert(!graph.get_is_reverse(handle));
            if (reversed) {
                handle = graph.flip(handle);
                assert(graph.get_is_reverse(handle));
            }
                   
            // wire up the edge to previous
            if (prevSeg != nullptr) {
#ifdef debug
                cerr << "creating edge from " << graph.get_id(prevHandle) << ":" << graph.get_is_reverse(prevHandle) << " -> "
                     << graph.get_id(handle) << ":" << graph.get_is_reverse(handle) << endl;
#endif
                graph.create_edge(prevHandle, handle);
            }

            // add the node to the path
            graph.append_step(pathHandle, handle);
            pathString += graph.get_sequence(handle);

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
        try {
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
            cerr << "successuflly converted " << seqName << endl;
        }
    }
}
