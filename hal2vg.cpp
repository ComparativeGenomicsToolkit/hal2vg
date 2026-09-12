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
#if defined(__GLIBC__) && !defined(HAVE_JEMALLOC)
#include <malloc.h>
#endif

#include "stPinchGraphs.h"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "hal.h"

using namespace std;
using namespace hal;
using namespace handlegraph;
using namespace bdsg;
using namespace handlegraph;

// std::toupper is locale-aware, and at one call per base it was about 4% of the runtime.
// hal2vg never touches the locale, so in the C locale this is equivalent for every byte value.
static inline char upper_base(char c) {
    return (c >= 'a' && c <= 'z') ? (char)(c - ('a' - 'A')) : c;
}

// The per-base pinch lists of one column.  This used to be a map<char, vector<...>> built and
// thrown away for every column; here the storage is reused between columns, and the entries are
// kept in ascending base order so that the pinches still happen in the order the map gave them.
struct BasePinches {
    typedef vector<tuple<stPinchThread*, hal_index_t, bool>> PinchList;

    vector<char> bases;
    vector<PinchList> lists;
    size_t used = 0;

    void clear() {
        for (size_t i = 0; i < used; ++i) {
            lists[i].clear();
        }
        used = 0;
    }

    PinchList& operator[](char base) {
        size_t i = 0;
        while (i < used && bases[i] < base) {
            ++i;
        }
        if (i < used && bases[i] == base) {
            return lists[i];
        }
        if (used == bases.size()) {
            bases.push_back(0);
            lists.push_back(PinchList());
        }
        // shift the later entries up to keep ascending order, moving the spare list into place
        // rather than copying it; a column holds a handful of distinct bases at most
        for (size_t j = used; j > i; --j) {
            std::swap(bases[j], bases[j - 1]);
            lists[j].swap(lists[j - 1]);
        }
        bases[i] = base;
        lists[i].clear();
        ++used;
        return lists[i];
    }
};

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("halFile", "input hal file");
    optionsParser->addOption("refGenomes",
                             "comma-separated (no spaces) genomes to treat as reference paths with all others as haplotype paths (default=all genomes)",
                             "\"\"");
    optionsParser->addOption("rootGenome", 
                             "process only genomes in clade with specified root"
                             " (HAL root if empty)", 
                             "\"\"");
    optionsParser->addOption("targetGenomes",
                             "comma-separated (no spaces) list of target genomes "
                             "(others are excluded) (all leaves if empty)",
                             "\"\"");
    optionsParser->addOptionFlag("noAncestors", 
                                 "don't write ancestral paths, nor sequence exclusive to ancestral genomes",
                                 false);
    optionsParser->addOption("ignoreGenomes",
                             "comma-separated (no spaces) list of genomes to ignore",
                             "\"\"");
    optionsParser->addOption("outputFormat",
                             "output graph format in {pg, hg} [default=pg]",
                             "pg");
    optionsParser->addOption("chop",
                             "chop up nodes in output graph so they are not longer than given length",
                             0);
    optionsParser->addOptionFlag("progress",
                                 "show progress",
                                 false);
    optionsParser->setDescription("Convert HAL alignment to handle graph");

}

static void add_genome_threads(const Genome* genome,
                               stPinchThreadSet* threads,
                               vector<string>& IDToName,
                               unordered_map<string, int64_t>& nameToID);

static void pinch_genome(const Genome* genome,
                         stPinchThreadSet* threads,
                         const unordered_map<string, int64_t>& nameToID,
                         const vector<string>& targetNames,
                         unordered_map<stPinchThread*, vector<bool>>& snp_cache);

static void pinch_snp(const Genome* genome,
                      stPinchThreadSet* threads,
                      const unordered_map<string, int64_t>& nameToID,
                      const TopSegmentIteratorPtr& topIt,
                      int64_t topOffset,
                      ColumnIteratorPtr& colIt,
                      char topBase,
                      stPinchThread* topThread,
                      unordered_map<stPinchThread*, vector<bool>>& snp_cache,
                      unordered_map<const Sequence*, stPinchThread*>& seqToThread,
                      BasePinches& base_pinches);

// Map a hal sequence to its pinch thread.  Doing this by name, as this used to, built a
// std::string and hashed it for every base of every column.  A sequence pointer is stable for as
// long as its genome is open, which covers a whole pinch_genome call.
static stPinchThread* thread_for_sequence(const Sequence* sequence,
                                          stPinchThreadSet* threads,
                                          const unordered_map<string, int64_t>& nameToID,
                                          unordered_map<const Sequence*, stPinchThread*>& seqToThread);

static void pinch_to_handle(const Genome* genome,
                            stPinchThreadSet* threadSet,
                            const vector<string>& IDToName,
                            const unordered_map<string, int64_t>& nameToID,
                            MutablePathMutableHandleGraph& graph,
                            const vector<string>& refNames);

static void chop_graph(MutablePathMutableHandleGraph& graph, size_t maxNodeLength);

static subrange_t resolve_subpath_naming(string& path_name);

static size_t resolve_haplotype_naming(string& genome_name);

int main(int argc, char** argv) {
    CLParser optionsParser;
    initParser(&optionsParser);
    string halPath;
    string refGenomes;
    string rootGenomeName;
    string targetGenomes;
    bool noAncestors;
    string ignoreGenomes;
    string outputFormat;
    size_t maxNodeLength;
    bool progress;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        refGenomes = optionsParser.getOption<string>("refGenomes");
        rootGenomeName = optionsParser.getOption<string>("rootGenome");
        targetGenomes = optionsParser.getOption<string>("targetGenomes");
        noAncestors = optionsParser.getFlag("noAncestors");
        ignoreGenomes = optionsParser.getOption<string>("ignoreGenomes");
        outputFormat = optionsParser.getOption<string>("outputFormat");
        if (outputFormat != "pg" && outputFormat != "hg") {
            throw hal_exception("--outputFormat must be one of {pg, hg}");
        }
        if (ignoreGenomes != "\"\"" && targetGenomes != "\"\"") {
            throw hal_exception("--ignoreGenomes and --targetGenomes options are "
                                "mutually exclusive");
        }
        
        maxNodeLength = optionsParser.getOption<size_t>("chop");
        progress = optionsParser.getFlag("progress");
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

        vector<string> refNames;
        if (refGenomes != "\"\"") {
            refNames = chopString(refGenomes, ",");
            std::sort(refNames.begin(), refNames.end());
        }

        // default to alignment root if none specified
        bool givenRoot = true;
        if (rootGenomeName == "\"\"") {
            givenRoot = false;
            rootGenomeName = alignment->getRootName();
            const Genome* rootGenome = alignment->openGenome(rootGenomeName);
            if (rootGenome == NULL) {
                throw hal_exception(string("Root genome, ") + rootGenomeName + 
                                    ", not found in alignment");
            }
            alignment->closeGenome(rootGenome);
        }

        vector<string> ignoreNames;
        if (ignoreGenomes != "\"\"") {
            ignoreNames = chopString(ignoreGenomes, ",");
            std::sort(ignoreNames.begin(), ignoreNames.end());
        }

        vector<string> targetNames;
        bool givenTargets;
        if (targetGenomes != "\"\"") {
            // if we're supplied targets, we use them
            targetNames = chopString(targetGenomes, ",");
            givenTargets = true;
        } else {
            // otherwise, we take all the leaves below the root, except any that are ignored
            vector<string> leafNames = alignment->getLeafNamesBelow(rootGenomeName);
            for (size_t i = 0; i < leafNames.size(); ++i) {
                if (!std::binary_search(ignoreNames.begin(), ignoreNames.end(), leafNames[i])) {
                    targetNames.push_back(leafNames[i]);
                }
            }
            givenTargets = false;
        }
        std::sort(targetNames.begin(), targetNames.end());

        // keep track of internal nodes needed to transitively align our targets
        vector<string> spanningNames;        
        set<const Genome*> targetSet;
        for (size_t i = 0; i < targetNames.size(); ++i) {
            const Genome* targetGenome = alignment->openGenome(targetNames[i]);
            if (targetGenome == NULL) {
                throw hal_exception(string("Target genome, ") + targetNames[i] + 
                                    ", not found in alignment");
            }
            targetSet.insert(targetGenome);
        }
        const Genome* rootGenome = getLowestCommonAncestor(targetSet);
        set<const Genome*> targetSetCpy = targetSet;
        getGenomesInSpanningTree(targetSetCpy, targetSet);
        if (!givenRoot) {
            // update our root if it wasn't user-specified
            rootGenomeName = rootGenome->getName();
        }
        for (set<const Genome*>::iterator i = targetSet.begin(); i != targetSet.end(); ++i) {
            if ((*i)->getNumChildren() > 0) {
                spanningNames.push_back((*i)->getName());
            }
            alignment->closeGenome(*i);            
        }
        std::sort(spanningNames.begin(), spanningNames.end());
        
        if (progress) {
            cerr << "Root: " << rootGenomeName << endl;
            if (!targetNames.empty()) {
                cerr << "Targets:";
                for (size_t i = 0; i < targetNames.size(); ++i) {
                    cerr << " " << targetNames[i];
                }
                cerr << endl;
            }
            if (!spanningNames.empty()) {
                cerr << "Spanning:";
                for (size_t i = 0; i < spanningNames.size(); ++i) {
                    cerr << " " << spanningNames[i];
                }
                cerr << endl;
            }
            if (!ignoreNames.empty()) {
                cerr << "Ignore:";
                for (size_t i = 0; i < ignoreNames.size(); ++i) {
                    cerr << " " << ignoreNames[i];
                }
                cerr << endl;
            }
        }

        // map Sequence pointers to integers (assumes sequence pointers stable within hal)
        vector<string> IDToName;
        unordered_map<string, int64_t> nameToID;
        
        // start up our pinch graph
        stPinchThreadSet* threadSet = stPinchThreadSet_construct();
        
        const Genome* parentGenome = nullptr;
        string parentName;

        deque<string> queue = {rootGenomeName};

        vector<string> pinchGenomes;
        
        while (!queue.empty()) {
            string genomeName = queue.front();
            queue.pop_front();

            // we have a target set, and this genome isn't in it, and this genome isn't needed to span it
            // so we can ignore it completely
            bool ignoreGenome = (!std::binary_search(targetNames.begin(), targetNames.end(), genomeName) &&
                                 !std::binary_search(spanningNames.begin(), spanningNames.end(), genomeName) &&
                                 genomeName != rootGenomeName);
            
            const Genome* genome = alignment->openGenome(genomeName);
            string curParent = alignment->getParentName(genomeName);

            // add the genome sequences as threads
            if (!ignoreGenome) {
                if (progress && !(!curParent.empty() && genomeName != rootGenomeName)) {
                    cerr << "adding threads from " << genome->getName() << endl;
                }
                add_genome_threads(genome, threadSet, IDToName, nameToID);
            }

            if (!ignoreGenome && !curParent.empty() && genomeName != rootGenomeName) {
                // load up the parent genome if it's not already open, taking care
                // to only ever have one parent open at a time
                if (curParent != parentName) {
                    if (parentGenome != nullptr) {
                        alignment->closeGenome(parentGenome);
                    }
                    parentName = curParent;
                    parentGenome = alignment->openGenome(parentName);
                }

                // pinching must now be done in second pass, so we queue up the genome here
                pinchGenomes.push_back(genome->getName());
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

        // do all the pinching
        unordered_map<stPinchThread*, vector<bool>> snp_cache;
        for (size_t i = 0; i < pinchGenomes.size(); ++i) {
            
            // pinch the child with its parent
            if (progress) {
                cerr << "pinching " << pinchGenomes[i] << endl;
            }
            pinch_genome(alignment->openGenome(pinchGenomes[i]), threadSet, nameToID, targetNames, snp_cache);
        }
        snp_cache.clear();

        // clean up the pinch graph
        if (progress) {
            cerr << "merging trivial segments and blocks in pinch graph" << endl;
        }
        stPinchThreadSet_joinTrivialBoundaries(threadSet);

        // building the pinch graph leaves the heap littered with small free chunks that the
        // much bigger allocations made below cannot reuse.  consolidate them and give whole
        // free pages back to the OS before switching over.  jemalloc does not suffer from
        // this and has no malloc_trim, so this is only for builds without it
#if defined(__GLIBC__) && !defined(HAVE_JEMALLOC)
        malloc_trim(0);
#endif

        // make a handle graph
        unique_ptr<MutablePathMutableHandleGraph> graph;
        if (outputFormat == "pg") {
            graph = unique_ptr<MutablePathMutableHandleGraph>(new PackedGraph());
        } else if (outputFormat == "hg") {
            graph = unique_ptr<MutablePathMutableHandleGraph>(new HashGraph());
        } else {
            assert(false);
        }

        // keep track of where blocks fit into the handle graph.  a block -> node id hash
        // table would cost some 48 bytes for every node in the output, so the id is parked
        // in the block's supporting-homology count instead, which nothing needs any more now
        // that pinching is done.  zero it first so that 0 means "not converted yet"
        stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
        for (stPinchBlock* block = stPinchThreadSetBlockIt_getNext(&blockIt); block != NULL;
             block = stPinchThreadSetBlockIt_getNext(&blockIt)) {
            stPinchBlock_setNumSupportingHomologies(block, 0);
        }

        // start iterating over the genomes again in order to export to handle graph
        queue = {rootGenomeName};
        while (!queue.empty()) {
            string genomeName = queue.front();
            queue.pop_front();

            // skip it if
            // it's an ancestor and we don't want ancestors or
            // if we have targets and it's not in it or
            // if it's on the ignore list
            bool ignoreGenome = ((noAncestors && !alignment->getChildNames(genomeName).empty()) ||
                                 (givenTargets && !std::binary_search(targetNames.begin(), targetNames.end(), genomeName)) ||
                                 (std::binary_search(ignoreNames.begin(), ignoreNames.end(), genomeName)));
            if (!ignoreGenome) {
                const Genome* genome = alignment->openGenome(genomeName);

                if (progress) {
                    cerr << "converting " << genomeName << " with " << genome->getNumSequences()
                         << " sequences and total length " << genome->getSequenceLength() << endl;
                }
                pinch_to_handle(genome, threadSet, IDToName, nameToID, *graph, refNames);

                alignment->closeGenome(genome);
            }
            
            vector<string> childs = alignment->getChildNames(genomeName);
            for (size_t i = 0; i < childs.size(); ++i) {
                queue.push_back(childs[i]);
            }
        }

        // free the pinch graph
        stPinchThreadSet_destruct(threadSet);

        // free the hal
        alignment = AlignmentConstPtr();

        // chop
        if (maxNodeLength > 0) {
            if (progress) {
                cerr << "chopping graph to max node size " << maxNodeLength << endl;
            }
            chop_graph(*graph, maxNodeLength);
        }

        // write out the graph
        if (progress) {
            cerr << "serializing graph" << endl;
        }
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
                        unordered_map<string, int64_t>& nameToID) {
    
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        hal_size_t seqLen = sequence->getSequenceLength();
        string name = sequence->getFullName();
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
                  const unordered_map<string, int64_t>& nameToID,
                  const vector<string>& targetNames,
                  unordered_map<stPinchThread*, vector<bool>>& snp_cache) {

    TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();
    BottomSegmentIteratorPtr botIt = genome->getParent()->getBottomSegmentIterator();

    // make a target set for column iterator pinching. unfortunately this means
    // opening every single genome
    const Alignment* alignment = genome->getAlignment();
    set<const Genome*> targets;
    for (size_t i = 0; i < targetNames.size(); ++i) {
        targets.insert(alignment->openGenome(targetNames[i]));
    }

    ColumnIteratorPtr colIt = genome->getColumnIterator(&targets);

    // avoid thread set lookups
    const Sequence* topSeq = nullptr;
    const Sequence* botSeq = nullptr;
    stPinchThread* topThread = nullptr;
    stPinchThread* botThread = nullptr;
    // sequence -> thread, so that no base costs a name lookup; and one reused set of per-base
    // pinch lists for the columns visited below
    unordered_map<const Sequence*, stPinchThread*> seqToThread;
    BasePinches base_pinches;
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

            if (topIt->tseg()->getSequence() != topSeq) {
                topSeq = topIt->tseg()->getSequence();
                topThread = thread_for_sequence(topSeq, threads, nameToID, seqToThread);
            }
            if (botIt->bseg()->getSequence() != botSeq) {
                botSeq = botIt->bseg()->getSequence();
                botThread = thread_for_sequence(botSeq, threads, nameToID, seqToThread);
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
                // these were four locale-aware toupper calls per base
                const char topUpper = upper_base(topString[i]);
                const char botUpper = upper_base(botString[i]);
                if (topUpper == botUpper) {
                    if (first_match == -1) {
                        first_match = i;
                    }
                    last_match = i;
                } else if (colIt.get() != NULL) {
                    pinch_snp(genome, threads, nameToID, topIt, i, colIt,
                              topUpper, topThread, snp_cache, seqToThread, base_pinches);
                }
                if (topUpper != botUpper || i == (int64_t)topString.length() - 1) {
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

// Use the column iterator to find all alignments of this snp and pinch accordingly
//
// Todo:  Worried this might be too slow to use at scale.  Also, it blows away all previous
// efforts in hal2vg to be cache-friendly by only loading 2 genomes at a time, so it may
// hog lots of memory.  On a star tree, it may just be better to manually scan the siblings
// before resorting to the column iterator.  Or perhaps do everything in the pinch graph
// by pinching snps then doing a pass over the graph to break them apart once its constructed.
void pinch_snp(const Genome* genome,
               stPinchThreadSet* threads,
               const unordered_map<string, int64_t>& nameToID,
               const TopSegmentIteratorPtr& topIt,
               int64_t topOffset,
               ColumnIteratorPtr& colIt,
               char topBase,
               stPinchThread* topThread,
               unordered_map<stPinchThread*, vector<bool>>& snp_cache,
               unordered_map<const Sequence*, stPinchThread*>& seqToThread,
               BasePinches& base_pinches) {

    const Sequence* topSeq = topIt->tseg()->getSequence();
    hal_index_t topStart = topIt->tseg()->getStartPosition() + topOffset - topSeq->getStartPosition();

    vector<bool>& cache_rec = snp_cache[topThread];
    if (!cache_rec.empty() && cache_rec[topStart] == true) {
        // we've already pinched this base
        return;
    }

    // move the column iterator into position
    colIt->toSite(topStart + topSeq->getStartPosition(), topStart + topSeq->getStartPosition() + 1);

    const ColumnIterator::ColumnMap* columnMap = colIt->getColumnMap();

    // remember all equivalence classes of pinches; the storage is the caller's and is reused
    base_pinches.clear();
    
    // scan through all the homologous bases, breaking them into lists for each possible nucleotide
    for (ColumnIterator::ColumnMap::const_iterator cmi = columnMap->begin(); cmi != columnMap->end(); ++cmi) {
        const Sequence* sequence = cmi->first;
        stPinchThread* otherThread = thread_for_sequence(sequence, threads, nameToID, seqToThread);
        for (ColumnIterator::DNASet::const_iterator dsi = cmi->second->begin(); dsi != cmi->second->end(); ++dsi) {
            char botBase = upper_base((*dsi)->getBase());
            
            hal_index_t otherStart = (*dsi)->getArrayIndex() - sequence->getStartPosition();

            base_pinches[botBase].push_back(make_tuple(otherThread, otherStart, !(*dsi)->getReversed()));
            
        }
    }

    // pinch through each nucleotde
    for (size_t bpi = 0; bpi < base_pinches.used; ++bpi) {
        vector<tuple<stPinchThread*, hal_index_t, bool>>& other_positions = base_pinches.lists[bpi];
        for (size_t i = 0; i < other_positions.size(); ++i) {
            if (i > 0) {
                stPinchThread_pinch(get<0>(other_positions[0]),
                                    get<0>(other_positions[i]),
                                    get<1>(other_positions[0]),
                                    get<1>(other_positions[i]),
                                    1,
                                    get<2>(other_positions[0]) == get<2>(other_positions[i]));
            }
            // update the cache
            vector<bool>& cache_vec = snp_cache[get<0>(other_positions[i])];
            if (cache_vec.empty()) {
                cache_vec.resize(stPinchThread_getLength(get<0>(other_positions[i])), false);
            }
            cache_vec[get<1>(other_positions[i])] = true;
        }
    }
}

static stPinchThread* thread_for_sequence(const Sequence* sequence,
                                          stPinchThreadSet* threads,
                                          const unordered_map<string, int64_t>& nameToID,
                                          unordered_map<const Sequence*, stPinchThread*>& seqToThread) {
    unordered_map<const Sequence*, stPinchThread*>::const_iterator cached = seqToThread.find(sequence);
    if (cached != seqToThread.end()) {
        return cached->second;
    }
    unordered_map<string, int64_t>::const_iterator found = nameToID.find(sequence->getFullName());
    if (found == nameToID.end()) {
        // operator[] used to insert a 0 here and carry on pinching the wrong thread
        throw runtime_error("[hal2vg] no pinch thread for sequence " + sequence->getFullName());
    }
    stPinchThread* thread = stPinchThreadSet_getThread(threads, found->second);
    seqToThread[sequence] = thread;
    return thread;
}

// create nodes and edges for a genome using the pinch graph
void pinch_to_handle(const Genome* genome,
                     stPinchThreadSet* threadSet,
                     const vector<string>& IDToName,
                     const unordered_map<string, int64_t>& nameToID,
                     MutablePathMutableHandleGraph& graph,
                     const vector<string>& refNames) {

    // iterate over the sequences of the genome
    for (SequenceIteratorPtr seqIt = genome->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
        const Sequence *sequence = seqIt->getSequence();
        PathSense sense = PathSense::REFERENCE;
        if (!refNames.empty() && !std::binary_search(refNames.begin(), refNames.end(), genome->getName())) {
            sense = PathSense::HAPLOTYPE;
        }
        int64_t seqID = nameToID.find(sequence->getFullName())->second;
        stPinchThread* thread = stPinchThreadSet_getThread(threadSet, seqID);

        // cactus_graphmap_split can make paths like contig_sub_1_3.  here we convert that
        // into a format vg can (sometimes) understand contig[1-3].
        // (the reason we go through this is that assembly hubs can't handle any special characters apparently)
        string parsed_name = sequence->getName();
        subrange_t subpath = resolve_subpath_naming(parsed_name);
        string parsed_genome_name = genome->getName();
        size_t haplotype = resolve_haplotype_naming(parsed_genome_name);
        if (haplotype == PathMetadata::NO_HAPLOTYPE) {
            haplotype = 0;
        }
        // create the path
        path_handle_t pathHandle = graph.create_path(sense,
                                                     parsed_genome_name,
                                                     parsed_name,
                                                     haplotype,
                                                     sense == PathSense::HAPLOTYPE ? 0 : PathMetadata::NO_PHASE_BLOCK,
                                                     subpath,
                                                     false);
        // the converted path gets checked against the hal a segment at a time (below), so we
        // never need to hold a whole chromosome's worth of sequence in memory to do it
        size_t pathLength = 0;
        size_t numMismatches = 0;
        vector<pair<size_t, pair<char, char>>> mismatches;

        // iterate over the segments of the sequence
        stPinchSegment* prevSeg = nullptr;
        handle_t prevHandle;
        stPinchSegment* lastSeg = stPinchThread_getLast(thread);
        hal_index_t segStart = 0;
        string seqString;
        string nodeString;
        for (stPinchSegment* seg = stPinchThread_getFirst(thread); ;
             seg = stPinchSegment_get3Prime(seg)) {

            // get the segment's block.  note that if it's not aligned to anything, it will have no block
            stPinchBlock* block = stPinchSegment_getBlock(seg);
            bool reversed = block != nullptr && stPinchSegment_getBlockOrientation(seg) == 0;
            handle_t handle;

            // get the segment's dna sequence from the hal.  seqString stays in path
            // orientation; nodeString is the block-relative orientation that gets stored
            sequence->getSubString(seqString, segStart, stPinchSegment_getLength(seg));

            // have we already converted this block?
            nid_t blockNode = block != nullptr ? (nid_t)stPinchBlock_getNumSupportingHomologies(block) : 0;
            if (blockNode == 0) {
                // no: it is a new block
                if (reversed) {
                    // we always work in block-relative orientation
                    nodeString = seqString;
                    reverseComplement(nodeString);
                    handle = graph.create_handle(nodeString);
                } else {
                    handle = graph.create_handle(seqString);
                }
                if (block != nullptr) {
                    assert(graph.get_id(handle) > 0);
                    stPinchBlock_setNumSupportingHomologies(block, (uint64_t)graph.get_id(handle));
                }
#ifdef debug
                cerr << "created node " << graph.get_id(handle) << " for block " << block << " from " << sequence->getFullName() << " at " << segStart
                     << " rev=" << reversed << " len=" << seqString.length()
                     << endl;
                cerr << "node seq " << graph.get_sequence(handle) << endl;
#endif
            } else {
                // yes: the id is stored on the block itself
                handle = graph.get_handle(blockNode);
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

            // make sure what we just appended agrees with the hal.  a node created above is
            // trivially identical to it, so only a block first converted from some other
            // sequence can actually disagree
            if (blockNode != 0) {
                nodeString = graph.get_sequence(handle);
                for (size_t i = 0; i < nodeString.size() && i < seqString.size(); ++i) {
                    if (upper_base(nodeString[i]) != upper_base(seqString[i])) {
                        if (mismatches.size() < 10) {
                            mismatches.push_back(make_pair(segStart + i, make_pair(nodeString[i], seqString[i])));
                        }
                        ++numMismatches;
                    }
                }
            }
            pathLength += seqString.length();

            prevSeg = seg;
            prevHandle = handle;
            
            segStart += stPinchSegment_getLength(seg);
            
            if (seg == lastSeg) {
                break;
            }
        }

        // make sure the path we added is the same as the hal
        if (pathLength != sequence->getSequenceLength()) {
            throw runtime_error("Incorrect length in coverted path for " + sequence->getFullName() + ": " + std::to_string(pathLength) +
                                ". Should be: " + std::to_string(sequence->getSequenceLength()));
        }
        if (numMismatches > 0) {
            stringstream msg;
            msg << numMismatches << " mismatches found in converted path for " << sequence->getFullName() << ":\n";
            for (size_t i = 0; i < mismatches.size(); ++i) {
                msg << " path[" << mismatches[i].first << "]=" << mismatches[i].second.first
                    << ". should be " << mismatches[i].second.second << "\n";
            }
            throw runtime_error(msg.str());
        }
    }
}

void chop_graph(MutablePathMutableHandleGraph& graph, size_t maxNodeLength) {
    // borrowed from https://github.com/vgteam/odgi/blob/master/src/subcommand/chop_main.cpp
    std::vector<handle_t> to_chop;
    graph.for_each_handle([&](const handle_t& handle) {
            if (graph.get_length(handle) > maxNodeLength) {
                to_chop.push_back(handle);
            }
        });

    for (auto& handle : to_chop) {
        // get divide points
        uint64_t length = graph.get_length(handle);
        std::vector<size_t> offsets;
        for (uint64_t i = maxNodeLength; i < length; i+=maxNodeLength) {
            offsets.push_back(i);
        }
        graph.divide_handle(handle, offsets);
    }
}

subrange_t resolve_subpath_naming(string& path_name) {
    size_t first_length = 0;
    size_t start_offset = 0;
    bool found_subpath = false;
    while (true) {
        size_t sp = path_name.rfind("_sub_");
        if (sp != string::npos) {
            size_t up = path_name.rfind("_");
            if (up != string::npos && up > sp + 1) {
                int64_t start;
                int64_t end;
                start = stol(path_name.substr(sp + 5, up - sp - 5));
                end = stol(path_name.substr(up + 1));
                stringstream new_name;
                start_offset += start; // final offset is sum of all nested offsets
                if (first_length == 0) {
                    first_length = end - start;
                    assert(first_length > 0);
                } else {
                    // in the case of nested subpaths, the end coordinate will always
                    // be derived from the start, plus the length of the "top" path
                    end = start_offset + first_length;
                }
                new_name << path_name.substr(0, sp);
                path_name = new_name.str();
                found_subpath = true;
            }
        } else {
            break;
        }
    }
    if (found_subpath) {
        return make_pair(start_offset, start_offset + first_length);
    } else {
        return PathMetadata::NO_SUBRANGE;
    }
}

size_t resolve_haplotype_naming(string& genome_name) {
    size_t haplotype = PathMetadata::NO_HAPLOTYPE;
    size_t dp = genome_name.rfind(".");
    if (dp != string::npos) {
        try {
            haplotype = stol(genome_name.substr(dp + 1));
            genome_name = genome_name.substr(0, dp);
        } catch(...) {
        }
    }
    return haplotype;
}
