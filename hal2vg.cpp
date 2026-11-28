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
#include <functional>
#include <unordered_map>
#include <sys/stat.h>
#include "low_memory_plan.hpp"
#if defined(__GLIBC__) && !defined(HAVE_JEMALLOC)
#include <malloc.h>
#endif

#include "stPinchGraphs.h"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "hal.h"
#include "hal3/hal3.h"

using namespace std;
using namespace hal;
using namespace handlegraph;
using namespace bdsg;
using namespace handlegraph;

struct LowMemoryBlock {
    vector<stPinchSegment*> members;
    hal2vg_low_memory::Plan plan;
};
using LowMemoryPlans = unordered_map<stPinchBlock*, unique_ptr<LowMemoryBlock>>;

// std::toupper is locale-aware, and at one call per base it was about 4% of the runtime.
// hal2vg never touches the locale, so in the C locale this is equivalent for every byte value.
static inline char upper_base(char c) {
    return (c >= 'a' && c <= 'z') ? (char)(c - ('a' - 'A')) : c;
}

// A block that has already been converted parks its node id in its supporting-homology
// count, which nothing needs once pinching is done.  A checksum of the node's sequence is
// parked in the same field, in the bits above the id: that is what lets the per-sequence
// check in pinch_to_handle confirm that the node really does spell what the hal says
// without reading the node's bases back out of the graph, which it otherwise does once per
// genome that touches the block.  The id gets the low 40 bits; a graph too big for that is
// refused rather than silently truncated.
static const unsigned BLOCK_ID_BITS = 40;
static const uint64_t BLOCK_ID_LIMIT = (uint64_t)1 << BLOCK_ID_BITS;
static const uint64_t BLOCK_ID_MASK = BLOCK_ID_LIMIT - 1;

static inline uint64_t pack_block_node(uint64_t node_id, uint32_t checksum) {
    return node_id | ((uint64_t)checksum << BLOCK_ID_BITS);
}
static inline uint64_t unpack_block_id(uint64_t packed) {
    return packed & BLOCK_ID_MASK;
}
static inline uint32_t unpack_block_checksum(uint64_t packed) {
    return (uint32_t)(packed >> BLOCK_ID_BITS);
}

// FNV-1a over the case-folded bases, truncated to the 22 bits that are free above the id.
// Case is folded because the check it stands in for compared bases case-insensitively.
static inline uint32_t sequence_checksum(const string& seq) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < seq.size(); ++i) {
        h = (h ^ (unsigned char)upper_base(seq[i])) * 16777619u;
    }
    // fold the discarded bits back in rather than dropping them
    return ((h >> 22) ^ h) & 0x3fffffu;
}

// One run of the branch above a genome: [childStart, childStart + length) of the genome's
// sequence childSequence aligned, gapless and colinear, to [parentStart, parentStart + length)
// of its parent's sequence parentSequence.  reversed pairs childStart + i with
// parentStart + length - 1 - i.  Positions are in the sequences, not the genomes.
struct AlignedRun {
    uint32_t childSequence;
    uint32_t parentSequence;
    int64_t childStart;
    int64_t parentStart;
    int64_t length;
    bool reversed;
};

// a sequence: its genome's number in the AlignmentReader, and its own index in that genome
struct SequenceRef {
    uint32_t genome;
    uint32_t sequence;
};

// What hal2vg reads of an alignment -- the tree, the genomes' sequences, the runs of each
// branch and the bases -- from either of two backends: HAL's own API, which reads HDF5 and
// mmap files, or HAL format 3 (a pack file, or a directory holding a manifest).  Genomes are
// numbered in the order the backend lists them, the root first and every genome after its
// parent; a sequence is known by its index in its genome.
class AlignmentReader {
public:
    virtual ~AlignmentReader() {}

    size_t num_genomes() const { return names.size(); }
    const string& genome_name(size_t genome) const { return names.at(genome); }
    // -1 if there is no such genome
    int64_t find_genome(const string& name) const {
        unordered_map<string, size_t>::const_iterator found = index.find(name);
        return found == index.end() ? -1 : (int64_t)found->second;
    }
    size_t root() const { return 0; }
    // -1 for the root
    int64_t parent(size_t genome) const { return parents.at(genome); }
    const vector<size_t>& children(size_t genome) const { return childLists.at(genome); }

    virtual size_t num_sequences(size_t genome) = 0;
    virtual string sequence_name(size_t genome, size_t sequence) = 0;
    virtual int64_t sequence_length(size_t genome, size_t sequence) = 0;
    // every run of the branch above a genome other than the root, in the order of the
    // genome's sequences and, within one, of position
    virtual void scan_runs(size_t genome, const function<void(const AlignedRun&)>& visit) = 0;
    // bases [start, start + length) of a sequence, forward strand, soft-masked ones in lower case
    virtual void fetch(size_t genome, size_t sequence, int64_t start, int64_t length, string& out) = 0;
    // nothing of the genome is wanted again soon: give back what is held for it
    virtual void release(size_t genome) {}

protected:
    // the root first, then every genome after its parent, a genome's children in their order
    void add_genome(const string& name, int64_t parentGenome) {
        if (parentGenome < 0 ? !names.empty() : (size_t)parentGenome >= names.size()) {
            throw runtime_error("[hal2vg] genome " + name + " is listed before its parent");
        }
        if (!index.insert(make_pair(name, names.size())).second) {
            throw runtime_error("[hal2vg] two genomes are named " + name);
        }
        if (parentGenome >= 0) {
            childLists[parentGenome].push_back(names.size());
        }
        names.push_back(name);
        parents.push_back(parentGenome);
        childLists.emplace_back();
    }

private:
    vector<string> names;
    vector<int64_t> parents;
    vector<vector<size_t>> childLists;
    unordered_map<string, size_t> index;
};

// HAL's API: HDF5 and mmap files.  A genome is opened when something of it is first asked for,
// and closed by release(), so that only the genomes in use hold their HDF5 caches.
class HalReader : public AlignmentReader {
public:
    HalReader(const string& path, CLParser* options) : alignment(openHalAlignment(path, options)) {
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("input hal alignmenet is empty");
        }
        // breadth first, which lists every genome after its parent
        deque<pair<string, int64_t>> queue = {make_pair(alignment->getRootName(), (int64_t)-1)};
        while (!queue.empty()) {
            pair<string, int64_t> next = queue.front();
            queue.pop_front();
            const int64_t genome = num_genomes();
            add_genome(next.first, next.second);
            for (const string& child : alignment->getChildNames(next.first)) {
                queue.push_back(make_pair(child, genome));
            }
        }
        open.assign(num_genomes(), nullptr);
        sequences.resize(num_genomes());
    }

    size_t num_sequences(size_t genome) override {
        return sequences_of(genome).size();
    }
    string sequence_name(size_t genome, size_t sequence) override {
        return sequences_of(genome).at(sequence)->getName();
    }
    int64_t sequence_length(size_t genome, size_t sequence) override {
        return sequences_of(genome).at(sequence)->getSequenceLength();
    }
    void fetch(size_t genome, size_t sequence, int64_t start, int64_t length, string& out) override {
        sequences_of(genome).at(sequence)->getSubString(out, start, length);
    }
    void release(size_t genome) override {
        if (open.at(genome) != nullptr) {
            alignment->closeGenome(open[genome]);
            open[genome] = nullptr;
            vector<const Sequence*>().swap(sequences[genome]);
        }
    }

    void scan_runs(size_t genome, const function<void(const AlignedRun&)>& visit) override {
        const int64_t parentGenome = parent(genome);
        if (parentGenome < 0) {
            return;
        }
        const Genome* child = genome_of(genome);
        // opened here rather than by getParent() below, so that release() knows to close it
        genome_of(parentGenome);
        TopSegmentIteratorPtr topIt = child->getTopSegmentIterator();
        BottomSegmentIteratorPtr botIt = child->getParent()->getBottomSegmentIterator();

        // consecutive segments are merged, for fewer pinches
        AlignedRun run;
        bool haveRun = false;
        for (; not topIt->atEnd(); topIt->toRight()) {
            if (!topIt->tseg()->hasParent()) {
                continue;
            }
            botIt->toParent(topIt);
            const Sequence* topSeq = topIt->tseg()->getSequence();
            const Sequence* botSeq = botIt->bseg()->getSequence();
            const bool reversed = botIt->getReversed();
            const int64_t length = topIt->getLength();
            const int64_t start1 = topIt->tseg()->getStartPosition() - topSeq->getStartPosition();
            const int64_t start2 = !reversed ? botIt->bseg()->getStartPosition() - botSeq->getStartPosition()
                : botIt->bseg()->getEndPosition() - length + 1 - botSeq->getStartPosition();
            const uint32_t childSequence = topSeq->getArrayIndex();
            const uint32_t parentSequence = botSeq->getArrayIndex();

            // are we dealing with two consecutive segments?
            bool canMerge = haveRun &&
                childSequence == run.childSequence &&
                parentSequence == run.parentSequence &&
                start1 == run.childStart + run.length &&
                reversed == run.reversed &&
                ((!reversed && start2 == run.parentStart + run.length) ||
                 (reversed && start2 + length == run.parentStart));
            if (canMerge) {
                run.length += length;
                if (reversed) {
                    run.parentStart = start2;
                }
            } else {
                if (haveRun) {
                    visit(run);
                }
                run = AlignedRun{childSequence, parentSequence, start1, start2, length, reversed};
                haveRun = true;
            }
        }
        if (haveRun) {
            visit(run);
        }
    }

private:
    const Genome* genome_of(size_t genome) {
        if (open.at(genome) == nullptr) {
            open[genome] = alignment->openGenome(genome_name(genome));
            if (open[genome] == nullptr) {
                throw hal_exception("genome " + genome_name(genome) + " not found in alignment");
            }
        }
        return open[genome];
    }
    // the genome's sequences, in order: a Sequence's index in its genome is its array index
    const vector<const Sequence*>& sequences_of(size_t genome) {
        const Genome* g = genome_of(genome);
        vector<const Sequence*>& list = sequences[genome];
        if (list.size() != g->getNumSequences()) {
            list.clear();
            for (SequenceIteratorPtr seqIt = g->getSequenceIterator(); not seqIt->atEnd(); seqIt->toNext()) {
                // not the iterator's pointer: an mmap hal's iterator hands out its own member,
                // which is gone once the iterator moves on.  the genome's own object lasts as
                // long as the genome is open, in either format
                const Sequence* sequence = g->getSequence(seqIt->getSequence()->getName());
                if (sequence == nullptr || sequence->getArrayIndex() != (hal_index_t)list.size()) {
                    throw runtime_error("[hal2vg] sequences of " + genome_name(genome) + " are not in array order");
                }
                list.push_back(sequence);
            }
        }
        return list;
    }

    AlignmentConstPtr alignment;
    vector<const Genome*> open;
    vector<vector<const Sequence*>> sequences;
};

// HAL format 3, through its API.  Nothing needs releasing but the DNA reader's block lists.
class Format3Reader : public AlignmentReader {
public:
    explicit Format3Reader(const string& path) :
        alignment(hal3::api::Alignment::open(path)), branches(alignment), dna(alignment) {
        // the root first, then a preorder, which lists every genome after its parent and a
        // genome's children in their order
        unordered_map<const hal3::api::Genome*, size_t> number;
        for (const hal3::api::Genome* genome : alignment->genomes()) {
            const int64_t parentGenome = genome->parent() == nullptr ? -1 : (int64_t)number.at(genome->parent());
            number[genome] = genomes.size();
            genomes.push_back(genome);
            add_genome(genome->name(), parentGenome);
        }
    }

    size_t num_sequences(size_t genome) override {
        return genomes.at(genome)->numSequences();
    }
    string sequence_name(size_t genome, size_t sequence) override {
        return genomes.at(genome)->sequenceName(sequence);
    }
    int64_t sequence_length(size_t genome, size_t sequence) override {
        return genomes.at(genome)->sequenceLength(sequence);
    }
    void fetch(size_t genome, size_t sequence, int64_t start, int64_t length, string& out) override {
        dna.fetch(genomes.at(genome), sequence, start, start + length, out);
    }
    void release(size_t genome) override {
        dna.recycle();
    }
    void scan_runs(size_t genome, const function<void(const AlignedRun&)>& visit) override {
        if (parent(genome) < 0) {
            return;
        }
        // runs are maximal already: there is nothing to merge
        branches.scan(genomes.at(genome), [&](const hal3::api::Run& run) {
            visit(AlignedRun{run.childSeq, run.parentSeq, run.childStart, run.parentStart, run.length, run.reversed});
        });
    }

private:
    shared_ptr<const hal3::api::Alignment> alignment;
    hal3::api::BranchReader branches;
    hal3::api::DnaReader dna;
    vector<const hal3::api::Genome*> genomes;
};

static bool is_url(const string& path) {
    return path.compare(0, 7, "http://") == 0 || path.compare(0, 8, "https://") == 0;
}

// A local alignment's format is told from the path: a directory is format 3, and a file by its
// first bytes, those of a format 3 pack or manifest or else HAL's.  A URL is tried as format 3
// first.
static unique_ptr<AlignmentReader> open_alignment(const string& path, CLParser* options) {
    if (is_url(path)) {
        string format3Error;
        try {
            return unique_ptr<AlignmentReader>(new Format3Reader(path));
        } catch (const exception& e) {
            format3Error = e.what();
        }
        try {
            return unique_ptr<AlignmentReader>(new HalReader(path, options));
        } catch (const exception& e) {
            throw runtime_error("cannot open " + path + " as HAL format 3 (" + format3Error +
                                ") nor with HAL's API (" + e.what() + ")");
        }
    }
    static const char packMagic[8] = {'\x89', 'H', 'A', 'L', '\r', '\n', '\x1a', '\n'};
    static const string manifestMagic = "#hal-manifest 3";
    bool format3 = false;
    struct stat info;
    if (stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode)) {
        format3 = true;
    } else {
        ifstream file(path, ios::binary);
        char head[16];
        file.read(head, sizeof(head));
        const string initial(head, (size_t)file.gcount());
        format3 = initial.compare(0, sizeof(packMagic), string(packMagic, sizeof(packMagic))) == 0 ||
            initial.compare(0, manifestMagic.size(), manifestMagic) == 0;
    }
    if (format3) {
        return unique_ptr<AlignmentReader>(new Format3Reader(path));
    }
    return unique_ptr<AlignmentReader>(new HalReader(path, options));
}

// the genomes that are leaves below a genome, in no particular order
static void leaves_below(const AlignmentReader& alignment, size_t genome, vector<string>& leaves) {
    if (alignment.children(genome).empty()) {
        leaves.push_back(alignment.genome_name(genome));
    }
    for (size_t child : alignment.children(genome)) {
        leaves_below(alignment, child, leaves);
    }
}

// the deepest genome that is an ancestor (or one) of every genome given, which must not be none
static size_t lowest_common_ancestor(const AlignmentReader& alignment, const vector<size_t>& genomes) {
    // the first genome's path to the root, from the genome up
    vector<size_t> path;
    unordered_map<size_t, size_t> onPath;
    for (int64_t genome = genomes.at(0); genome >= 0; genome = alignment.parent(genome)) {
        onPath[genome] = path.size();
        path.push_back(genome);
    }
    // each other genome meets that path somewhere: the highest meeting is the answer
    size_t highest = 0;
    for (size_t i = 1; i < genomes.size(); ++i) {
        int64_t genome = genomes[i];
        while (!onPath.count(genome)) {
            genome = alignment.parent(genome);
        }
        highest = max(highest, onPath[genome]);
    }
    return path[highest];
}

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
    optionsParser->addOptionFlag("legacyConversion",
                                 "split the pinch graph wherever a block's members disagree, as hal2vg used to. "
                                 "needs much more memory, and makes the same graph up to node ids, node "
                                 "orientations and unchopping",
                                 false);
    optionsParser->addOption("baseBuffer",
                             "bytes of DNA read at a time to find where a block's members disagree "
                             "(not a limit on memory use)",
                             size_t(64 * 1024 * 1024));
    optionsParser->setDescription("Convert HAL alignment to handle graph");

}

// Every pinch thread is a sequence of a genome in the pinch graph.  The threads are numbered
// genome by genome, so genome g's sequence s is thread threadOffset[g] + s (threadOffset[g] is -1
// for a genome with no threads), and that thread is threads[threadOffset[g] + s].
// threadToSequence says which sequence a thread number is.
struct PinchThreads {
    vector<int64_t> threadOffset;
    vector<SequenceRef> threadToSequence;
    vector<stPinchThread*> threads;

    stPinchThread* thread(size_t genome, size_t sequence) const {
        return threads.at(threadOffset.at(genome) + sequence);
    }
    const SequenceRef& sequence(int64_t thread) const {
        if (thread < 0 || (uint64_t)thread >= threadToSequence.size()) {
            throw runtime_error("[hal2vg] no hal sequence for pinch thread " + std::to_string(thread));
        }
        return threadToSequence[thread];
    }
};

static void add_genome_threads(AlignmentReader& alignment, size_t genome,
                               stPinchThreadSet* threadSet,
                               PinchThreads& threads);

static void pinch_genome(AlignmentReader& alignment, size_t genome,
                         const PinchThreads& threads);

// Pinching runs straight through mismatching bases, so a block can hold more than one base at
// a position.  This separates them again: every block whose segments disagree anywhere is cut
// at the positions where the grouping changes, and each such piece is repinched into one block
// per base.  It is what the column iterator used to do one SNP at a time, without needing a
// column iterator, a set of every genome open at once, or a bit per base of every sequence.
static void split_blocks_by_base(stPinchThreadSet* threadSet,
                                 AlignmentReader& alignment,
                                 const PinchThreads& threads,
                                 bool progress);

// both conversions make their paths here, so that they cannot disagree on how one is named
static path_handle_t create_sequence_path(const string& genomeName,
                                          const string& sequenceName,
                                          const vector<string>& refNames,
                                          MutablePathMutableHandleGraph& graph);

static void pinch_to_handle(AlignmentReader& alignment, size_t genome,
                            const PinchThreads& threads,
                            MutablePathMutableHandleGraph& graph,
                            const vector<string>& refNames);

static void build_low_memory_plans(stPinchThreadSet* threadSet,
                                   AlignmentReader& alignment,
                                   const PinchThreads& threads,
                                   LowMemoryPlans& plans, bool progress, size_t bufferBytes);

static void pinch_to_handle_low_memory(AlignmentReader& alignment, size_t genome,
                                        const PinchThreads& threads,
                                        LowMemoryPlans& plans,
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
    bool legacyConversion;
    size_t baseBuffer;
    bool halInMemory;
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
        legacyConversion = optionsParser.getFlag("legacyConversion");
        baseBuffer = optionsParser.getOption<size_t>("baseBuffer");
        if (baseBuffer == 0) {
            throw hal_exception("--baseBuffer must be positive");
        }
        // hal reads these itself; getFlagAlt would warn about --inMemory a second time
        halInMemory = optionsParser.getFlag("hdf5InMemory") || optionsParser.getFlag("inMemory");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
    try {
        unique_ptr<AlignmentReader> alignment = open_alignment(halPath, &optionsParser);

        vector<string> refNames;
        if (refGenomes != "\"\"") {
            refNames = chopString(refGenomes, ",");
            std::sort(refNames.begin(), refNames.end());
        }

        // default to alignment root if none specified
        bool givenRoot = true;
        if (rootGenomeName == "\"\"") {
            givenRoot = false;
            rootGenomeName = alignment->genome_name(alignment->root());
        } else if (alignment->find_genome(rootGenomeName) < 0) {
            throw hal_exception(string("Root genome, ") + rootGenomeName +
                                ", not found in alignment");
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
            vector<string> leafNames;
            leaves_below(*alignment, alignment->find_genome(rootGenomeName), leafNames);
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
        vector<size_t> targets;
        for (size_t i = 0; i < targetNames.size(); ++i) {
            const int64_t targetGenome = alignment->find_genome(targetNames[i]);
            if (targetGenome < 0) {
                throw hal_exception(string("Target genome, ") + targetNames[i] +
                                    ", not found in alignment");
            }
            targets.push_back(targetGenome);
        }
        if (targets.empty()) {
            throw hal_exception("no target genomes: every leaf below " + rootGenomeName + " is ignored");
        }
        const size_t lca = lowest_common_ancestor(*alignment, targets);
        if (!givenRoot) {
            // update our root if it wasn't user-specified
            rootGenomeName = alignment->genome_name(lca);
        }
        // the spanning tree: every genome on the way from a target up to their common ancestor
        set<size_t> spanning;
        for (size_t target : targets) {
            for (size_t genome = target; spanning.insert(genome).second && genome != lca;
                 genome = alignment->parent(genome)) {
            }
        }
        for (size_t genome : spanning) {
            if (!alignment->children(genome).empty()) {
                spanningNames.push_back(alignment->genome_name(genome));
            }
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

        // start up our pinch graph
        stPinchThreadSet* threadSet = stPinchThreadSet_construct();
        PinchThreads threads;
        threads.threadOffset.assign(alignment->num_genomes(), -1);

        const size_t rootGenome = alignment->find_genome(rootGenomeName);
        deque<size_t> queue = {rootGenome};

        vector<size_t> pinchGenomes;
        vector<size_t> threadGenomes;

        while (!queue.empty()) {
            const size_t genome = queue.front();
            queue.pop_front();
            const string& genomeName = alignment->genome_name(genome);
            const bool hasParent = alignment->parent(genome) >= 0;

            // we have a target set, and this genome isn't in it, and this genome isn't needed to span it
            // so we can ignore it completely
            bool ignoreGenome = (!std::binary_search(targetNames.begin(), targetNames.end(), genomeName) &&
                                 !std::binary_search(spanningNames.begin(), spanningNames.end(), genomeName) &&
                                 genome != rootGenome);

            // add the genome sequences as threads
            if (!ignoreGenome) {
                if (progress && !(hasParent && genome != rootGenome)) {
                    cerr << "adding threads from " << genomeName << endl;
                }
                add_genome_threads(*alignment, genome, threadSet, threads);
                threadGenomes.push_back(genome);
            }

            if (!ignoreGenome && hasParent && genome != rootGenome) {
                // pinching must now be done in second pass, so we queue up the genome here
                pinchGenomes.push_back(genome);
            }

            // recurse on children
            for (size_t child : alignment->children(genome)) {
                queue.push_back(child);
            }

            alignment->release(genome);
        }

        // a thread number is looked up once here, and never while pinching
        threads.threads.resize(threads.threadToSequence.size());
        for (size_t i = 0; i < threads.threads.size(); ++i) {
            threads.threads[i] = stPinchThreadSet_getThread(threadSet, i);
        }

        // do all the pinching
        for (size_t genome : pinchGenomes) {

            // pinch the child with its parent
            if (progress) {
                cerr << "pinching " << alignment->genome_name(genome) << endl;
            }
            pinch_genome(*alignment, genome, threads);
            if (!legacyConversion && !halInMemory) {
                // nothing reads it again until the splits are planned, which reopens it.
                // closing it gives back its hdf5 chunk caches, but a genome held in memory
                // would only be read in again whole, for no saving: its arrays are back by
                // the time the splits are planned, which is where memory peaks
                alignment->release(genome);
            }
        }

        // clean up the pinch graph
        if (progress) {
            cerr << "merging trivial segments and blocks in pinch graph" << endl;
        }
        stPinchThreadSet_joinTrivialBoundaries(threadSet);

        // the pinching above ran straight through mismatching bases, so blocks can hold more
        // than one base at a position, and the graph needs a node for each.  work out where
        // every such block divides, and how, without touching the pinch graph: the nodes are
        // made from that plan as the graph is built.  --legacyConversion instead separates
        // them in the pinch graph and then merges again, which costs a pinch segment for
        // every member of every piece and so far more memory
        LowMemoryPlans lowMemoryPlans;
        if (legacyConversion) {
            split_blocks_by_base(threadSet, *alignment, threads, progress);
            stPinchThreadSet_joinTrivialBoundaries(threadSet);
        } else {
            build_low_memory_plans(threadSet, *alignment, threads, lowMemoryPlans, progress, baseBuffer);
            // Plans own no HAL pointers. Release the cohort's arrays now; conversion
            // reopens only its current genome and reads one sequence at a time.
            for (size_t genome : threadGenomes) {
                alignment->release(genome);
            }
        }

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
        queue = {rootGenome};
        while (!queue.empty()) {
            const size_t genome = queue.front();
            queue.pop_front();
            const string& genomeName = alignment->genome_name(genome);

            // skip it if
            // it's an ancestor and we don't want ancestors or
            // if we have targets and it's not in it or
            // if it's on the ignore list or
            // it has no threads, like an ancestor whose descendants are all ignored, which is
            // not needed to span the targets: there is nothing in the pinch graph to convert
            bool ignoreGenome = ((noAncestors && !alignment->children(genome).empty()) ||
                                 (givenTargets && !std::binary_search(targetNames.begin(), targetNames.end(), genomeName)) ||
                                 (std::binary_search(ignoreNames.begin(), ignoreNames.end(), genomeName)) ||
                                 threads.threadOffset[genome] < 0);
            if (!ignoreGenome) {
                if (progress) {
                    const size_t numSequences = alignment->num_sequences(genome);
                    int64_t totalLength = 0;
                    for (size_t i = 0; i < numSequences; ++i) {
                        totalLength += alignment->sequence_length(genome, i);
                    }
                    cerr << "converting " << genomeName << " with " << numSequences
                         << " sequences and total length " << totalLength << endl;
                }
                if (legacyConversion) {
                    pinch_to_handle(*alignment, genome, threads, *graph, refNames);
                } else {
                    pinch_to_handle_low_memory(*alignment, genome, threads, lowMemoryPlans, *graph, refNames);
                }

                alignment->release(genome);
            }

            for (size_t child : alignment->children(genome)) {
                queue.push_back(child);
            }
        }

        lowMemoryPlans.clear();

        // free the pinch graph
        stPinchThreadSet_destruct(threadSet);

        // free the hal
        alignment.reset();

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

// Add every sequence from the genome into the pinch graph, numbered on from the threads there
void add_genome_threads(AlignmentReader& alignment, size_t genome,
                        stPinchThreadSet* threadSet,
                        PinchThreads& threads) {
    threads.threadOffset.at(genome) = threads.threadToSequence.size();
    const size_t numSequences = alignment.num_sequences(genome);
    for (size_t i = 0; i < numSequences; ++i) {
        int64_t seqID = threads.threadToSequence.size();
        int64_t seqLen = alignment.sequence_length(genome, i);
        threads.threadToSequence.push_back(SequenceRef{(uint32_t)genome, (uint32_t)i});
        // add to thread set
#ifdef debug
        cerr << "Adding sequence " << alignment.sequence_name(genome, i) << " as thread " << seqID << " with length " << seqLen << endl;
#endif
        stPinchThreadSet_addThread(threadSet, seqID, 0, seqLen);
    }
}

// Use exact pairwise alginments from genome to its parent to make the pinch graph
void pinch_genome(AlignmentReader& alignment, size_t genome,
                  const PinchThreads& threads) {
    const int64_t parent = alignment.parent(genome);
    if (parent < 0 || threads.threadOffset.at(genome) < 0 || threads.threadOffset.at(parent) < 0) {
        throw runtime_error("[hal2vg] no pinch threads for " + alignment.genome_name(genome) +
                            (parent < 0 ? string() : " or its parent " + alignment.genome_name(parent)));
    }
    stPinchThread* const* childThreads = threads.threads.data() + threads.threadOffset[genome];
    stPinchThread* const* parentThreads = threads.threads.data() + threads.threadOffset[parent];
    const size_t numChild = alignment.num_sequences(genome);
    const size_t numParent = alignment.num_sequences(parent);

    alignment.scan_runs(genome, [&](const AlignedRun& run) {
        if (run.childSequence >= numChild || run.parentSequence >= numParent) {
            throw runtime_error("[hal2vg] a run of " + alignment.genome_name(genome) + " is on no sequence");
        }
        // the whole aligned run is pinched, mismatching bases included; blocks holding more
        // than one base are separated afterwards.  not comparing the bases means neither
        // sequence has to be read here at all
        stPinchThread_pinch(childThreads[run.childSequence],
                            parentThreads[run.parentSequence],
                            run.childStart,
                            run.parentStart,
                            run.length,
                            !run.reversed);
    });
}

// one segment of the block being examined
namespace {
struct BlockMember {
    stPinchThread* thread;
    int64_t start;    // of the segment, in thread coordinates
    bool forward;     // the segment's orientation within the block
    string bases;     // block-relative and upper case, so the members line up position by position
};

// Number the members of a block at one position by which base they carry, first base seen
// getting group 0.  Two positions belong in the same block exactly when this is the same at
// both, whatever the bases themselves are.
inline void grouping_at(const vector<BlockMember>& members, int64_t position, vector<uint32_t>& key,
                        uint32_t* base_stamp, uint32_t* base_group, uint32_t& stamp) {
    ++stamp;
    key.resize(members.size());
    uint32_t next_group = 0;
    for (size_t i = 0; i < members.size(); ++i) {
        unsigned char base = (unsigned char)members[i].bases[position];
        if (base_stamp[base] != stamp) {
            base_stamp[base] = stamp;
            base_group[base] = next_group++;
        }
        key[i] = base_group[base];
    }
}
}

void split_blocks_by_base(stPinchThreadSet* threadSet,
                          AlignmentReader& alignment,
                          const PinchThreads& threads,
                          bool progress) {

    vector<BlockMember> members;
    vector<uint32_t> key, next_key;
    vector<int64_t> group_first;
    string buffer;
    uint32_t base_stamp[256] = {0};
    uint32_t base_group[256] = {0};
    uint32_t stamp = 0;
    size_t blocks_split = 0;

    // the block of a segment is visited when its first segment is reached, which is how the
    // block iterator does it.  the repinching below splits segments to the right of the one
    // in hand and makes new blocks, all of which agree by construction: reaching one of those
    // later costs a second read of its bases and nothing else
    stPinchThreadSetSegmentIt segIt = stPinchThreadSet_getSegmentIt(threadSet);
    for (stPinchSegment* seg = stPinchThreadSetSegmentIt_getNext(&segIt); seg != NULL;
         seg = stPinchThreadSetSegmentIt_getNext(&segIt)) {

        stPinchBlock* block = stPinchSegment_getBlock(seg);
        if (block == NULL || stPinchBlock_getFirst(block) != seg || stPinchBlock_getDegree(block) < 2) {
            continue;
        }
        const int64_t length = stPinchBlock_getLength(block);

        // read what every member of the block says, in the block's orientation
        members.clear();
        stPinchBlockIt blockIt = stPinchBlock_getSegmentIterator(block);
        for (stPinchSegment* member = stPinchBlockIt_getNext(&blockIt); member != NULL;
             member = stPinchBlockIt_getNext(&blockIt)) {
            const SequenceRef& sequence = threads.sequence(stPinchSegment_getName(member));
            BlockMember entry;
            entry.thread = stPinchSegment_getThread(member);
            entry.start = stPinchSegment_getStart(member);
            entry.forward = stPinchSegment_getBlockOrientation(member) != 0;
            alignment.fetch(sequence.genome, sequence.sequence, entry.start, length, buffer);
            for (size_t i = 0; i < buffer.size(); ++i) {
                buffer[i] = upper_base(buffer[i]);
            }
            if (!entry.forward) {
                reverseComplement(buffer);
            }
            entry.bases = buffer;
            members.push_back(entry);
        }

        // the common case by far: every member says the same thing everywhere
        bool disagrees = false;
        for (int64_t position = 0; position < length && !disagrees; ++position) {
            for (size_t i = 1; i < members.size(); ++i) {
                if (members[i].bases[position] != members[0].bases[position]) {
                    disagrees = true;
                    break;
                }
            }
        }
        if (!disagrees) {
            continue;
        }
        ++blocks_split;

        // rebuild the block as one block per run of positions that group the same way.  the
        // segments keep their coordinates through this, so the members stay usable, and
        // repinching splits them wherever a run ends
        stPinchBlock_destruct(block);
        block = NULL;

        grouping_at(members, 0, key, base_stamp, base_group, stamp);
        int64_t run_start = 0;
        for (int64_t position = 1; position <= length; ++position) {
            bool end_of_block = (position == length);
            if (!end_of_block) {
                grouping_at(members, position, next_key, base_stamp, base_group, stamp);
            }
            if (!end_of_block && next_key == key) {
                continue;
            }
            // close the run [run_start, position): pinch each member onto the first member
            // that carries its base, and leave a member whose base is unique unpinched
            const int64_t run_length = position - run_start;
            group_first.assign(members.size(), -1);
            for (size_t i = 0; i < members.size(); ++i) {
                int64_t& first = group_first[key[i]];
                if (first < 0) {
                    first = (int64_t)i;
                    continue;
                }
                const BlockMember& to = members[first];
                const BlockMember& from = members[i];
                int64_t to_start = to.forward ? to.start + run_start : to.start + length - position;
                int64_t from_start = from.forward ? from.start + run_start : from.start + length - position;
                stPinchThread_pinch(to.thread, from.thread, to_start, from_start, run_length,
                                    to.forward == from.forward);
            }
            run_start = position;
            if (!end_of_block) {
                key.swap(next_key);
            }
        }
    }

    if (progress) {
        cerr << "separated " << blocks_split << " blocks whose members disagreed" << endl;
    }
}

// make the (still empty) path for a hal sequence: a reference path if its genome is one of
// refNames, or if there are no refNames at all, and a haplotype path otherwise
path_handle_t create_sequence_path(const string& genomeName,
                                   const string& sequenceName,
                                   const vector<string>& refNames,
                                   MutablePathMutableHandleGraph& graph) {
    PathSense sense = PathSense::REFERENCE;
    if (!refNames.empty() && !std::binary_search(refNames.begin(), refNames.end(), genomeName)) {
        sense = PathSense::HAPLOTYPE;
    }

    // cactus_graphmap_split can make paths like contig_sub_1_3.  here we convert that
    // into a format vg can (sometimes) understand contig[1-3].
    // (the reason we go through this is that assembly hubs can't handle any special characters apparently)
    string parsed_name = sequenceName;
    subrange_t subpath = resolve_subpath_naming(parsed_name);
    string parsed_genome_name = genomeName;
    size_t haplotype = resolve_haplotype_naming(parsed_genome_name);
    if (haplotype == PathMetadata::NO_HAPLOTYPE) {
        haplotype = 0;
    }
    return graph.create_path(sense,
                             parsed_genome_name,
                             parsed_name,
                             haplotype,
                             sense == PathSense::HAPLOTYPE ? 0 : PathMetadata::NO_PHASE_BLOCK,
                             subpath,
                             false);
}

// create nodes and edges for a genome using the pinch graph
void pinch_to_handle(AlignmentReader& alignment, size_t genome,
                     const PinchThreads& threads,
                     MutablePathMutableHandleGraph& graph,
                     const vector<string>& refNames) {

    // iterate over the sequences of the genome
    const string& genomeName = alignment.genome_name(genome);
    const size_t numSequences = alignment.num_sequences(genome);
    for (size_t sequence = 0; sequence < numSequences; ++sequence) {
        const string sequenceName = alignment.sequence_name(genome, sequence);
        // as HAL names it
        const string fullName = genomeName + "." + sequenceName;
        stPinchThread* thread = threads.thread(genome, sequence);

        // create the path
        path_handle_t pathHandle = create_sequence_path(genomeName, sequenceName, refNames, graph);
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
            alignment.fetch(genome, sequence, segStart, stPinchSegment_getLength(seg), seqString);

            // have we already converted this block?
            uint64_t blockPacked = block != nullptr ? stPinchBlock_getNumSupportingHomologies(block) : 0;
            nid_t blockNode = (nid_t)unpack_block_id(blockPacked);
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
                    if ((uint64_t)graph.get_id(handle) >= BLOCK_ID_LIMIT) {
                        // refuse rather than truncate the id into the checksum bits.  a node
                        // needs at least one base, so this wants a graph of a trillion bases
                        throw runtime_error("node id " + std::to_string(graph.get_id(handle)) +
                                            " does not fit the per-block node index");
                    }
                    // nodeString is only set in the reversed branch above; the node spells
                    // seqString as it stands otherwise
                    stPinchBlock_setNumSupportingHomologies(
                        block, pack_block_node((uint64_t)graph.get_id(handle),
                                               sequence_checksum(reversed ? nodeString : seqString)));
                }
#ifdef debug
                cerr << "created node " << graph.get_id(handle) << " for block " << block << " from " << fullName << " at " << segStart
                     << " rev=" << reversed << " len=" << seqString.length()
                     << endl;
                cerr << "node seq " << graph.get_sequence(handle) << endl;
#endif
            } else {
                // yes: the id is stored on the block itself
                handle = graph.get_handle(blockNode);
#ifdef debug
                cerr << "found node " << graph.get_id(handle) << " for block " << block << " from " << fullName << " at " << segStart
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
                // the node's own bases are not read back: its checksum was stored on the
                // block when it was created, so checksumming what the hal says here is
                // enough to agree.  the length is compared too, which the old base-by-base
                // loop did not do -- it stopped at the shorter of the two.
                if (reversed) {
                    nodeString = seqString;
                    reverseComplement(nodeString);
                }
                bool lengthAgrees = graph.get_length(handle) == seqString.size();
                bool checksumAgrees = sequence_checksum(reversed ? nodeString : seqString) ==
                    unpack_block_checksum(blockPacked);
                if (!lengthAgrees || !checksumAgrees) {
                    // something is wrong: read the node out and say exactly what, which is
                    // the only place that pays for the slow comparison
                    nodeString = graph.get_sequence(handle);
                    if (!lengthAgrees) {
                        throw runtime_error("node " + std::to_string(blockNode) + " has length " +
                                            std::to_string(nodeString.size()) + " but " +
                                            fullName + " covers " +
                                            std::to_string(seqString.size()) + " bases of it at " +
                                            std::to_string(segStart));
                    }
                    for (size_t i = 0; i < nodeString.size() && i < seqString.size(); ++i) {
                        if (upper_base(nodeString[i]) != upper_base(seqString[i])) {
                            if (mismatches.size() < 10) {
                                mismatches.push_back(make_pair(segStart + i, make_pair(nodeString[i], seqString[i])));
                            }
                            ++numMismatches;
                        }
                    }
                    if (numMismatches == 0) {
                        throw runtime_error("checksum mismatch on node " + std::to_string(blockNode) +
                                            " for " + fullName + " at " +
                                            std::to_string(segStart) + ", but its bases agree");
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
        if ((int64_t)pathLength != alignment.sequence_length(genome, sequence)) {
            throw runtime_error("Incorrect length in coverted path for " + fullName + ": " + std::to_string(pathLength) +
                                ". Should be: " + std::to_string(alignment.sequence_length(genome, sequence)));
        }
        if (numMismatches > 0) {
            stringstream msg;
            msg << numMismatches << " mismatches found in converted path for " << fullName << ":\n";
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

// a genome name is SAMPLE#HAPLOTYPE (PanSN), or just SAMPLE.  a '.' is never a haplotype: Cactus
// used to spell one SAMPLE.HAPLOTYPE, but now only reads that spelling in its input seqfile and
// writes the genome as SAMPLE#HAPLOTYPE, so KOLF2.1J and GCA_054883195.1 are just samples.  (a HAL
// from an older Cactus that names its genomes HG002.1 needs them renamed, with halRenameGenomes)
size_t resolve_haplotype_naming(string& genome_name) {
    size_t sep = genome_name.find('#');
    if (sep == string::npos) {
        return PathMetadata::NO_HAPLOTYPE;
    }
    if (sep == 0 || sep + 1 == genome_name.length() ||
        genome_name.find_first_not_of("0123456789", sep + 1) != string::npos) {
        throw runtime_error("[hal2vg] genome name " + genome_name + " has a # but is not SAMPLE#HAPLOTYPE "
                            "(with HAPLOTYPE an integer), so it cannot be made into PanSN path names");
    }
    size_t haplotype = stoul(genome_name.substr(sep + 1));
    genome_name = genome_name.substr(0, sep);
    return haplotype;
}

static void build_low_memory_plans(stPinchThreadSet* threadSet,
                                   AlignmentReader& alignment,
                                   const PinchThreads& threads,
                                   LowMemoryPlans& plans, bool progress, size_t bufferBytes) {
    uint64_t visited = 0, virtualNodes = 0, labelBytes = 0, storageBytes = 0;
    stPinchThreadSetBlockIt blockIt = stPinchThreadSet_getBlockIt(threadSet);
    for (stPinchBlock* block = stPinchThreadSetBlockIt_getNext(&blockIt); block != nullptr;
         block = stPinchThreadSetBlockIt_getNext(&blockIt)) {
        ++visited;
        if (stPinchBlock_getDegree(block) < 2) continue;
        auto result = unique_ptr<LowMemoryBlock>(new LowMemoryBlock());
        // the plans are kept until the graph is built, so none of their space should be spare
        result->members.reserve(stPinchBlock_getDegree(block));
        stPinchBlockIt memberIt = stPinchBlock_getSegmentIterator(block);
        for (stPinchSegment* member = stPinchBlockIt_getNext(&memberIt); member != nullptr;
             member = stPinchBlockIt_getNext(&memberIt)) {
            result->members.push_back(member);
        }
        // Stable pointer order gives one compact member lookup per ORIGINAL segment,
        // rather than a hash entry for every expanded segment.
        sort(result->members.begin(), result->members.end(), less<stPinchSegment*>());
        const size_t degree = result->members.size();
        const uint64_t length = stPinchBlock_getLength(block);
        const size_t window = max(size_t(1), bufferBytes / degree);
        vector<SequenceRef> sequences(degree);
        vector<int64_t> starts(degree);
        vector<bool> forward(degree);
        vector<string> bases(degree);
        for (size_t i = 0; i < degree; ++i) {
            sequences[i] = threads.sequence(stPinchSegment_getName(result->members[i]));
            starts[i] = stPinchSegment_getStart(result->members[i]);
            forward[i] = stPinchSegment_getBlockOrientation(result->members[i]) != 0;
        }
        hal2vg_low_memory::Builder builder(degree);
        vector<uint8_t> key(degree), previous;
        uint64_t stamp = 0, seen[256] = {0};
        uint8_t group[256] = {0};
        for (uint64_t begin = 0; begin < length;) {
            const size_t count = min<uint64_t>(window, length - begin);
            for (size_t i = 0; i < degree; ++i) {
                const uint64_t start = starts[i] + (forward[i] ? begin : length - begin - count);
                alignment.fetch(sequences[i].genome, sequences[i].sequence, start, count, bases[i]);
                for (char& base : bases[i]) base = upper_base(base);
                if (!forward[i]) reverseComplement(bases[i]);
            }
            for (size_t position = 0; position < count; ++position) {
                ++stamp;
                unsigned groups = 0;
                for (size_t i = 0; i < degree; ++i) {
                    const unsigned char base = bases[i][position];
                    if (seen[base] != stamp) {
                        seen[base] = stamp;
                        group[base] = groups++;
                    }
                    key[i] = group[base];
                }
                if (previous.empty()) {
                    previous = key;
                } else if (key != previous) {
                    builder.append_run(previous, begin + position);
                    previous.swap(key);
                }
            }
            begin += count;
        }
        builder.append_run(previous, length);
        if (builder.has_disagreement()) {
            result->plan = builder.finish();
            virtualNodes += result->plan.nodes.size();
            labelBytes += result->plan.labels.size();
            // the block itself and its entry in the plans table are about as big as its
            // vectors are on a small block, so they are counted too
            storageBytes += result->plan.storage_bytes() + result->members.capacity() * sizeof(stPinchSegment*) +
                sizeof(LowMemoryBlock) + sizeof(LowMemoryPlans::value_type) + 2 * sizeof(void*);
            plans.emplace(block, std::move(result));
        }
        if (progress && visited % 1000000 == 0) {
            cerr << "examined " << visited << " blocks: " << plans.size() << " to separate into "
                 << virtualNodes << " nodes" << endl;
        }
    }
    if (progress) {
        cerr << "found " << plans.size() << " of " << visited << " blocks whose members disagreed,"
             << " to be separated into " << virtualNodes << " nodes (" << labelBytes
             << " bytes of labels, about " << storageBytes << " bytes of plans in all)" << endl;
    }
}

static void pinch_to_handle_low_memory(AlignmentReader& alignment, size_t genome,
                                        const PinchThreads& threads,
                                        LowMemoryPlans& plans,
                                        MutablePathMutableHandleGraph& graph,
                                        const vector<string>& refNames) {
    const string& genomeName = alignment.genome_name(genome);
    const size_t numSequences = alignment.num_sequences(genome);
    for (size_t sequence = 0; sequence < numSequences; ++sequence) {
        const string sequenceName = alignment.sequence_name(genome, sequence);
        // as HAL names it
        const string fullName = genomeName + "." + sequenceName;
        path_handle_t path = create_sequence_path(genomeName, sequenceName, refNames, graph);
        stPinchThread* thread = threads.thread(genome, sequence);
        // One chromosome string avoids an HDF5 iterator setup for every output step.
        string dna;
        alignment.fetch(genome, sequence, 0, alignment.sequence_length(genome, sequence), dna);
        uint64_t pathLength = 0;
        bool havePrevious = false;
        handle_t previous;
        auto emit = [&](uint64_t start, uint64_t length, bool reversed, uint64_t packed) {
            // the one empty interval is an empty sequence's: its thread is a single segment of
            // length 0, which --legacyConversion makes into an empty node, and so does this
            if (start > dna.size() || length > dna.size() - start || (length == 0 && !dna.empty())) {
                throw runtime_error("invalid path interval for " + fullName);
            }
            string bases = dna.substr(start, length);
            if (reversed) reverseComplement(bases);
            const uint32_t checksum = sequence_checksum(bases);
            nid_t id = unpack_block_id(packed);
            handle_t handle;
            if (id == 0) {
                handle = graph.create_handle(bases);
                id = graph.get_id(handle);
                if (id <= 0 || static_cast<uint64_t>(id) >= BLOCK_ID_LIMIT) {
                    throw runtime_error("node id " + std::to_string(id) + " does not fit the per-block node index");
                }
                packed = pack_block_node(id, checksum);
            } else {
                handle = graph.get_handle(id);
                if (graph.get_length(handle) != length || unpack_block_checksum(packed) != checksum) {
                    throw runtime_error("converted path disagrees with the hal for " + fullName +
                                        " at " + std::to_string(start));
                }
            }
            if (reversed) handle = graph.flip(handle);
            if (havePrevious) graph.create_edge(previous, handle);
            graph.append_step(path, handle);
            previous = handle;
            havePrevious = true;
            pathLength += length;
            return packed;
        };
        uint64_t segmentStart = 0;
        for (stPinchSegment* segment = stPinchThread_getFirst(thread); segment != nullptr;
             segment = stPinchSegment_get3Prime(segment)) {
            stPinchBlock* block = stPinchSegment_getBlock(segment);
            const uint64_t length = stPinchSegment_getLength(segment);
            auto found = plans.find(block);
            if (found == plans.end()) {
                const bool reversed = block != nullptr && stPinchSegment_getBlockOrientation(segment) == 0;
                uint64_t packed = block != nullptr ? stPinchBlock_getNumSupportingHomologies(block) : 0;
                packed = emit(segmentStart, length, reversed, packed);
                if (block != nullptr) stPinchBlock_setNumSupportingHomologies(block, packed);
            } else {
                LowMemoryBlock& split = *found->second;
                const auto member = lower_bound(split.members.begin(), split.members.end(), segment,
                                                less<stPinchSegment*>());
                if (member == split.members.end() || *member != segment) {
                    throw runtime_error("original pinch member missing from split plan");
                }
                const size_t index = member - split.members.begin();
                const bool forward = stPinchSegment_getBlockOrientation(segment) != 0;
                uint64_t lastNode = std::numeric_limits<uint64_t>::max();
                for (size_t step = 0; step < split.plan.runs; ++step) {
                    const size_t run = forward ? step : split.plan.runs - 1 - step;
                    const uint64_t nodeId = split.plan.node_for(index, run);
                    if (nodeId == lastNode) continue;
                    lastNode = nodeId;
                    auto& node = split.plan.nodes[nodeId];
                    const uint64_t start = segmentStart + (forward ? node.begin : length - node.end);
                    node.packed = emit(start, node.end - node.begin, !forward, node.packed);
                }
            }
            segmentStart += length;
        }
        if (pathLength != dna.size() || segmentStart != dna.size()) {
            throw runtime_error("Incorrect length in converted path for " + fullName);
        }
    }
}
