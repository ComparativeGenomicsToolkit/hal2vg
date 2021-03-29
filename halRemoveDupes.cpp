/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

//#define debug

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>
#include <unordered_map>

#include "hal.h"

using namespace std;
using namespace hal;

static void initParser(CLParser* optionsParser) {
    optionsParser->addArgument("halFile", "input hal file");
    optionsParser->addArgument("genome", "remove all paralogy edges from this genome");
    optionsParser->setDescription("Remove paralogy edges from given genome (in place)");
}

int main(int argc, char** argv) {
    CLParser optionsParser(WRITE_ACCESS);
    initParser(&optionsParser);
    string halPath;
    string genomeName;
    try {
        optionsParser.parseOptions(argc, argv);
        halPath = optionsParser.getArgument<string>("halFile");
        genomeName = optionsParser.getArgument<string>("genome");
    }
    catch(exception& e) {
        cerr << e.what() << endl;
        optionsParser.printUsage(cerr);
        exit(1);
    }
    try {
        AlignmentPtr alignment(openHalAlignment(halPath, &optionsParser, READ_ACCESS | WRITE_ACCESS));
        if (alignment->getNumGenomes() == 0) {
            throw hal_exception("input hal alignmenet is empty");
        }

        Genome* genome = alignment->openGenome(genomeName);
        if (genome == NULL) {
            throw hal_exception("Genome " + genomeName + " not found in alignment");
        }

        if (genomeName == alignment->getRootName()) {
            throw hal_exception("Cannot run on root");
        }

        TopSegmentIteratorPtr topIt = genome->getTopSegmentIterator();

        size_t total_length = 0;
        size_t total_edges = 0;
        for (; not topIt->atEnd(); topIt->toRight()) {
            TopSegment* topSeg = topIt->tseg();
            if (topSeg->hasNextParalogy()) {
                topSeg->setNextParalogyIndex(NULL_INDEX);
                if (!topSeg->isCanonicalParalog()) {
                    topSeg->setParentIndex(NULL_INDEX);
                    total_length += topSeg->getLength();
                    ++total_edges;
                }
            }
        }

        if (total_length > 0) {
            cerr << "[halRemoveDupes]: " << total_edges << " paralogy edges removed from " << genomeName
                 << " with total length " << total_length << endl;
        } else {
            cerr << "[halRemoveDupes] : No paralogy edges found in " << genomeName << endl;
        }
    }

    catch(exception& e) {
        cerr << e.what() << endl;
        exit(1);
    }
     
    return 0;
}
