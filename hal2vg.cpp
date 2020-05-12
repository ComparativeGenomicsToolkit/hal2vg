/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.txt
 */

// This file was created by merging hal2sg.cpp and sg2vg.cpp with
// a small amount of glue for the interface. 

#include <cstdlib>
#include <iostream>
#include <cassert>
#include <fstream>
#include <deque>

#include "sgbuilder.h"
#include "side2seq.h"
#include "sg2vghandle.h"
#include "bdsg/packed_graph.hpp"
#include "bdsg/hash_graph.hpp"
#include "bdsg/odgi.hpp"

using namespace std;
using namespace hal;
using namespace handlegraph;
using namespace bdsg;

static bool isCamelHal(AlignmentConstPtr aligment);
static void breadthFirstGenomeSearch(const Genome* reference,
                                     const vector<const Genome*>& targets,
                                     vector<const Genome*>& outTraversal);
static set<const Sequence*> parseRefSequences(
  const Genome* refGenome, const string& refSequenceFilePath);


static void initParser(CLParser* optionsParser)
{
  optionsParser->addArgument("halFile", "input hal file");
  optionsParser->addOption("refGenome", 
                           "name of reference genome (HAL root if empty)", 
                           "\"\"");
  optionsParser->addOption("rootGenome", 
                           "process only genomes in clade with specified root"
                           " (HAL root if empty)", 
                           "\"\"");
  optionsParser->addOption("targetGenomes",
                           "comma-separated (no spaces) list of target genomes "
                           "(others are excluded) (vist all if empty)",
                           "\"\"");
  optionsParser->addOptionFlag("noAncestors", 
                               "don't write ancestral paths. Note that "
                               "ancestral *sequences* may still get written"
                               " as they can be necessary for expressing"
                               " some alignments. IMPORTANT: "
                               "Must be used in conjunction with --refGenome"
                               " to set a non-ancestral genome as the reference"
                               " because the default reference is the root.", 
                               false);
  optionsParser->addOptionFlag("refDupes",
                               "process duplications in reference genome, which will"
                               " be written as vg cycles.  This is off by default.",
                               false);
  optionsParser->addOptionFlag("onlySequenceNames",
                               "use only sequence names for output names.  By "
                               "default, the UCSC convention of "
                               "Genome.Sequence is used",
                               false);
  optionsParser->addOptionFlag("keepCase",
                               "don't convert all nucleotides to upper case",
                               false);
  optionsParser->addOption("refSequenceFile",
                           "white-space delimited list of sequence names in the "
                           "reference genome which will *not* be collapsed by duplications."
                           "  Overrides --refDupes", "\"\"");
  optionsParser->addOption("outputFormat",
                           "output graph format in {pg, hg, odgi} [default=pg]",
                           "pg");

  optionsParser->setDescription("Convert HAL alignment to handle graph");

}

int main(int argc, char** argv)
{
  CLParser optionsParser;
  initParser(&optionsParser);
  string halPath;
  string refGenomeName;
  string rootGenomeName;
  string targetGenomes;
  bool noAncestors;
  bool refDupes;
  bool onlySequenceNames;
  bool keepCase;
  // Hardcode a large chop value here.  The logic as implemented
  // is not efficient enough to chop into more standard node size (<1000) for
  // larger graphs.  So we only use to make sure we don't overflow protobuf.
  // Todo: tune down?
  const int chop = 1000000;
  string refSequenceFile;
  string outputFormat;
  try
  {
    optionsParser.parseOptions(argc, argv);
    halPath = optionsParser.getArgument<string>("halFile");
    refGenomeName = optionsParser.getOption<string>("refGenome");
    rootGenomeName = optionsParser.getOption<string>("rootGenome");
    targetGenomes = optionsParser.getOption<string>("targetGenomes");
    noAncestors = optionsParser.getFlag("noAncestors");
    refDupes = optionsParser.getFlag("refDupes");    
    onlySequenceNames = optionsParser.getFlag("onlySequenceNames");
    keepCase = optionsParser.getFlag("keepCase");
    outputFormat = optionsParser.getOption<string>("outputFormat");
    refSequenceFile = optionsParser.getOption<string>("refSequenceFile");
    if (rootGenomeName != "\"\"" && targetGenomes != "\"\"")
    {
      throw hal_exception("--rootGenome and --targetGenomes options are "
                          "mutually exclusive");
    }
    if (refSequenceFile != "\"\"" && refGenomeName == "\"\"")
    {
      throw hal_exception("--refSequenceFile must be used in conjunction "
                          " with --refGenome");
    }
    if (outputFormat != "pg" && outputFormat != "hg" && outputFormat != "odgi") {
      throw hal_exception("--outputFormat must be one of {pg, hg, odgi}");
    }
  }
  catch(exception& e)
  {
    cerr << e.what() << endl;
    optionsParser.printUsage(cerr);
    exit(1);
  }
  try
  {
    AlignmentConstPtr alignment(openHalAlignment(halPath, 
                                                 &optionsParser,
                                                 hal::READ_ACCESS));
    if (alignment->getNumGenomes() == 0)
    {
      throw hal_exception("hal alignmenet is empty");
    }

    // root is specified either by the parameter or as the alignment root
    // by default
    const Genome* rootGenome = NULL;
    if (rootGenomeName != "\"\"")
    {
      rootGenome = alignment->openGenome(rootGenomeName);
    }
    else
    {
      rootGenome = alignment->openGenome(alignment->getRootName());
    }
    if (rootGenome == NULL)
    {
      throw hal_exception(string("Root genome, ") + rootGenomeName + 
                          ", not found in alignment");
    }

    // if target set not specified we default to all leaves under the
    // given root. 
    vector<const Genome*> targetVec;
    if (targetGenomes == "\"\"")
    {
      set<const Genome*> allGenomes;
      getGenomesInSubTree(rootGenome, allGenomes);
      for (set<const Genome*>::iterator i = allGenomes.begin();
           i != allGenomes.end(); ++i)
      {
        if ((*i)->getNumChildren() == 0)
        {
          targetVec.push_back(*i);
        }
      }
      // Throw in the root if empty (ie case where reference is a leaf
      // and there are no other leaves)
      if (targetVec.size() == 1 && targetVec[0] != rootGenome)
      {
        targetVec.push_back(rootGenome);
      }
    }
    // target genomes pulled from list.  
    else
    {
      vector<string> targetNames = chopString(targetGenomes, ",");
      for (size_t i = 0; i < targetNames.size(); ++i)
      {
        const Genome* tgtGenome = alignment->openGenome(targetNames[i]);
        if (tgtGenome == NULL)
        {
          throw hal_exception(string("Target genome, ") + targetNames[i] + 
                              ", not found in alignment");
        }
        targetVec.push_back(tgtGenome);
      }
    }

    // open the reference genome (root genome if unspecified)
    const Genome* refGenome = NULL;
    set<const Sequence*> refSequences;
    if (refGenomeName != "\"\"")
    {
      refGenome = alignment->openGenome(refGenomeName);
      if (refGenome == NULL)
      {
        throw hal_exception(string("Reference genome, ") + refGenomeName + 
                            ", not found in alignment");
      }
      set<const Genome*> genomeSet;
      genomeSet.insert(refGenome);
      genomeSet.insert(rootGenome);
      if (getLowestCommonAncestor(genomeSet) != rootGenome)
      {
        throw hal_exception(string("reference genome must be under root"));
      }

      if (refSequenceFile != "\"\"")
      {
        refSequences = parseRefSequences(refGenome, refSequenceFile);
        refDupes = true;
      }
    }
    else
    {
      refGenome = rootGenome;
    }

    // make sure refGenome not in target genomes
    for (vector<const Genome*>::iterator i = targetVec.begin();
         i != targetVec.end(); ++i)
    {
      if (*i == refGenome)
      {
        targetVec.erase(i);
        break;
      }
    }

    // get a breadth-first ordering of all target genomes
    // starting at the reference and including any ancestors
    // that need to get walked across
    vector<const Genome*> breadthFirstOrdering;
    breadthFirstGenomeSearch(refGenome, targetVec, breadthFirstOrdering);

    bool camelMode = isCamelHal(alignment);
    if (camelMode)
    {
      cout << "CAMEL output detected.  Will infer root sequence "
           << "from children" << endl;
    }
    SGBuilder sgbuild;
    sgbuild.init(alignment, rootGenome, refDupes, isCamelHal(alignment),
                 onlySequenceNames, true);
    
    // add the genomes in the breadth first order
    for (size_t i = 0; i < breadthFirstOrdering.size(); ++i)
    {
      sgbuild.addGenome(breadthFirstOrdering[i], NULL,
                        i == 0 ? &refSequences : NULL);
    }

    // compute all the joins in second pass (and do sanity check
    // on every path in graph)
    sgbuild.computeJoins(!noAncestors);

    // write out the sequence strings in order
    const SideGraph* sg = sgbuild.getSideGraph();
    vector<string> sequenceStrings(sg->getNumSequences());
    for (sg_int_t i = 0; i < sg->getNumSequences(); ++i)
    {
      const SGSequence* seq = sg->getSequence(i);
      sgbuild.getSequenceString(seq, sequenceStrings[i]);
    }

    // get the paths
    vector<const Sequence*> halSequences = sgbuild.getHalSequences();
    vector<SGNamedPath> namedPaths;
    for (size_t i = 0; i < halSequences.size(); ++i)
    {
      if (halSequences[i]->getGenome()->getNumChildren() == 0 ||
          !noAncestors)
      {
        string pathName = sgbuild.getHalSeqName(halSequences[i]);
      
        vector<SGSegment> path;
        sgbuild.getHalSequencePath(halSequences[i], path);
        namedPaths.push_back(SGNamedPath(pathName, path));
      }
    }

    // destroy the hal alignment to save a bit on memory
    alignment = AlignmentConstPtr();
    sgbuild.clear_except_sg();
    
    // convert side graph into sequence graph 
    cerr << "Converting Side Graph to VG Sequence Graph" << endl;
    Side2Seq converter;
    converter.init(sg, &sequenceStrings, &namedPaths, !keepCase,
                   false, "", chop, true);
    converter.convert();

    // free up the sidegraph from sgbuild because it's no longer
    // needed (and will leak if we don't do so)
    delete sg;
    sg = NULL;
    
    const SideGraph* outGraph = converter.getOutGraph();
    const vector<string>& outBases = converter.getOutBases();
    const vector<SGNamedPath>& outPaths = converter.getOutPaths();

    
    // convert to vg handle
    cerr << "Converting SideGraph to HandleGraph" << endl;
    
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

    SG2VGHandle vgConverter;    
    vgConverter.convert(outGraph, outBases, outPaths, graph.get());

    // write tot stdout
    cerr << "Writing HandleGraph to stdout" << endl;
    dynamic_cast<SerializableHandleGraph*>(graph.get())->serialize(cout);

    //cout << *sgbuild.getSideGraph() << endl;

  }
/*  catch(hal_exception& e)
  {
    cerr << "hal exception caught: " << e.what() << endl;
    return 1;
  }
  catch(exception& e)
  {
    cerr << "Exception caught: " << e.what() << endl;
    return 1;
  }
*/
  catch(int e) {}
  return 0;
}

/** CAMEL writes the root's DNA sequence as N's.  This screws up SNP detection
 * in the conversion when the root is used as the first anchor.  So we run a
 * check here to see if the root's all N's, if it is we set the camel flag
 * in sgbuilder to tell it to infer the root from its children (possible as
 * there are no substitutions in CAMEL HAL files */
bool isCamelHal(AlignmentConstPtr alignment)
{
  const Genome* rootGenome = alignment->openGenome(alignment->getRootName());
  if (rootGenome->getSequenceLength() == 0)
  {
    return false;
  }
  DnaIteratorPtr dnaIt = rootGenome->getDnaIterator();
  size_t length = rootGenome->getSequenceLength();
  for (size_t i = 0; i < length; ++i, dnaIt->toRight())
  {
    if (toupper(dnaIt->getBase()) != 'N')
    {
      return false;
    }
  }
  return true;
}

void breadthFirstGenomeSearch(const Genome* reference,
                              const vector<const Genome*>& targets,
                              vector<const Genome*>& outTraversal)
{
  // find all genomes we need to visit
  set<const Genome*> inputSet;
  inputSet.insert(reference);
  for (size_t i = 0; i < targets.size(); ++i)
  {
    inputSet.insert(targets[i]);
  }
  set<const Genome*> visitSet;
  getGenomesInSpanningTree(inputSet, visitSet);

  // find our breadth first order through the visit set, starting at
  // reference.
  set<const Genome*> flagged;
  deque<const Genome*> bfsQueue;
  bfsQueue.push_back(reference);
  while (!bfsQueue.empty())
  {
    const Genome* genome = bfsQueue.front();
    bfsQueue.pop_front();
    outTraversal.push_back(genome);
    vector<const Genome*> neighbours;
    for (hal_size_t i = 0; i < genome->getNumChildren(); ++i)
    {
      neighbours.push_back(genome->getChild(i));
    }
    if (genome->getParent() != NULL)
    {
      neighbours.push_back(genome->getParent());
    }
    for (hal_size_t i = 0; i < neighbours.size(); ++i)
    {
      const Genome* neighbour = neighbours[i];
      if (visitSet.find(neighbour) != visitSet.end() &&
          flagged.find(neighbour) == flagged.end())
      {
        bfsQueue.push_back(neighbour);
      }
    }
    flagged.insert(genome);
  }
}

set<const Sequence*> parseRefSequences(const Genome* refGenome,
                                       const string& refSequenceFilePath)
{
  ifstream refSequenceFile(refSequenceFilePath.c_str());
  if (!refSequenceFile)
  {
    throw hal_exception("Unable to load reference sequences file: " +
                        refSequenceFilePath);
  }
  set<const Sequence*> refSet;
  while (refSequenceFile)
  {
    string buf;
    refSequenceFile >> buf;
    if (!buf.empty())
    {
      const Sequence* seq = refGenome->getSequence(buf);
      if (seq == NULL)
      {
        throw hal_exception("Sequence " + buf + ", specified in " +
                            refSequenceFilePath + ", not found in reference genome " +
                            refGenome->getName());
      }
      refSet.insert(seq);
    }
  }
  return refSet;
}
