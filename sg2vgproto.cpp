/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <iostream>
#include <sstream>

#include "sg2vgproto.h"
#include "stream.hpp"

using namespace std;
using namespace vg;

SG2VGProto::SG2VGProto() : _os(0), _sg(0), _bases(0), _paths(0)
{
}

SG2VGProto::~SG2VGProto()
{

}

void SG2VGProto::init(ostream* os)
{
  _os = os;
  _graph = Graph();
}

void SG2VGProto::writeGraph(const SideGraph* sg,
                           const vector<string>& bases,
                           const vector<SGNamedPath>& paths)
{
  _sg = sg;
  _bases = &bases;
  _paths = &paths;

  // add every node to proto doc
  for (int i = 0; i < _sg->getNumSequences(); ++i)
  {
    addNode(_sg->getSequence(i));
  }

  // add every edge to proto doc
  const SideGraph::JoinSet* joinSet = _sg->getJoinSet();
  for (SideGraph::JoinSet::const_iterator i = joinSet->begin();
       i != joinSet->end(); ++i)
  {
    addEdge(*i);
  }

  // add every path to proto doc
  for (int i = 0; i < paths.size(); ++i)
  {
    addPath(paths[i].first, paths[i].second);
  }

  function<Graph(uint64_t)> lambda = [this](uint64_t i) -> Graph {
    return _graph;
  };
  
  stream::write(*_os, 1, lambda);
}

void SG2VGProto::writeChunkedGraph(const SideGraph* sg,
                                   const std::vector<std::string>& bases,
                                   const std::vector<SGNamedPath>& paths,
                                   int sequencesPerChunk,
                                   int joinsPerChunk,
                                   int pathSegsPerChunk)
{
  _sg = sg;
  _bases = &bases;
  _paths = &paths;

  // we begin by counting up how many chunks we'll write, using
  // the three difference chunking parameters
  const SideGraph::JoinSet* joinSet = _sg->getJoinSet();
  int totalPathSegments = 0;
  for (int i = 0; i < paths.size(); ++i)
  {
    totalPathSegments += paths[i].second.size();
  }
  int seqChunks = (int)ceil((double)_sg->getNumSequences() / sequencesPerChunk);
  int joinChunks = (int)ceil((double)joinSet->size() / joinsPerChunk);
  int segmentChunks = (int)ceil((double)totalPathSegments / pathSegsPerChunk);
  int totalChunks = seqChunks + joinChunks + segmentChunks;

  // need to keep track of where we're at in the input
  int curSequence = 0;
  SideGraph::JoinSet::const_iterator curJoin = joinSet->begin();
  int curPath = 0;
  int curSegment = 0;

  // fill in the graph and write it to the stream
  function<Graph(uint64_t)> lambda = [&](uint64_t i) -> Graph {
    init(_os);

    // write a chunk's worth of nodes if we're still in the node chunks
    if (i < seqChunks)
    {
      if (i == 0)
      {
        cerr << "Writing " << seqChunks << " chunks of up to "
             << sequencesPerChunk << " nodes" << endl;
      }
      for (int count = 0; count < sequencesPerChunk &&
              curSequence < _sg->getNumSequences();
           ++count, ++curSequence)
      {
        addNode(_sg->getSequence(curSequence));
      }
    }

    // write a chunk's worth of edges  if we're still in the edge chunks
    else if (i < joinChunks + seqChunks)
    {
      if (i == seqChunks)
      {
        cerr << "Writing " << joinChunks << " chunks of up to "
             << joinsPerChunk << " edges" << endl;
      }
      for (int count = 0; count < joinsPerChunk && curJoin != joinSet->end();
           ++count, ++curJoin)
      {
        addEdge(*curJoin);
      }
    }

    // write a chunk's worth of path segments if we're done nodes and edges
    else
    {
      assert(i < totalChunks);
      if (i == joinChunks + seqChunks)
      {
        cerr << "Writing " << segmentChunks << " chunks of up to "
             << pathSegsPerChunk << " path segments" << endl;
      }

      // iterate paths
      for (int count = 0; count < pathSegsPerChunk && curPath < paths.size();
           ++curPath)
      {
        const SGNamedPath& path = paths[curPath];
        if (curSegment == 0 && pathSegsPerChunk - count >= path.second.size())
        {
          // we're adding entire path
          addPath(path.first, path.second);
          count += path.second.size();
        }
        else
        {
          // we need to chop path
          while (curSegment < path.second.size() && count < pathSegsPerChunk)
          {
            int pathChunkSize = min(pathSegsPerChunk,
                                    (int)path.second.size() - curSegment);
            vector<SGSegment>::const_iterator a =
               path.second.begin() + curSegment;
            vector<SGSegment>::const_iterator b = a + pathChunkSize;
            addPath(path.first, vector<SGSegment>(a, b), curSegment);
            count += pathChunkSize;
            curSegment += pathChunkSize;
          }
        }
        curSegment = 0;
      }
    }    
    return _graph;
  };

  stream::write(*_os, totalChunks, lambda);
}


void SG2VGProto::addNode(const SGSequence* seq)
{
  Node* node = _graph.add_node();
  // node id's are 1-based in VG! 
  node->set_id(seq->getID() + 1);
  node->set_sequence(_bases->at(seq->getID()));
}

void SG2VGProto::addEdge(const SGJoin* join)
{
  Edge* edge = _graph.add_edge();
  // node id's are 1-based in VG!
  edge->set_from(join->getSide1().getBase().getSeqID() + 1);
  edge->set_to(join->getSide2().getBase().getSeqID() + 1);
  edge->set_from_start(join->getSide1().getForward() == true);
  edge->set_to_end(join->getSide2().getForward() == false);
}

void SG2VGProto::addPath(const string& name, const vector<SGSegment>& path,
                        int rank)
{
  Path* vgPath = _graph.add_path();
  vgPath->set_name(name);
  
  int inputPathLength = 0;
  int outputPathLength = 0;
  for (int i = 0; i < path.size(); ++i)
  {
    sg_int_t sgSeqID = path[i].getSide().getBase().getSeqID();
    
    if (path[i].getLength() != _sg->getSequence(sgSeqID)->getLength())
    {
      stringstream ss;
      ss << "Sanity check fail for Mapping " << i << " of path " << name
         << ": Segment size " << path[i].getLength() << " does not span "
         << "all of node " << (sgSeqID + 1) << " which has length "
         << _sg->getSequence(sgSeqID)->getLength();
      throw runtime_error(ss.str());
    }
    inputPathLength += path[i].getLength();

    Mapping* mapping = vgPath->add_mapping();
    mapping->set_rank(rank + i + 1);
    Position* position = mapping->mutable_position();
    // node id's are 1-based in VG!
    position->set_node_id(sgSeqID + 1);

    // Offsets are along the strand of the node that is being visited.
    // We always use the whole node.
    position->set_offset(0);
    position->set_is_reverse(!path[i].getSide().getForward());
    
    outputPathLength += _sg->getSequence(sgSeqID)->getLength();    
  }
  if (inputPathLength != outputPathLength)
  {
    stringstream ss;
    ss << "Sanity check fail for path " << name << ": input length ("
       << inputPathLength << ") != output length (" << outputPathLength << ")";
    throw runtime_error(ss.str());
  }
}
