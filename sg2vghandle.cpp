/*
 * Copyright (C) 2015 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#include <iostream>
#include <sstream>
#include <cmath>

#include "sg2vghandle.h"

using namespace std;
using namespace handlegraph;

SG2VGHandle::SG2VGHandle() : _graph(0),  _sg(0), _bases(0), _paths(0)
{
}

SG2VGHandle::~SG2VGHandle()
{

}

void SG2VGHandle::convert(const SideGraph* sg,
                          const vector<string>& bases,
                          const vector<SGNamedPath>& paths,
                          MutablePathMutableHandleGraph* graph)
{
  _graph = graph;
  _sg = sg;
  _bases = &bases;
  _paths = &paths;

  // add every node to handle doc
  for (int i = 0; i < _sg->getNumSequences(); ++i)
  {
    addNode(_sg->getSequence(i));
  }

  // add every edge to handle doc
  const SideGraph::JoinSet* joinSet = _sg->getJoinSet();
  for (SideGraph::JoinSet::const_iterator i = joinSet->begin();
       i != joinSet->end(); ++i)
  {
    addEdge(*i);
  }

  // add every path to handle doc
  for (int i = 0; i < paths.size(); ++i)
  {
    addPath(paths[i].first, paths[i].second);
  }
}


void SG2VGHandle::addNode(const SGSequence* seq)
{
  _graph->create_handle(_bases->at(seq->getID()),
                        // node id's are 1-based in VG!
                        seq->getID() + 1);
}

void SG2VGHandle::addEdge(const SGJoin* join)
{
  handle_t node1 = _graph->get_handle(join->getSide1().getBase().getSeqID() + 1,
                                      join->getSide1().getForward() == true);
  handle_t node2 = _graph->get_handle(join->getSide2().getBase().getSeqID() + 1,
                                      join->getSide2().getForward() == false);
  _graph->create_edge(node1, node2);
}

void SG2VGHandle::addPath(const string& name, const vector<SGSegment>& path,
                        int rank)
{
  path_handle_t path_handle = _graph->create_path_handle(name);
  
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

    handle_t handle = _graph->get_handle(sgSeqID + 1,
                                         !path[i].getSide().getForward());    
    step_handle_t step_handle = _graph->append_step(path_handle, handle);
    
    size_t step_len = _sg->getSequence(sgSeqID)->getLength();
    assert(step_len == _graph->get_length(handle));
    
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
