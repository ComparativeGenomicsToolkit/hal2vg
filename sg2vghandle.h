/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SG2VGHANDLE_H
#define _SG2VGHANDLE_H

#include <vector>
#include <string>
#include <stdexcept>
#include "handlegraph/mutable_path_mutable_handle_graph.hpp"
#include "sidegraph.h"

/** This class replaces the old sg2vgproto.h which itself replaced

https://github.com/glennhickey/sg2vg/blob/master/sg2vgjson.h

changing the JSON/protobuf output to any handle graph implementation in libbdsg

It will require more memory than the old streaming implementation, but will
write something that's more useful and compatible. 
*/

class SG2VGHandle
{
public:
   SG2VGHandle();
   ~SG2VGHandle();

   /** write nodes and edges and paths*/
   void convert(const SideGraph* sg,
                const std::vector<std::string>& bases,
                const std::vector<SGNamedPath>& paths,
                handlegraph::MutablePathMutableHandleGraph* graph);
   
protected:

   handlegraph::MutablePathMutableHandleGraph* _graph;

   // add to handle 
   void addNode(const SGSequence* seq);
   void addEdge(const SGJoin* join);
   void addPath(const std::string& name, const std::vector<SGSegment>& path,
                int rank = 0);

   const SideGraph* _sg;
   const std::vector<std::string>* _bases;
   const std::vector<std::pair<std::string, std::vector<SGSegment> > >* _paths;
};


#endif
