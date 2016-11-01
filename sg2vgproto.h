/*
 * Copyright (C) 2016 by Glenn Hickey (hickey@soe.ucsc.edu)
 *
 * Released under the MIT license, see LICENSE.cactus
 */

#ifndef _SG2VGPROTO_H
#define _SG2VGPROTO_H

#include <vector>
#include <string>
#include <stdexcept>

#include "vg.pb.h"

#include "sidegraph.h"

/** This class replaces 

https://github.com/glennhickey/sg2vg/blob/master/sg2vgjson.h

changing the JSON output to protobuf.  We write protobuf 
directly instead of going through the VG class to save memory
and preserve some of the old (naive) chunking logic from sg2vgjson. 
 
We are writing a SideGraph object, but one that was created with 
side2seq -- ie all joins are to ends of sequences, so can be translated
directly to vg sequence graph... 
*/

class SG2VGProto
{
public:
   SG2VGProto();
   ~SG2VGProto();

   /** init output stream and proto document */
   void init(std::ostream* os);

   /** write nodes and edges and paths*/
   void writeGraph(const SideGraph* sg,
                   const std::vector<std::string>& bases,
                   const std::vector<SGNamedPath>& paths);

   /** write a graph chunk by chunk. chunks are not subgraphs, they are just
    * bags of nodes, then bags of edges, then bags of paths (streamed out in
    * that order) */
   void writeChunkedGraph(const SideGraph* sg,
                          const std::vector<std::string>& bases,
                          const std::vector<SGNamedPath>& paths,
                          int sequencesPerChunk = 5000,
                          int joinsPerChunk = 100000,
                          int pathSegsPerChunk = 10000);
   
protected:

   vg::Graph _graph;

   // add to proto 
   void addNode(const SGSequence* seq);
   void addEdge(const SGJoin* join);
   void addPath(const std::string& name, const std::vector<SGSegment>& path,
                int rank = 0);

   std::ostream* _os;
   const SideGraph* _sg;
   const std::vector<std::string>* _bases;
   const std::vector<std::pair<std::string, std::vector<SGSegment> > >* _paths;
};


#endif
