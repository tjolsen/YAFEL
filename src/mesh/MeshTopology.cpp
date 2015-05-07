#include "mesh/MeshTopology.hpp"

YAFEL_NAMESPACE_OPEN

MeshTopology::MeshTopology(const Mesh &M) :
  points(), lines(), faces()
{
  
  unsigned nPoints = M.get_n_nodes();

  // create vector of TopoPoints;
  for(unsigned i=0; i<nPoints; ++i) {
    TopoPoint tp(i);
    points.push_back(tp);
  }

  // Loop over elements, creating lines and cells as needed
  unsigned nElems = M.get_n_elems();
  for(unsigned e=0; e<nElems; ++e) {
    unsigned etype = M.element_type[e];

    if( etype==1 || etype==8 || etype==26 || etype==27 || etype==28 ||
	etype==62 || etype==63 || etype==64 || etype==65 || etype==66) {
      //handle 2-node line element
      unsigned id1, id2;
      id1 = M.elements[e][0];
      id2 = M.elements[e][1];
      bool line_exists = false;
      unsigned lineid = get_line_id(id1, id2, line_exists);
      if(!line_exists) {
	std::cout << "Creating line from line elem" << std::endl;
	unsigned tail = (id1<id2) ? id1 : id2;
	unsigned head = (id1<id2) ? id2 : id1;
	lineid = lines.size();
	lines.emplace_back(new TopoLine(lineid, &points[head], &points[tail]));
	points[head].incoming.push_back(lines[lineid]);
	points[tail].outgoing.push_back(lines[lineid]);
      }
    }
    else if(etype == 2 || etype == 3) {
      
      // handle 3-node triangle or 4-node quad
      unsigned Nverts = M.elements[e].size();
      faces.emplace(e,new TopoFace(e));
      TopoFace *TFp = faces[e];
      
      //get vertices
      for(unsigned i=0; i<Nverts; ++i) {
	TFp->vertices.push_back(&points[M.elements[e][i]]);
      }
      
      //get edges
      for(unsigned edge=0; edge<Nverts; ++edge) {
	unsigned id1,id2;
	id1 = M.elements[e][edge];
	id2 = M.elements[e][(edge+1)%Nverts];
	
	bool line_exists = false;
	unsigned lineid = get_line_id(id1, id2, line_exists);
	if(!line_exists) {
	  unsigned tail = (id1<id2) ? id1 : id2;
	  unsigned head = (id1<id2) ? id2 : id1;
	  lineid = lines.size();
	  TopoLine *TL = new TopoLine();//(lineid, &(points[head]), &(points[tail]));
	  TL->id = lineid;
	  TL->head = &points[head];
	  TL->tail = &points[tail];
	  TL->right = nullptr;
	  TL->left = nullptr;
	  lines.push_back(TL);
	  
	  points[head].incoming.push_back( lines[lineid] );
	  points[tail].outgoing.push_back( lines[lineid] );
	}
	
	TFp->boundary.push_back(lines[lineid]);
	
	if(id1 == lines[lineid]->tail->id) {
	  lines[lineid]->left = TFp;
	}
	else {
	  lines[lineid]->right = TFp;
	}
      }
    }
    else {
      //emit warning/error?
      std::cerr << "MeshTopology: Warning, unsupported element type " 
		<< etype << std::endl;
    }
    
  } // end e-loop
  
  
}

MeshTopology::~MeshTopology() {
  
  for(auto l : lines) {
    delete l;
  }
  
  for(auto f : faces) {
    delete f.second;
  }
  
}



unsigned MeshTopology::get_line_id(unsigned vid1, unsigned vid2, bool &exists) {
  
  unsigned retid = 0;
  unsigned minid = (vid1<vid2 ? vid1 : vid2);
  unsigned maxid = (vid1<vid2 ? vid2 : vid1);
  TopoPoint & TP = points[minid];

  for(unsigned i=0; i<TP.outgoing.size(); ++i) {
    if(TP.outgoing[i]->head->id == maxid) {
      exists = true;
      retid = TP.outgoing[i]->id;
      return retid;
    }
  }

  retid = 0;
  exists = false;
  
  return retid;
}

void MeshTopology::print(std::ostream &out) {
  
  // print point info
  for(unsigned i=0; i<points.size(); ++i) {
    out << points[i] << "\n";
  }
  out << "\n";
  
  //print line info
  out << "===== " << lines.size() << " Lines =====\n";
  
  for(auto l=lines.begin(); l != lines.end(); ++l) {
    out << *(*l) << "\n";
  }
  out << "\n";
  
  //print face info
  for(auto f = faces.begin(); f!=faces.end(); ++f) {
    out << "Face:\n\tid = "<< f->second->id << "\n\tVertices:";
    for(auto v=f->second->vertices.begin(); v != f->second->vertices.end(); ++v) {
      out << " " << (*v)->id;
    }
    out << "\n\tBoundary:";
    for(auto e=f->second->boundary.begin(); e != f->second->boundary.end(); ++e) {
      out << " " << (*e)->id;
    }
    out << "\n";
  }
  out << "\n";
}

unsigned MeshTopology::getCellNeighbor(unsigned faceNum, unsigned edgenum) {

  TopoFace &TF = *(faces[faceNum]);
  
  //check number of neighbors (= number boundary edges), return facenum if edgenum too large
  if(edgenum >= TF.boundary.size()) {
    return faceNum;
  }

  TopoFace *neighbor;
  TopoLine *edge = TF.boundary[edgenum];
  if(edge->right!=nullptr && (faceNum == edge->right->id) ) {
    neighbor = edge->left;
  }
  else {
    neighbor = edge->right;
  }

  return (neighbor != nullptr) ? neighbor->id : faceNum;
}

YAFEL_NAMESPACE_CLOSE
