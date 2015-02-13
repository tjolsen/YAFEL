#include "yafel.hpp"
#include "heap.hpp"
#include <iostream>
#include <vector>
#include <cmath>

using namespace yafel;

inline unsigned index_of(unsigned x, unsigned y, unsigned N) {
  return x + y*N;
}

std::ostream& operator<<(std::ostream&os, const std::pair<double,unsigned> &p) {
  os << p.first << "," << p.second;
  return os;
}

double calcMinT(unsigned x, unsigned y, unsigned N, 
		const std::vector<double> &Tfield,
		const std::vector<bool> &visited) {

  bool xneighbor = false;
  bool yneighbor = false;
  double Tx=0, Ty=0;
  //not on edge/corner
  if(x>0 && x<N-1 && y>0 && y<N-1) {
    if(visited[index_of(x-1,y,N)]) {
      xneighbor=true;
      Tx = Tfield[index_of(x-1,y,N)];
    }
    else if(visited[index_of(x+1,y,N)]) {
      xneighbor=true;
      Tx = Tfield[index_of(x+1,y,N)];
    }
    if(visited[index_of(x,y-1,N)]) {
      yneighbor=true;
      Ty = Tfield[index_of(x,y-1,N)];
    }
    else if(visited[index_of(x,y+1,N)]) {
      yneighbor=true;
      Ty = Tfield[index_of(x,y+1,N)];
    }

    double T;
    if(xneighbor && yneighbor) {
      double T1 = (2*(Tx+Ty) + std::sqrt(4*(Tx+Ty)*(Tx+Ty) - 8*(Tx*Tx + Ty*Ty - 1)))/4.0;
      //double T2 = (2*(Tx+Ty) - std::sqrt(4*(Tx+Ty)*(Tx+Ty) - 8*(Tx*Tx + Ty*Ty - 1)))/4.0;
      
      T = T1;
    }
    else if(xneighbor) {
      T = 1 + Tx;
    }
    else if(yneighbor){
      T = 1 + Ty;
    }
    else {
      T = 1.0;
    }
    
    return T;
  }
  
  //on vertical edge
  
  
  //in corner
  return 0;
}

int main(int argc, char **argv) {

  unsigned N;
  if(argc<2)
    N = 101;
  else
    N = atoi(argv[1]);
  
  std::vector<double> Tfield(N*N,0.0);
  std::vector<bool> visited(N*N,false);
  
  typedef std::pair<double, unsigned> node_t;
  double Tinit = -4.0;
  Heap<node_t> PQ;
  PQ.insert(node_t{Tinit,index_of(N/2, N/2, N)});
  //starting node: 0 -> (x,y) = (0,0)


  for(unsigned i=1; i<N-1; ++i) {
    PQ.insert(node_t{Tinit, index_of(i,1,N)});
    PQ.insert(node_t{Tinit, index_of(i,N-2,N)});
    PQ.insert(node_t{Tinit, index_of(1,i,N)});
    PQ.insert(node_t{Tinit, index_of(N-2,i,N)});
    Tfield[index_of(i,1,N)] = Tinit;
    Tfield[index_of(i,N-2,N)] = Tinit;
    Tfield[index_of(1,i,N)] = Tinit;
    Tfield[index_of(N-2,i,N)] = Tinit;
    visited[index_of(i,1,N)] = true;
    visited[index_of(i,N-2,N)] = true;
    visited[index_of(1,i,N)] = true;
    visited[index_of(N-2,i,N)] = true;
  }


  //iterate until all nodes have been handled
  while(PQ.size() > 0) {

    node_t n = PQ.extract();
    unsigned index = n.second;
    double t = n.first;
    if(visited[index]) 
      continue;
    
    //std::cout << index << std::endl;
    visited[index] = true;
    Tfield[index] = t;
    unsigned x = index%N;
    unsigned y = index/N;
    
    //add un-visited neighbors to priority queue
    //up
    if(y < N-1 && !visited[index_of(x,y+1,N)]) {
      double T = calcMinT(x,y+1,N,Tfield, visited);
      PQ.insert(node_t{T, index_of(x,y+1,N)});
    }
    //right
    if(x < N-1 && !visited[index_of(x+1,y,N)]) {
      double T = calcMinT(x+1,y,N,Tfield, visited);
      PQ.insert(node_t{T, index_of(x+1,y,N)});      
    }
    //down
    if(y > 0 && !visited[index_of(x,y-1,N)]) {
      double T = calcMinT(x,y-1,N,Tfield, visited);
      PQ.insert(node_t{T, index_of(x,y-1,N)});      
    }
    //left
    if(x > 0 && !visited[index_of(x-1,y,N)]) {
      double T = calcMinT(x-1,y,N,Tfield, visited);
      PQ.insert(node_t{T, index_of(x-1,y,N)});      
    }

    //PQ.print();
  }

  FullMatrix Tmat(N,N,0.0);
  for(unsigned i=0; i<N; ++i) {
    for(unsigned j=0; j<N; ++j) {
      Tmat(j,i) = Tfield[index_of(j,i,N)];
    }
  }

  Tmat.print();
  
  MatrixVisualization::contour(Tmat);

  return 0;
}
