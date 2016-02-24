#include "utils/DirBC.hpp"

YAFEL_NAMESPACE_OPEN

void DirBC::apply(sparse_csr<double> &Ksys, Vector<double> &Fsys) {
  std::cout << "Applying bcs" << std::endl;

  for(std::size_t i=0; i<bcvals.size(); ++i) {
    ubc(bcdofs[i]) = bcvals[i];
  }

  Fsys -= Ksys*ubc;

  for(std::size_t i=0; i<bcvals.size(); ++i) {
    Fsys(bcdofs[i]) = bcvals[i];
  }
  
  //set Ksys entries to 0 or 1 to apply BC's
  for(std::size_t r=0; r<Ksys.rows(); ++r) {
    for(std::size_t i=Ksys.row_ptr[r]; i<Ksys.row_ptr[r+1]; ++i) {
      std::size_t c = Ksys.col_index[i];
      if(bcmask[r] || bcmask[c]) {
        std::cout << "\tzeroing element (" << r << "," << c <<")\n";
	Ksys._data[i] = 0.0;
      }
      if(r==c && bcmask[r]) {
        std::cout << "\tzSetting element ("<<r<<","<<c <<") to 1.0\n";
	Ksys._data[i] = 1.0;
      }
    }
  }

}

YAFEL_NAMESPACE_CLOSE
