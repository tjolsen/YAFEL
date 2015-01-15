#include "lin_alg/sparse_coo.hpp"

YAFEL_NAMESPACE_OPEN

sparse_coo::sparse_coo() {
  init();
}

sparse_coo::sparse_coo(unsigned n) {
  init();
  rows = n;
  cols = n;
}

sparse_coo::sparse_coo(unsigned m, unsigned n) {
  init();
  rows = m;
  cols = n;
}

sparse_coo::sparse_coo(const sparse_coo & src) {
  
  row_index = new unsigned[src.getSize()];
  col_index = new unsigned[src.getSize()];
  data = new double[src.getSize()];
  
  for(unsigned i=0; i<src.getSize(); i++) {
    row_index[i] = src.row_index[i];
    col_index[i] = src.col_index[i];
    data[i] = src.data[i];
  }
  
  size = src.getSize();
  capacity = src.getSize();
  rows = src.getRows();
  cols = src.getCols();
  oldNNZ = src.oldNNZ;
  isSorted = src.isSorted;
  consistentNNZ = src.consistentNNZ;
}

sparse_coo & sparse_coo::operator=(const sparse_coo & rhs) {

  if(this == &rhs) {
    return *this;
  }
  
  unsigned *tmp_ri, *tmp_ci;
  double *tmp_data;
  unsigned sz = rhs.getSize();
  
  tmp_ri = new unsigned[sz];
  tmp_ci = new unsigned[sz];
  tmp_data = new double[sz];
  
  for(unsigned i=0; i<sz; ++i) {
    tmp_ri[i] = rhs.row_index[i];
    tmp_ci[i] = rhs.col_index[i];
    tmp_data[i] = rhs.data[i];
  }
  
  delete[] row_index;
  delete[] col_index;
  delete[] data;
  
  rows = rhs.getRows();
  cols = rhs.getCols();
  size = sz;
  capacity = sz;
  oldNNZ = rhs.oldNNZ;
  isSorted = rhs.isSorted;
  consistentNNZ = rhs.consistentNNZ;
  row_index = tmp_ri;
  col_index = tmp_ci;
  data = tmp_data;
  
  return *this;
}


sparse_coo::sparse_coo(const sparse_csr & csr) {
  
  init();
  rows = csr.getRows();
  cols = csr.getCols();

  for(unsigned row=0; row<rows; ++row) {
    for(unsigned k=csr.row_ptr[row]; k<csr.row_ptr[row+1]; ++k) {
      unsigned col = csr.col_index[k];
      double val = csr.data[k];
      
      this->add(row, col, val);
    }
  }
  
}


sparse_coo::~sparse_coo() {
  delete[] row_index;
  delete[] col_index;
  delete[] data;
}

void sparse_coo::init() {
  row_index = new unsigned[default_capacity];
  col_index = new unsigned[default_capacity];
  data  = new double[default_capacity];
  rows = 0;
  cols = 0;
  size = 0;
  capacity = default_capacity;
  isSorted = true;
  consistentNNZ = true;
  oldNNZ = 0;
}

unsigned sparse_coo::nnz() {
  
  if(consistentNNZ) {
    return oldNNZ;
  }
  
  consistentNNZ = true;
  sort();

  unsigned NNZ = (size!=0) ? 1 : 0;
  for(unsigned i = 1; i<size; i++) {
    if (!( (row_index[i]==row_index[i-1]) &&
	   (col_index[i]==col_index[i-1])) ) {
	++NNZ;
    }
  }
  
  oldNNZ = NNZ;
  return NNZ;
}

void sparse_coo::add(unsigned row, unsigned col, double value) {
  //track consistency of matrix state
  isSorted = false;
  consistentNNZ = false;

  //track correct size of matrix
  if(row+1 > rows) {
    rows = row+1;
  }
  if(col+1 > cols) {
    cols = col+1;
  }
  
  //memory management
  if(size==capacity) {
    unsigned *tmp_row = row_index;
    unsigned *tmp_col = col_index;
    double *tmp_data = data;
    
    row_index = new unsigned[2*capacity];
    col_index = new unsigned[2*capacity];
    data = new double[2*capacity];
    for(unsigned i=0; i<size; i++) {
      row_index[i] = tmp_row[i];
      col_index[i] = tmp_col[i];
      data[i] = tmp_data[i];
    }
    delete[] tmp_row;
    delete[] tmp_col;
    delete[] tmp_data;
    capacity = capacity*2;
  }
  
  row_index[size] = row;
  col_index[size] = col;
  data[size] = value;
  
  size++;
}

void sparse_coo::sort() {
  if(isSorted) {
    return;
  }
  isSorted = true;
  unsigned *Brow = new unsigned[size];
  unsigned *Bcol = new unsigned[size];
  double *Bval = new double[size];
  
  mergeSort(0, size-1, Brow, Bcol, Bval); 
  
  delete[] Brow;
  delete[] Bcol;
  delete[] Bval;
  
}

void sparse_coo::mergeSort(unsigned start, unsigned end, unsigned *Brow, unsigned *Bcol, double *Bval) {
  
  if(end-start < INSERTION_SORT_THRESHOLD) {
    insertionSort(start, end);
    return;
  }
  
  unsigned middle = start + (end-start)/2;
  
  mergeSort(start, middle, Brow, Bcol, Bval);
  mergeSort(middle+1, end, Brow, Bcol, Bval);
  merge(start, middle, end, Brow, Bcol, Bval);
  
  for(unsigned i=start; i<=end; i++) {
    row_index[i] = Brow[i];
    col_index[i] = Bcol[i];
    data[i] = Bval[i];
  }
}

void sparse_coo::merge(unsigned start, unsigned middle, unsigned end, unsigned *Brow, unsigned *Bcol, double *Bval) {
  unsigned i0 = start;
  unsigned i1 = middle+1;

  for(unsigned j=start; j<=end; j++) {
    if(i0 <= middle &&
       (i1>end || row_index[i0]<row_index[i1] ||
	(row_index[i0]==row_index[i1] && col_index[i0]<col_index[i1]))) {
      Brow[j] = row_index[i0];
      Bcol[j] = col_index[i0];
      Bval[j] = data[i0];
      i0++;
    }
    else {
      Brow[j] = row_index[i1];
      Bcol[j] = col_index[i1];
      Bval[j] = data[i1];
      i1++;
    }
    
  }

}

void sparse_coo::insertionSort(unsigned start, unsigned end) { 

  for(unsigned i=start; i<=end; i++) {
    unsigned newRow = row_index[i];
    unsigned newCol = col_index[i];
    double newVal = data[i];
    unsigned holePos = start;

    for(unsigned j=i; j> start; j--) {
      if(newRow < row_index[j-1] && j > start) {
	continue;
      }
      else {
	if(newRow==row_index[j-1] && newCol < col_index[j-1]) {
	  continue;
	}
	holePos = j;
	break;
      }
    }
    
    for(unsigned j=i; j>holePos; j--) {
      row_index[j] = row_index[j-1];
      col_index[j] = col_index[j-1];
      data[j] = data[j-1];
    }
    row_index[holePos] = newRow;
    col_index[holePos] = newCol;
    data[holePos] = newVal;
  }
     
}

void sparse_coo::compress() {
  
  unsigned compressed_size = nnz();
  
  unsigned *temp_row_index = row_index;
  unsigned *temp_col_index = col_index;
  double *temp_data = data;
  
  row_index = new unsigned[compressed_size];
  col_index = new unsigned[compressed_size];
  data = new double[compressed_size];
  
  unsigned comp_i = 0;
  row_index[0] = temp_row_index[0];
  col_index[0] = temp_col_index[0];
  data[0] = temp_data[0];
  for(unsigned i=1; i<size; ++i) {
    if( (temp_row_index[i]==temp_row_index[i-1]) &&
	(temp_col_index[i]==temp_col_index[i-1]) ) {
      data[comp_i] += temp_data[i];      
    }
    else {
      //next element in new arrays
      comp_i++;
      
      row_index[comp_i] = temp_row_index[i];
      col_index[comp_i] = temp_col_index[i];
      data[comp_i] = temp_data[i];
    }
    
  }
  
  delete[] temp_row_index;
  delete[] temp_col_index;
  delete[] temp_data;
  
  size = compressed_size;
  capacity = compressed_size;
  
}

void sparse_coo::print_sparse() {
  for(unsigned i=0; i<size; i++) {
    printf("(%d,\t%d,\t%f)\n", row_index[i], col_index[i], data[i]);
  }
}

void sparse_coo::load_from_file(const char *fname) {
  
  delete[] row_index;
  delete[] col_index;
  delete[] data;

  row_index = new unsigned[10];
  col_index = new unsigned[10];
  data = new double[10];
  
  size = 0;
  capacity = 10;
  
  FILE *f;
  f = fopen(fname,"r");

  if(f==NULL){
    perror("Failed to open file");
    exit(1);
  }
  
  //fscanf(f, "%d, %d", &rows, &cols);

  while (true) {
    unsigned nextRow, nextCol;
    double nextVal;
    unsigned n = fscanf(f, "%d, %d, %lf", &nextRow, &nextCol, &nextVal);
    if(n != 3 || feof(f)) {
      break;
    }
    add(nextRow, nextCol, nextVal);
    
  }
  
}

YAFEL_NAMESPACE_CLOSE
