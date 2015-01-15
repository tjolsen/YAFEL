#include "output/VTKTimeOutput.hpp"
#include <cstdio>

VTKTimeOutput::VTKTimeOutput() {
  fnames.clear();
  times.clear();
}

void VTKTimeOutput::addDataFile(const std::string &fname, double time) {
  fnames.push_back(fname);
  times.push_back(time);
}

void VTKTimeOutput::write(const std::string &fname) {
  
  if(fnames.size() != times.size()) {
    perror("VTKTimeOutput::write() : Inconsistent data member lengths. Aborting.");
    return;
  }
  
  
  if(fnames.size() == 0) {
    perror("VTKTimeOutput::wite() : No data to write. Aborting.");
    return;
  }
  
  
  FILE *fp = fopen(fname.c_str(),"w");
  if(fp == NULL) {
    perror("VTKTimeOutput::write() : Could not open file.");
    exit(1);
  }
  
  fprintf(fp, "<?xml version=\"1.0\"?>\n");
  fprintf(fp, "<VTKFile type=\"Collection\" version=\"0.1\">\n");
  
  fprintf(fp, "<Collection>\n");
  
  for(unsigned i=0; i<fnames.size(); ++i) {
    fprintf(fp, "<DataSet timestep=\"%f\" part=\"0\" file=\"%s\"/>\n",
	    times[i], fnames[i].c_str());
  }

  fprintf(fp, "</Collection>\n");
  fprintf(fp, "</VTKFile>");
  
  fclose(fp);

}
