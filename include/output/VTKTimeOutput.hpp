#ifndef _YAFEL_VTKTIMEOUTPUT_HPP
#define _YAFEL_VTKTIMEOUTPUT_HPP

#include "yafel_globals.hpp"
#include <string>
#include <vector>


class VTKTimeOutput {
  
private:
  std::vector<std::string> fnames;
  std::vector<double> times;
  
public:
  VTKTimeOutput();
  void addDataFile(const std::string &fname, double time);
  void write(const std::string &fname);
  
};

#endif
