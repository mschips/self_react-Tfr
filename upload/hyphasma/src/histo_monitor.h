#ifndef i_histo_monitor
#define i_histo_monitor
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
/** Allows to monitor the evolution of some variable in the code as histogram.
 ** The variable is projected onto regular bins.
 ** The histogram has two dimensions x[d][n], where d is a time variable and
 ** n denotes the bin.
 ** A new histogram is initialized with the same-named routine.
 ** As parameters, are needed:
 ** (maximal possible value of the variable, number of bins,
 **  minimum and maximum time in days,
 **  name of the output filename, name of the variable in the output file)
 ** Note that entry of 0,0 for the min and max times generates a one-dimensional histogram.
 **/

class histo_monitor {
 public:
  histo_monitor();
  ~histo_monitor();
  void initialize(double, int, double, double, std::string, std::string);
  void add(double value);
  void rm(double value);
  double get_mean(int day, long& totaln);
  double get_SD(int day, double mean, long totaln);
  void reset(int day);
  void update_day();
  void write2file();
  void write_mean2file(double time);
  void show_all();
 private:
  std::vector <std::vector<int> > histo_day;
  std::vector<int> histo_total;
  unsigned int bin_max;
  double dvalue;
  unsigned int thisday, ndays;
  int get_bin(double value);
  double get_value(int bin);
  double get_mean(long& totaln);
  double get_SD(double mean, long totaln);
  void reset();
  std::string filename, quantity_name;
};

#endif
