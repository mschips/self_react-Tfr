#include "histo_monitor.h"
#include <math.h>
histo_monitor::histo_monitor() { };
histo_monitor::~histo_monitor() { }
void histo_monitor::initialize(double max_value, int set_bin_max, 
			       double tmin, double tmax, 
			       std::string set_filename,
			       std::string set_quantity_name) {
  //  const unsigned int bin_max_tmp = set_bin_max;
  ndays = int((tmax - tmin) / (24.0)) + 1;
  thisday = int((tmin / 24.0));
  bin_max = set_bin_max;
  dvalue = max_value / double(bin_max);
  std::vector<int> vtottmp(bin_max + 1);
  histo_total = vtottmp;
  vtottmp.clear();
  std::vector <std::vector<int> > vtmp(ndays, std::vector<int>(bin_max + 1));
  histo_day = vtmp;
  vtmp.clear();
  for (unsigned int day = 0; day < histo_day.size(); day++) {
    for (unsigned int nc = 0; nc < histo_day[day].size(); nc++) {
      histo_day[day][nc] = 0;
    }
  }
  for (unsigned int nc = 0; nc < histo_total.size(); nc++) {
    histo_total[nc] = 0;
  }
  filename = set_filename;
  quantity_name = set_quantity_name;
  /* When tmin==tmax, a file is initialized for further appending later on.
   * In this case the first columns is meant as an externally provided time stamp.
   * Writing is then performed by calling write_mean_SD_to_file(double time).
   * All operations are performed on histo[day=0][bin].
   */
  if (tmin == tmax) {
    if (tmin > 0) { 
      std::cerr << "tmin>0 but tmax==tmin not allowed in histo_monitor --> Exit\n";
      exit(1);
    }
    std::ofstream f(filename);
    f << "! time[hr] : " << quantity_name << " (mean : SD : N)\n";
    f.close();
  }
}
int histo_monitor::get_bin(double value) {
  int bin = int(value / dvalue + 0.5);
  if (bin > int(bin_max)) {
    std::cout << "\n                               "
              << "In int histo_monitor::get_bin(value="
              << value << ", " << quantity_name << "): bin = "
              << bin << " larger than max = "
              << bin_max << " was generated.\n";
    bin = bin_max;
  }
  return bin;
}
double histo_monitor::get_value(int bin) {
  double value = double(bin) * dvalue;
  return value;
}
void histo_monitor::add(double value) {
  ++histo_day[thisday][get_bin(value)];
}
void histo_monitor::rm(double value) {
  --histo_day[thisday][get_bin(value)];
}
double histo_monitor::get_mean(int day, long& totaln) {
  double mean = 0.;
  totaln = 0; 
  for (unsigned n = 0; n <= bin_max; ++n) {
    mean += get_value(n) * histo_day[day][n];
    totaln += histo_day[day][n];
  }
  if (totaln > 0) { mean /= double(totaln); }
  return mean;
}
double histo_monitor::get_mean(long& totaln) {
  // Same as above for 1D-histograms
  return get_mean(0, totaln);
}
double histo_monitor::get_SD(int day, double mean, long totaln) {
  double sd = 0., v = 0.;
  for (unsigned n = 0; n <= bin_max; ++n) { 
    v = ( get_value(n) - mean );
    sd += histo_day[day][n] * v * v;
  }
  if ( totaln > 1 ) { 
    sd /= double(totaln - 1); 
    sd = sqrt(sd);
  } else { sd = 0.; }
  return sd;
}
double histo_monitor::get_SD(double mean, long totaln) {
  // Same as above for 1D-histograms
  return get_SD(0, mean, totaln);
}
void histo_monitor::reset(int day) {
  for (unsigned n = 0; n <= bin_max; ++n) { histo_day[day][n] = 0; }
  // Note that histo_total is not updated in this routine.
}
void histo_monitor::reset() {
  // Same as above for 1D-histograms
  reset(0);
}
void histo_monitor::show_all() {
  std::cout << "bin : value : total_frey : day 0 : ... : day last\n";
  for (unsigned ni = 0; ni <= bin_max; ++ni) {
    std::cout << ni << "  " << get_value(ni) << "    " << histo_total[ni] << "    ";
    for (unsigned ti = 0; ti < ndays; ++ti)
      { std::cout << histo_day[ti][ni] << "  "; }
    std::cout << "\n";
  }
}  
void histo_monitor::update_day() {
  // save this day in histo_total
  for (unsigned nc = 0; nc < histo_total.size(); nc++) 
    { histo_total[nc] += histo_day[thisday][nc]; }
  // switch to next day
  ++thisday;
}
void histo_monitor::write2file() {
  std::ofstream outfile(filename.c_str());
  outfile << "! bin : " << quantity_name << " : total freq : day 0 : ... : day last\n";
  for (unsigned ni = 0; ni <= bin_max; ++ni) {
    outfile << ni << "  " << get_value(ni) << "    " << histo_total[ni] << "    ";
    for (unsigned ti = 0; ti < ndays; ++ti) 
      { outfile << histo_day[ti][ni] << "  "; }
    outfile << "\n";
  }
  outfile.close();
}
void histo_monitor::write_mean2file(double time) {
  // This routine is used with an external time stamp and a 1-dimensional histogram.
  long totalN = 0;
  double mean = get_mean(totalN);
  double SD = get_SD(mean, totalN);
  std::ofstream f(filename, std::ofstream::app);
  f << time << "   " << mean << "   " << SD << "   " << totalN << "\n";
  f.close();
  // The 1D-histogram is set to zero after writing!
  reset();
}
