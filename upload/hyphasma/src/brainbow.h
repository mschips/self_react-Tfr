/* August 12, 2015 Michael Meyer-Hermann:
 * Class for analysis of cell mutations.
 *
 * Initialisation:
 * Define an brainbow object to make the settings.
 *
 * Usage:
 * Write every birth or death event by calling:
 * long write_brainbow(double time, bool founder, bool birth, long mother_index, long ss_position);
 * this returns the index in the brainbow vector at which the event was saved.
 *
 */

#ifndef i_BRAINBOW
#define i_BRAINBOW

#include <vector>
#include <iostream>
#include <fstream>
#include <math.h>
#include <algorithm>
using namespace std;

class brainbow {
public:
  // define different possible modes of how to attribute colours to cells:
  enum staining_modes {
    gabriel_4colours,
    gabriel_10colours,
    gabriel_10colours_2founder,
    Ncolours_random,
    N_staining_modes
  };
  
  // Starting time of track
  double trackfrom;
  // track until
  double trackuntil;
  
  // construct by default that events are collected all the time
  // and analysis is done with all written events:
  brainbow();
  // As in experiment, analysis might be restricted to particular periods.
  // Then use this one and read_from_file();
  // write(..) events will also check, whether the event lies within the period.
  brainbow(double from, double until);
  ~brainbow();
  
  // Sets the private variable Ncolour to a value
  void set_Ncolour(int c);
  void set_stain_time(double t);
  void set_stain_fraction(double f);
  void set_stain_dt(double t);
  void set_lineage_time(double t);
  void set_tamoxifen_decay(double t);
  void set_include_seclargest_colour(bool t);
  void set_include_late_founders_in_lineages(bool t);
  void set_tamoxifen_stop_time(double t);
  void set_staining_mode(staining_modes b);
  void set_allow_restaining(bool b);
  void set_origin_same_colour(bool b);
  void set_merge_identical_origin(bool b);
  void set_cell_subset_size(int x);
  void set_do(bool x, bool y, bool z); // set the function of make_brainbow here
  void set_dt_Muller(double dt);
  void set_Muller_maxtime(double maxt);
  
  // write events, called from outside, returns the event index (-1 if not written)
  long write(double xtime, bool xfounder, bool xbirth, long xmother_index, long xss_position);
  void write_number_of_different_cell_types(double deltat, int Nevals);
  
  // save all events in a file that can be read with read_from_file() for further analysis:
  void save_to_file();
  void save_Muller_diagram();
  
  // used to read the data from file ./brainbow.out for analysis:
  void read_from_file();
  
  // generates the file clonality.out, clonality_real.out, clone_fractions(_rand).out
  // set the detailed options with set_do and other set_xxx functions
  void make_brainbow(double * evaluation_times, int Nevaluations);
  // perform analysis (default version)
  void analysis();
  
private:
  struct brainbow_data {
    // Absolute time of the position:
    double time;
    // flag for founder cells
    bool founder;
    // flag for lineage start
    bool lineage;
    // flag for birth (true=birth, false=death)
    bool birth;
    // pointer to the index of the daughter cells
    long daughter1_index, daughter2_index;
    // mother index
    long mother_index;
    // position in shape space
    long ss_position;
    // brainbow colour
    int colour;
  };
  
  // flag to overwrite clonality.out and other files
  // which serves to reduce the analysis to parts:
  bool do_clonality, do_clone_fractions, do_clone_fractions_rand;
  
  // define a vector of brainbow_data
  vector<brainbow_data> rb_list;
  // initial length of rb_list:
  long rb_dim_ini;
  // infinitesimal time for avoidance of rounding errors
  static constexpr double infinitesimal_time = 1.e-08;
  
  long get_founder_cell_index(long &cell_index);
  long get_lineage_cell_index(long &cell_index);
  long get_clone_of_founder(long &cell_index);
  long get_index_for_time(double &time);
  int cell_subset_size;
  vector<long> get_cell_subset(vector<long> all_cells);
  void get_cell_indices_at_time(double &time, vector<long> &cellind);
  int get_number_of_different_cell_types(double time, ofstream& diverse);
  long flist_index(long i, vector<long>& v, long max);
  double dt_Muller, Muller_maxtime;
  
  // staining procedures
  int Ncolours;
  double stain_fraction;
  staining_modes how2do_staining;
  int get_random_colour(int &Ncols);
  void transfer_colour2daughters(long cell_index);
  double tamoxifen_decay, tamoxifen_stop_time;
  bool include_late_founders_in_lineages, include_seclargest_colour;
  double get_stain_amplitude();
  double get_stain_probability(double &time, double &amplitude);
  // void stain(double time, double& amplitude, int& Ncols);
  void stain_all_founder(double stain_prob, int &Ncolours);
  void stain(double time, double &stain_prob, int &Ncols, int * Nfc,
	     vector<long> &founder_list,
	     vector<int> &projection);
  double stain_time, stain_dt, lineage_time;
  bool merge_identical_origin, origin_same_colour, allow_restaining;
  void stain(vector<long> &founder_list,
              vector<int> &projection);
  
  unsigned long evaluate_colours(long * colour_frequency, int &Ncols, vector<long> &cellind);
  void get_fractions(long Ncells, long * frequency, double * fraction, int &Norigins);
  void put_fraction_header(ofstream &xf);
  void save_fractions(double time, long cell_num, int &Ncolumns, int &maxNcolumns,
		      double * fraction, ofstream &fx);
  void save_clonality(double &evaluation_time, ofstream &outfile,
		      long unsigned &Nstained_cells, double &largest_colour_fraction,
		      int &dominant_origin_in_colour, long &colour_ss_position,
                       long &Nreal_cells, double &largest_origin_fraction,
		      int &dominant_origin_in_real, long &real_ss_position,
		      bool clonal);
  int find_largest(long * frequency, int &Nentries);
  int find_seclargest(long * frequency, int ilargest, int &Nentries);
  
  int get_all_lineages(double time, vector<long> &linlist);
  int add_late_founders_to_lineages(double fromtime,  
				     vector<long> &foundlist, vector<long> &linlist);
  int get_all_founders(double time, vector<long> &founder_list);
  int merge_origins(int Norigin,
		    vector<long> &origin_list,
		    vector<int> &origin_projection,
		    vector<long> &origin_short_list);
  void project_origin(int &listpos, vector<int> &projection);
  int get_origin_arraypos(long origin_index, int Norigins, vector<long> &origin_list);
  int get_dominant_origin(double time, vector<long> &cellind, // time and event indices
			  int Norigins, vector<long> &origin_list, // list of founder/lineage and
			  // length
			  vector<int> &projection, // map from origin_list to its short_list
			  int largest_colour, bool use_founder, // mode of action
			  long * freq); // result
};

#endif
