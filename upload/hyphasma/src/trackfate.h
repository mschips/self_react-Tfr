/* 03.10.2017 MMH:
 * Class for analysis of centrocyte fate development.
 * Initialisation: Define an instance of the fateTRACK-class for each CC.
 * Analysis: Call global analysis files which recollect all saved CCs.
 */

#ifndef i_fatetrack
#define i_fatetrack

#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <string.h>
using namespace std;

/*
enum fate_action_types {FATECCstart, FATEcontactFDC, FATEleaveFDC, 
			FATEsignalprogress, FATEcontactTFH, FATEleaveTFH, 
			FATEselected, FATEcelldeath};
*/

class fateTRACK {
 public:
  fateTRACK();
  ~fateTRACK();
  void mk_fate_track_legend();
  bool mk_fate_track_file();
  void mk_new_fate_track();
  void clear_fate_track();
  void show_track();
  void write_fate(short, int, double, double, int, double, double, double, double, double, 
		  double, double, double, int, double, double, double, double, double);
  void write2file();
  
 private:
  struct fatetrack_data {
    // state (= CC state from enum list in cellthis.h)
    // unselected,contact,FDCselected,TCcontact,selected,apoptosis
    // MMH2Marta unselected,contact,FDCselected,TCcontact,TFRcontact,selected,apoptosis
    short cellstate;
    // position on shape space
    int pos_ss;
    // affinity to nearest antigen
    double affinity;
    // absolute time of saving the fate:
    double t;
    // number of contacts with FDCs
    int nFDCcontacts;
    // pMHC level (amount of collected and processed antigen)
    double pMHC;
    // duration of being in state FDCselected
    double FDCselected_clock;
    // tc_search_time
    double tc_search_period;
    // t-b-interaction time
    double t_b_interaction_time;
    // FoxO level
    double FoxO;
    // FoxO upregulation rate
    double FoxO_upreg;
    // mTORC1 level
    double mTORC1;
    // ICOSL level
    double ICOSL;
    // number of contacts with Tfh
    int nTCcontacts;
    // Tfh signal at time of selection
    double TFHsignal;
    // AUC of Tfh signal
    double TFHsignalAUC;
    // Induced DND
    double DND;
    // BCRsignal level (same as nFDCcontacts without retaining antigen)
    double BCRsignal;
    // cMyc level
    double cMyc;
  };
  long trackindex;
  // Standard size of cell fate tracks:
  static int TRACKDIM;
  // absolute cell index
  static long abs_cell_index;
  static long abs_track_index;
  // Filename to save the data:
  static string filename;
  // For every tracked object a dynamic array of positions and times is defined
  vector<fatetrack_data> fate_list;
};

#endif
