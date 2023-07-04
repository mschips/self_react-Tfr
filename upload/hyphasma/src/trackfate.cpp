#include "trackfate.h"
#include <math.h>

long fateTRACK::abs_cell_index = 0;
long fateTRACK::abs_track_index = 0;
string fateTRACK::filename = "trackfate.out";
int fateTRACK::TRACKDIM = 2000;

fateTRACK::fateTRACK() {
}
fateTRACK::~fateTRACK() {
}
void fateTRACK::mk_fate_track_legend() {
  cout << "Generate new trackfate{,_legend}.out files ... \n";
  ofstream cols;
  cols.open("trackfate_legend.out");
  cols << "01 = absolute cell index\n" 
       << "02 = cellstate:\n"
    // MMH2Marta << "     0=unselected,1=contact,2=FDCselected,3=TCcontact,4=selected,5=apoptosis\n"
    // This will have impact on the trackfate-analysis which will not be functional anymore
    // We might add the autoreactivity flag (selfMutation) here as an additiona column.
          ///MS (TODO!)
       << "     0=unselected,1=contact,2=FDCselected,3=TCcontact,4=TRFcontact,5=selected,6=apoptosis\n"
       << "03 = position in shape space\n"
       << "04 = best affinity\n"
       << "05 = time in hours\n"
       << "06 = number of FDC contacts\n"
       << "07 = amount of pMHC\n"
       << "08 = clock since state FDCselected\n"
       << "09 = duration of the Tfh search period in hours\n"
       << "10 = duration of Tfh-B-interaction\n"
       << "11 = FoxO level\n"
       << "12 = FoxO upregulation rate\n"
       << "13 = mTORC1 level\n"
       << "14 = ICOSL level\n"
       << "15 = number of successful contacts to Tfh\n"
       << "16 = amount of collected TFH signals\n"
       << "17 = AUC of Tfh signal\n"
       << "18 = DND\n"
       << "19 = BCRsignal\n"
       << "20 = cMyc\n";
  cols.close();
}
bool fateTRACK::mk_fate_track_file() {
  // Check whether the file to save those objects exists:
  ofstream f(filename.c_str());
  bool fileexists = f.good();
  // if not, open it for further append of data later on
  if (fileexists == false) {
    f.open(filename.c_str());
    f.close();
  }
  return not(fileexists); // true when the file was newly generated, false when it existed before
}
void fateTRACK::mk_new_fate_track() {
  // Initialise a fate_list for a single cell object
  if (fate_list.size() > 0) { 
    //cerr << "Cleared fate_list with index " << trackindex << " in mk_new_fate_track()\n";
    fate_list.clear(); 
  }
  fate_list.reserve(TRACKDIM); 
  // this does not allocate memory if the vector has sufficient capacity
  trackindex = abs_track_index;
  //cerr << "Made " << trackindex << " mk_new_fate_track()\n";
  ++abs_track_index;
}
void fateTRACK::clear_fate_track() {
  ++abs_cell_index;
  fate_list.clear();
}
void fateTRACK::show_track() {
  for (unsigned int i = 0; i < fate_list.size(); i++) {
    cout << abs_cell_index << "  "
	 << trackindex << "  "
	 << fate_list[i].cellstate << "  "
	 << fate_list[i].pos_ss << "  "
	 << fate_list[i].affinity << "  "
	 << fate_list[i].t << "  "
	 << fate_list[i].nFDCcontacts << "  "
	 << fate_list[i].pMHC << "  "
	 << fate_list[i].FDCselected_clock << "  "
	 << fate_list[i].tc_search_period << "  "
	 << fate_list[i].t_b_interaction_time << "  "
	 << fate_list[i].FoxO << "  "
	 << fate_list[i].FoxO_upreg << "  "
	 << fate_list[i].mTORC1 << "  "
	 << fate_list[i].ICOSL << "  "
	 << fate_list[i].nTCcontacts << "  "
	 << fate_list[i].TFHsignal << "  "
	 << fate_list[i].TFHsignalAUC << "  "
	 << fate_list[i].DND << "  "
	 << fate_list[i].BCRsignal << "  "
	 << fate_list[i].cMyc << "\n";
  }
}
void fateTRACK::write_fate(short c, int pos_ss, double affinity,
			   double t, int nFDCcontacts, double pMHC, double FDCclock, 
			   double tc_search_period, double tbinter_time,
			   double foxo, double foxo_upreg, double mtor, 
			   double ICOSL, int nTCcontacts, 
			   double TFHsignal, double TFHsignalAUC, double DND,
			   double BCRsignal, double cMyc) {
  /* Adds an entry to the vector with cell fate development */
  fatetrack_data tmp;
  tmp.cellstate = c;
  tmp.pos_ss = pos_ss;
  tmp.affinity = affinity;
  tmp.t = t;
  tmp.nFDCcontacts = nFDCcontacts;
  tmp.pMHC = pMHC;
  tmp.FDCselected_clock = FDCclock;
  tmp.tc_search_period = tc_search_period;
  tmp.t_b_interaction_time = tbinter_time;
  tmp.FoxO = foxo;
  tmp.FoxO_upreg = foxo_upreg;
  tmp.mTORC1 = mtor;
  tmp.ICOSL = ICOSL;
  tmp.nTCcontacts = nTCcontacts;
  tmp.TFHsignal = TFHsignal;
  tmp.TFHsignalAUC = TFHsignalAUC;
  tmp.DND = DND;
  tmp.BCRsignal = BCRsignal;
  tmp.cMyc = cMyc;
  fate_list.push_back(tmp);
}
void fateTRACK::write2file() {
  /* Writes a finished fate track to the file <filename> */
  ofstream f;
  f.open(filename.c_str(), ofstream::app);
  for (unsigned long i=0; i < fate_list.size(); i++) {
    f << abs_cell_index << "  "
      //      << trackindex << "  "
      << fate_list[i].cellstate << "  "
      << fate_list[i].pos_ss << "  "
      << fate_list[i].affinity << "  "
      << fate_list[i].t << "  "
      << fate_list[i].nFDCcontacts << "  "
      << fate_list[i].pMHC << "  "
      << fate_list[i].FDCselected_clock << "  "
      << fate_list[i].tc_search_period << "  "
      << fate_list[i].t_b_interaction_time << "  "
      << fate_list[i].FoxO << "  "
      << fate_list[i].FoxO_upreg << "  "
      << fate_list[i].mTORC1 << "  "
      << fate_list[i].ICOSL << "  "
      << fate_list[i].nTCcontacts << "  "
      << fate_list[i].TFHsignal << "  "
      << fate_list[i].TFHsignalAUC << "  "
      << fate_list[i].DND << "  "
      << fate_list[i].BCRsignal << "  "
      << fate_list[i].cMyc << "\n";
  }
  f.close();
}

