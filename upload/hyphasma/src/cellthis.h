/*
 * cellthis.h
 * --------------------------------------------------
 * Deklaration von spezifischen Zelltypen unter Verwendung der Routinen  in cell.h
 * --------------------------------------------------------------------------------
 *
 * Struktur der Deklaration:
 *
 * 1) Eintragen des neuen Typs (XX) in den states-enumerator in point.h
 *
 * 2) Erstelle neue Klasse cellXX in dieser Datei:
 *
 * class cellXX : public cell {...} // Zelle aus einem Fragment
 *
 * oder
 *
 * class cellXX : public frag_cell { // Zelle mit mehr Fragmenten
 * public:
 *  cellXX() { }               // notwendiger Konstruktur
 *  cellXX(const cellXX& x);   // notwendiger Selbst-Konstruktur mit Zuordnung von Werten
 *  cellXX(cellYY& x);         // optional fuer Zuordnung zwischen verwandten Zelltypen
 *  ~cellXX();                 // notwendiger Destruktur
 *  void ini(const long& i, const long& li, lattice& l, AffinitySpace& shape);
 *                             // schreibt eine neue Zelle auf lattice und AffinitySpace ein.
 *                 // Dies koennte in den Konstruktur wenn jede erzeugte
 *                 // Zelle auch auf lattice und AffinitySpace erscheinen soll.
 * }
 *
 * 3) Optional Erstellung eines zellspezifischen Zell-enumerators zur Unterscheidung
 * verschiedener Zellzustaende:
 *
 * enum XXstates{state1, state2};
 *
 * und in die cellXX-Klasse Zustands variable einfuehren:
 *
 * XXstates state;
 *
 * 4) Einfuehrung der statischen Variablen, die die Prozesse
 * quantitativ steuern (Wahrscheinlichkeiten p_):
 *
 * private:
 * static double p_pro,             // Proliferation
 *              max_pro_distance,  // fuer 1-Fragment: pushing distance
 *              p_mut,             // Mutation (verwendet AffinitySpace)
 *      p_apo,             // Apoptose
 *      p_difu,            // Diffusion/Motility
 *      p_difu_width,
 *      distance_tolerace, // Formparameter bei Diffusion
 *      p_grow,            // Wachstum
 *      p_shrink,          // Schrumpfung
 *      p_dif;             // Differenzierung
 * public:
 * static int target_volume; // Zielzellvolumen
 * static void set_statics(const Parameter& par, lattice& l, ofstream& ana);
 * static void set_statics(const double& time, const Parameter& par);
 *
 * Es ist zu beachten, dass nur die fuer den Zelltyp relevanten
 * p_-Variablen zu verwenden sind.
 *
 * target_volume ist nicht unbedingt von aussen sichtbar zu definieren.
 *
 * set_statics() muss von aussen zu steuern sein. Die Argumente haengen
 * vom neuen Typ ab. Die erste Variante dient der Initialisierung, wobei
 * die Resultate in ana rausgeschrieben werden. Mit der zweiten Variante
 * koennen die Parameter waehrend der Laufzeit geaendert werden.
 *
 * Will man die p_ fuer einige Prozesse dynamisch gestalten, sind die
 * entsprechenden p_ nicht static zu deklarieren.
 *
 * 4a) Die Werte p_ und andere werden aus Parameter eingelesen. Die Parameter
 * muessen dort deklariert und definiert sein. Dies ist in setparam.x zu tun.
 *
 * 5) Einfuegung der erwuenschten Prozesse (Bsp nur aus cell und frag_cell).
 * Bei einigen Prozessen ist die Abfrage, ob diese durchgefuehrt werden
 * getrennt zu programmieren, etwa wie durch ask_mitosis().
 *
 * public:
 * // Diffusion fuer mehrfragmentige Zellen
 * double diffuse(const long& li, lattice& l)
 *  { return fragdiffuse(XX,distance_tolerance,li,p_difu,l); }
 * // Diffusion fuer einfragmentige Zellen
 * short diffuse(const long& li, lattice& l)
 *  { return do_diffuse(XX,li,l); }
 *
 * // Proliferation
 * long ask_mitosis(long* pp, lattice& l) {
 *  if (volume==1 && target_volume==1) {
 *    return find_mitosis_place(p_pro,max_pro_distance,pp,l);
 *  } // im Fall der Klasse cell genuegt find_mitosis_place(...)
 *  else { // case of more fragment object
 *    // probabilistic decision if proliferation is done
 *    if ( volume>0.9*target_volume && drandom(1.) < p_pro ) return 0;
 *    // Proliferation is allowed if the total volume is near the target-volume of a cell!
 *    return 1;
 *  }
 * }
 * short mitosis(const long& i, const long& li,
 *      const long& newli, frag_cell& newCell,
 *      lattice& l)
 *  { return do_mitosis(XX,i,li,newli,newCell,l); }
 *
 * // Mutation
 * short mutate(AffinitySpace& shape)
 *  { return do_mutate(shape,p_mut); }
 *
 * // Wachstum und Schrumpfung (### letzteres nicht programmiert)
 * short grow(const long& li, lattice& l)
 *  { if (target_volume==1) return 1; // d.h. kein Wachstum!
 *    else return do_grow(XX,li,target_volume,p_grow,p_shrink,l); }
 *
 * // Kontaktfragen
 * short find_contact(states celltype, const long& i, lattice& l); //cell
 * short check_for_border(const long& i, lattice& l); //cell
 * short contact(states celltype, lattice& l); //frag_cell
 *
 * // Signal Produktion
 * void signal_production(const long& i, lattice& l) {
 *  // Produce signal for differentiation of CB to CC
 *  signal_secretion(i,sig_differ2CC,p_mksignal,l);
 *  // other signals may be included here
 * }
 *
 * 5a) Initialisierung einer neuen Zellliste fuer XX.
 *
 * 5b) Ansteuerung der Prozesse in cellman.x durch eine neue Routine
 * calc_XX() und eventuelle weitere Routinen, die Zelltypumwandlung
 * oder Zellerzeugung -- also Zell-Listen-relevante Prozesse betreffen.
 *
 * 6) Das sollte es sein.
 */

#ifndef i_cellthis
#define i_cellthis
#include "setparam.h"
#include "antibody.h"
#include "cell.h"
#include "ss.h"
#include "sequencespace.h"
#include "space.h"
#include "ode.h"
#include "trackfate.h"
#include "histo_monitor.h"
// <sb>
#include <array>
// </sb>

// ===============================================================
// ===============================================================
// ===============================================================
// =================== Specific cell types =======================
// ===============================================================
// ===============================================================

// ===============================================================
// ======================= cellCB ================================
// ===============================================================

class cellCC;
class cellFDC;
class cellTC;
class cellTFR;

class immunoglobulin_class {
  private:
   static const int matrix_dimension = nIg_classes * nIg_classes;
   static int get_matrix_index(const Ig_classes &i, const Ig_classes &j);
   static double switch_matrix[matrix_dimension];
  public:
   immunoglobulin_class();
   ~immunoglobulin_class();
   static short do_switch_classes;
   static void load_matrix(const Parameter &par, ofstream &ana);
   Ig_classes Ig_class;
   void class_switch();
   void set_class(const immunoglobulin_class &c);
};

// ### dies in cellCB integrieren und dann ueberall cellCB:: davor schreiben
// +++ muss enum static gemacht werden? Bekommt jedes Objekt ein neues enum?
// +++ nein, laut blog im Netz ist das nur eine Typendeklaration und jede Instanz
// +++ der Klasse belegt dafuer keinen Speicher.
enum centroblasts {
   cb_normal,cb_differentiate,
   cb_G1,cb_G0,cb_S,cb_G2,cb_M,        // keep these together!
   cb_divide,cb_stop_dividing,
   cb_statenumber
};

class cellCB: public frag_cell {
  private:
   static double
     limit_volume,
     p_mut,
     p_mut_after_selection,
     p_mut_after_dec_selection,
     p_mut_affinity_exponent,
     p_difu,
     p_difu_width,
     p_tension,
     p_grow,
     p_shrink,
     p_dif,
     p_dif_target,
     start_differentiate,
     max_pro_distance,
     diffusion_tolerance_min,
     diffusion_tolerance_steepness,
   // for receptors:
     receptors,
     receptor_activation,
     receptor_binding,
     receptor_dissociation,
   // movement and shape
     tension,
     elongation,
     K_elongation,
     smoothmove,
     persistence,
     v_slow_factor,
     p_switch_v,
     max_adhesion;
   static double dtphase[cb_statenumber];
   static const short phase_number = 5;
   static const int dtphase_resolution = 21;
   static int dtphase_frequency[phase_number][dtphase_resolution];
   static double fraction_of_phase_width;
   static short v_modi,n_v_states;
   // static long cell_cycle_delay,DEC205_cell_cycle_delay;
   static bool p_mut_affinity_dependent;
   static double total_n_of_DEC205_divisions;
   static double total_n_of_divisions;
   static bool transmit_CCdelay2cellcycle;
   centroblasts progress_cycle_phase();

   double cycle_state_time;
   double receptor_ligand;               // receptors with ligand bound

 public:
   static const int max_n_of_divisions = 12;
   static long cummulative_attributed_n_of_divisions[max_n_of_divisions + 1];
 private:
   static long attributed_n_of_divisions[max_n_of_divisions + 1];
   //MSchips
   static long selfBC_attributed_n_of_divisions[max_n_of_divisions + 1];

   static const int mutation_bins = 20;
   static long attributed_mutation_prob[mutation_bins];
   static long cummulative_attributed_mutation_prob[mutation_bins];
   static long attributed_mutation_prob_SELF[mutation_bins];
   static long cummulative_attributed_mutation_prob_SELF[mutation_bins];

  public:
   static int target_volume;
   static short receptor_use;
   static bool ag_loaded_CB_stop_mutation;
   centroblasts state;                   // Differentiation state of CB
   void transmit_CCdelay2cycle(double waited_time);
   bool DEC205,DEC205_ova,MHCdeficient,diff2output;
   // +++++++++++++++ OPTION pMHCdeficient ++++++++++++++++++++++++++++++++
   static constexpr double p_pMHCdeficient = 0.8;
   // +++++++++++++++ OPTION pMHCdeficient ++++++++++++++++++++++++++++++++
   void attribute_DEC205(double fracofpos);
   double retained_ag;
   vector<int> collected_ag_portions;
   bool iamhighag;
   //MSchips --> self could be moved to cell property -- actually all of them
   bool selfMutation, redeemed, inheritance, ccl3KO;
   static short TFR_CC_interaction_mode; ///modify this!!!
   static bool ag_loaded_CB_diff2output;
   static double asymmetric_polarity_index,smooth_PI;
   static double p_pro,delta_p_pro,average_seeder_affinity,p_CXCR4down;

   cellCB();
   cellCB(const cellCB &x);
   cellCB(const cellCC &x);
   ~cellCB();
   void destruct() { delete[] fragments; } // Destruktor des Fragment-Feldes
   void ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape);
   void make_CB_new();

   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   static void set_statics(const double &time, const Parameter &par);
   static void set_differentiation(const double &time);
   static bool SMOOTH_DIFFERENTIATION;
   static double smooth_differentiation_time;

   static double total_cell_cycle_duration();
   static bool fixed_number_of_divisions();
   double time_of_cycle_state_switch;
   int n_divisions2do;
   static void show_number_of_divisions(double&, ofstream&);
   static void show_cummulative_number_of_divisions();
   static void show_mutation_prob(double&, double&, double&, ofstream&);
   static void show_cummulative_mutation_prob();
   static void show_cell_cycle_phase_duration();
   //MSchips
   static void selfBC_number_of_divisions(double&, ofstream&);
   static void show_mutation_prob_SELF(double&, double&, double&, ofstream&);
   int nTFRcontacts;

   immunoglobulin_class IgX;
   static double IgE_factor_cellcycle, IgE_factor_divisions;

   // Aktionen:
   // =========
   void set_p_move();
   double move(const long &li, space &l, sigs &s, TRACK &td, double &time);

   void set_CXCR4expression();
   void resensitise4CXCL12(sigs &s);
   void adapt_specific_ag(double factor);

   static double ag_preloaded;
   void preload_with_ag();

   void set_remaining_divisions();
   /*  double inverse_erf(double x);
    * double get_sample_from_normal(const double& mean, const double& width);
    * double get_positive_sample_from_normal(const double& mean, const double& width);
    */
   double set_cycle_state_duration(centroblasts &s);
   static centroblasts get_virtual_cell_cycle_phase(double waited);
   static bool shiftCCdelay2CBcycle();
   long ask_mitosis(long * pp, space &l);
   short mitosis(const long &i, const long &li,
                 const long &newli, frag_cell &newCell,
                 space &l)
   { return do_mitosis(CB,i,li,newli,newCell,l); }

   void set_mutation_after_TC(AffinitySpace &shape);
   short mutate(AffinitySpace &shape)
   { return do_mutate(shape); }

   short grow(const long &li, space &l) {
      if (target_volume == 1) {
         return 1;                    // d.h. kein Wachstum!
      } else { return do_grow(CB,li,target_volume,p_grow,p_shrink,l); }
   }
   void get_new_state(const long &i, double &dt, space &l, sigs &s);
   centroblasts set2differentiation();
   void set_adhesion() { get_adhesion(max_adhesion); }
   short ask_differentiate();

   cellCB&operator =(const cellCB &x);
   cellCB&operator =(const cellCC &x);
};

// ====================================================
// ================== cellCC ==============================
// ====================================================

//MSchips
//I need an additional state 'TFRcontact'
//double check whether this needs to be added at the end!
/* MMH2Marta: there is a comparison of centrocytes state < selected in cellthis (only for tracking).
   But TFRcontact should definitely be before selection, so you have to put TFRcontact before selected,
   which I did.
*/
enum centrocytes {
   unselected,contact,FDCselected,TCcontact,TFRcontact,selected,apoptosis
};

class cellCC: public cell {
private:
  static double p_apo,
    p_apo4FDCselected,
    p_mph,
    p_FDCdetach,
    TCselect_prob,
    p_FDCsignalling,
    p_dif,
  // p_dif_DEC,
    p_dif2out,
    p_dif2out_target,
    p_dif2out_DEC,
    p_dif2out_DEC_target,
    p_final_differentiation,
    start_differentiate,
    p_difu,
    p_difu_width,
    persistence,
    v_slow_factor,
    p_switch_v;
  static short v_modi,n_v_states;
  static short apoptotic_motility_mode;
  static double p_apo_randomwalk;
  static short TC_CC_selection;
  static bool force_selection_at_TFHthreshold, negativeTCselection, BCstaysonTCbyTCtime;
  static short TFHsignal_delivery_mode;
  static double tc_time, tc_dec205ova_binding_time, tc_time_width,
    BTtime_K, BTtime_min, BTtime_max, BTtime_n,
    TFHsignal_delivery_min, TFHsignal_delivery_max, 
    TFHsignal_delivery_n, tfr_time;//Msc

public:
  static double TFHsignal_delivery_KpMHC;
  //MSchips
  static short TFR_CC_interaction_mode;
  static double TfhHelp_reduction_factor;
  static double time_firstTFRinteraction;
  static short TFR_bind_onlyself;
  static bool exp_CCL3KO;
  //MSchips
  long ccInd_TFRbound;
  long ccInd_TFHbound;
  double tc_signal, tfr_signal, time_before_tfr; //MS --> this was private!
  double individual_dif_delay; //MS --> this was private!
private:
  static bool TFHadapt2pMHCseen;
  static double TFHsignal_decay, cMyc_decay_factor;
  static double get_TFHsignal_delivery_factor(double, bool);
  static double get_TFHsignal_delivery_factor(double, cellTC&, bool);
  static const int TFHsignal_delivery_max_pMHC = 200;
  static double TFHsignal_delivery_factor[TFHsignal_delivery_max_pMHC];
  static short mode_of_setting_tc_time;
  static double max_tc_signal_portion, rescue_signal_threshold, SST_tc_signal;
  static bool ag_deleted_in_fresh_CC;
  static bool inhibitFDCapoptosis;
  static double prob2kill_noFDCcontactBCs;
  static bool no_rebinding_to_moved_TCs;
  short get_tc_selected(AffinitySpace &shape);
  void progress_selection_state(AffinitySpace &shape);
  void make_apoptosis(const double& time, AffinitySpace &shape, bool fromTfh);
  //MS-- CHANGE!!
  void make_apoptosis2(const double& time, AffinitySpace &shape, bool fromTfh);
  void return2unselected(AffinitySpace &shape);
  int bound_ag_index;
  void add_collected_ag();
  void add_tc_signal(double&);
  
  static bool cMyc_FoxO_break;
  static double cMyc_FoxO_break_K, cMyc_FoxO_break_n;
  double get_cMyc_FoxO_break();
  void add_cMyc_signal(double&);
  void add_mTOR_signal(double&);
  
  double FDCselected_clock;
  double tc_clock,
//            tc_signal,
            TFHsignalAUC, tc_interaction_time;
  double tfr_interaction_time;
  bool SSTactive;
  static double AgThreshold4Selection;
  static short tc_search_duration_mode, tc_search_time_determination_mode;
  static bool FCRB;
  static double tc_search_duration_per_FDCcontact, tc_search_duration_fixed;
  double tc_search_duration;
  double get_tc_search_duration();
  double get_pMHC_presentation();
  double get_tc_interaction_time();
  double set_selected_CC_delay();
//  double individual_dif_delay;
  //MSchips
  double get_tfr_interaction_time();
  double tfr_clock;
  
  static bool pMHC_dependent_division, SIND, use_cMyc4DND;
  static short mode_of_DND_of_Tfh;
  static double DND_P_standard;
  static double DND_P_min;
  static double DND_P_max;
  static double pMHC_dependent_K, TFHsignal_dependent_K, cMyc_dependent_K;
  static double DND_nHill;
  static double pMHC_dependent_pMHC_of_2divisions, TFHsignal_of_P0divisions;
  static double TFHgradient_dependent_K, TFHgradient_of_P0divisions;
  static double pMHC_of_DEC205_ova;
  double get_pMHC_dependent_division();
  double get_cMyc_dependent_DND();
  double get_signal_induced_number_of_divisions();
  static const int max_n_of_ag_portions = 100;
  static short present_specific_ag2TC;
  int get_max_collected_ag(bool returnindex);
  static short outputfiles;
  static double write_trackfate;
  
  /* Signal molecules */
  static short ICOSL_upregulation_mode; 
  static bool ICOSL_dependent_Tfh_signals, ICOSL_memory;
  static double ICOSL_upregulation_time;
  static short FoxO_mode, dT_FoxO_start, dT_FoxO_reg;
  static bool stopFoxOonTFH;
  static double FoxO_ini, FoxO_production, 
    FoxOup_min, FoxOup_max, FoxOup_K, FoxOup_n, KFoxO, nFoxO;
  static short mTOR_mode;
  static double mTORC1_production;
  static double mTOR_Tfh_production;
  static double mTOR_BCR_K;
  static double mTOR_BCR_n;
  double mTOR_BCRcontrol;
  void inhibitFoxO(double& pMHClevel);
  double get_FoxO_rate(double& pMHClevel);
  double ICOSL;
  double get_ICOSL_expression();
  double mTORC1, FoxO, FoxO_upregulation;
  bool hadcontact2Tfh;
  fateTRACK fatetracker;
  static int fatetracker_ndt;
  int fatetracker_n;
  bool came_from_FDCselected;
  double TimeOfLastTbinding;
  
  static void antigen_collection_statistics(int, bool);
  // count the frequency of the duration of the period of pure search for Ag on FDC
  static histo_monitor monitor_pureFDCsearch;
  // count the frequency of different numbers of CC-TC interactions
  static histo_monitor monitor_nTCcontacts_selected;
  static histo_monitor monitor_nTCcontacts_deleted;
  //MSchips
  // count the frequency of different numbers of CC-TC interactions
  static histo_monitor monitor_nTFRcontacts_SelfCCselected;
  static histo_monitor monitor_nTFRcontacts_SelfCCdeleted;
  static histo_monitor monitor_nTFRcontacts_CCselected;
  static histo_monitor monitor_nTFRcontacts_CCdeleted;
  static histo_monitor monitor_time_first_TFRcontact_Selfselected;
  static histo_monitor monitor_time_first_TFRcontact_NonSelfselected;
  static histo_monitor monitor_TFHsig_postTFR_CCselected;
  static histo_monitor monitor_TFHsig_postTFR_SelfCCselected;
  static histo_monitor monitor_BinteractTFRtime_self;
  static histo_monitor monitor_BinteractTFRtime_nonself;
  static histo_monitor monitor_TFRsig_select;
  static histo_monitor monitor_TFRsig_apopt;
  //CD138+ cells
  static histo_monitor monitor_nMut_nonSelfCD138;
  static histo_monitor monitor_nMut_SelfCD138;
  static histo_monitor monitor_nMut_RedeemedCD138;
  static histo_monitor monitor_nTFRcont_SelfCD138;
  static histo_monitor monitor_nTFRcont_NonSelfCD138;

  // count the frequency of acquired BC-Tfh search times
  static histo_monitor monitor_BsearchTtime_selected;
  static histo_monitor monitor_BsearchTtime_deleted;
  // count the frequency of acquired BC-Tfh interaction times
  static histo_monitor monitor_BinteractTtime;
  // count the frequency of times between two Tfh bindings
  static histo_monitor monitor_TimeBetweenTbindings;
  // count the frequency of the amount of integrated Tfh signal at time of selection:
  static histo_monitor monitor_TfhSignalAtSelection;
  static histo_monitor monitor_TfhSignalAtSelection_selected;
  static histo_monitor monitor_TfhSignalAtSelection_deleted;
  // count the frequency of pMHC levels at the time of fate decision
  static histo_monitor monitor_pMHCatFateDecision_selected;
  static histo_monitor monitor_pMHCatFateDecision_deleted;
  // count the frequency of Tfh-signal-speed in BCs
  double calc_TfhSignalSpeed();
  static double get_decay_factor(double, double);
  static histo_monitor monitor_TfhSignalSpeed_selected;
  static histo_monitor monitor_TfhSignalSpeed_deleted;
  // count the frequency of the AUC of Tfh signal at time of selection:
  static histo_monitor monitor_TFHsignalAUC;
  static histo_monitor monitor_TFHsignalAUC_selected;
  static histo_monitor monitor_TFHsignalAUC_deleted;
  // cMyc and mTOR
  static histo_monitor monitor_cMycAtSelection;
  static histo_monitor monitor_cMycAtSelection_selected;
  static histo_monitor monitor_cMycAtSelection_deleted;
  static histo_monitor monitor_mTORatSelection;
  static histo_monitor monitor_mTORatSelection_selected;
  static histo_monitor monitor_mTORatSelection_deleted;
  
public:
  cellCC();
  cellCC(const cellCC &x);
  cellCC(cellCB &x);
  ~cellCC();
  static void set_statics(const Parameter &par, space &l, ofstream &ana);
  static void set_statics(const double &time, const Parameter &par);
  void trackfate_initialize();
  void trackfate(double, bool);
  void trackfate_show();
  void writetrackfate();
  static const int fdc_max_encounters = 50;
  static int fdc_encounters[fdc_max_encounters];
  static unsigned short CC_FDC_selection;
  static double p_CXCR5down;
  static void set_differentiation(const double &time);
  static bool SMOOTH_DIFFERENTIATION;
  static double smooth_differentiation_time;
  static bool collectFDCsignals;
  static bool multipleTFHcontacts;
  static int reset_antigen_after_collection;
  static double ignore_affinity;
  static double dif_delay,dif_delay_DEC;
  static bool ag_loaded_CC_directly2TFH;
  static long test_delay, ICAM_delay;
  static int AgThreshold4Diff;
  static double collectFDCperiod;
  static double IgE_BCRlevel, IgE_prob_CXCR5down;
  static void show_cummulative_ag_collection();
  static long cummulative_ag_collection_all[max_n_of_ag_portions + 1];
  static long cummulative_ag_collection_selected[max_n_of_ag_portions + 1];

  //MS
  /** @brief: #contacts(Tfr:xCC)/#xCC --> x:type self/nonself
   * */
  static void show_freq_of_TfrCont(double&, ofstream&, ofstream&);
  static void show_freq_of_TfrCont_nonself(double&, ofstream&);
  static void show_freq_of_TfrCont_self(double&, ofstream&);
  static const int ntfr = 50;
  static long num_tfr_contacts_nonself[ntfr + 1];
  static long num_tfr_contacts_self[ntfr + 1];
  
  static bool simultaneousTfhFDC;
  static double prob2bindFDCfirst;
  double get_prob2bindFDCfirst();
  
  short selectable,mobile;
  centrocytes state;
  double affinity;
  double selected_clock;
  bool selected4output;
  bool DEC205,DEC205_ova,MHCdeficient;
  long tc_index, last_tc_index, last_tc_id,
  tfr_index, last_tfr_index, last_tfr_id; //MSchips
  long tfr_tmp_index;
  double fdc_clock;
  int nFDCcontacts,nTCcontacts,
  nTFRcontacts, alterednFDCcontacts; //MSchips
  double BCRsignal, cMyc;
  vector<int> collected_ag_portions;
  immunoglobulin_class IgX;
  double BCRexpression;
  short CXCR5failure;
  double DND; // pMHC-, cMyc-, Tfh-signal-dependent number of divisions
  double get_FoxO();
  //MSchips
  bool selfMutation, redeemed, came_from_selected, ccl3KO;
  double tfhHelpReduction;
  bool CD138;
  
  void set_p_move();
  short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
  void set_selectable();
  void set_CXCR5expression();
  void resensitise4CXCL13(sigs &s);
  void set_CXCR4expression();
  void resensitise4CXCL12(sigs &s);
   bool set_apoptotic_motility(sigs &s);
  void go2TCselection(AffinitySpace &shape);
  void delete_antigen();
  void process_antigen();
  long contact2FDC(space &l);
  short bind_antigen(cellFDC &fdc, int frag_index, int &ag_index,
		     AffinitySpace &shape, double &threshold);
  short FDCdetach(AffinitySpace &shape);
  short stop_collecting_FDCsignals(AffinitySpace &shape, const double& time, double &dt);
  short dif2OUTorCB(double&);
  bool final_differentiation();
  void attribute_tc_search_duration();
  bool same_TC_as_before(long new_tc_index, long new_tc_id);
  bool try2findTFH();
  void bind_TC(const double& time, cellTC &tcell, space &l, AffinitySpace &shape);
  //MSchips
  void bind_TFR(const double& time, cellTFR &tcell, space &l, AffinitySpace &shape);
  void progress_signalling(); 
  void decay_TFHsignal();
  void decay_cMyc();
  void addTFHsignalAUC(double& dt);
  short got_tc_signals(const double &time, const double &dt, cellTC &tcell, 
		       space &l, AffinitySpace &shape);
  bool same_TFR_as_before(long new_tfr_index, long new_tfr_id);

  //MSchips
  short got_tfr_interaction(const double &time, const double &dt, cellTFR &tcell,
                       space &l, AffinitySpace &shape);
  void tfr_fate(cellTFR &tcell);
  bool canbindTFR(const double &time);
  short tfh_pointing(const double &time, const double &dt, cellTC &tcell,
                       space &l, AffinitySpace &shape);

  bool do_selection();
  short try_selection(const double &time, const double &dt, AffinitySpace &shape);
  static short mode_of_apoptosis;
  short apoptose(const double &time, AffinitySpace &shape);
  short macrophagocyte(AffinitySpace &shape);
  // monitor TFH signalling intensity
  static histo_monitor monitor_TFHintensity;
  // update and write histograms
  static void update_histograms();
  static void write_histograms();
  
  cellCC&operator =(const cellCB &x);
  cellCC&operator =(const cellCC &x);
};

//MSchips -- new class for TFRcells
// for the time being this is a pure copy of the cellTC class
// one difference would be: active search of BC ???
// ====================================================
// ================== cellTFR ==========================
// ====================================================

enum tfrstates {
   TFRnormal,TFR_CCcontact,TFRdivide
};
class cellTFR: public cell {
  private:
   static double p_difu, p_difu_width, persistence;
   static double proliferation, meancycle, cyclewidth, max_distance2divide;
   static int Ndivisions;

   // <sb>
   //long * CC_nn_TFR;
   //double * CC_affinity;
   std::array<long, 6> CC_nn_TFR;
   std::array<double, 6> CC_affinity;
   // </sb>
   short n_CC_nn_TFR;
   short changed_nn;
   double cell_cycle_clock,cell_cycle_duration;

   static double dTrepolarise, dTsignal;
   double TsinceBCpolarisation;
   bool able2repolarise;
   long CurrentSignalTarget;
   void ini_polarisation_state();

  public:
   static bool do_division;
   // state of TFR
   tfrstates state;
   bool able2signal;
   bool tmp_blocked, to_delete;
   long cc_ind, last_cc_ind;

   cellTFR();
   cellTFR(const cellTFR &x);
   ~cellTFR();
   void ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape);
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   //  static void set_statics(const double& time, const Parameter& par);
   void evolve_polarisation_state(double&);
   void set_changed_nn(short x);

   // Aktionen:
   // =========
   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void make_tc_cc_link(const long &index, const double &nFDCcontacts);
   void make_tc_cc_link(const long &index,
                        const long &CCpos,
                        int ag_index,
                        AffinitySpace &shape,
                        const bool &highag);
   //MSchips
   void TFR_bind_CC(const double& time, cellCC &cccell, space &l, AffinitySpace &shape);
   void make_tfr_cc_link(const long &index);
   void liberateCC(const long &index);
   void set_polarity(space &l);
   static double pMHCsignal_adaptation_weight;
   static double K_signal_intensity_pMHC_default;
   double K_signal_intensity_pMHC;
   double adapt_K_signal_intensity(int);
   short get_n_boundCC();
   void set_CXCR5expression();
   void resensitise4CXCL13(sigs &s);
   void set_CXCR4expression();
   void resensitise4CXCL12(sigs &s);

   // Division:
   double set_cell_cycle_duration();
   void reset_cycle_times();
   void ask_enter_cell_cycle();
   bool progress_cell_cycle(double &dt);
   long ask_mitosis(long * pp, space &l);

   cellTFR&operator =(const cellTFR &x);
};


// ====================================================
// ================== cellTC ==========================
// ====================================================

enum tcells {
   TCnormal,TC_CCcontact,TCdivide
};

// class cellTC: public frag_cell {
/* This can be activated together with a change in cellTC::cellTC(..) in cellthis.cpp.
 * It generates a normal run (tested with bcinflow09rand0.par).
 * But I want to check what exactly happens to all the additional parameters
 * introduced by frag_cell. Probably they are just there and do nothing. Better to check.
 * ### no frag_cell yet for TC!!!
 */
class cellTC: public cell {
  private:
   static double p_difu, p_difu_width, persistence;
   static double proliferation, meancycle, cyclewidth, max_distance2divide;
   static int Ndivisions;

   long * CC_nn;
   short n_CC_nn;
   double * CC_affinity;
   short changed_nn;
   double cell_cycle_clock,cell_cycle_duration;

   static double dTrepolarise, dTsignal;
   double TsinceBCpolarisation;
   bool able2repolarise;
   long CurrentSignalTarget;
   void ini_polarisation_state();

  public:
   static bool do_division;
   // state of TC
   tcells state;
   bool able2signal;

   cellTC();
   cellTC(const cellTC &x);
   ~cellTC();
   void ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape);
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   //  static void set_statics(const double& time, const Parameter& par);
   void evolve_polarisation_state(double&);
   void set_changed_nn(short x);

   // Aktionen:
   // =========
   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void make_tc_cc_link(const long &index, const double &nFDCcontacts);
   void make_tc_cc_link(const long &index,
                        const long &CCpos,
                        int ag_index,
                        AffinitySpace &shape,
                        const bool &highag);
   void liberateCC(const long &index);
   void set_polarity(space &l);
   static double pMHCsignal_adaptation_weight;
   static double K_signal_intensity_pMHC_default;
   double K_signal_intensity_pMHC;
   double adapt_K_signal_intensity(int);
   short get_n_boundCC();
   //MSchips
   void set_CXCR5expression();
   void resensitise4CXCL13(sigs &s);


   // Division:
   double set_cell_cycle_duration();
   void reset_cycle_times();
   void ask_enter_cell_cycle();
   bool progress_cell_cycle(double &dt);
   long ask_mitosis(long * pp, space &l);

   cellTC&operator =(const cellTC &x);
};

// ===============================================================
// ======================= cellFDC ===============================
// ===============================================================

enum FDCstates {
   none,soma,dendrite
};

class cellFDC: public frag_cell {
  private:
   static double p_mksignal,p_mkCXCL13,p_mkSEMA4D;
   static short vesicle;
   static double ic_k_on,ic_k_off;
   static int FDCmaxFrags,n_Antigen_max,n_Antigen_dim_factor;
   static vector<double> ag_fraction;
   static short ag_distribution_mode, ag_detection_mode;
   static bool BCreduceFDCag;
   int get_highest_amount_ag(const int &frag_index, AffinitySpace &AS);
   int get_highest_affinity_ag(const int &frag_index, const long &BCRposAS, AffinitySpace &AS);
   double get_total_antigen_at_site(int frag);
  public:
   cellFDC();
   cellFDC(const cellFDC&);
   ~cellFDC();
   static vector<double> ini_ag_fraction();
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   static void set_statics(const double &time, const Parameter &par);
   void signal_production(const long &i, sigs &l);
   void add_antigen(int n_Antigen);
   double get_total_antigen();
   double get_total_antigen(int ag_index);
   int get_voidsites();
   int get_fragment_index(const long &fragpos);
   short consume_ag(const long &frag_index, const int &ag_index);
   int local_interaction_with_ag(const int &frag_index, const long &BCRpos_ss,
                                 AffinitySpace &shape);
   double get_immune_complex(int fragment);
   double get_total_immune_complex();
   double get_total_immune_complex(const int &ag_index);
   void add2histo_ag_per_site(int* free_ag_site, int* tot_ag_site, int resolution);

   FDCstates state;
   /*multiAg: evtl. do double** antigen_amount, 
    * where the 2nd dim runs over the different antigens,
    * i.e. antigen_amount[f][a] with f running over FDC fragments and a running over antigens.
    * As a vector one could add antigen later on. But not necessary, as it will be known how many
    * antigens
    * will be presented in the time of the GC, so dimension of array can be set accordingly.
    */
   double * * antigen_amount;
   /*multiAg: 
    * Now: ic_amount[f][b] with f running over FDC fragments and b running over affinity bins.
    * For each FDC fragment we have a fixed Ag and a set of Ab affinity bins associated with it.
    * This recollects all Abs in the GC.
    * For multiple Ags, each Ag requires its own Ab affinity bin set, 
    * in order to calculate and save the ic_amount,
    * thus, we need ic_amount[f][b][a] with <a> running over the antigens.
    */
   double * * * ic_amount;
   static unsigned short use_antigen;
   static double antigen_amount_per_FDC,antigen_saturation;
   short antigen_presence(const long &fragpos, int &ag_index);
   void mk_immune_complex(sigs &l);
   void mk_immune_complex(const double &d_t, sigs &l);
   void mk_immune_complex(const double &d_t, AntibodyDyn &ABS, AffinitySpace &AS);
   static double ag_threshold;
   static long ab_sign_errors,ag_sign_errors,ic_sign_errors,ic_calculations;
   void set_antigen_amount(double time, double dt);

   static int DendriteLength;

   cellFDC&operator =(const cellFDC &x);
};

// ===============================================================
// ======================= cellOUT ===============================
// ===============================================================

enum OutState {
    Outfree, OutTFRcontact
};

class cellOUT: public frag_cell {
  private:
   static double 
     p_difu,
     p_difu_width,
     persistence,
     v_slow_factor,
     p_switch_v;
   static short v_modi,n_v_states;
   static short vesicle;
   static bool exit2tz;

  public:
   cellOUT();
   cellOUT(const cellOUT &x);
   ~cellOUT();
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   static double average_affinity,max_affinity;
   static short use_threshold;
   static double initial_ab_affinity;
   static double p_mk_ab;

   immunoglobulin_class IgX;
   bool DEC205;
   short try_eject(long i, space &l, AffinitySpace &shape);
   void set_p_move();
   short move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   void signal_production(const long &i, sigs &l);
   //MSchips
   bool canbindTFR(const double &time);
   void bind_TFR(const double& time, cellTFR &tcell, space &l, AffinitySpace &shape);
   short got_tfr_interaction(const double &time, const double &dt, cellTFR &tcell,
                        space &l, AffinitySpace &shape);
   bool selfMutation, out_to_die;
   int nTFRcontacts;
   static short TFR_CC_interaction_mode;
   long tfr_index;
   double tfr_clock, tfr_interaction_time;
   OutState state;

   cellOUT&operator =(const cellOUT &x);
   cellOUT&operator =(const cellCC &x);
   cellOUT&operator =(const cellCB &x);
};

// ===============================================================
// ======================= cellbeta ==============================
// ===============================================================

// enum betacells{};

class cellbeta: public frag_cell {
  private:
   // constants:
   static constexpr double Avogadro = 6.02205e+23; // mol^-1
   static constexpr double MFaraday = 9.6485309e-02; // Faraday constant in C/(micromol)
   // 9.6e+04 C/mol = 9.6e+04*1e-06 C / 1e-06mol = 9.6e-02 C/micromol
   static constexpr double Faraday = 9.6485309e+04; // Faraday constant in C/(mol)
   static constexpr double Rydberg = 8.315; // in J/(K*mol)
   static constexpr double pi = 3.141592654;

   enum beta_currents {
      NaK,K_ATP,K_V,K_Ca,sK_Ca,Na_V,fNa_V,NCX,PMCA,Ca_L,Ca_T,
      SERCA,cIP3,N_currents
   };
   static double I[N_currents];
   static void get_currents(double * y,
                            double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER,
                            double &C_IP3_inh);
   static void get_current_factors(double t);

   static double
      max_pro_distance,
      p_tension,
      p_grow,
      p_shrink,
      p_difu,
      diffusion_tolerance_min,
      diffusion_tolerance_steepness,
   // movement and shape
      tension,
      elongation,
      K_elongation,
      smoothmove,
      persistence,
      v_slow_factor,
      p_switch_v,
      max_adhesion;
   static short v_modi,n_v_states;

   // needed for electrophysiology calculation:
   static double xi,xi_ER,xi_ERC,xi_S_ERC;
   static double J_K,J_Na,J_Ca,J_ions;
   static double V_ER_0;

   static double get_2sigmoidal(double &t, double t_a, double t_b,
                                double rest, double factor,
                                double kappa_a, double kappa_b);
   double get_glucose(const long &i, sigs &s);
   static double get_IP3(double &t);
   static double get_tau_IP3(double &ip3);
   static double get_sigmoidal(double &half, double &x, double &kappa);
   static double get_inactivation(double &half, double &x, double &kappa);
   static double get_ca_buffer(double &b_0, double &c, double &dissociation);
   static double get_Ca_free_fraction(double &c);
   static double get_Ca_ER_free_fraction(double &c);
   static double get_dynamic_half(double &x, double a, double b);
   static double get_tau_V(double &x, double a, double b, double c, double Vx);
   static double get_tau_K_V_sherman88(double &x, double &VbarK);
   static double get_K_ext(double &t, double &k_0);
   static void get_all_Nernst(double &t, double &K, double &Na, double &Ca, double &Ca_ER,
                              double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER);
   static void insert(double &ion, double &voltage, double I_load, double valence_sign);

   void show(double t, double * y_n, ofstream &file, ofstream &file2);
   void show_rev(double t, double * v, ofstream &file);
   void show_c(double t, ofstream &file);
   void show_cr(double t, ofstream &file);
   void show_gap(double t, double * y, ofstream &file);

   // Declare cell specific variables
   // betacell specific variables for each cell
   double * y_n1; // try these to be static! ###
   double * y_n_old;
   static double betadt;
   static long betandt;
   void set_initial_values();

   // solver related functions:
   // static const ode_method method=RungeKutta_4th;     // available methods are listed in ode.h
   ode_method method;
   static ode solver;
   void (*prhs)(double, double*, double*); // Pointer on rhs(...)
   void get_gap_junction(dynarray<cellbeta> &bl, space &l);

   // each cell gets its own files
   ofstream result,result2,currents,currents_rho,reversal,gapfile;
   void open_files();
   static long n_write;
   long step_count;

  public:
   // ++++++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++
   // if the number of betacells is larger than 200 this flag shall be true
   static const bool LOCAL_FILES = false;
   // choose true here, if gap-junctions shall also be treated with Runge-Kutta
   static const bool FULL_RUNGE = true;
   // set whether the sequence of cells in cellman::calc_CB shall be randomised
   static const bool RANDOMISE_SEQUENCE = false;
   // ++++++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++
   static constexpr double z_K = 1.0; // valence of potassium ions
   static constexpr double z_Na = 1.0; // valence of sodium ions
   static constexpr double z_Ca = 2.0; // valence of calcium ions
   static betaWerte p;
   static double p_pro;
   static double glucose_rest;
   static int target_volume;
   enum beta_quantities {
      V,K,Na,Ca,                                                  // 0-3
      g_K_ATP,g_K_V,g_Na_V,g_Ca_L,g_Ca_T,                 // 4-8
      h_K_V,h_Na_V,h_Ca_L,h_Ca_T,                         // 9-12
      g_K_Ca,C_K_Ca,                                      // 13-14
      Ca_ER,V_ER,IP3,g_IP3,h_IP3,                         // 15-19
      g_sK_Ca,glu,                                        // 20-21
      gap_K,gap_K_0,gap_K_1,gap_K_2,gap_K_3,gap_K_4,gap_K_5,                 // 22-28
      gap_Na,gap_Na_0,gap_Na_1,gap_Na_2,gap_Na_3,gap_Na_4,gap_Na_5,          // 29-35
      gap_Ca,gap_Ca_0,gap_Ca_1,gap_Ca_2,gap_Ca_3,gap_Ca_4,gap_Ca_5,          // 36-42
      g_fNa_V,                                                               // 43
      N_equations
   };                                                                        // 44

   void randomise_protein_expression();

   // The actual quantities for each cell
   double * y_n;
   // double rho_gap,rho_K_ATP;
   double * rho;

   cellbeta();
   cellbeta(const cellbeta &x);
   ~cellbeta();
   void destruct() { delete[] fragments; } // Destruktor des Fragment-Feldes
   void ini(const long &i, const long &li, const double &t, space &l);
   static void set_statics(const Parameter &par, space &l, ofstream &ana);
   void synchronise();

   // Aktionen:
   // =========

   static double get_Nernst(double &x, double &x_ext, double valence);
   static void rhs(double t, double * y, double * derivative);
   void electrophysiology(double thr, double dthr, sigs &s,
                          dynarray<cellbeta> &bl, space &l);

   void set_p_move();
   double move(const long &li, space &l, sigs &s, TRACK &td, double &time);
   long ask_mitosis(long * pp, space &l);
   short mitosis(const long &i, const long &li,
                 const long &newli, frag_cell &newCell,
                 space &l)
   { return do_mitosis(CB,i,li,newli,newCell,l); }
   short grow(const long &li, space &l) {
      if (target_volume == 1) {
         return 1;                    // d.h. kein Wachstum!
      } else { return do_grow(CB,li,target_volume,p_grow,p_shrink,l); }
   }
   void get_new_state(const long &i, sigs &s);
   void set_adhesion() { get_adhesion(max_adhesion); }
   void show_all(double t);
   void show_all(double t,
                 ofstream &f_a,ofstream &f_b,ofstream &f_i,
                 ofstream &f_r,ofstream &f_n,ofstream &f_g);
   static suffix beta_index;

   cellbeta&operator =(const cellbeta &x);
};

// ===========================================
// (Vergleichs)-Operatoren
char operator ==(const cellCB &a, const cellCB &b);
char operator !=(const cellCB &a, const cellCB &b);

#endif
