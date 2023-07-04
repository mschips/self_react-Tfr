#ifndef i_cellman
#define i_cellman
#include "cellthis.h"
#include "kinetics.h"
#include "track.h"
#include "brainbow.h"
#include "odelist.h"
#include "tamoxifen.h"

const short int xlogs = N_cells;

class cellman {

public:
   cellman();
   cellman(Parameter &par, space &l, sigs &s, AffinitySpace &shape, ofstream &ana);
   ~cellman();
   void close_files();

   // move this to private sector when deform-file is removed
   dynarray<cellCB> CB_list;

   representation show_mode;
   short int outputfiles;
   void time_step(short int ss_save, short int day_save,
                  space &l, sigs &s,
                  AntibodyDyn &Ab,
                  AffinitySpace &shape,
                  ofstream &ana);
   void write_cell_specific();
   void set_pars(Parameter &par, space &l);
   // Space picture of all cells:
   void set_pov_sphere(int,int,int,double,double,double,double,ofstream&);
   void xfiles(suffix tnr, space &l);
   // Spacial Distribution of CBs and CCs relative to the FDCs
   void zone_files(suffix tnr, space &l);
   // Ki67 injection
   void inject_Ki67();
   void label_BrdU(bool);
   // void preload_CB_with_ag(const double&);
   // perform all checks
   short check_all(space&, AffinitySpace&);
   // write betacell files:
   void show_BETA();
   // write files at the very end
   void write_final(const double &deltax);
   void show_ag_distribution_on_FDC();
   ofstream movie;
   void show_BrdU();
   void show_cell_in_BrdU();
   void show_ag_loaded();
   void show_ag_BrdU(bool);
   void write_mutations(double time);

   double time;
   double t_dark_end,CC2CBratio,t_1st_under100;
   double first_record, time_window;

   // Average of centroblasts on last 2 days
   double CB_end,CB_average,CB_variance;

   GCkinetics kinetics;

   TRACK trackdata;

   // antibody handling
   static bool use_antibody_bins;

   bool track_mutations; // true if brainbow class shall be used
   brainbow trackmutations; // uses default constructor

   void mk_single_beta_files(space &l);
   
   void write_log_bc_aff(AffinitySpace &AS);

  private:
   // *************************************
   // Cell tracking ***********************
   // *************************************
   // Zahl der Punkte auf der v-Achse fuer die v-Verteilung der diffusiven Bewegung
   static int v_resolution;
   // intervall of velocities in microns/min (used for time intervall averaged velocities)
   static double delta_v;

   // Number of cells of different celltypes to be tracked
   // If tALL<=0 the individual values are used
   // If tALL>0 only the total number of tracked cells is fixed irrespective of their type
   static int tALL;
   static int tCB;
   static int tCC;
   static int tOUT;
   static int tTC;
   static int tBETA;
   // *************************************

   // ++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++++
   // Chose the resolution of the data output for the affinity
   static const unsigned short affinity_resolution = 10;
   // +++++++++++++++end OPTION +++++++++++++++++++++++++++++++++
   // returns the bin attributed to value assuming max_value of 1 and affinity_resolution
   short get_bin(double value);

   void set_tracking(Parameter &p);
   void get_average_and_sd(int * values, int &n, double &average, double &sd);
   short checkit;
   short int norm;
   long GCpointnumber,FDCpointnumber;
   long BC_integral;
   double dt;
   double last2days;
   short int tcounter;

   void make_fluorescent(space &l);
   void stop_fluorescent();
   tamoxifen tmx;
   void recombine_cells(double t);

   // movement handling
   long * velocity;
   long * delta_velocity, * delta_velocity_euklid;
   bool allow_exchange;

   bool photoactivation;
   long * photoactivation_rmin, * photoactivation_rmax;
   double photoactivation_t0;
   void make_photoactivation(space &l);
   bool def_DEC205,inject_antiDEC205OVA,DEC205_induce_differentiation,retain_DEC205_ag;
   double def_DEC205_t0,inject_antiDEC205OVA_t0,p_DEC205,DEC205_ova_activity,
          TC_factor_dec205ova,DEC205_forces_output,p_CB2OUT;
   void attribute_DEC205();
   void do_inject_antiDEC205OVA(space &l, AffinitySpace &shape, ofstream &ana);
   static double dec205_max_antigen;

   // BrdU staining:
   bool do_inject_BrdU;
   double tfirst_inject_BrdU, tnow_inject_BrdU, tnext_inject_BrdU, 
     deltat_inject_BrdU, BrdU_detection_threshold;
   int n_inject_BrdU;
   void show_BCstate_in_BrdUpositive();
   void show_BrdU_at_time();

   // Some checks of cell-fragments in the cell-lists
   void fragment_consistency(space &l);
   //   void check_TC_neighbours(space& l);
   void check_listi(space &l);

   // Definiere Listen fuer die verschiedenen Zelltypen
   dynarray<cellCC> CC_list;
   dynarray<cellTC> TC_list;
   dynarray<cellFDC> FDC_list;
   dynarray<cellOUT> OUT_list;
   dynarray<long> STROMA_list;
   dynarray<cellbeta> BETA_list;
   dynarray<cellTFR> TFR_list; //MSchips

   // type-independent cell-handling
   long try2exchange_cells(long &i, long &li, states my_type, double * pol,
                           dynarray<long> &redolist, space &l);
   long retry_movement(dynarray<long> &redolist, space &l); // calls try2exchange

   // CB handling on the CB-lists:
   // put/del a centroblast at point i and position pos_ss in shapespace
   long put_CB(long int i, cellCB &newCB, space &l, AffinitySpace &shape);
   // long int find_CB(long i);
   short int del_CB(long i, long li, space &l, AffinitySpace &shape);
   double p_macrophage;
   short macrophagocyte_CB(long i, long li, space &l, AffinitySpace &shape);
   short int differ2CC_CB(long i, long li, space &l, AffinitySpace &shape);
   short int CB_differ2OUT(long i, long li, space &l, AffinitySpace &shape);
   // cell division
   short proliferate_CB(long i, long li, space &l, AffinitySpace &shape);
   static double mutation_start_time;
   bool retain_ag;
   double divide_ag_asymmetric;

   // founder cell division
   int founder_ndiv,n_founder;
   double p_founder_ndiv_plus,p_BCinflux,t_stopBCinflux;
   void determine_founder_times_of_division(Werte &p);
   int get_founder_times_of_division();
   static double smooth_stopBCinflux;
   static bool do_smooth_stopBCinflux;
   short BCinflux(space &l, AffinitySpace &shape, sigs &s);

   // do everything:
   void calc_CB(long * m,long mlang, short int ss_save,
                space &l, sigs &s, AffinitySpace &shape, dynarray<long> &redo);

   // Probabilities for proliferation etc...
   // number of missed proliferation trials
   long pp[4];

   // CC handling on the CC-lists
   // put/del a centrocyte at point i, position pos_ss in ss, and in state ccs
   short int put_CC(long i, cellCC &newCC, space &l, AffinitySpace &shape);
   // long int find_CC(long i);
   short int del_CC(long i, long li, space &l);
   short int get_into_FDC_contact(long li, AffinitySpace &shape, AntibodyDyn &Ab,
				  space &l, TRACK &td, double &time);
   long TestedFDCsiteVoidOfAg; // event counter
   short CC_differ2CB(long i, long li, space &l, AffinitySpace &shape);
   short CC_differ2out(long i, long li, space &l, AffinitySpace &shape);
   short int macrophagocyte_CC(long i, long li, space &l, AffinitySpace &shape);
   long CC_total; // Zaehle alle jemals erzeugten CCs
   void calc_CC(long * m,
                long mlang,
                space &l,
                sigs &s,
                AntibodyDyn &Ab,
                AffinitySpace &shape,
                dynarray<long> &redo);

   // TC handling on the TC-lists
   // put/del a TC at point i, position pos_ss in ss, and in state ccs
   long put_TC(long i, cellTC &newTC, space &l, AffinitySpace &shape);
   void addTC(long nTC, space &l, double restrict, AffinitySpace &shape, ofstream &anafile);
   short int del_TC(long i, long li, space &l);
   void divide_TC(long i, long li, long j, space &l, AffinitySpace &shape);
   void calc_TC(long * m, long mlang, space &l, sigs &s,
                AffinitySpace &shape, dynarray<long> &redo);

   //MSchips
   // TFR handling on the TFR-lists
   // put/del a TFR at point i, position pos_ss in ss, and in state ccs
   long put_TFR(long i, cellTFR &newTFR, space &l, AffinitySpace &shape);
   void addTFR(short tfr_mode, long nTFR, space &l, double restrict, AffinitySpace &shape, ofstream &anafile);
   short int del_TFR(long i, long li, space &l);
   void divide_TFR(long i, long li, long j, space &l, AffinitySpace &shape);
   void calc_TFR(long * m, long mlang, space &l, sigs &s,
                AffinitySpace &shape, dynarray<long> &redo);
   void free_tfr(long cc_ind, long tfr_ind);
   //This is meant for CB
   double p_selfMut, p_redeem;
   short TFR_mode, SelfReact_mode;
   double tfr_depletion_time;
   bool exp_stop_apo, exp_CCL3KO;
   double CCL3KO;
   int asc_got_contact;
   int newly_generated_sCBs, newly_generated_redeemed;


   // FDC handling on the lists
   // put a FDC with soma at i and with 6 arms of length armlength
   short FDCtransparent;
   short int put_FDC(long i, space &l, AffinitySpace &shape, short int transparence);
   void calc_FDC(space &l, sigs &s, AntibodyDyn &Ab, AffinitySpace &AS);
   long getFDCnetwork_volume();
   long getNofVoidFDCsites();
   void showFDCfailure(double& t);

   double p_mkCXCL12;
   void calc_stroma(sigs &s);

   // OUT handling
   // Diffusion/Ejection of out cells
   // long int find_OUT(long i);
   void add_extrafollicular_PC(AffinitySpace& shape, long pos_ss);
   short int del_OUT(long i, long li, space &l);
   void calc_OUT(space &l, sigs &s, AffinitySpace &shape, dynarray<long> &redo);

   // BETA-CELLS
   long put_BETA(long int i, cellbeta &newBETA, space &l);
   short del_BETA(long i, long li, space &l);
   double p_macrophage_BETA;
   short macrophagocyte_BETA(long i, long li, space &l);
   short proliferate_BETA(long i, long li, space &l);
   void calc_BETA(long * m, long mlang, space &l, sigs &s, dynarray<long> &redo);

   ode_method method;
   odelist solver;
   static long get_global_index(long &li);
   static void get_gap_junction_static(double * y, long startat,
                                       const long &celli, dynarray<cellbeta> &list,
                                       long xypos, space &l);
   static void get_gap_junction_dynamic(double * y, long startat,
                                        const long &celli, dynarray<cellbeta> &list,
                                        long xypos, space &l);
   //   static void get_gap_junction(double* y, long startat, space& l);
   static void get_gap_rhs(double * y, double * derivative,
                           const long &startat, const long &celli,
                           dynarray<cellbeta> &list, space &l);
   static void rhs(double t, double * y, double * derivative, dynarray<cellbeta> &list, space &l);
   void (*prhs)(double, double*, double*, dynarray<cellbeta>&, space&); // Pointer on rhs(...)
   void call_solver(dynarray<cellbeta> &list, space &l);

   // Ki67
   short show_Ki67;

   // file handling
   ofstream dark_zone;
   ofstream light_zone;
   void zone_add(long i, long &n_CB, long &n_CBnr, long * n_CC,
                 long &n_all, space &l);
   void zone_put(double t, long n_CB, long n_CBnr, long * n_CC,
                 long n_all, ofstream &datfile);
   void mk_cell_sum(space &l, sigs &s, AntibodyDyn &Ab, AffinitySpace &shape);
   void save_foxo(ofstream& foxofile, short& foxores, int totalCC, int* nfoxo);
   void show_foxo(double time, long* nCC);
   void count_dec205_ova_positive();
   void show_cycle_phases(long zone_separator);
   void show_K_signal_intensity_pMHC();
   bool ignore_apoptotic_CC;
   void get_average_mutation_prob(double&,double&);
   // void get_signal_name(const unsigned short& sx, char* name);

   /* @brief writes affinities of BC and accumulated output cells to files */
   
   void write_log_bc_aff_this(short celltyp, long &ASpos, AffinitySpace &AS);
   void write_log_out_aff(double bestaff, long ASpos, AffinitySpace &shape);
   void write_Ig_class_output();

   //MSchips
   void BcAnalysis(double time, space &l, AffinitySpace &shape);
   void BcAnalysis2(double time, space &l, AffinitySpace &shape);
   void writeAffinity(double time, space &l, AffinitySpace &shape);
   void writeAffinity2(double time, space &l, AffinitySpace &shape);
   void contactscheck(double time, space &l);
   void writeZPositions(double time, space &l);
   void write_info(double time, double tfhsig, double pmhc, double affinity, short status, bool self);
   void get_average_mutation_prob_SELF(double&, double&);
  
   ofstream xsums[xlogs];
   ofstream xvolume;
   ofstream xsumBCig;
   ofstream xsumBCOig;
   ofstream xapoig,xapoigdz,xapoiglz;
   ofstream xaffig;
   ofstream log_out_aff, log_bc_aff;
   ofstream ag_presentation;
   ofstream cum_ag;
   ofstream prolog;
   // ofstream cbquality;
   ofstream cbhighaff;
   ofstream movements;
   ofstream move_vdt;
   ofstream axis;
   ofstream polarities;
   ofstream corrs;
   ofstream xtc;
   ofstream xdec205;
   ofstream aff_cb_lz;
   ofstream aff_cc_lz;
   ofstream aff_out_lz;
   ofstream aff_cb_dz;
   ofstream aff_cc_dz;
   ofstream aff_out_dz;
   //MSchips
   ofstream xtfr;
   ofstream freq_self; //nTfr:cont
   ofstream freq_nonself;

   // Totale Zahl der recycling events
   long int n_recycling_events,n_recycling_events_last;
   // Totale Zahl der Mutationen von output-Zellen
   long int n_outs,n_outs_with_recycling,n_muts,n_recmuts;
   ofstream mutation_out;
   static const int max_mutation_bin = 30;
   static int mutation_frequency[max_mutation_bin]; // write this to file in write_final(..)
   ofstream cum_vol;
   // count the number of generated output cells of different Ig-classes
   int integrate_out[nIg_classes+1];
   int integrate_out_old[nIg_classes+1];

   ofstream ag_collected;
   void show_ag_collected();
   ofstream ndivtime;
   ofstream mutation_time;
   //MSchips
   ofstream selfBC_ndivtime;
   ofstream Selfmutation_time;

   // about betacells:
   ofstream beta_a;
   ofstream beta_b;
   ofstream beta_i;
   ofstream beta_r;
   ofstream beta_n;
   ofstream beta_g;
   void mk_beta1file(const long &ind, const char filename[30], const int columns);

   long get_GCvolume();
};

#endif
