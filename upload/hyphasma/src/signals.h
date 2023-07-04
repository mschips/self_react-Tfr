#ifndef i_signals
#define i_signals
// #include "random.h"
#include "setparam.h"
#include "grid.h"
#include "space.h"

enum signal_molecule {
   sig_differ2CC,    // hypothetic signal inducing CB differentiation to CC
   sig_proliferate,  // hypothetic signal inducing mitosis
   CXCL12,           // Ligand for CXCR4 expressed on B cells
   CXCL13,           // Ligand for CXCR5 expressed on B cells
   antibody,         // soluble antibody specific for antigen
   antigen,          // soluble antigen
   SEMA4D,           // soluble SEMA4D semaphorin
   glucose,
   oxygen,
   signals           // number of signals
};
/*
 * Doku:
 * Um ein neues Molekuel einzufuehren, sind die folgenden Schritte durchzufuehren:
 * 1) Molekuel XXX (fuer new signal) oben in die Liste eintragen.
 * 2) Molekuel-Produktion in parameter-file einfuehren. Das heisst, eine neue
 *   Variable mkXXX in setparam.x einzufuehren, die die Produktionsrate des Molekuels
 *   durch eine Zelle (z.B. FDC) bestimmt.
 * 3) Boundary-conditions am Rand des Reaktionsvolumens angeben. Dazu ist ebenfalls
 *   eine neue Variable boundXXX in setparam.x einzufuehren.
 * 3a)Ausserdem ist eine Variable fuer die Diffusionskonstante notwendig (neu).
 * 4) Die obigen neuen Variablen muessen im Konstruktor der Klasse lattice (in
 *   lattice.C) eingelesen werden. Dazu ist:
 *     signal_use[XXX]=1
 *     boundary.signal[XXX]=par.Value.bound_XXX; (Wert aus parameter-file)
 *     D[XXX]=par.Value.D_XXX;
 *   zu setzen.
 * 5) Die Produktionsrate ist als Wahrscheinlichkeit in die Klasse der produzierenden
 *   Zelle einzutragen (Bsp. FDC). Dazu ist:
 *     static double p_mkXXX; in cellthis.h
 *     double cellFDC::p_mkXXX=0.; in cellthis.C
 *     p_mkXXX=par.Value.mkXXX*par.Value.deltat; in beide cellFDC::set_statics(...)
 *       (die eine ist zur Initialisierung, die andere zur zeitabhaengigen Aktualisierung)
 *   einzutragen.
 * 6) Die Zeile
 *     signal_secretion(i,XXX,vesicle,p_mkXXX,l);
 *   ist in die Routine cellFDC::signal_production(...) in cellthis.C einzutragen.
 *   Note: Fuer andere produzierende Zellen muss eventuell eine neue signal_production-Routine
 *         eingefuehrt werden (in der entsprechenden Zellklasse). Diese muss dann in
 *     cellthis.h deklariert und von cellman.C aufgerufen werden (wie in calc_FDC(...)).
 * 7) Den Dateinamen fuer den Output in cellman.C in die Routine get_signal_name(..) eintragen.
 * 8) Fertig!
 *
 */

class spoint {
  public:
   // Konstruktoren
   spoint();
   ~spoint();

   // save the new-variables as standard
   void actualize();
   void actualize(const signal_molecule &a);

   double signal[signals];
   double signal_new[signals];
   double signal_tmp[signals];
};

// (Vergleichs)-Operatoren
// char operator==(const point& a, const point& b);
// char operator!=(const point& a, const point& b);
// (Vergleichs)-Operatoren
// char operator==(const spoint& a, const spoint& b);
// char operator!=(const spoint& a, const spoint& b);

class sigs: public grid {
  public:
   static constexpr double SIGNAL_MAX = 1.0e+25;

   static bool TEST_MODE;
   static double TEST_AMPLITUDE;
   static double TEST_SOURCE;
   static suffix TEST_SUF;

   sigs();
   sigs(space &xyz, Parameter &par, ofstream &ana);
   ~sigs();

   spoint * sigsknot;
   double mk_concentration;

   short diffusion_mode;
   short signal_use[signals];
   double signal_total[signals];  // Zaehle alle jemals sekretierten Signale
   double signal_diffout[signals]; // Zaehle alle jemals durch Diffusion
   // in dem Volumen verlorene Signale
   double signal_used[signals];   // Alle durch Rezeptoren gebundenen Signale
   double signal_produced[signals]; // Zaehle die zwischen zwei Ausgaben produzierten Quanta

   static double get_const_signal_value(double &t, double &ampl);
   static double get_const_signal_value(double &t, double &ampl, long * wo);
   void mk_const_signal(const signal_molecule sig_type, double t, double ampl);
   void mk_const_signal(const signal_molecule sig_type, double &value);
   void mk_grad1d_signal(const signal_molecule sig_type,
                         double min, double max, short direction);
   void load_signal_file(const signal_molecule sig_type);
   void signal_put(const long &i, const signal_molecule &sig_type, const double &zahl);
   void signal_set(const long &i, const signal_molecule &sig_type, const double &zahl);
   void boundary_set(const long &i, const signal_molecule &sig_type, const double &zahl);
   short signal_get(const long &i, const signal_molecule &sig_type, const double &howmuch);
   void signal_diffuse(space &xyz);
   bool overcritical_signal(const long &index, const signal_molecule &s, const double &crit);
   bool undercritical_signal(const long &index, const signal_molecule &s, const double &crit);
   void get_gradient(const signal_molecule &s, const long &index, double * gradient);
   double get_signal_total(const signal_molecule &sig_type);

   void get_signal_name(const unsigned short &sx, char * name);
   void write_siglog(double &t);
   void write_TEST(double &t);
   void write_files(suffix tnr, bool forceit);
   void write_all_signals(double &chemo_max, double &chemo_steep,double &chemo_half);

  private:
   // ++++++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++
   static const int SIGNAL_TIME_RESOLUTION = 1; // recommend 2
   static const int N_PRE_WORDS = 20;
   static double vonNEUMANNinwards;
   // ++++++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++

   // GLUCOSE HANDLING:
   // Initialise glucose distribution as a gradient field (allows for later changes)
   // static bool FIX_GLU_GRADIENT; // not needed anymore as shifted to parameter file
   // Use a constant but dynamic in time glucose field (determine it in mk_const_signal)
   static double CONST_DYN_GLU_FIELD;
   // Note that if glucose is not treated as reaction-diffusion-soluble
   // at least one of the latter two options has to be set true!

   // double dt; // global time step (info for diffusion)
   spoint boundary;

   double D[signals];
   bool fix_signal[signals];
   bool dirichlet[signals];
   double p_diffuse_signal[signals];
   unsigned int diff_steps[signals]; // number of internal time steps for diffusion

   static double get_2sigmoidal(double &t, double &t_a, double &t_b,
                                double &rest, double &factor,
                                double &kappa_a, double &kappa_b);

   short ignore_objects;
   void set_initial_condition(const signal_molecule &bb);
   double reference(long * v, double &t, const signal_molecule &s);
   double signal_anregung(const long &n, const signal_molecule &sig);
   void signal_diffuse_QUANTA(const long &i, const signal_molecule &sig_type);
   void signal_diffuse_EULER(const signal_molecule &sig_type);
   void signal_diffuse_EULER(const long &i, const signal_molecule &sig_type, space &xyz);
   void signal_diffuse_CN(const signal_molecule &sig_type);
   void signal_diffuse_ADI(const signal_molecule &sig_type);
   double get_rhs_x_CN(const long &r,const double &alpha,
                       const signal_molecule &sig_type);
   double get_rhs_y_CN(const long &r,const double &alpha,
                       const signal_molecule &sig_type);
   double get_rhs_x_ADI(const long &r, const double &alpha,
                        const signal_molecule &sig_type);
   double get_rhs_y_ADI(const long &r, const double &alpha,
                        const signal_molecule &sig_type);
   double get_rhs_z_ADI(const long &r, const double &alpha,
                        const signal_molecule &sig_type);
   void actualize();
   ofstream siglog;
   // help-field for matrix inversion
   double * gam;
   long maxprodim;
};

#endif
