/**
 * @file cell.h
 * @brief Cell class providing cell representation on grid, cell
 *        movement, nutrient consumption, chemotactic properties etc.
 *
 * @author Michael Meyer-Hermann
 */

#ifndef i_cell
#define i_cell
#include "ss.h"
#include "sequencespace.h"
#include "space.h"
#include "signals.h"
#include "track.h"
enum cell_status {
   proliferate,quiescent,necrotic
};

/**
 * @brief General cell class for representation on a grid, using
 *        the subelement method.
 */
class cell {
  friend class cellTest;
public:
  /**
   * @brief Default constructor
   * @return cell instance
   *
   * Sets all cell properties to default values
   */
  cell();
  
  /**
   * @brief Copy constructor
   * @param x cell instance to copy values from
   * @return cell instance
   *
   * Create a copy of cell x
   */
  cell(const cell &x);
  
  /**
   * @brief Destructor
   */
  ~cell();
  
  /**
   * @brief Set static cell attributes based on parameter values in `par`
   * @param par Parameter instance
   * @param ana ofstream for writing out the calculated attributes.
   */
  static void set_statics(const Parameter &par, ofstream &ana);
  
  /// Lattice dimension
  static short lattice_dim;
  
  /**
   * @brief Returns the lattice index of a neighbour of type celltype.
   * @param celltype states
   * @param l space
   * @return int index of a contact cell or -1
   *
   * Checks if there is a cell of type celltype in contact to self-lattice
   * point index. Returns the index of the cell or -1 if none is found. If more
   * than one contact exists, it is chosen randomly.
   */
  long get_contact(const states &celltype, space &l);
  
  // Index on lattice
  long index, born_index;
  // Index for the brainbow analysis
  long brainbow_index;
  // global identifier (1..\infty) for all cell types
  long id;
private:
  // global identifier last attributed (starts with zero before first attribution)
  static long id_last;
public:
  // Internal cell clock in time steps
  double age;
  void aging(const double &dt);
  long clock;
  void set_clock();
  
  // Time of cell production
  double born_time;
  
  // Status of the cell
  cell_status status;

  /// Marker for cell-tracking
  bool trackit;
  
  action_types writethis2track;
  
  /// Absolute number of tracked cell
  int trackno;
  
  /// Nutrient consumption
  void use_nutrient(sigs &l, long &sigindex);
  
  /// saves on short term whether initiated movements are suppressed for lack
  /// of space
  bool contact_inhibited;
  
  /// Total number of cell fragments = cell volume (1 within class cell)
  int volume;
  
  // Probability of cell displacement defined in derived classes of cell
  double p_move;
  
  /// Normed vector pointing into direction of a cell polarity.
  /// This may be used for persistence of cell movement for example
  double polarity[3];
  
  /// Changes of polarity (1 for changed, 0 for unchanged)
   short changed_polarity;
  
  // Signalling state
  bool responsive2signal[signals];
  
  /// Position in shape space
  long int pos_ss;
  
  /// @brief Mutation probability
  /// (has to be on class `cell` level even though only used in derived class
  /// `cellCB` because the value has to be saved and set when the cell becomes
  /// a member of class `cellCC`)
  /// mimehe: I would keep it here, as I will also eventually use it for TC in the future!
  double p_mutation;
  
  /// Proliferative markers: Ki67 marking: =1 else =0 and BrdU
  short int Ki67;
  double BrdU;
  
  /* Counter of actions for each cell 
   * ### this is actually a BC property rather than a cell property
   * mimehe: yes, but again (as for double p_mutation) these counters have to be 
   * kept through the transformation to the cellCC class and back to cellCB class. 
   * that's why they are here. 
   * It might be useful to introduce an intermediate class cellBC where such things 
   * are declared and cellCC/CB would both be derived from cellBC.
   */
  int n_recycling, n_mutation, n_recandmute,n_fdc_encounters;
  
  /// Added movements
  double add_moves;
  
  /// Position for evaluation of v
  double last_position[3];
  
  /// Number of time steps without movement
  int n_immobile;
  
  /// Time resolution with which the velocity of the cell is analyzed
  static double deltat_v;
  
  /// Chemotaxis parameters
  static double chemo_max,chemo_steep,chemo_half;
  static double north_weight;
  static double TFR_north_weight;
  static bool exp_stop_apo;
  
  /// = operator
  cell&operator =(const cell &x);
  
protected:
  /// Expression level of adhesion molecules in %
  double adhesive;
  
  /// Time the cell stays in place upon adhesive contact
  static long adhesion_time;
  
  /// @brief Calculate expression of adhesion molecules
  /// @param max double maximal expression level
  ///
  /// Currently simoky sets expression level to function argument
  /// `max`. A more sophisticated calculation of the expression of adhesion
  /// molecules can be introduced here.
  void get_adhesion(const double &max);
  
  static int ab_resolution;
  
  // velocity states
  double v_state;
  
  // cell pressure in units of external pressure
  double pressure;
  
  short do_diffuse(states celltype, const long &li, space &l);
  double get_v_factor(const short &N_v, const double &p_slow_factor);
  
  void set_south_polarity(space &l);
  //MSc
  /*set_north_polarity fct modified to take the north weight as input*/
  void set_ext_polarity(space &l, double N_weight);
  void set_north_polarity(space &l, double N_weight);
  void set_southTfr_polarity(space &l, double N_weight, const signal_molecule &chemokine, double *pol, sigs &s);
  void set_Tfrchemotaxis_polarity(const signal_molecule &chemokine, double * pol, sigs &s);
  void set_polarity_velocity(const double &persistence,
			     const short &v_modus,
			     const double &p_switch_v,
			     const short &N_v,
			     const double &p_slow_factor,
			     space &l,
			     sigs &s);
  // This one is to be used if the cell don't cares of velocity states:
  void set_polarity_velocity(const double &persistence,space &l, sigs &s);
   
  short do_mutate(AffinitySpace &shape);
  // short do_mutate(AffinitySpace& shape, const double& p);
  short find_contact(states celltype, const long &i, space &l);
  //  short check_for_border(const long& i, space& l); // now no_border(const long&) in grid.h
  long find_mitosis_place(const double &p, bool forceit, const double &dx_max,
			  long * pp, space &l);
  double get_ss_receptor_ligand(const double &K, const double &r0, const double &s);
  double get_receptor_ligand(double &kplus, double &kminus,
			     double &s,
			     double &r0, double &rold);
  void signal_secretion(const long &i, const signal_molecule &sig_type, const short &mode,
			double &p, sigs &l);
  
  static double CXCL12recrit,CXCL13recrit;

  static double CXCL12crit,CXCL13crit;
  
  // tamoxifen sensitivity
  static bool tmx_MHC, tmx_MHC_noAgPresentation, tmx_MHC_noDivision, tmx_MHC_noTfhSignal;
  
private:
  void set_chemotaxis_polarity(const signal_molecule &chemokine, double * pol, sigs &s);
  
  static double use_glucose,use_oxygen,use_glucose_pro,use_oxygen_pro,critical_nutrient;
  
  static short use_specific_turning_angles;
};

/* derived cell objects with more than one fragment,
 * i.e. occupying more than one lattice point: */

class cellCB;

class frag_cell: public cell {

public:
   frag_cell();
   frag_cell(const frag_cell &x);
   ~frag_cell();
   void deliberate_memory();

   // Observable cell properties:
   // ===========================
   // Indices of other parts of the same cell
   long * fragments;
   // barycenter
   double barycenter[3];
   // Calculates and saves the barycenter of the cell
   long get_barycenter(space &l);
   // calculates the ratio of the long and the short axis of the cell
   // double longaxis,shortaxis;
   double get_long2short_axis(const long &li, const states &celltype, space &l);
   double get_long2short_axis(const long &li, const states &celltype, space &l,
                              double &longaxis, double &shortaxis);
   double get_elongation(const long &li, const states &celltype, space &l);
   // calculates and saves cell radius
   double get_radius(const short unsigned &d);
   // time setting
   void set_clock();
   // You can do set-up work for each test here.

   // Help and check routines:
   // ========================
   // Returns the position of the lattice point i on the fragment list:
   int where_fragment(const long &i);
   int check_connection(states celltype, const long &li, space &l);
   // long get_barycenter(const states& celltype, const long& li, space& l);
   void check_barycenter(space &l);
   // void check_barycenter(const states& celltype, const long& li, space& l);

   // Statistics and cell analysis:
   // =============================
   // for deformation properties:
   /* fraction of free neighbour points of the cell that are used for diffusion
    * of fragments of the cell */
   static double frac_average; // alpha_mean averaged over all cells
   double alpha_mean;
   // count the number of calls of get_free_nn(); for each cell
   long n_moves;
   // average displacements by directed moves
   long n_directed_moves;
   double performed2aimed_move,performed2aimed_move_now;

   // OPTIONAL: look for key MOVE_ANALYSIS for activation (in files cell.C and hyphasma.C)
   // counts the total number of tries and successes of movement (summed over all cells)
   /*
    * static long n_try_move,n_move_done,n_move_removed,n_move_forbidden,n_move_self,
    * n_move_forbidden_back,n_move_self_back;
    */
   // delete up to here.
   // Is it possible to keep it for more than one cell? Why not.
   // But what does it mean? It is the same number with better statistics.

   // OPERATORS
   // =========
   frag_cell&operator =(const frag_cell &x);

  protected:
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // +++ OPTIONS ...
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // maximal size of a cell in number of fragments
   static const int FRAGMENT_STEP = 77;

   /* switch between fragment diffusion with corrected D constant (=1)
    * and the use of force-based cell movement by shift of the bary center (=0)*/
   static const short use_D_correction = 0;

   /* Measure the performed movement and use a dynamic correction factor.
    * =1 use the dynamic correction factor only
    * =0 use the same with an additional cutoff criterion that stops the movement of subunits
    *    when the barycentre has been moved by one lattice constant
    */
   static const short use_dynamic_correction = 0;
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // +++ end OPTIONS.
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   // immobility clock for every fragment:
   long * t_immobile;

   // ACTIONS on the frag_cell
   // ========================
   short do_grow(states celltype, const long &li,
                 const int &target_vol,
                 const double &p_grow, const double &p_shrink,
                 space &l);

   short do_mitosis(states celltype,
                    const long &i, const long &li,
                    const long &newli, frag_cell &newCell,
                    space &l);

   // movemet of fragments (sum of diffusion and other forces:
   double fragmove(states celltype,
                   const long &li,
                   const double &tolerance_min,
                   const double &tolerance_steepness,
                   const double &elongation,
                   const double &eta_max,
                   const double &p_tension,
                   const double &K_elongation,
                   space &l);

   short int contact(states celltype, space &l);

   // Receptor handling
   // =================
   double bind_ss_receptor_ligand(const signal_molecule &sig,
                                  double &K, double &r0, double &old,
                                  space &l, sigs &s);
   double bind_receptor_ligand(const signal_molecule &sig,
                               double &kplus, double &kminus,
                               double &r0, double &rold,
                               space &l, sigs &s);

  private:
   // Geometry:
   // =========
   // number of border points (as found in the last call of fragdiffuse(...) )
   int borderpoints;
   // cell radius
   double radius;
   // elongation (deltax/radius)
   double elongation;
   // Routines for calculation of cell shape index
   double get_axis_length(const long &li, const states &celltype, space &l,
                          double * direction);
   long find_last_fragment(const long &li, const states &celltype, space &l,
                           double * direction);

   // time reset:
   // ===========
   void reset_clock(const int &i);

   // some flags and factors:
   // =======================
   /* fraction of free neighbour points of the cell that are used for diffusion
    * of fragments of the cell */
   double alpha,alpha_mean10;
   // fraction of fragments that (i) had to move and (ii) didn't move for some reason
   // and (iii) are not to be considered in the correction factor
   double flag_no_correction;

   // Fragment handling and analysis:
   // ===============================
   // deletes one fragment of the fragment-list of a cell without changing the lattice
   void del_fragment(const long &nr);
   void del_fragment(const long &nr, long * frags, int &max);

   void attribute_neighbors(states celltype, const long &li,
                            const long &j, const long &lj, frag_cell &targetcell,
                            space &l);

   // checks, if the barycenter belongs to the cell
   short check_for_ring(const states &celltype, const long &li, space &l);
   // Checks recursively the connectivity of the cell fragments
   int check_connection(states celltype,
                        const long &i, const long &li,
                        long * unconnected, int &n_unconnected,
                        space &l);

   /* starts from lattice point k
    * and find the nearast neighbor belonging to cell li */
   long find_nn2cell(space &l, long * k);
   int get_free_nn(space &l,
                   const long &li,
                   long * n_list,
                   const double &tolerance_min,
                   const double &tolerance_steepness,
                   double * target_point,
                   const long &excluded,
                   const long &include);

   short enclosed_object(const long &source, const long &target,
                         const long &li, const states &celltype, space &l);
   short find_path(const long &source, const long &target,
                   space &l, const long &li, const states &celltype,
                   double &max_dist, long &max_dist_index);
   short find_path(const long &source, const long &target,
                   space &l, const long &li, const states &celltype);
   short narrow_neck(const long &source, long &target,
                     const long &li, const states &celltype,
                     space &l);
   // some correction factors for the diffusion constant
   double volume_D_factor();
   double deformable_correction();

   double get_reshaping_force(const double &p_max_tension,
                              const double &K_elongation);
   double get_distances(char * which, double * ref_point, space &l, double * distance);
   void set_forbidden_counts(const long &j, double * target, double * old, space &l);
   void fragdiffuse(states celltype,
                    const double &tolerance_min,
                    const double &tolerance_steepness,
                    const long &li, double p,
                    const short &nocutoff,
                    const short &ordered,
                    double * target_point,
                    space &l);

   // FRAGMENT-HANDLING
   void show_fragments();
   void show_fragment_immobility();
   void show_barycenter(short unsigned&);
   // returns the position of i on the fragment list frags
   int where_fragment(const long &i, long * frags, const int &max);
};

// ===========================================
// (Vergleichs)-Operatoren
char operator ==(const cell &a, const cell &b);
char operator !=(const cell &a, const cell &b);

#endif
