/* 05.02.2007 MMH:
 * Class for analysis of cell tracking data.
 * The routines can cope with single- or multi-node objects.
 *
 * Initialisation: Define an object and call init(..)
 * Best usage is achieved if all movements are sent to
 * the data-structure during simulation.
 * Evaluation is meant to be performed at the end of the
 * simulations only. The needed analysis has to be chosen.
 * Parameters of evaluation are set after the keyword OPTION.
 *
 * Bei allen Analysen wird vorausgesetzt, dass in allen N_OBJECT-Zellen
 * mindestens ein Eintrag extern geschrieben wurde. Werden extern weniger
 * Zellen verfolgt sollte trotzdem ein Eintrag mit status=nocell geschrieben
 * werden.
 *
 * Note that the files trace_endXXXX.out contain the end-points of the
 * cells of type XXXX irrespective of the cell-state at time TRACKUNTIL.
 *
 * 15.02.2007 MMH:
 * Write_tracks() fuer absolute und relative Koordinaten erweitert.
 *
 *
 #### Write more gle-file generation and try to generate N_OBJECTS dynamically
 */

#ifndef i_track
#define i_track

#include "dynarray.h"
#include "space.h"
#include <fstream>
#include <math.h>
using namespace std;

//MSchips --> TFR contact addd in enum -- track needs to be adapted (TODO!)
enum action_types {
  trackini,movement,polarisation,fdcdetachment,deathORend,incontact,offcontact,
    inTFRcontact,offTFRcontact
};

class TRACK {
  public:
   // Starting time of track
   static double TRACKFROM;
   // track until
   static double TRACKUNTIL;
   // time interval between measurements
   static double DELTA_T;
   // Number of velocities considered on the v-axis
   static int V_RESOLUTION,ALPHA_RESOLUTION;
   // interval of velocities in microns/min (for histograms)
   static double DELTA_V,DELTA_ALPHA;
   // Number of shape-indices considered on the shape-index-axis
   static int S_RESOLUTION;
   // interval of shape-indices (for histograms)
   static double DELTA_S;
   // how to treat no movement during contact to FDC or TC
   static bool INCLUDE_INCONTACT;

   TRACK();
   ~TRACK();

   void init(const double &xr_ext, const double &dt_ext,
             const short &d, const int * prodim, const int &objno,
             const double &from, const double &until, const double &t_interval,
             const int &v_resol, const double &delta_v,
             const int &alpha_resol, const double &delta_alpha,
             const int &s_resol, const double &delta_s);
   void get_trackraw();
   void get_trackraw(int Npertype);

   // A movement of the object with global index i is inserted into the motility data:
   // The new position at time t is r=(x,y,z).
   // Eventual shape-values are given optionally:
   void Write_movement(long i, states s, double t, double * r, double * pol,
                       action_types action_type);
   void Write_movement(long i, states s, double t, double * r, double * pol,
                       double elong, double l2s, action_types action_type);
   long add_new_track(states s, double t, double * r, double * pol,
                      double elong, double l2s);
   void Write_movement(long i, states s, double t, double * r, double * pol,
                       double elong, double l2s, double fct,
                       action_types action_type);
   void Stop_movement(long i, states s, double t, double * r, double * pol);
   void Write_files();
   void Write_files(bool generate_raw);

   // hess vector for zone analysis
   static double nhess[3];

   // allows to set the number of objects externally
   void set_N_OBJECTS(int n);

  private:
   // OPTION ++++++++++++++++++++++++++++++++++++++++++++++++
   // By setting this factor different from 1
   // the maximum number of tracked objects is increased.
   // This allows to increase the number of tracked objects
   // later on after initialistion of tracking data.
   const static int factor_N_OBJECTS = 10;
   // end OPTION ++++++++++++++++++++++++++++++++++++++++++++

   // Maximum number of tracked objects saved per file
   const static int MAX_COLUMNS = 10;

   void Show_tracks();

   typedef char fileno[5];
   void increment_fileno(fileno &tmp);
   int get_first_index(int obj);
   int get_last_index(int obj);
   int get_min_a(int obj);
   int get_max_a(int obj);
   double get_deltar(double * a, double * b);
   double get_deltar(const int &i, const int &ab, const int &bis);

   double get_scalarproduct(const double * a, const double * b);
   double get_2norm(const double * k, const double * l);

   void set_nhess(short &d);
   void get_color(int i, char color[10]);
   void get_marker(short c, char marker[10]);
   void get_cellname(short c, char cc[10]);
   void make_traces_r();
   void make_trace_types_a();
   void make_traces_a();
   void make_reach_dist_x();
   void make_gle();

   // dimension of the simulation
   short dim;
   // time and time step
   double time,dt;
   // lattice constant in microns
   double x_resolution;
   // Points per dimension
   int PPDIM[3];
   // Number of objects tracked
   int N_OBJECTS,max_N_OBJECTS;
   // Flag to save that data are written
   bool wrote_data;

   struct track_data {
      // state or type of the objects
      states status;
      // position of the object (r_z=0 in 2D)
      double r[3];
      // polarisation of the object
      double pol[3];
      // Shape values are saved for composed objects:
      double elongation,l2s_axis;
      // duration of last contact to FDC
      double fdc_contact_time;
      // Absolute time of the position:
      double t;
      // index for the action that induced the writing to TRACK
      action_types action;
   };

   bool used_cell_type[N_cells];

   // For every tracked object a dynamic array of positions and times is defined
   dynarray<track_data> * movement_list;

   // Write into files:
   // -----------------
   void Write_trace_begin(char prename[20]);
   void Write_trace_end(char prename[20]);
   void Write_trace_dead(char prename[20]);
   /* These files contain the start- end- and dead-positions of tracked cells
    */

   void Write_tracks(char prename[20], int from, bool relative);
   void Write_track_types(char prename[20], int from, bool relative);
   /* The tracks with the initial position of every object in the centre:
    * The time between two plotted points is DELTA_T.
    */

   // Analyse trans-migration between zones
   void Write_trans_zone(char prename[20], double * hess, double thick);

   void Write_cellnumber_in_zones(char prename[20], double * hess, double zpos0);

   void Write_reached_distance();
   /* The reached distance of the objects is plotted agains t and sqrt(t):
    * The time between two plotted distances is DELTA_T.
    */

   void Write_turning_angle();
   /* Histogram of turning angles derived from the polarity data.
    * Every polarity change is counted (thus no restriction by DELTA_T).
    */

   void Write_speed();
   /* The time course of the object speed:
    * Write in N_OBJECT columns (t speed).
    * The speed is calculated on the basis of time intervals of DELTA_T
    * and of the initial and final position in these intervals.
    *
    * Frequency of occurence of objects speeds:
    * Use V_RESOLUTION different speeds with DELTA_V interval.
    * The smallest speed is 0.
    * The speed is calculated on the basis of time intervals of DELTA_T
    * and of the initial and final position in these intervals.
    */

   void Write_fdc_contact_times();

   void Write_shape();
   /* The time course of the object shape index:
    * Long to short axis.
    * Measurements are done at intervals DELTA_T.
    *
    * Frequency of occurence of object shape:
    * Minimum shape index is 1.
    * Use S_RESOLUTION different speeds with DELTA_S interval.
    */
};

#endif
