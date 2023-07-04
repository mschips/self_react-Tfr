#ifndef i_grid
#define i_grid
#include "random.h"
#include "setparam.h"
#include "gridpoint.h"

enum grid_types {
   cellspace,molecules,alltypes
};

class grid {
  public:
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // +++ OPTIONS ...
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   /*
    * // The following is set in the parameter files:
    * // Border conditions in form of obstacles:
    * static const short obstacles = 3;
    * // 0: nothing; 1: wall; 2: wall with slit; 3: collagen matrix
    * // Height (level) of a wall in the reaction volume (% of # of lattice points from top)
    * static const double wall_level = 0.3;
    * // Width of the wall in lattice constants
    * static const int wall_width = 2;
    * // Number of slits to insert
    * static const int slit_number = 5;
    * // Width of slit (# of lattice points)
    * static const int slit_width = 1;
    * // density of collagen matrix (% of occupied points)
    * static const double collagen_density = 0.3;
    * // clustering of obstacles: percent of collagen points per cluster (0: no cluster)
    * static const double collagen_cluster = 0.00333;
    * // end of obstacles.
    */
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   // +++ end OPTIONS.
   // +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   grid();
   grid(const grid_types&, Parameter&, ofstream&);
   grid(const grid_types &t,
        const short &par_system,
        const short &par_dim,
        const double &par_dx,
        const double &par_radius,
        const double * size,
        const short &vol_shape,
        const short &obstacles,
        const double &wall_level,
        const double &collagen_density,
        const double &collagen_cluster,
        const int &wall_width,
        const int &slit_number,
        const int &slit_width,
        short theta_type,
        ofstream &ana);
   ~grid();
   grid_types grid_type;

   gridpoint * knot;

   static constexpr double pi = 3.141592654;

   // Initialisiere die Permutationen der Nachbarn
   pre_randomize nn_permuts;
   // save them as set of possible neighbors in (global variable)
   short * nn;

   short int system;

   // dimension of grid space and two times of it
   unsigned short dim,dim2;
   // total number of points on the grid
   long pointnumber;
   // points per dimension
   // int prodim,prodim2;
   int prodim2;
   int prodimvec[3];
   // grid constant
   double dx;

   // calculates a grid index from a more component vector
   long Index(const long * k);
   long Index(const double * k);
   long Index(const long * k, const int * per_dim);
   long Index(const double * k, const int * per_dim);
   long Index(const int &i, const int &j);
   long Index(const int &i, const int &j, const int &k);

   // calculates a more component vector k from a grid index wo
   void get_koord(long int wo, long * k);
   void get_koord(long int wo, double * k);

   // get index of neighbor in direction of v
   short get_nn_directed2(const double * v, short * nnindex);
   short get_nn_directed2(const double * v, const long &i, long * indices);
   long get_nn_directed2(const long &i, const double * v);

   // get vector between point i and j:
   void get_diff_vector(const long &i, const long &j, double * v);

   // Berechnet den (euklidischen) Abstand der Gitterpunkte n und n2
   double Abstand(const long int &n, const long int &n2);
   double get_distance2line(const long * c, const long * a, const long * b);

   double get_scalarproduct(const double * a, const double * b);
   double get_scalarproduct(const long * a, const double * b);
   double get_scalarproduct(const long * a, const long * b);
   double get_1norm(const long * k, const long * l);
   double get_2norm(const long * k, const long * l);
   double get_2norm(const long * n);
   double get_2norm(const double * n);
   double get_1norm(const double * k, const long * l);
   double get_2norm(const double * k, const long * l);
   double get_1norm(const double * k, const double * l);
   double get_2norm(const double * k, const double * l);

   // defines a normed vector with random direction:
   short get_random_direction(double * n);
   short get_distribution_direction(double * n);

   // Relation to grids with different resolution:
   long point_project(long &which, const double &dq, const int * nq);
   void guest_grid_project(const double &dq, const int * nq);

   short no_external(const long &i);
   short no_border(const long &i);

  private:
   // double dt; // global time step (info for diffusion)
   void get_diags(short int &n, short int &i, short int &j, short int &di,
                  short int &dj, short int &setno, long * k,
                  gridpoint &pwert);
   void fix_neighbors(long int wo);

   gauss_randomize thetas;

   int set_wall(const double &wall_level, const int &wall_width);
   void makehole(const long * start, const int &wall_width, const int &slit_width);
   void set_slit(const double &wall_level, const int &wall_width,
                 const int &slit_number, const int &slit_width);
   void set_collagen(const double &collagen_density, const double &collagen_cluster);

   void out_sphere();
   // short int sphere(long ind, int radius, states stat, long pos_ss);
};

#endif
