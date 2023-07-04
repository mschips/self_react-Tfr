#include "grid.h"
#include <math.h>
#include <string.h>

// ### flag Neumann oder Moore einfuehren

grid::grid()
   : nn_permuts(4),thetas(0) {
   //  :nn_permuts(4),thetas(0.8,0.7,181,0,pi)
   nn = new short[4];
   for (unsigned short i = 0; i < 2; i++) {
      prodimvec[i] = 50;
   }
   prodim2 = prodimvec[0] * prodimvec[1];
   pointnumber = 2500;
   knot = new gridpoint[pointnumber];
}
grid::grid(const grid_types &t,
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
           ofstream &ana)
   : nn_permuts(2 * par_dim),thetas(theta_type) {
   // :nn_permuts(2*par_dim),thetas(0.8,0.7,181,0,3.141592654)
   cout << "Start grid constructor ...\n";
   grid_type = t;
   if (grid_type == cellspace) {
      ana << "Type of grid is cellspace.\n";
      cout << "Type of grid is cellspace.\n";
   } else {
      ana << "Type of grid is molecules.\n";
      cout << "Type of grid is molecules.\n";
   }

   system = par_system;
   if (system == 1) { ana << "Run under DOS/Windows!\n"; } else {
      ana << "Run under Unix/Linux!\n";
   }

   ana << " Use euclidean norm.\n";

   // ====================================
   ana << " Get grid parameters ...\n";
   // ====================================
   // dimension
   dim = par_dim;
   ana << "  grid dimension = " << dim << "\n";
   dim2 = 2 * dim;
   dx = par_dx;
   ana << "  grid constant = " << dx << "\n";
   if ((vol_shape == 0) || (vol_shape == 1)) {
      int prodimtmp = int ((2.0 + 1.e-08) * par_radius / dx) + 1;
      for (short ii = 0; ii < dim; ii++) {
         prodimvec[ii] = prodimtmp;
      }
      prodim2 = prodimtmp * prodimtmp;
      // Total number of points
      pointnumber = long (pow(double (prodimtmp),dim));
   } else if (vol_shape == 2) {
      // build a non-regular cube with double size[3]
      pointnumber = 1;
      for (short ii = 0; ii < dim; ii++) {
         prodimvec[ii] = int ((size[ii] + 1.e-08) / dx) + 1;
         pointnumber *= prodimvec[ii];
      }
      prodim2 = prodimvec[0] * prodimvec[1];
   }
   ana << "  total number of grid points = " << pointnumber << "\n";
   knot = new gridpoint[pointnumber];
   // Dimension of global neighbor field (containing sequences and no index)
   nn = new short[dim2];
   ana << "... done\n";
   cout << "... grid parameters read\n";

   // =====================================
   ana << " Initialize space-grid ...";
   // =====================================
   for (long n = 0; n < pointnumber; n++) {
      // Index anpassen
      knot[n].index = n;
      // Determine the corresponding coordinates and save them
      get_koord(n,knot[n].x);
      // Dann die Nachbarn bestimmen
      fix_neighbors(n);
   }

   cout << "... lattice initialized\n";
   ana << " done.\n";

   // Input of initial conditions
   // Initialize(logit);

   // ====================================================
   ana << " Set border conditions on the grid ...\n";
   // ====================================================
   if (vol_shape == 0) {
      ana << "     define a sphere ...\n";
      cout << "... sphere on the grid defined.\n";
      out_sphere();
   } else if (vol_shape == 1) {
      ana << "     define a cube ...\n";
      cout << "... cube on the grid defined.\n";
   } else if (vol_shape == 2) {
      ana << "     define a rectangular space ...\n";
      cout << "... rectangle on the grid defined.\n";
   } else { cerr << "Undefined space structure!\n"; exit(1); }

   if (obstacles == 1) {
      ana << "     put a wall on level " << wall_level << " ...\n";
      set_wall(wall_level,wall_width);
   }
   if (obstacles == 2) {
      ana << "     put a wall on level " << wall_level << " ...\n";
      ana << "     put " << slit_number
          << " slits of width " << slit_width << " into the wall ...\n";
      set_slit(wall_level,wall_width,slit_number,slit_width);
   }
   if (obstacles == 3) {
      ana << "     put collagen matrix of density " << collagen_density
          << " and clustering " << collagen_cluster
          << " ...\n";
      set_collagen(collagen_density,collagen_cluster);
   }
   ana << " ... done\n";
}
grid::grid(const grid_types &t, Parameter &par, ofstream &ana)
   : nn_permuts(2 * par.Value.DimSpace) {
   // ### This one is not used for the moment: delete?
   cout << "Start grid constructor ...\n";
   grid_type = t;
   if (grid_type == cellspace) {
      ana << "Type of grid is cellspace.\n";
      cout << "Type of grid is cellspace.\n";
   } else {
      ana << "Type of grid is molecules.\n";
      cout << "Type of grid is molecules.\n";
   }

   system = par.Value.system;
   if (system == 1) { ana << "Run under DOS/Windows!\n"; } else {
      ana << "Run under Unix/Linux!\n";
   }

   ana << " Use euclidean norm.\n";

   // ====================================
   ana << " Get grid parameters ...\n";
   // ====================================
   // dimension
   dim = par.Value.DimSpace;
   ana << "  grid dimension = " << dim << "\n";
   dim2 = 2 * dim;
   dx = par.Value.dx;
   ana << "  grid constant = " << dx << "\n";
   if ((par.Value.vol_shape == 0) || (par.Value.vol_shape == 1)) {
      int prodimtmp = int ((2.0 + 1.e-08) * par.Value.GC_radius / par.Value.dx) + 1;
      for (short ii = 0; ii < dim; ii++) {
         prodimvec[ii] = prodimtmp;
      }
      prodim2 = prodimtmp * prodimtmp;
      // Total number of points
      pointnumber = long (pow(double (prodimtmp),dim));
   } else if (par.Value.vol_shape == 2) {
      // build a non-regular cube with double size[3]
      pointnumber = 1;
      for (short ii = 0; ii < dim; ii++) {
         prodimvec[ii] = int ((par.Value.gridsize[ii] + 1.e-08) / dx) + 1;
         pointnumber *= prodimvec[ii];
      }
      prodim2 = prodimvec[0] * prodimvec[1];
   }
   ana << "  total number of grid points = " << pointnumber << "\n";
   knot = new gridpoint[pointnumber];
   // Dimension of global neighbor field (containing sequences and no index)
   nn = new short[dim2];
   ana << "... done\n";
   cout << "... grid parameters read\n";

   // =====================================
   ana << " Initialize space-grid ...";
   // =====================================
   for (long n = 0; n < pointnumber; n++) {
      // Index anpassen
      knot[n].index = n;
      // Determine the corresponding coordinates and save them
      get_koord(n,knot[n].x);
      // Dann die Nachbarn bestimmen
      fix_neighbors(n);
   }

   cout << "... lattice initialized\n";
   ana << " done.\n";

   // Input of initial conditions
   // Initialize(logit);

   // ====================================================
   ana << " Set border conditions on the grid ...\n";
   // ====================================================
   if (par.Value.vol_shape == 0) {
      ana << "     define a sphere ...\n";
      cout << "... sphere on the grid defined.\n";
      out_sphere();
   } else if (par.Value.vol_shape == 1) {
      ana << "     define a cube ...\n";
      cout << "... cube on the grid defined.\n";
   } else if (par.Value.vol_shape == 2) {
      ana << "     define a rectangular space ...\n";
      cout << "... rectangle on the grid defined.\n";
   } else { cerr << "Undefined space structure!\n"; exit(1); }

   if (par.Value.obstacles == 1) {
      ana << "     put a wall on level " << par.Value.wall_level << " ...\n";
      set_wall(par.Value.wall_level,par.Value.wall_width);
   }
   if (par.Value.obstacles == 2) {
      ana << "     put a wall on level " << par.Value.wall_level << " ...\n";
      ana << "     put " << par.Value.slit_number
          << " slits of width " << par.Value.slit_width << " into the wall ...\n";
      set_slit(par.Value.wall_level,par.Value.wall_width,
               par.Value.slit_number,par.Value.slit_width);
   }
   if (par.Value.obstacles == 3) {
      ana << "     put collagen matrix of density " << par.Value.collagen_density
          << " and clustering " << par.Value.collagen_cluster
          << " ...\n";
      set_collagen(par.Value.collagen_density,par.Value.collagen_cluster);
   }
   ana << " ... done\n";
}
grid::~grid() {
   // cout<<"in ~grid()\n";
   delete[] nn;
   delete[] knot;
}
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// Relation of Index and Vector ==========================
// =======================================================
// =======================================================
// =======================================================
// =======================================================

long int grid::Index(const long * k) {
   long int i = 0;
   /*
    * short int out=0;
    * if (k[0]>=0 && k[0]<prodim) i+=k[0]; else out=1;
    * if (k[1]>=0 && k[1]<prodim) i+=k[1]*prodim; else out=1;
    * if (dim>2) {
    * if (k[2]>=0 && k[2]<prodim) i+=k[2]*prodim2;
    * else out=1;
    * }
    * if (out==0) return i;
    * else return -1;
    */
   /* this one ist tested and works for equal size in each dimension!
    * if (k[0]>=0 && k[0]<prodim) {
    * i+=k[0];
    * if (k[1]>=0 && k[1]<prodim) {
    *  i+=k[1]*prodim;
    *  if (dim>2) {
    * if (k[2]>=0 && k[2]<prodim) i+=k[2]*prodim2;
    * else i=-1;
    *  }
    * }
    * else i=-1;
    * }
    * else i=-1;
    */
   if ((k[0] >= 0) && (k[0] < prodimvec[0])) {
      i += k[0];
      if ((k[1] >= 0) && (k[1] < prodimvec[1])) {
         i += k[1] * prodimvec[0];
         if (dim > 2) {
            if ((k[2] >= 0) && (k[2] < prodimvec[2])) { i += k[2] * prodim2; } else { i = -1; }
         }
      } else { i = -1; }
   } else { i = -1; }
   return i;
}
long int grid::Index(const long * k, const int * per_dim) {
   long int i = 0;
   /*
    * short int out=0;
    * if (k[0]>=0 && k[0]<per_dim) i+=k[0]; else out=1;
    * if (k[1]>=0 && k[1]<per_dim) i+=k[1]*per_dim; else out=1;
    * if (dim>2) {
    * if (k[2]>=0 && k[2]<per_dim) i+=k[2]*per_dim*per_dim;
    * else out=1;
    * }
    * if (out==0) return i;
    * else return -1;
    */
   if ((k[0] >= 0) && (k[0] < per_dim[0])) {
      i += k[0];
      if ((k[1] >= 0) && (k[1] < per_dim[1])) {
         i += k[1] * per_dim[0];
         if (dim > 2) {
            if ((k[2] >= 0) && (k[2] < per_dim[2])) {
               i += k[2] * per_dim[0] * per_dim[1];
            } else { i = -1; }
         }
      } else { i = -1; }
   } else { i = -1; }
   return i;
}
long grid::Index(const double * k) {
   long kl[dim];
   for (short i = 0; i < dim; i++) {
      kl[i] = long (k[i] + 0.5);
   }
   return Index(kl);
}
long grid::Index(const double * k, const int * per_dim) {
   long kl[dim];
   for (short i = 0; i < dim; i++) {
      kl[i] = long (k[i] + 0.5);
   }
   return Index(kl,per_dim);
}
long grid::Index(const int &i, const int &j, const int &k) {
   /*
    * long res=0;
    * short out=0;
    * if (i>=0 && i<prodim) res+=i; else out=1;
    * if (j>=0 && j<prodim) res+=j*prodim; else out=1;
    * if (k>=0 && k<prodim) res+=k*prodim2; else out=1;
    * if (out==1) { res=-1; cout<<"\nFehler in Index!\n"; exit(1); }
    * return res;
    */
   // return i+j*prodim+k*prodim2;
   return i + j * prodimvec[0] + k * prodim2;
}
long grid::Index(const int &i, const int &j) {
   /*
    * long res=0;
    * short out=0;
    * if (i>=0 && i<prodim) res+=i; else out=1;
    * if (j>=0 && j<prodim) res+=j*prodim; else out=1;
    * if (out==1) { res=-1; cout<<"\nFehler in Index!\n"; exit(1); }
    * return res;
    */
   // return i+j*prodim;
   return i + j * prodimvec[0];
}
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// Geometry ==============================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================

void grid::get_koord(long int wo, long * k) {
   long int reduce,p;
   for (int n = dim - 1; n >= 0; n--) {
      p = 1;
      for (short i = 0; i < n; i++) {
         p *= prodimvec[i];
      }
      // p=long(pow(double(prodim),n));
      reduce = wo / p;
      k[n] = reduce;
      wo -= reduce * p;
   }
   // for (short i=0; i<dim; i++) k[i]=knot[wo].x[i];
}
void grid::get_koord(long int wo, double * k) {
   long int reduce,p;
   for (int n = dim - 1; n >= 0; n--) {
      p = 1;
      for (short i = 0; i < n; i++) {
         p *= prodimvec[i];
      }
      // p=long(pow(double(prodim),n));
      reduce = wo / p;
      k[n] = reduce;
      wo -= reduce * p;
   }
   // for (short i=0; i<dim; i++) k[i]=double(knot[wo].x[i]);
}
double grid::get_scalarproduct(const double * a, const double * b) {
   double sum = 0.0;
   for (short d = 0; d < dim; d++) {
      sum += (a[d] * b[d]);
   }
   return sum;
}
double grid::get_scalarproduct(const long * a, const double * b) {
   double sum = 0.0;
   for (short d = 0; d < dim; d++) {
      sum += (double (a[d]) * b[d]);
   }
   return sum;
}
double grid::get_scalarproduct(const long * a, const long * b) {
   double sum = 0.0;
   for (short d = 0; d < dim; d++) {
      sum += (double (a[d]) * double (b[d]));
   }
   return sum;
}
double grid::get_distance2line(const long * c, const long * a, const long * b) {
   /* returns the distance of gridpoint c to the line defined by the gridpoints a and b */
   double c_a[dim];
   double b_a[dim];
   for (short i = 0; i < dim; i++) {
      c_a[i] = double (c[i] - a[i]);
      b_a[i] = double (b[i] - a[i]);
   }
   double whatever = get_scalarproduct(c_a,b_a) / get_2norm(b_a);
   return sqrt(get_scalarproduct(c_a,c_a) - whatever * whatever);
}
short grid::get_nn_directed2(const double * w, short * nnindex) {
   /* Returns the index on the neighbor list corresponding
    * mostly to the direction of w.
    * The lattice-index of this neighbour is then called
    * with knot[bla].near_n[get_nn_directed2()].
    * nnindex[1++] contains further indices with equal projection on w.
    * The number of equal projections is returned as result of the routine.
    ### Note that diagonal neighbors are not considered here!
    */
   // Neighbors are saved in the rank: 3 left neighbors then 3 right neighbors
   // normalisation of w:
   double v[dim];
   double normfaktor = sqrt(get_scalarproduct(w,w));
   for (short d = 0; d < dim; d++) {
      v[d] = w[d] / normfaktor;
   }
   double u[dim];
   double max = -10.0;  // The norm of v is one therefore this is smaller than every value!
   double actual;
   short count = 0;
   for (short d = 0; d < dim; d++) {
      u[d] = 0.0;
      nnindex[d] = -1;
   }
   for (short d = 0; d < dim; d++) {
      // neighbor dim+d:
      u[d] = 1.0;
      actual = get_scalarproduct(u,v);
      if (actual > max) {
         count = 0;
         nnindex[count] = dim + d;
         max = actual;
      } else if (actual == max) {
         ++count;
         nnindex[count] = dim + d;
      }
      // neighbor d:
      u[d] = -1.0;
      actual = get_scalarproduct(u,v);
      if (actual > max) {
         count = 0;
         nnindex[count] = d;
         max = actual;
      } else if (actual == max) {
         ++count;
         nnindex[count] = d;
      }
      // reset:
      u[d] = 0.0;
   }
   // cout<<"v=("<<v[0]<<","<<v[1]<<") count="<<count+1
   //    <<" nnindex=("<<nnindex[0]<<","<<nnindex[1]<<")\n";
   return count + 1;
}
short grid::get_nn_directed2(const double * v, const long &i, long * indices) {
   // as get_nn_directed(const double *) but defining the indices on the lattice in "indices"
   // and returning the number of values.
   short tmp[dim];
   short howmany = get_nn_directed2(v,tmp);
   for (short d = 0; d < howmany; d++) {
      indices[d] = knot[i].near_n[tmp[d]];
   }
   return howmany;
}
long grid::get_nn_directed2(const long &i, const double * v) {
   /* as get_nn_directed2(const double*, const long&, long*) but
    * directly returning the grid index of the neighbour in direction of v
    * and chosing one randomly if more than one result is found.
    */
   short tmp[dim2];
   short nn_howmany = get_nn_directed2(v,tmp);
   long nn_index = knot[i].near_n[tmp[0]];
   if (nn_howmany > 1) { nn_index = knot[i].near_n[tmp[irandom(nn_howmany)]]; }
   return nn_index;
}
void grid::get_diff_vector(const long &i, const long &j, double * v) {
   for (short a = 0; a < dim; a++) {
      v[a] = double (knot[j].x[a] - knot[i].x[a]);
   }
}
double grid::get_1norm(const double * k, const double * l) {
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      distance += fabs(k[i] - l[i]);
   }
   return distance;
}
double grid::get_2norm(const double * k, const double * l) {
   double distance = 0.;
   for (short i = 0; i < dim; i++) {
      double x = k[i] - l[i];
      distance += (x * x);
   }
   return sqrt(distance);
}
double grid::get_2norm(const double * n) {
   double distance = 0.;
   for (short i = 0; i < dim; i++) {
      distance += (n[i] * n[i]);
   }
   return sqrt(distance);
}
double grid::get_2norm(const long * n) {
   double distance = 0.;
   for (short i = 0; i < dim; i++) {
      distance += (n[i] * n[i]);
   }
   return sqrt(distance);
}
double grid::get_1norm(const double * k, const long * l) {
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      distance += fabs(k[i] - l[i]);
   }
   return distance;
}
double grid::get_2norm(const double * k, const long * l) {
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      double x = k[i] - double (l[i]);
      distance += (x * x);
   }
   return sqrt(distance);
}
double grid::get_1norm(const long * k, const long * l) {
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      distance += fabs(double (k[i] - l[i]));
   }
   return distance;
}
double grid::get_2norm(const long * k, const long * l) {
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      double x = k[i] - l[i];
      distance += (x * x);
   }
   return sqrt(distance);
}
double grid::Abstand(const long int &n, const long int &n2) {
   return get_2norm(knot[n].x,knot[n2].x);
}
short grid::get_random_direction(double * n) {
   double phi = drandom(2.0 * pi);
   n[0] = cos(phi);
   n[1] = sin(phi);
   if (dim == 3) {
      double theta = drandom(pi);
      n[0] *= sin(theta);
      n[1] *= sin(theta);
      n[2] = cos(theta);
   }
   // cout<<"n=("<<n[0]<<","<<n[1]<<")\n";
   return 1;
}
/*
 * short grid::get_gauss_random_direction(double * n) {
 * double phi;
 * if (dim==3) phi=drandom(2.0*pi);
 * else {
 *  phi=thetas.get_gauss_distributed_value(); // returns a value between [0,pi]
 *  // Project this on 2pi while keeping the turning angle distribution Gaussian
 *  if (irandom(2)==1) {
 *    phi=2.0*pi-phi;
 *  }
 * }
 * n[0]=cos(phi);
 * n[1]=sin(phi);
 * if (dim==3) {
 *  double theta=thetas.get_gauss_distributed_value();
 *  n[0]*=sin(theta);
 *  n[1]*=sin(theta);
 *  n[2]=cos(theta);
 * }
 * //cout<<"n=("<<n[0]<<","<<n[1]<<")\n";
 * return 1;
 * }
 */

short grid::get_distribution_direction(double * n) {
   // TURNS vector n by an angle from a gauss-weighted random distribution!
   if (dim == 2) {
      double phi = thetas.get_distribution_value();   // returns a value between [0,pi]
      // randomly chose to turn either to the right or to the left
      if (irandom(2) == 1) { phi *= -1.0; }
      // turn the original vector n
      double tmp[2];
      tmp[0] = cos(phi) * n[0] - sin(phi) * n[1];
      tmp[1] = sin(phi) * n[0] + cos(phi) * n[1];
      for (short a = 0; a < dim; a++) {
         n[a] = tmp[a];
      }
   } else {
      // if (dim==3)
      // double nnorm=sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
      // cout<<"polarity=("<<n[0]<<","<<n[1]<<","<<n[2]<<") norm="<<nnorm<<"\n";

      double theta = thetas.get_distribution_value();
      double phi = drandom(2.0 * pi);
      double tmp[3],ttmp[3];

      // find the phi of the old polarity:
      double nphi = 0.;
      if ((n[0] == 0.) && (n[1] == 0.)) {
         nphi = 0.;
      } else { nphi = acos(n[0] / sqrt(n[0] * n[0] + n[1] * n[1])); }
      if (n[1] < 0) { nphi = 2.0 * pi - nphi; }
      // cout<<"nphi="<<nphi<<"\n";

      // turn n onto the x-z-plane by rotation of -nphi around the z-axis:
      nphi *= -1.0;
      tmp[0] = cos(nphi) * n[0] - sin(nphi) * n[1];
      tmp[1] = sin(nphi) * n[0] + cos(nphi) * n[1];
      tmp[2] = n[2];
      // cout<<"in x-z=("<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<")\n";

      // turn the vector in the z-plane by theta around the y-axis
      ttmp[0] = cos(theta) * tmp[0] + sin(theta) * tmp[2];
      ttmp[1] = tmp[1];
      ttmp[2] = -1.0 * sin(theta) * tmp[0] + cos(theta) * tmp[2];
      // cout<<"turned by theta=("<<ttmp[0]<<","<<ttmp[1]<<","<<ttmp[2]<<")\n";

      // turn back by nphi around z-axis
      nphi *= -1.0;
      tmp[0] = cos(nphi) * ttmp[0] - sin(nphi) * ttmp[1];
      tmp[1] = sin(nphi) * ttmp[0] + cos(nphi) * ttmp[1];
      tmp[2] = ttmp[2];
      // cout<<"turned back by nphi=("<<tmp[0]<<","<<tmp[1]<<","<<tmp[2]<<")\n";
      // double ntheta=acos(tmp[0]*n[0]+tmp[1]*n[1]+tmp[2]*n[2]);
      // cout<<"theta="<<theta<<"; newtheta="<<ntheta<<"\n";

      // turn the new vector tmp by phi around the old vector n
      ttmp[0] = (cos(phi) + n[0] * n[0] * (1.0 - cos(phi))) * tmp[0]
                + (n[0] * n[1] * (1.0 - cos(phi)) - n[2] * sin(phi)) * tmp[1]
                + (n[0] * n[2] * (1.0 - cos(phi)) + n[1] * sin(phi)) * tmp[2];
      ttmp[1] = (n[0] * n[1] * (1.0 - cos(phi)) + n[2] * sin(phi)) * tmp[0]
                + (cos(phi) + n[1] * n[1] * (1.0 - cos(phi))) * tmp[1]
                + (n[1] * n[2] * (1.0 - cos(phi)) - n[0] * sin(phi)) * tmp[2];
      ttmp[2] = (n[0] * n[2] * (1.0 - cos(phi)) - n[1] * sin(phi)) * tmp[0]
                + (n[1] * n[2] * (1.0 - cos(phi)) + n[0] * sin(phi)) * tmp[1]
                + (cos(phi) + n[2] * n[2] * (1.0 - cos(phi))) * tmp[2];

      // save the result as new polarity
      // ntheta=acos(ttmp[0]*n[0]+ttmp[1]*n[1]+ttmp[2]*n[2]);
      // cout<<"theta="<<theta<<"; nnewtheta="<<ntheta<<"\n";
      double newnorm = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
      for (short a = 0; a < dim; a++) {
         n[a] = ttmp[a] / newnorm;
      }
      // cout<<"new-polarity=("<<n[0]<<","<<n[1]<<","<<n[2]<<") norm="<<nnorm<<"\n\n";
   }
   return 1;
}
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// Neighbours ============================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================

void grid::get_diags(short int &n,
                     short int &i, short int &j,
                     short int &di,short int &dj,
                     short int &setno, long * k,
                     gridpoint &pwert) {
   /* Berechnet den Index eines diagonalen Nachbarn zum Punkt k.
    * Es werden die Nachbarn in einer Ebene berechnet.
    * Dabei werden die Koordinaten i und j um di bzw. dj veraendert.
    * Im Fall des ueberschreitens des Gitterrandes erhaelt man -1.
    * Die Speicherung wird an Position n + setno * dim in pwert vorgenommen.
    * n bezeichnet eine (nummerierte) Dimension und
    * setno durchlaeuft die Werte 0..3
    */
   long tmp[dim];
   // Change of coordinates
   long int zahli = k[i] + di;
   long int zahlj = k[j] + dj;
   // save a copy of k[] in tmp[]
   for (short int ii = 0; ii < dim; ii++) {
      tmp[ii] = k[ii];
   }
   // boundary reached?
   // if (zahli>=prodim[i] || zahli<0) zahli=-1;
   // if (zahlj>=prodim || zahlj<0) zahlj=-1;
   if ((zahli >= prodimvec[i]) || (zahli < 0)) { zahli = -1; }
   if ((zahlj >= prodimvec[j]) || (zahlj < 0)) { zahlj = -1; }
   // save new coordinates
   /* cout<<"zahli="<<zahli<<"  zahlj="<<zahlj<<"\n";
    * cout<<"n="<<n<<"  setno="<<setno<<"  dim="<<dim<<"\n"; */
   if ((zahli >= 0) && (zahlj >= 0)) {
      tmp[i] = zahli;
      tmp[j] = zahlj;
      pwert.diag_n[n + setno * dim] = Index(tmp);
   } else { pwert.diag_n[n + setno * dim] = -1; }
}
void grid::fix_neighbors(long int wo) {
   /* Calculates the indizes of the next neighbors and
    * calls the routine for calculation of diagonal neighbors.
    * Only working for 2,3 dimensions!!! */

   // Load the gridpoint wo on the lattice
   gridpoint pwert = knot[wo];
   // Calculate the coordinates on the lattice corresponding to wo
   long k[dim];
   get_koord(wo,k);
   // cout<<"  k=("<<k[0]<<" "<<k[1]<<" "<<k[2]<<")\n";
   // Temporary variables
   long int zahl;
   short int i,j,di,dj,setno,m;
   long tmp[dim];
   // Go through the dimensions
   short int n;
   for (n = 0; n < dim; n++) {
      // Left neighbors
      zahl = k[n] - 1;
      // transmit coordinates
      for (m = 0; m < dim; m++) {
         tmp[m] = k[m];
      }
      // boundary condition
      if (zahl < 0) { zahl = -1; }
      // save coordinate of the left neighbor
      tmp[n] = zahl;
      // recalculate corresponding Index for the neighbor tmp
      if (zahl >= 0) { pwert.near_n[n] = Index(tmp); } else { pwert.near_n[n] = -1; }
      // The same for right neighbors
      zahl = k[n] + 1;
      for (m = 0; m < dim; m++) {
         tmp[m] = k[m];
      }
      // if (zahl>=prodim) zahl=-1;
      if (zahl >= prodimvec[n]) { zahl = -1; }
      tmp[n] = zahl;
      if (zahl >= 0) { pwert.near_n[n + dim] = Index(tmp); } else { pwert.near_n[n + dim] = -1; }

      // Now going to the diagonal neighbors:
      // Here, two coordinates are changed for each neighbor.
      if (dim == 3) {
         switch (n) {
            case 0:
               i = 1;
               j = 2;
               break;

            case 1:
               i = 0;
               j = 2;
               break;

            case 2:
               i = 0;
               j = 1;
               break;
         }
         // Parameters diag_n[n]=...
         di = -1;
         dj = -1;
         setno = 0;
         // calc diagonal neighbors
         get_diags(n,i,j,di,dj,setno,k,pwert);
         // Parameters diag_n[n+3]=...
         di = -1;
         dj = 1;
         setno = 1;
         // calc diagonal neighbors
         get_diags(n,i,j,di,dj,setno,k,pwert);
         // Parameters diag_n[n+6]=...
         di = 1;
         dj = -1;
         setno = 2;
         // calc diagonal neighbors
         get_diags(n,i,j,di,dj,setno,k,pwert);
         // Parameters diag_n[n+9]=...
         di = 1;
         dj = 1;
         setno = 3;
         // calc diagonal neighbors
         get_diags(n,i,j,di,dj,setno,k,pwert);
      }
   }
   if (dim == 2) {
      i = 0;
      j = 1;
      // Parameters diag_n[0]=... links unten
      di = -1;
      dj = -1;
      n = 0;
      setno = 0;
      // calc diagonal neighbors
      get_diags(n,i,j,di,dj,setno,k,pwert);
      // Parameters diag_n[1]=... links oben
      di = -1;
      dj = 1;
      n = 1;
      setno = 0;
      // calc diagonal neighbors
      get_diags(n,i,j,di,dj,setno,k,pwert);
      // Parameters diag_n[1]=... rechts unten
      di = 1;
      dj = -1;
      n = 2;
      setno = 0;
      // calc diagonal neighbors
      get_diags(n,i,j,di,dj,setno,k,pwert);
      // Parameters diag_n[1]=... rechts oben
      di = 1;
      dj = 1;
      n = 3;
      setno = 0;
      // calc diagonal neighbors
      get_diags(n,i,j,di,dj,setno,k,pwert);
   }
   /* cout<<"Neighbors: ";
    * for (i=0; i<2*dim; i++) cout<<pwert.near_n[i]<<" ";
    * cout<<"\nDiags: ";
    * for (i=0; i<4*dim; i++) cout<<pwert.diag_n[i]<<" ";
    * cout<<"\n";
    */
   // Write back on the lattice
   knot[wo] = pwert;
}
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// Guest grids ===========================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================

long grid::point_project(long &which, const double &dq, const int * nq) {
   // Returns the index of point q on a grid with lattice constant dq
   // that is nearest to the point p on this lattice.
   // Returns -1 if the nearest point is outside the q-grid.
   double q[dim];
   double p[dim];
   get_koord(which,p);
   for (short i = 0; i < dim; i++) {
      q[i] = p[i] * dx / dq;
      if (q[i] > nq[i] - 1) { return -1; }
   }
   return Index(q,nq);
}
void grid::guest_grid_project(const double &dq, const int * nq) {
   for (long i = 0; i < pointnumber; i++) {
      knot[i].guest_grid_pointer = point_project(i,dq,nq);
   }
}
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// Border conditions =====================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================
// =======================================================

int grid::set_wall(const double &wall_level, const int &wall_width) {
   // int height=int(wall_level*(double(prodim-1)+0.5));
   int height = int (wall_level * (double (prodimvec[dim - 1] - 1) + 0.5));
   // always set the wall in the y-direction, coordinate 1
   // determine first and last level of the wall
   int a = -wall_width / 2;
   int b = wall_width / 2;
   if (((wall_width / 2) * 2) == wall_width) {
      if (irandom(2) == 0) { ++a; } else { --b; }
   }
   // go through all slices of the wall
   for (int w = a; w <= b; w++) {
      // go through the wall
      // for (int n=0; n<prodim; n++) {
      for (int n = 0; n < prodimvec[0]; n++) {
         if (dim == 2) { knot[Index(n,height + w)].status = external; } else {
            // dim==3
            // for (int m=0; m<prodim; m++)
            for (int m = 0; m < prodimvec[1]; m++) {
               knot[Index(n,m,height + w)].status = external;
            }
         }
      }
   }
   return height + a;
}
void grid::makehole(const long * start, const int &wall_width, const int &slit_width) {
   long v[dim];
   for (short d = 0; d < dim; d++) {
      v[d] = start[d];
   }
   // v[dim-1] enthaelt die Position der Wand
   // diese hat eine Tiefe von "wall_width"
   // get the number of neighbours to change per direction
   int n = slit_width / 2;

   // find if width even or odd
   short even = ((slit_width / 2) * 2 == slit_width);
   // if it is even, n has to be reduced by 1 on one of the directions in each dimension!

   // Change all gridpoints but not the borders
   int x = v[0];
   int y = v[1];
   for (int i = x - n + even; i <= x + n - even; i++) {
      v[0] = i;
      if (dim == 3) {
         for (int j = y - n + 1; j < x + n; j++) {
            v[1] = j;
            for (int k = 0; k < wall_width; k++) {
               v[2] = start[2] + k;
               if (knot[Index(v)].status
                   != external) { cout << "Error in grid::makehole(...)!\n"; }
               knot[Index(v)].status = nothing;
            }
         }
      } else {
         for (int j = 0; j < wall_width; j++) {
            v[1] = start[1] + j;
            if (knot[Index(v)].status != external) { cout << "Error in grid::makehole(...)!\n"; }
            knot[Index(v)].status = nothing;
         }
      }
   }
   // In the case even==0 we are done.

   if (even == 1) {
      // Now treat border gridpoints of the hole (x+-n, y+-n):
      // decide right or left for each dimension
      if (dim == 2) {
         if (irandom(2) == 0) { v[0] = x - n; } else { v[0] = x + n; }
         for (int j = 0; j < wall_width; j++) {
            v[1] = start[1] + j;
            if (knot[Index(v)].status != external) { cout << "Error in grid::makehole(...)!\n"; }
            knot[Index(v)].status = nothing;
         }
      } else {
         // if dim==3
         short leftx = 0;
         short lefty = 0;
         if (irandom(2) == 0) { leftx = 1; }
         if (irandom(2) == 0) { lefty = 1; }
         if (leftx == 1) { x -= n; } else { x += n; }
         if (lefty == 1) { y -= n; } else { y += n; }
         // (x,y) zeigt auf die Ecke, die "external" wird
         long count = 0;    // correct n-1 values from there
         while (count < n - 1) {
            if (leftx == 1) { v[0] = x + count; } else { v[0] = x - count; }
            v[1] = y;
            for (int k = 0; k < wall_width; k++) {
               v[2] = start[2] + k;
               if (knot[Index(v)].status
                   != external) { cout << "Error in grid::makehole(...)!\n"; }
               knot[Index(v)].status = nothing;
            }
            if (count > 0) {
               if (lefty == 1) { v[1] = y + count; } else { v[1] = y - count; }
               v[0] = x;
               for (int k = 0; k < wall_width; k++) {
                  v[2] = start[2] + k;
                  if (knot[Index(v)].status
                      != external) { cout << "Error in grid::makehole(...)!\n"; }
                  knot[Index(v)].status = nothing;
               }
            }
            ++count;
         }
      }
   }  // end if (even==1)
}
void grid::set_slit(const double &wall_level, const int &wall_width,
                    const int &slit_number, const int &slit_width) {
   int height = set_wall(wall_level,wall_width);
   /* Find the length of the wall!
    * This may be different of "prodim" when a non-cubic reaction volume has been defined.
    * A criterion is to look on the grid points on the level "height": If these
    * grid points have at least one neighbour with .status!=external the wall is within
    * the reaction volume. In 3D this is a measure for the surface of the wall.
    * Count the number of these points:
    */
   int wallsurface = 0;
   // for (int n=0; n<prodim; n++)
   for (int n = 0; n < prodimvec[0]; n++) {
      for (int nn = 0; nn < dim2; nn++) {
         if (dim == 2) {
            if (knot[knot[Index(n,height)].near_n[nn]].status != external) {
               ++wallsurface;
               nn = dim2;
            }
         } else {
            // dim==3
            // for (int m=0; m<prodim; m++)
            for (int m = 0; m < prodimvec[1]; m++) {
               if (knot[knot[Index(n,m,height)].near_n[nn]].status != external) {
                  ++wallsurface;
                  nn = dim2;
               }
            }
         }
      }
   }
   // wall surface defined!

   // Get the surface of the holes
   long holesurface = slit_number * slit_width;
   if (dim == 3) { holesurface *= slit_width; }
   // Stop if the holesurface is too large:
   if (holesurface > wallsurface) {
      cout << "There are more holes than wall: holesurface=" << holesurface
           << " and wallsurface=" << wallsurface << "!\n";
      exit(1);
   }

   /* Distribute the holes regularly on the surface.
   * This assumes a compact (not fragmented) wall.*/
   /* The distance between two holes in each dimension
    * is just sqrt(wallsurface/number), sqrt only for 3D:*/
   int dist;
   if (dim == 2) {
      dist = int (double (wallsurface) / double (slit_number) + 0.5);
   } else { dist = int (sqrt(double (wallsurface) / double (slit_number)) + 0.5); }
   /* Where to set the first hole?
    * Start at the center of the reaction volume?
    * But what if the reaction volume is not symmetric to the center?
    * Then it is a starting point as any other. So take it.
    * Or better, take it and calculate the point most up and left that
    * is some "dist-steps" apart. This enables to go in one way only
    * through the lattice.
    */
   // Find the centre:
   long v[dim];
   // for (short d=0; d<dim-1; d++) v[d]=prodim/2;
   for (short d = 0; d < dim - 1; d++) {
      v[d] = prodimvec[d] / 2;
   }
   v[dim - 1] = height;
   // Find the left and up starting point
   for (short d = 0; d < dim - 1; d++) {
      while (v[d] - dist >= 0) {
         v[d] -= dist;
      }
   }
   // Go through all points and make a hole there
   // while (v[0]<prodim) {
   while (v[0] < prodimvec[0]) {
      if (dim == 3) {
         // while (v[1]<prodim) {
         while (v[1] < prodimvec[1]) {
            makehole(v,wall_width,slit_width);
            v[1] += dist;
         }
         v[1] -= (v[1] / dist) * dist;
      } else { makehole(v,wall_width,slit_width); }
      v[0] += dist;
   }
   // that's it.
}
void grid::set_collagen(const double &collagen_density, const double &collagen_cluster) {
   // Use "collagen_density" and "collagen_cluster"
   /*
    * "collagen_density" is interpreted as the fraction of points in the
    *    reaction volume to be occupied by collagen. The number of
    * points in the reaction volume is counted and includes all
    * points that are not "external" at the moment of the call. All
    * objects within the reaction volume will be ignored and possibly
    * overwritten.
    * !!! Thus call this routine before the definition of objects !!!
    * "collagen_cluster" is the fraction of collagen points that are to
    *    be inserted in each cluster. Thus after having inserted a first collagen
    * point, a number of further points (corresponding to this fraction)
    * is inserted as direct neighbours. Note, that all clusters will have
    * the same size but random shape.
    * For "collagen_cluster==0.0" each collagen point is inserted randomly.
    */

   if (dim == 2) {
      /*
       * dim==2: A collagen matrix is represented by a random set of
       *        single points. A fibre lying within the plane represented
       *    in the simulation has to be ignored, as no possibility to
       *    walk around the fibre exist, i.e. to deviate into the 3rd
       *    dimension. What about connected points. For the moment
       *    they may come up by chance.
       */
      // get the total number of points
      long rea_vol = 0;
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int j = 0; j < prodimvec[1]; j++) {
            if (knot[Index(i,j)].status != external) { ++rea_vol; }
         }
      }
      // get the number of collagen points to insert:
      long n = long (double (rea_vol) * collagen_density);
      // get the number of points per cluster:
      long c = 1;
      if (collagen_cluster > 0) { c = long (collagen_cluster * double (n) + 0.5); }
      // calculate cluster radius:
      long r_cluster = long (sqrt(double (c) / pi) + 0.5);
      //    cout<<"(points per fragment="<<c<<")... ";
      long actual_index = -1;
      // insert the points:
      long i = 0;
      while (i < n) {
         short found = 0;
         while (found == 0) {
            long x = irandom(prodimvec[0]);
            long y = irandom(prodimvec[1]);
            if (knot[Index(x,y)].status != external) {
               // remember the successful index (for clusters)
               actual_index = Index(x,y);
               knot[actual_index].status = external;
               ++i;
               // remember success in order to switch to the next cluster to insert
               found = 1;
               // if clusters are wanted, initiate further points to be clustered with it
               if (c > 1) {
                  for (long xx = x - r_cluster; xx < x + r_cluster; xx++) {
                     for (long yy = y - r_cluster; yy < y + r_cluster; yy++) {
                        if ((xx >= 0) && (xx < prodimvec[0]) && (yy >= 0)
                            && (yy < prodimvec[1])) {
                           long this_index = Index(xx,yy);
                           if ((Abstand(this_index,
                                        actual_index) <= double (r_cluster)) && (i < n)) {
                              if (knot[this_index].status != external) {
                                 knot[this_index].status = external;
                                 ++i;
                              }
                           }
                        }
                     }
                  }
               }      // end if (c>1)
            }
         }
      }
   }  // end if (dim==2)
   else {
      // Im Prinzip mit dim==2 identisch bis auf 3 Befehle, die markiert sind !!!
      /*
       * dim==3: For the moment points or extended obstacles are randomly
       *        distributed and set to status=external.
       *        It may be interesting to construct real fibres here. ###
       */
      // get the total number of points
      long rea_vol = 0;
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int j = 0; j < prodimvec[1]; j++) {
            for (int k = 0; k < prodimvec[2]; k++) {
               // additional 3D line
               if (knot[Index(i,j,k)].status != external) { ++rea_vol; }
            }
         }
      }
      // get the number of collagen points to insert:
      long n = long (double (rea_vol) * collagen_density);
      // get the number of points per cluster:
      long c = 1;
      if (collagen_cluster > 0) { c = long (collagen_cluster * double (n) + 0.5); }
      //    cout<<"(points per fragment="<<c<<")... ";
      long r_cluster = long (pow(3. * double (c) / (4. * pi),1. / 3.) + 0.5);
      long actual_index = -1;
      // insert the points:
      long i = 0;
      while (i < n) {
         short found = 0;
         while (found == 0) {
            long x = irandom(prodimvec[0]);
            long y = irandom(prodimvec[1]);
            long z = irandom(prodimvec[2]);
            if (knot[Index(x,y,z)].status != external) {
               // remember the successful index (for clusters)
               actual_index = Index(x,y,z);
               knot[actual_index].status = external;
               ++i;
               // remember success in order to switch to the next cluster to insert
               found = 1;
               // if clusters are wanted, initiate further points to be clustered with it
               if (c > 1) {
                  for (long xx = x - r_cluster; xx <= x + r_cluster; xx++) {
                     for (long yy = y - r_cluster; yy <= y + r_cluster; yy++) {
                        for (long zz = z - r_cluster; zz <= z + r_cluster; zz++) {
                           if ((xx >= 0) && (xx < prodimvec[0]) && (yy >= 0)
                               && (yy < prodimvec[1])
                               && (zz >= 0) && (zz < prodimvec[2])) {
                              long this_index = Index(xx,yy,zz);
                              if ((Abstand(this_index,
                                           actual_index) <= double (r_cluster))
                                  && (i < n)) {
                                 if (knot[this_index].status != external) {
                                    knot[this_index].status = external;
                                    ++i;
                                 }
                              }
                           }
                        }
                     }
                  }
               }      // end if (c>1)
            }
         }
      }
   }  // end if (dim==3)
   cout << "... collagen matrix defined.\n";
}
void grid::out_sphere() {
   // das Zentrum:
   long k[dim];
   for (unsigned short i = 0; i < dim; i++) {
      k[i] = long (prodimvec[i] / 2);
   }
   long zero = Index(k);
   // find axis of minimum extension
   int minprodim = prodimvec[0];
   for (unsigned short i = 1; i < dim; i++) {
      if (prodimvec[i] < minprodim) {
         minprodim = prodimvec[i];
      }
   }
   int radius = int ((minprodim - 1) / 2);
   // Save external outside the sphere
   // Test euclidean distance for all grid points
   for (long i = 0; i < pointnumber; i++) {
      if (Abstand(i,zero) > radius) { knot[i].status = external; }
   }
}
short grid::no_external(const long &i) {
   short err = 1;
   long j;
   // Pruefe Zelltyp der Nachbarn
   for (int n = 0; n < dim2; n++) {
      j = knot[i].near_n[n];
      if (j != -1) {
         if (knot[j].status == external) { err = 0; }
      }
   }
   return err;
}
short grid::no_border(const long &i) {
   /* Checks if there is a neighbour of knot i which belongs to the
    * border of the reaction volume (i.e. is in state external).
    * Returns 0, if a border point was found (1 otherwise).
    */
   short err = 1;
   long j;
   // Pruefe Zelltyp der Nachbarn
   for (int n = 0; n < dim2; n++) {
      j = knot[i].near_n[n];
      if (j != -1) {
         if (knot[j].status == external) { err = 0; }
      }
   }
   // #### This can be unified !!!
   if (err > 0) {
      // der Fall: kein ext als Nachbarn gefunden
      for (int n = 0; n < dim2; n++) {
         if (knot[i].near_n[n] == -1) { err = 0; }
      }
      // d.h. finde Randpunkte des ganzen Gitters
   }
   return err;
}
/*
 * short int lattice::sphere(long ind, int radius, states stat, long pos_ss) {
 * // ###
 * // This routine is not used for the moment !!! Actualize if it is activated.
 * // geht nur in 3 D
 * // Achtung! Die Einfuehrung von FDC_pos_ss ist hier nicht beruecksichtigt
 *
 * long center[dim];
 * long ctmp[dim];
 * //cout<<"dimsavepos="<<long(4.2*pow(radius,dim)+1)<<"\n";
 * dynarray<long> savepos(long(4.2*pow(radius,dim)+1),1,0);
 * long i;
 * long count=0;
 * short int error=0;
 * get_koord(ind,center);
 * ctmp=center;
 * // Test euclidean distance for all points
 * for (int delta0=-1*radius; delta0<=radius; delta0++) {
 *   ctmp[0]=center[0]+delta0;
 *   if (ctmp[0]>=0 && ctmp[0]<prodim) {
 *     for (int delta1=-1*radius; delta1<=radius; delta1++) {
 *       ctmp[1]=center[1]+delta1;
 *       if (ctmp[1]>=0 && ctmp[1]<prodim) {
 *         for (int delta2=-1*radius; delta2<=radius; delta2++) {
 *           ctmp[2]=center[2]+delta2;
 *           if (ctmp[2]>=0 && ctmp[2]<prodim) {
 *             i=Index(ctmp);
 *             // Check if all the points are nothing
 *             if (knot[i].status==nothing || stat==nothing) {
 *               savepos.write(count,i); // savepos[count]=i macht Fehler!
 *               count++;
 *             }
 *             // if not all points nothing do not write!
 *             else error=1;
 *           }
 *         }
 *       }
 *     }
 *   }
 * }
 * if (error==0) {
 *     // save all fragments on the lattice and point to new cell on list
 *     for (i=0; i<savepos.benutzt(); i++) {
 *     knot[savepos[i]].status=stat;
 *     knot[savepos[i]].listi=CB_list.benutzt();
 *     }
 *     // Write new cell on cell-list (prepared for CB only ##)
 *     cellCB newCB;
 *     newCB.pos_ss=pos_ss;
 *     newCB.index=ind;
 *     //newCB.center=ind;
 *     newCB.volume=savepos.benutzt();
 *     // save index-list in newCB.fragments also ##
 *     CB_list.add(newCB);
 * }
 * return error;
 * }
 */
