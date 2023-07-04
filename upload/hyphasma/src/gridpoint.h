#ifndef i_gridpoint
#define i_gridpoint

enum grid_states {
   nothing,object,object_border,external,N_grid_states
};

class gridpoint {
  public:
   // Konstruktoren
   gridpoint();
   gridpoint(long int&, grid_states&);
   ~gridpoint();

   // Index im Gitter
   long int index;
   // Index on guest grids
   long guest_grid_pointer;
   // coordinates of this point in the lattice
   long x[3];
   // Neighbors in the rank: 3 left neighbors then 3 right neighbors
   long int near_n[6];
   // Neighbors in the rank:
   // 3 2xleft, 3 1xleft 1x right, 3 1xright 1xleft, 3 2xright
   long int diag_n[12];
   // Art der Belegung
   grid_states status;
};

#endif
