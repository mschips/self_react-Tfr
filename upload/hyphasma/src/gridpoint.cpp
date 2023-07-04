#include "gridpoint.h"

gridpoint::gridpoint() {
   index = 0;
   guest_grid_pointer = -1;
   short int i;
   for (i = 0; i < 3; i++) {
      x[i] = 0;
   }
   for (i = 0; i < 6; i++) {
      near_n[i] = 0;
   }
   for (i = 0; i < 12; i++) {
      diag_n[i] = 0;
   }
   status = nothing;
}
gridpoint::gridpoint(long int &i, grid_states &s) {
   index = i;
   guest_grid_pointer = -1;
   short int ii;
   for (ii = 0; ii < 3; ii++) {
      x[ii] = 0;
   }
   for (ii = 0; ii < 6; ii++) {
      near_n[ii] = 0;
   }
   for (ii = 0; ii < 12; ii++) {
      diag_n[ii] = 0;
   }
   status = s;
}
gridpoint::~gridpoint() { }
char operator ==(const gridpoint &a, const gridpoint &b) {
   if ((a.index != b.index) || (a.status != b.status)) { return 0; }
   return 1;
}
char operator !=(const gridpoint &a, const gridpoint &b) {
   return !((a == b) == 1);
}
