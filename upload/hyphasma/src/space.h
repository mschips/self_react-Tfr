#ifndef i_space
#define i_space
#include "setparam.h"
#include "grid.h"

enum states {
   // state of knots, attribution to cell types
   nocell,
   CB,
   CC,
   FDC,
   TC,
    TFR, //MSchips --> check whether it is fine to add it here or should be at the end of enum //MMH2Marta: should be fine
   out,
   blast1,
   blast2,
   BETA,
   N_cells
};

// states operator++(const states& a);

class spacepoint {
  public:
   // Konstruktoren
   spacepoint();
   spacepoint(long int &i, states &c, long &listindex, long &FDClistindex);
   ~spacepoint();

   // Art der Belegung
   states cell;
   // Index on cell list
   long listi;
   // Index of possibly underlying FDCs (its center)
   long FDClisti;
};

class space: public grid {
  public:
   space();
   space(Parameter &par, ofstream &ana);
   ~space();

   spacepoint * cellknot;

   void set_knot(const long &i, const states &s, const long &li);
   void clear_knot(const long &i);
   long get_random_empty_position();
   
   // Zones of GCs
   long zone_separator;

   // Objekt tests
   short object_border(const long &i);
   long next_border(const long &i);

   // Routines for cell-consistencies ### better in cell classes
   short self(const long &i, const long &li, const states &status);
   short get_n_self(const long &i, const long &li,
                    const states &status, const long excluded);
   long end_of_tail(const long &start, long back,
                    const long &object, const states &status);
   long decompose(const long &i, const long &li, const states &status);

  private:
   // Objekt tests
   long check_neighbors_for_border(const long i, const long &ref, double &min_distance,
                                   dynarray<long> &done);
   // calculate zone separation index for GCs:
   long get_zone_separation_index();

};

#endif
