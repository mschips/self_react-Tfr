#ifndef i_random
#define i_random

#include <stdlib.h>
#include "dynarray.h"
#include <fstream>
using namespace std;

// const double RAND_MAXp1=double(RAND_MAX)+double(1);
// Liefert Pseudozufalls-int 0 <= x < bis
int irandom(int bis);
// Liefert Pseudozufalls-long int 0 <= x <= bis
long int lrandom(long int bis);
// Liefert Pseudozufalls-double 0 <= x < bis
double drandom(double bis);
double drandom(); // bis=1. vorausgesetzt
// Liefert eine Sequenz von long int Zahlen zwischen [0,max[
// in zufaelliger Reihenfolge
void random_sequence(dynarray<long> &m, long int max);
void random2_sequence(long * m, long int max);

double inverse_erf(double x);
double get_sample_from_normal(const double &mean, const double &width);
double get_positive_sample_from_normal(const double &mean, const double &width);

class pre_randomize {
  public:
   pre_randomize();
   pre_randomize(short dim2);
   ~pre_randomize();
   void get_rand_set(dynarray<short> &xx);
   void get_set(long n, dynarray<short> &xx);
   short* get_rand_set();
   short* get_set(long n);
  private:
   long permut;
   short sl;
   short * * field;
   void add_number(short &n);
   ofstream out;
};

class gauss_randomize {
  public:
   gauss_randomize();
   gauss_randomize(short dataset);
   ~gauss_randomize();
   double get_distribution_value();
  private:
   void gauss_initialize(double, double, int, double, double);
   double gauss(double&,double&,double&);
   void cyster_initialize();
   double cyster07angle_wt(int angle);
   double * field;
   int arraydim;
   static constexpr double Ny = 10;
};

#endif
