#ifndef i_odelist
#define i_odelist

#include <iostream>
#include "cellthis.h"
using namespace std;

// enum odelist_method{Euler,RungeKutta_2nd,RungeKutta_4th,all_ode_methods};

class odelist {
  public:
   odelist() { dydt = new double[10]; }
   odelist(const int &N) { dydt = new double[N]; }
   ~odelist() { }
   //  ~odelist() { delete[] dydt; };
   void ode_step(double * y_n,
                 double * y_n1,
                 dynarray<cellbeta> &list,
                 space &l,
                 int N,
                 double &t_n,
                 double _dt,
                 double max_error,
                 void (*callrhs)(double, double*, double*, dynarray<cellbeta>&, space &l),
                 const ode_method &num_method);
  private:
   double * dydt;
   void _Euler(double * y_n,
               double * y_n1,
               dynarray<cellbeta> &list,
               space &l,
               int N,
               double t_n,
               double _dt,
               void (*callrhs)(double, double*, double*, dynarray<cellbeta>&, space &l));
   void _RungeKutta_2nd(double * y_n,
                        double * y_n1,
                        dynarray<cellbeta> &list,
                        space &l,
                        int N,
                        double t_n,
                        double _dt,
                        void (*getrhs)(double, double*, double*, dynarray<cellbeta>&, space &l));
   void _RungeKutta_4th(double * y_n,
                        double * y_n1,
                        dynarray<cellbeta> &list,
                        space &l,
                        int N,
                        double t_n,
                        double _dt,
                        void (*getrhs)(double, double*, double*, dynarray<cellbeta>&, space &l));
   void run_step(double * y1,
                 double * y2,
                 dynarray<cellbeta> &list,
                 space &l,
                 int N,
                 double t_n,
                 double _dt,
                 void (*callrhs)(double, double*, double*, dynarray<cellbeta>&, space &l),
                 const ode_method &num_method);
};

#endif
