#ifndef i_ode
#define i_ode

#include <iostream>
using namespace std;

enum ode_method {
   Euler,RungeKutta_2nd,RungeKutta_4th,all_ode_methods
};

class ode {
  public:
   ode() { dydt = new double[10]; }
   ode(const int &N) { dydt = new double[N]; }
   ~ode() { }
   //  ~ode() { delete[] dydt; };
   void ode_step(double * y_n,
                 double * y_n1,
                 int N,
                 double &t_n,
                 double _dt,
                 double max_error,
                 void (*callrhs)(double, double*, double*),
                 const ode_method &num_method);

  private:
   double * dydt;
   void _Euler(double * y_n,
               double * y_n1,
               int N,
               double t_n,
               double _dt,
               void (*callrhs)(double, double*, double*));
   void _RungeKutta_2nd(double * y_n,
                        double * y_n1,
                        int N,
                        double t_n,
                        double _dt,
                        void (*getrhs)(double, double*, double*));
   void _RungeKutta_4th(double * y_n,
                        double * y_n1,
                        int N,
                        double t_n,
                        double _dt,
                        void (*getrhs)(double, double*, double*));
   void run_step(double * y1,
                 double * y2,
                 int N,
                 double t_n,
                 double _dt,
                 void (*callrhs)(double, double*, double*),
                 const ode_method &num_method);
};

#endif
