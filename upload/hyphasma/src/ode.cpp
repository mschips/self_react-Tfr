// #include <iostream>
// #include <math>
// #include <stdio>
// #include "random.h"

#include "ode.h"
// #include <stdlib.h>

void ode::_Euler(double * y_n,
                 double * y_n1,
                 int N,
                 double t_n,
                 double _dt,
                 void (*getrhs)(double, double*, double*)) {
   //  cout<<"Run Euler ...";
   getrhs(t_n,y_n,dydt);
   for (int i = 0; i < N; ++i) {
      y_n1[i] = y_n[i] + _dt * dydt[i];
   }
}
void ode::_RungeKutta_2nd(double * y_n,
                          double * y_n1,
                          int N,
                          double t_n,
                          double _dt,
                          void (*getrhs)(double, double*, double*)) {
   double k[N];
   getrhs(t_n,y_n,k);
   double t_middle = t_n + _dt / 0.5;
   for (int i = 0; i < N; ++i) {
      k[i] = 0.5 * _dt * k[i] + y_n[i];
   }
   getrhs(t_middle,k,dydt);
   for (int i = 0; i < N; ++i) {
      y_n1[i] = y_n[i] + _dt * dydt[i];
   }
}
void ode::_RungeKutta_4th(double * y_n,
                          double * y_n1,
                          int N,
                          double t_n,
                          double _dt,
                          void (*getrhs)(double, double*, double*)) {
   double k[N];
   double k1[N];
   double k2[N];
   double k3[N];
   double k4[N];
   getrhs(t_n,y_n,k);
   for (int i = 0; i < N; ++i) {
      k1[i] = _dt * k[i];
   }

   double t_x = t_n + _dt / 0.5;
   for (int i = 0; i < N; ++i) {
      k[i] = 0.5 * k1[i] + y_n[i];
   }
   getrhs(t_x,k,k2);
   for (int i = 0; i < N; ++i) {
      k2[i] *= _dt;
   }

   for (int i = 0; i < N; ++i) {
      k[i] = 0.5 * k2[i] + y_n[i];
   }
   getrhs(t_x,k,k3);
   for (int i = 0; i < N; ++i) {
      k3[i] *= _dt;
   }

   t_x = t_n + _dt;
   for (int i = 0; i < N; ++i) {
      k[i] = k3[i] + y_n[i];
   }
   getrhs(t_x,k,k4);
   for (int i = 0; i < N; ++i) {
      k4[i] *= _dt;
   }

   for (int i = 0; i < N; ++i) {
      y_n1[i] = y_n[i] + k1[i] / 6. + k2[i] / 3. + k3[i] / 3. + k4[i] / 6.;
   }
}
void ode::run_step(double * y1,
                   double * y2,
                   int N,
                   double t_n,
                   double _dt,
                   void (*callrhs)(double, double*, double*),
                   const ode_method &num_method) {
   if (num_method == Euler) { _Euler(y1, y2, N, t_n, _dt, callrhs); }
   if (num_method == RungeKutta_2nd) { _RungeKutta_2nd(y1, y2, N, t_n, _dt, callrhs); }
   if (num_method == RungeKutta_4th) { _RungeKutta_4th(y1, y2, N, t_n, _dt, callrhs); }
}
void ode::ode_step(double * y_n,
                   double * y_n1,
                   int N,
                   double &t_n,
                   double _dt,
                   double max_error,
                   void (*callrhs)(double, double*, double*),
                   const ode_method &num_method) {
   double y_1[N];
   double y_2[N];
   short redo = 1;
   int n_steps = 1;
   while (redo == 1) {
      redo = 0;
      run_step(y_n, y_2,  N, t_n,        _dt,    callrhs, num_method);
      run_step(y_n, y_1,  N, t_n,        _dt / 2., callrhs, num_method);
      run_step(y_1, y_n1, N, t_n + _dt / 2., _dt / 2., callrhs, num_method);
      int i = 0;
      while (i < N && redo == 0) {
         if ((y_n1[i] - y_2[i] > max_error) || (y_2[i] - y_n1[i] > max_error)) {
            _dt /= 2.;
            cout << "i=" << i << ": y_n1[" << i << "]=" << y_n1[i] << " and y_2[" << i
                 << "]=" << y_2[i] << "\n"
                 << "    ... the difference is " << y_n1[i] - y_2[i] << "\n"
                 << "    ... which is bigger than max_error " << max_error << "\n"
                 << "    ==> time step reduced to " << _dt << " at t=" << t_n << "\n";
            n_steps *= 2;
            redo = 1;
         }
         ++i;
      }
   }
   // Save the new values as old ones:
   t_n += _dt;
   for (int i = 0; i < N; ++i) {
      y_n[i] = y_n1[i];
   }

   // Go through the additional sub-time-steps if necessary:
   for (int k = 1; k < n_steps; ++k) {
      ode_step(y_n, y_n1, N, t_n, _dt, max_error, callrhs, num_method);
   }
}
