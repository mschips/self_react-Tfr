#include "signals.h"
#include <math.h>
#include <string.h>

spoint::spoint() {
   for (short i = 0; i < signals; i++) {
      signal[i] = 0.0;
      signal_new[i] = 0.0;
      signal_tmp[i] = 0.0;
   }
}
spoint::~spoint() { }
void spoint::actualize() {
   for (short i = 0; i < signals; i++) {
      signal[i] = signal_new[i];
   }
   // wenn hier nur die benutzten aktualisiert werden ginge das schneller ####
}
void spoint::actualize(const signal_molecule &a) {
   signal[a] = signal_new[a];
}
char operator ==(const spoint &a, const spoint &b) {
   char back = 1;
   for (short i = 0; i < signals; i++) {
      if (a.signal[i] != b.signal[i]) { back = 0; }
   }
   return back;
}
char operator !=(const spoint &a, const spoint &b) {
   return !((a == b) == 1);
}
// ==================================================================
// ### flag Neumann oder Moore einfuehren
// ==================================================================
// ignore_object==0 fuer Euler umprogrammieren und fuer ADI neu programmieren!!!
// ==================================================================

bool sigs::TEST_MODE = false;
suffix sigs::TEST_SUF = "0000";
double sigs::TEST_SOURCE = 0;
double sigs::TEST_AMPLITUDE = 0.;
double sigs::vonNEUMANNinwards = -1.0;
double sigs::CONST_DYN_GLU_FIELD = 10.;

sigs::sigs()
   : grid() {
   sigsknot = new spoint[pointnumber];
   maxprodim = 0;
   for (short i = 0; i < dim; i++) {
      if (prodimvec[i] > maxprodim) {
         maxprodim = prodimvec[i];
      }
   }
   for (short i = 0; i < signals; i++) {
      dirichlet[i] = true;
   }
   gam = new double[maxprodim + 2];
}
sigs::sigs(space &xyz, Parameter &par, ofstream &ana)
   : grid(molecules,
          par.Value.system,
          par.Value.DimSpace,
          par.Value.dx_signal,
          par.Value.GC_radius,
          par.Value.gridsize,
          par.Value.vol_shape,
          par.Value.obstacles,
          par.Value.wall_level,
          par.Value.collagen_density,
          par.Value.collagen_cluster,
          par.Value.wall_width,
          par.Value.slit_number,
          par.Value.slit_width,
          par.Value.use_specific_turning_angles,
          ana) {
   //  :  grid(molecules,par,ana)
   // Made new grid for signals:

   // ++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++
   // For the TESTMODE:
   // Switch TEST_MODE on or off:
   // TEST_MODE=true;
   // AMPLITUDE of the initial configuration:
   // TEST_AMPLITUDE=50.;  // molecules per sqrt(unit volume)
   // Additional SOURCE for the reference solution only:
   // Use the same value as calculated in cellthis.
   // Note that the lattice constant of the space lattice has to be used (as in cellthis.C):
   // TEST_SOURCE=par.Value.mkCXCL13*par.Value.dx*par.Value.dx*par.Value.dx*1.e-15*par.N_A;
   // Alternatively:
   // The value per FDC and for one time step is given as a value in ana_ini.out
   // and is to be divided by dt to get the production per hour.
   // This allows to put a value by hand!
   // TEST_SOURCE=75275.6;
   if (TEST_MODE) {
      cout << "TEST_SOURCE=" << TEST_SOURCE << " is used in the signals.C constructor!\n";
      ana << "TEST_SOURCE=" << TEST_SOURCE << " is used in the signals.C constructor!\n";
   }
   // ++++++++++++++++ end OPTION +++++++++++++++++++++++++++++

   // Define the field of sigspoint
   sigsknot = new spoint[pointnumber];

   // Project space-grid on signal-grid
   xyz.guest_grid_project(dx,prodimvec);

   // ==========================
   ana << "Initialize signals ...\n";
   // ==========================
   if (par.Value.mksignal < 1.e-08) {
      signal_use[sig_differ2CC] = 0;
      if (par.Value.mksignal < -0.5) { signal_use[sig_differ2CC] = -1; }
      if (signal_use[sig_differ2CC] == 0) {
         ana << "CB rate differentiation independent of FDC.\n";
      } else {
         ana << "CB differentiate in response to FDC-contact (mk_signal<-0.5).\n";
      }
   } else {
      signal_use[sig_differ2CC] = 1;
      ana << "CB differentiate in response to FDC-derived signal.\n";
   }

   if ((par.Value.mkCXCL12 <= 0.) && (par.Value.bound_CXCL12 <= 0.)) {
      signal_use[CXCL12] = 0;
   } else { signal_use[CXCL12] = 1; }
   if ((par.Value.mkCXCL13 <= 0.) && (par.Value.bound_CXCL13 <= 0.)) {
      signal_use[CXCL13] = 0;
   } else { signal_use[CXCL13] = 1; }
   if (par.Value.mk_ab <= 0.) { signal_use[antibody] = 0; } else { signal_use[antibody] = 1; }
   if (par.Value.bound_ag <= 0.) { signal_use[antigen] = 0; } else { signal_use[antigen] = 1; }
   if ((par.Value.mk_SEMA4D <= 0.) && (par.Value.bound_SEMA4D <= 0.)) {
      signal_use[SEMA4D] = 0;
   } else { signal_use[SEMA4D] = 1; }
   if (par.Value.bound_glucose <= 0.) { signal_use[glucose] = 0; } else { signal_use[glucose] = 1; }
   if (par.Value.bound_oxygen <= 0.) { signal_use[oxygen] = 0; } else { signal_use[oxygen] = 1; }

   // Diffusion of signals:
   D[sig_differ2CC] = par.Value.D_differ2CC;
   D[CXCL12] = par.Value.D_CXCL12;
   D[CXCL13] = par.Value.D_CXCL13;
   D[antibody] = par.Value.D_antibody;
   D[antigen] = par.Value.D_antigen;
   D[SEMA4D] = par.Value.D_SEMA4D;
   D[glucose] = par.Value.D_glucose * 60.;  // transform in per hour
   D[oxygen] = par.Value.D_oxygen * 60.;  // transform in per hour
   // #### D[nutrient] in H2O is not used here!!!

   // dt=par.Value.deltat;
   for (short b = 0; b < signals; b++) {
      if (signal_use[b] == 1) {
         p_diffuse_signal[b] = D[b] * par.Value.deltat
                               / (par.Value.dx_signal * par.Value.dx_signal);
         if (par.Value.signal_mode != 2) {
            // not ADI
            diff_steps[b] = (int (SIGNAL_TIME_RESOLUTION * dim2 * D[b] * par.Value.deltat
                                  / (par.Value.dx_signal * par.Value.dx_signal)) + 1);
         } else {
            // case ADI:
            // determine internal timestep delta_t such that alpha are approx 0.5
            double delta_t = par.Value.deltat
                             / (2.0 * double (SIGNAL_TIME_RESOLUTION) * p_diffuse_signal[b]);
            /* Typischerweise steht hier nur ein Faktor 2.0. Mit SIGNAL_TIME_RESOLUTION=2 oder 4
             * wird die Zeitaufloesung verfeinert.
             */
            diff_steps[b] = static_cast<unsigned int> (par.Value.deltat / delta_t);
            if (diff_steps[b] < 1) { diff_steps[b] = 1; }
            /* diff_steps erlaubt es mit einer groben aeusseren Zeitaufloesung die schnelle
             * Signal-Diffusion dennoch mit einer adaequaten Zeitaufloesung zu berechnen.
             */
         }
         // cout<<"p_dif="<<p_diffuse_signal<<"\n";
         // cout<<"correct it with diffsteps="<<diff_steps<<"\n";
         // dt/=double(diff_steps);
         p_diffuse_signal[b] /= double (diff_steps[b]);
         // cout<<"p_dif="<<p_diffuse_signal<<"\n";
         ana << "  signal " << b << " diffuse = " << p_diffuse_signal[b]
             << " using " << diff_steps[b] << " sub-timesteps.\n";
      }
   }

   if (signal_use[glucose] == 1) {
      cout << "diffsteps[glucose]=" << diff_steps[glucose] << "\n";
   }
   CONST_DYN_GLU_FIELD = par.Value.const_dynamic_glucose_field;

   if (signal_use[oxygen] == 1) {
      cout << "diffsteps[oxygen]=" << diff_steps[oxygen] << "\n";
   }
   // memory for the help-field in ADI matrix inversion
   maxprodim = 0;
   for (short i = 0; i < dim; i++) {
      if (prodimvec[i] > maxprodim) {
         maxprodim = prodimvec[i];
      }
   }
   gam = new double[maxprodim + 2];

   // mode of diffusion
   diffusion_mode = par.Value.signal_mode;
   ana << "  signal diffusion mode = " << diffusion_mode << "\n";
   // do signals ignore the presence of objects?
   ignore_objects = par.Value.objects_transparent;
   ana << "  object transparency for signals = " << ignore_objects << "\n";
   if (ignore_objects == 0) {
      /* ### Objects have to be ignored for the moment.
       *     Make routines for objects treated as obstacles
       * with different diffusion algorithms,
       * with variable resolution in cell and signal grid.
       */
      cout << "Objects have to be transparent for signals!\n";
      exit(1);
   }
   // Zaehler setzen
   for (unsigned short n = 0; n < signals; n++) {
      signal_total[n] = 0.0;
      signal_diffout[n] = 0.0;
      signal_produced[n] = 0.0;
      signal_used[n] = 0.0;
   }

   ana << "... done.\n";

   // =====================================
   ana << "Set boundary for signals ...";
   // =====================================

   // Boundary conditions:

   // Set Dirichlet boundaries as standard
   for (short i = 0; i < signals; i++) {
      dirichlet[i] = true;
   }
   char signame[30] = "";
   // ++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++
   // Add deviations if needed
   // dirichlet[CXCL13]=false;
   // Only one boundary value is given in case of von Neumann boundary conditions.
   // If it is not no-flow but constant flow consider the following:
   // For vonNEUMANNinwards==1.0  a positive value for the derivative at the boundary
   // corresponds to an inflow at the left and an outflow at the right border.
   // For vonNEUMANNinwards==-1.0 a positive value for the derivative at the boundary
   // corresponds to an inflow at every border.
   vonNEUMANNinwards = -1.0;
   // ++++++++++++++++ end OPTION +++++++++++++++++++++++++++++
   for (short i = 0; i < signals; i++) {
      if ((dirichlet[i] == false) && (signal_use[i] == 1)) {
         get_signal_name(CXCL13,signame);
         cout << "WARNING!!! Use von Neumann boundary for " << signame
              << ". End of WARNING!!!\n";
         ana << "WARNING!!! Use von Neumann boundary for " << signame << ". End of WARNING!!!\n";
      }
   }

   boundary.signal[sig_differ2CC] = par.Value.bound_differ2CC;
   // Im folgenden wird auf die Einheit [Zahl der Molekuele] umgeformt:
   /* Der Faktor "mk_concentration" bringt die Einheit Molekuele auf Mol=mol/dm^3.
    * Dabei ist dx in micrometern angeben.
    * Falls 2D simuliert wird sind zwei dx vom Signalgitter und ein dx vom Zellgitter zu waehlen,
    * sonst alle 3 dx vom Signalgitter (in 3D).
    */
   if (dim
       == 2) {
      mk_concentration = 1. / (pow(par.Value.dx_signal,2.) * par.Value.dx * 1.e-15 * par.N_A);
   } else { mk_concentration = 1. / (pow(par.Value.dx_signal,3.) * 1.e-15 * par.N_A); }

   boundary.signal[CXCL12] = par.Value.bound_CXCL12 / mk_concentration;
   boundary.signal[CXCL13] = par.Value.bound_CXCL13 / mk_concentration;
   boundary.signal[antibody] = par.Value.bound_ab / mk_concentration;
   boundary.signal[antigen] = par.Value.bound_ag / mk_concentration;
   boundary.signal[SEMA4D] = par.Value.bound_SEMA4D / mk_concentration;
   boundary.signal[glucose] = par.Value.bound_glucose / mk_concentration;
   boundary.signal[oxygen] = par.Value.bound_oxygen / mk_concentration;
   cout << "bound_glu=" << boundary.signal[glucose] << ", bound_ox=" << boundary.signal[oxygen]
        << "\n";

   for (long n = 0; n < pointnumber; n++) {
      // +++ OPTION For nutrients in vitro (make as external option) OPTION +++
      if (signal_use[glucose] == 1) {
         sigsknot[n].signal_new[glucose] = boundary.signal[glucose];
         sigsknot[n].signal[glucose] = boundary.signal[glucose];
      }
      if (signal_use[oxygen] == 1) {
         sigsknot[n].signal_new[oxygen] = boundary.signal[oxygen];
         sigsknot[n].signal[oxygen] = boundary.signal[oxygen];
      }
      // end OPTION
   }

   // For diffusion: set the boundary values on "external" points:
   for (long n = 0; n < pointnumber; n++) {
      if (knot[n].status == external) {
         for (short a = 0; a < signals; a++) {
            sigsknot[n].signal[a] = boundary.signal[a];
            sigsknot[n].signal_new[a] = boundary.signal[a];
         }
      }
   }

   // If no reaction-diffusion-system but a fixed distribution is used
   // load from a file:
   if (int (signals) > MAXDIMSMALL) {
      cout << "ERROR!\n"
           << "Please increase MAXDIMSMALL in setparam.h to capture all signals!\n";
      exit(1);
   }
   for (short bb = 0; bb < short (signals); bb++) {
      fix_signal[bb] = 0;
      if (par.Value.fix_signals[bb]) {
         fix_signal[bb] = 1;
         if (signal_use[bb] == 1) {
            load_signal_file(signal_molecule(bb));
         } else { signal_use[bb] = 1; }
         /* If fix_signal==1 and signal_use==0, signal_use==1 is set!
          * The signal is not read from a file but calculated according to some algorithm.
          * The explicit initialisation has to be done in the calling routine.
          * The signal is updated by recalculating the distribution (called from outside).
          * Diffusion is switched off and boundaries are ignored.
          */
      }
   }

   if (par.Value.fix_glucose_gradient) {
      mk_grad1d_signal(glucose,par.Value.fix_glucose_gradient_min,
                       par.Value.fix_glucose_gradient_max,0);
   }

   siglog.open("siglog.out");
   siglog << "! Signal types: ";
   short n_used_signals = 0;
   for (short bb = 0; bb < short (signals); bb++) {
      if (signal_use[bb] == 1) {
         char signalname[30] = "";
         get_signal_name(bb,signalname);
         siglog << signalname;
         if (bb < short (signals) - 1) { siglog << ", "; }
         ++n_used_signals;
      }
   }
   siglog << "\n";
   siglog << "! time & " << n_used_signals << " signals";
   siglog << ": [current & diffused out & used & produced since last write]\n";
   write_siglog(par.Value.tmin);

   // If the class is run in TESTMODE set the initial signal distribution
   if (TEST_MODE) {
      cout << "WARNING!!! Class signals runs in TESTMODE. End of WARNING!!!\n";
      ana << "WARNING!!! Class signals runs in TESTMODE. End of WARNING!!!\n";
      for (short bb = 0; bb < signals; bb++) {
         if (signal_use[bb] == 1) { set_initial_condition(signal_molecule(bb)); }
      }
   }

   ana << " done\n";
}
sigs::~sigs() {
   // cout<<"in ~sigs()\n";
   // delete[] nn;
   siglog.close();
   delete[] sigsknot;
   delete gam;
}
void sigs::signal_put(const long &i, const signal_molecule &sig_type, const double &zahl) {
   if (fix_signal[sig_type] == 0) { sigsknot[i].signal_new[sig_type] += zahl; }
}
void sigs::signal_set(const long &i, const signal_molecule &sig_type, const double &zahl) {
   sigsknot[i].signal_new[sig_type] = zahl;
}
void sigs::boundary_set(const long &i, const signal_molecule &sig_type, const double &zahl) {
   boundary.signal[sig_type] = zahl;
}
void sigs::load_signal_file(const signal_molecule sig_type) {
   signal_use[sig_type] = 1;
   char fname[30] = "";
   get_signal_name(sig_type,fname);
   // Chose the file according to the grid properties
   if ((sig_type == CXCL12) || (sig_type == CXCL13)) {
      if (dim == 3) { strcat(fname,"_3d"); }
      else { strcat(fname,"_2d"); }
      if ((dx < 5.0 + 1.e-08) && (dx > 5.0 - 1.e-08)) { strcat(fname,"_5micron"); }
      if ((dx < 6.0 + 1.e-08) && (dx > 6.0 - 1.e-08)) { strcat(fname,"_6micron"); }
      if ((dx < 6.5 + 1.e-08) && (dx > 6.5 - 1.e-08)) { strcat(fname,"_6_5micron"); }
      if ((dx < 7.5 + 1.e-08) && (dx > 7.5 - 1.e-08)) { strcat(fname,"_7_5micron"); }
      if ((dx < 10.0 + 1.e-08) && (dx > 10.0 - 1.e-08)) { strcat(fname,"_10micron"); }
   }
   strcat(fname,".sig");
   cout << "Load " << fname << " from file ... ";
   ifstream ff(fname);
   string tmp = "";
   for (short a = 0; a < N_PRE_WORDS; a++) {
      ff >> tmp;
   }
   long d1,d2,d3,d3a;
   double d4,d5,d6,d7,d8,d9,d9a;
   for (long i = 0; i < pointnumber; i++) {
      ff >> d1 >> d2 >> d3;
      if (dim == 3) { ff >> d3a; }
      ff >> d4 >> d5 >> d6 >> d7 >> d8 >> d9;
      if (dim == 3) { ff >> d9a; }
      if (ff.fail()) { exit(1); } else {
         sigsknot[i].signal[sig_type] = d4;
         sigsknot[i].signal_new[sig_type] = d4;
      }
   }
   ff.close();
   cout << " done.\n";
}
void sigs::mk_const_signal(const signal_molecule sig_type, double &value) {
   for (long i = 0; i < pointnumber; i++) {
      sigsknot[i].signal[sig_type] = value;
      sigsknot[i].signal_new[sig_type] = value;
   }
   boundary.signal[sig_type] = value;
}
double sigs::get_2sigmoidal(double &t, double &t_a, double &t_b,
                            double &rest, double &factor,
                            double &kappa_a, double &kappa_b) {
   double f = 1. / ((1 + exp((t_a - t) / kappa_a)) * (1 + exp((t - t_b) / kappa_b)));  // weight
                                                                                       // function
   double g = rest * (1 - f) + factor * rest * f;  // mMol
   return g;
}
double sigs::get_const_signal_value(double &t, double &ampl) {
   // ++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++
   double t_a = 3.;  // 20 s, time of external potassium increase
   double t_b = 120.;  // 110 s, time of external potassium decrease
   // double t_b=53.; // 110 s, time of external potassium decrease
   // factor of glucose increase with respect to the resting level:
   double factor = CONST_DYN_GLU_FIELD;

   // EITHER use the step-function ...
   // if (t>t_a && t<t_b) return factor*p.glu_0; else return p.glu_0;

   // ... OR use a linear increase in time ...
   /*
    * double factor_0=6.0;
    * factor=7.0; // factor of glucose increase with respect to the resting level
    * if (t>t_a && t<t_b) return factor_0+(t-t_a)*(factor-factor_0)/(t_b-t_a);
    * else return p.glu_0;
    */

   // ... OR use the smooth double sigmoidal function ...
   ///*
   double kappa_a = 0.10;  // s, width of external potassium increase
   double kappa_b = 0.10;  // s, width of external potassium decrease
   return get_2sigmoidal(t,t_a,t_b,ampl,factor,kappa_a,kappa_b);
   // cerr<< "t="<<t<<"s;  glucose="<<value<<" mMol;  ";
   // */
   // +++++++++++++++++++++ end OPTION ++++++++++++++++++++++++++
}
double sigs::get_const_signal_value(double &t, double &ampl, long * wo) {
   // +++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++
   // const long BLOCK0=4;
   double value = get_const_signal_value(t,ampl);
   // wo contains the coordinates on the cellular space grid
   //
   // WARNING: This is incompatible for cell- and signal-grid with different lattice constant!
   //
   //  if (wo[0]>=BLOCK0) value=ampl+0.8*(value-ampl);
   // +++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++
   return value;
}
void sigs::mk_const_signal(const signal_molecule sig_type, double t, double ampl) {
   if (CONST_DYN_GLU_FIELD > 0) {
      // cerr<<"I am in mk_const_signal(...)!\n";
      double value = get_const_signal_value(t,ampl);

      // +++++++++++++++++++++++ OPTION ++++++++++++++++++++++++
      mk_const_signal(sig_type,value);

      // If single positions shall be changed do this here
      // 1) inactivate the previous mk_const_signal()
      // 2) replace it by the one here
      // 3) Chose the position where to change the value
      /*
       * double reduced_value=ampl+0.8*(value-ampl);
       * mk_const_signal(sig_type,reduced_value);
       * long i,j,ii;
       *
       * for (j=0; j<prodimvec[1]; j++)
       * for (i=0; i<4; i++) {
       * ii=Index(i,j);
       * sigsknot[ii].signal[sig_type]=value;
       * sigsknot[ii].signal_new[sig_type]=value;
       * }
       */
      /*
       * i=1; j=0; ii=Index(i,j);
       * //i=4; j=1; ii=Index(i,j);
       * //cerr<<"Index(i,j)="<<ii<<"\n";
       * sigsknot[ii].signal[sig_type]=value;
       * sigsknot[ii].signal_new[sig_type]=value;
       */
      // +++++++++++++++++++++++ end OPTION ++++++++++++++++++++++++
   }
}
void sigs::mk_grad1d_signal(const signal_molecule sig_type,
                            double min, double max, short direction) {
   double value = 0;
   if (dim == 3) {
      for (long i = 0; i < prodimvec[0]; i++) {
         for (long j = 0; j < prodimvec[1]; j++) {
            for (long k = 0; k < prodimvec[2]; k++) {
               long ind = Index(i,j,k);
               if (direction
                   == 0) {
                  value = min + double (i) * (max - min) / double (prodimvec[0] - 1);
               } else if (direction == 1) {
                  value = min + double (j) * (max - min) / double (prodimvec[1] - 1);
               } else if (direction == 2) {
                  value = min + double (k) * (max - min) / double (prodimvec[2] - 1);
               } else {
                  cerr << "Wrong value of direction=" << direction
                       << " in sigs::mk_grad1d_signal(...)!\n";
                  exit(1);
               }
               sigsknot[ind].signal[sig_type] = value;
               sigsknot[ind].signal_new[sig_type] = value;
            }
         }
      }
   }
   if (dim == 2) {
      for (long i = 0; i < prodimvec[0]; i++) {
         for (long j = 0; j < prodimvec[1]; j++) {
            long ind = Index(i,j);
            if (direction
                == 0) {
               value = min + double (i) * (max - min) / double (prodimvec[0] - 1);
            } else if (direction == 1) {
               value = min + double (j) * (max - min) / double (prodimvec[1] - 1);
            } else {
               cerr << "Wrong value of direction=" << direction
                    << " in sigs::mk_grad1d_signal(...)!\n";
               exit(1);
            }
            sigsknot[ind].signal[sig_type] = value;
            sigsknot[ind].signal_new[sig_type] = value;
         }
      }
   }
}
short sigs::signal_get(const long &i, const signal_molecule &sig_type, const double &howmuch) {
   if (sigsknot[i].signal_new[sig_type] > howmuch) {
      sigsknot[i].signal_new[sig_type] -= howmuch;
      return 0;
   } else {
      sigsknot[i].signal_new[sig_type] = 0.;
      // cout<<"signal drops to zero!!\n";
      return 1;   // return 1 if signal is completely used, 0 otherwise.
   }
}
double sigs::signal_anregung(const long &n, const signal_molecule &sig) {
   /* 2/2008: It seems this routine is not used!
    * prodim is not updated exactly to prodimvec[] #####
    * just maxprodim is used.
    * the routine makes only sense if prodimvec[] is equal in all dimensions!
    */
   /*
    * long k[dim];
    * get_koord(n,k);
    */
   double k_sum = 0.0;
   //  for (short i=0; i<dim; i++) k_sum+=k[i];
   for (short i = 0; i < dim; i++) {
      k_sum += knot[n].x[i];
   }
   double R = ((maxprodim - 1) * dx);
   return double (dim) * pi * pi * D[sig] * sin(pi * k_sum / R);
}
void sigs::set_initial_condition(const signal_molecule &s) {
   if (dim == 3) {
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int j = 0; j < prodimvec[1]; j++) {
            for (int k = 0; k < prodimvec[2]; k++) {
               long m = Index(i,j,k);
               if (dirichlet[s]) {
                  sigsknot[m].signal[s] = sqrt(8. / (double (prodimvec[0] + 1)
                                                     * double (prodimvec[1] + 1)
                                                     * double (prodimvec[2] + 1)
                                                     * dx * dx * dx));
                  sigsknot[m].signal[s]
                     *= sin(pi * double (i + 1) / double (prodimvec[0] + 1));
                  sigsknot[m].signal[s]
                     *= sin(pi * double (j + 1) / double (prodimvec[1] + 1));
                  sigsknot[m].signal[s]
                     *= sin(pi * double (k + 1) / double (prodimvec[2] + 1));
                  sigsknot[m].signal[s] *= TEST_AMPLITUDE;
               } else {
                  // von Neumann: use a flat distribution
                  sigsknot[m].signal[s] = TEST_AMPLITUDE;
               }
               sigsknot[m].signal_new[s] = sigsknot[m].signal[s];
            }
         }
      }
   }
   if (dim == 2) {
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int j = 0; j < prodimvec[1]; j++) {
            long m = Index(i,j);
            if (dirichlet[s]) {
               sigsknot[m].signal[s] = sqrt(4. / (double (prodimvec[0] + 1)
                                                  * double (prodimvec[1] + 1)
                                                  * dx * dx));
               sigsknot[m].signal[s] *= sin(pi * double (i + 1) / double (prodimvec[0] + 1));
               sigsknot[m].signal[s] *= sin(pi * double (j + 1) / double (prodimvec[1] + 1));
               sigsknot[m].signal[s] *= TEST_AMPLITUDE;
            } else {
               // von Neumann: use a flat distribution
               sigsknot[m].signal[s] = TEST_AMPLITUDE;
            }
            sigsknot[m].signal_new[s] = sigsknot[m].signal[s];
         }
      }
   }
}
double sigs::reference(long * v, double &t, const signal_molecule &s) {
   if (dim == 3) {
      if (dirichlet[s]) {
         double lambda = pi * pi
                         * (1.0 / (dx * dx * double ((prodimvec[0] + 1) * (prodimvec[0] + 1)))
                            + 1.0
                            / (dx * dx * double ((prodimvec[1] + 1) * (prodimvec[1] + 1)))
                            + 1.0
                            / (dx * dx * double ((prodimvec[2] + 1) * (prodimvec[2] + 1))));
         double expo = exp(-1.0 * lambda * D[s] * t);
         double value = TEST_AMPLITUDE * expo;
         value += TEST_SOURCE * (1.0 - expo) / (lambda * D[s]);
         // cerr<<prodimvec[0]<<"\n";
         /*    cerr<<"la="<<lambda<<"  ex="<<expo<<"  AM="<<TEST_AMPLITUDE
          * <<"  SO="<<TEST_SOURCE<<"  VA="<<value<<"\n"; */
         value *= sqrt(8. / (double (prodimvec[0] + 1)
                             * double (prodimvec[1] + 1)
                             * double (prodimvec[2] + 1)
                             * dx * dx * dx));
         value *= sin(pi * double (v[0] + 1) / double (prodimvec[0] + 1));
         value *= sin(pi * double (v[1] + 1) / double (prodimvec[1] + 1));
         value *= sin(pi * double (v[2] + 1) / double (prodimvec[2] + 1));
         // cerr<<"value="<<value<<"\n";
         return value;
      } else {
         // von Neumann
         return TEST_AMPLITUDE;
      }
   }
   if (dim == 2) {
      if (dirichlet[s]) {
         double lambda = pi * pi
                         * (1.0 / (dx * dx * double ((prodimvec[0] + 1) * (prodimvec[0] + 1)))
                            + 1.0
                            / (dx * dx * double ((prodimvec[1] + 1) * (prodimvec[1] + 1))));
         double expo = exp(-1.0 * lambda * D[s] * t);
         double value = TEST_AMPLITUDE * expo;
         value += TEST_SOURCE * (1.0 - expo) / (lambda * D[s]);
         // cerr<<prodimvec[0]<<"\n";
         /*    cerr<<"la="<<lambda<<"  ex="<<expo<<"  AM="<<TEST_AMPLITUDE
          * <<"  SO="<<TEST_SOURCE<<"  VA="<<value<<"\n"; */
         value *= sqrt(4. / (double (prodimvec[0] + 1)
                             * double (prodimvec[1] + 1)
                             * dx * dx));
         value *= sin(pi * double (v[0] + 1) / double (prodimvec[0] + 1));
         value *= sin(pi * double (v[1] + 1) / double (prodimvec[1] + 1));
         // cerr<<"value="<<value<<"\n";
         return value;
      } else {
         // von Neumann
         return TEST_AMPLITUDE;
      }
   }
   return 0;
}
bool sigs::overcritical_signal(const long &index, const signal_molecule &s, const double &crit) {
   if (sigsknot[index].signal[s] > crit) { return true; }
   return false;
}
bool sigs::undercritical_signal(const long &index, const signal_molecule &s, const double &crit) {
   if (sigsknot[index].signal[s] < crit) { return true; }
   return false;
}
void sigs::get_gradient(const signal_molecule &s, const long &index, double * gradient) {
   for (unsigned short i = 0; i < dim; i++) {
      if (knot[index].near_n[dim + i] == -1) {
         gradient[i] = boundary.signal[s];
      } else { gradient[i] = sigsknot[knot[index].near_n[dim + i]].signal[s]; }
      if (knot[index].near_n[i] == -1) {
         gradient[i] -= boundary.signal[s];
      } else { gradient[i] -= sigsknot[knot[index].near_n[i]].signal[s]; }
      /* Note that it is not necessary to check for cell=external as at these points
       * the value of signal[s] is identical to the value in boundary.signal[s].
       * Therefore it is sufficient to exclude non-existing points from the
       * calculation.
       */
   }
   /* Note that it is not necessary to divide by the lattice constant
    * because the gradient is used as difference only in the calling routine.
    * !!! Take care of this, if other routines use this get_gradient(...)!!!
    */
}
double sigs::get_signal_total(const signal_molecule &sig_type) {
   // This routine is used for ADI, where signal_total is not actualized at every time step!
   signal_total[sig_type] = 0.0;
   for (long r = 0; r < pointnumber; r++) {
      if (knot[r].status != external) {
         signal_total[sig_type] += sigsknot[r].signal_new[sig_type];
      }
   }
   //  cout<<"total["<<sig_type<<"]="<<signal_total[sig_type]<<"\n";
   return signal_total[sig_type];
}
void sigs::signal_diffuse(space &xyz) {
   /* Situation: signal[] enthaelt die alte Signalverteilung
    *            signal_new[] die neue inklusive neuer Produktion
    *    signal_tmp[] wird hier die Produktion pro internen
    *                 Diffusions-Zeitschritt erhalten.
    * Um die Signal-Diffusion zu starten ist die Erzeugung der neuen
    * Signale auf diff_steps Unterzeitschritte zu verteilen!
    */
   /*
    * // eventuelle zusaetzliche raeumlich differenzierte Produktion von Signalen:
    * for (long n=0; n<pointnumber; n++) {
    * sigsknot[n].signal_new[sig_type]+=anregung(n);
    * }
    */
   for (short sig = 0; sig < signals; sig++) {
      if ((signal_use[sig] == 1) && (fix_signal[sig] == 0)) {
         if (diffusion_mode != 0) {
            // =============================================
            // This for diffuse_ADI _CN _Euler:
            // =============================================
            // berechne die Signalproduktion pro Unterzeitschritt
            for (long n = 0; n < pointnumber; n++) {
               sigsknot[n].signal_tmp[sig]
                  = ((sigsknot[n].signal_new[sig] - sigsknot[n].signal[sig])
                     / diff_steps[sig]);
               /*
                * // Falls Signal entfernt wird, ist dies vor der Diffusion zu tun,
                * // denn sonst koennte das Signal, das entfernt werden soll, vor der
                * // Entfernung aus dem Punkt herausdiffundieren!!! Dies muss sogar fuer
                * // alle Unterzeitschritte mitgemacht werden!
                * if (sigsknot[n].signal_tmp[sig]<0) {
                * sigsknot[n].signal[sig]=sigsknot[n].signal_new[sig];
                * sigsknot[n].signal_tmp[sig]=0.0;
                * }
                * // Bei positivem _tmp wird verzoegert produziert!
                * else sigsknot[n].signal_new[sig]=sigsknot[n].signal[sig];
                */
               // ### Es wird jetzt Verbrauch und Produktion verteilt --> testen
               sigsknot[n].signal_new[sig] = sigsknot[n].signal[sig];
            }

            if (diffusion_mode == 2) {
               if (dim == 3) {
                  signal_diffuse_ADI(signal_molecule(sig));
               } else { signal_diffuse_CN(signal_molecule(sig)); }
            } else {
               signal_diffuse_EULER(signal_molecule(sig));
            }
         } else {
            // =============================================	// This is for diffuse_QUANTA:
            // =============================================
            /* Hier wird die Produktion nicht auf Zeitschritte verteilt.
             * Es stellt sich die Frage, ob die Produktion vor oder nach
             * dem ersten Diffusionszeitschritt einbezogen werden sollte.
             * Um die im gleichen Zeitschritt produzierten Molekuele
             * bereits jetzt mit zu diffundieren ist die folgende Zeile zu aktivieren: */
            // for (long n=0; n<pointnumber; n++)
            // sigsknot[n].signal[sig]=sigsknot[n].signal_new[sig];
            /* Negative Signaldichten treten nicht auf, da QUANTA-Diffusion immer
             * mit receptor_use=0 laeuft und dort Signalverbrauch sofort umgesetzt wird. */
            // Dann gehe durch die Zeitschritte
            for (unsigned int nn = 0; nn < diff_steps[sig]; nn++) {
               for (long n = 0; n < pointnumber; n++) {
                  signal_diffuse_QUANTA(n,signal_molecule(sig));
               }
               // save the new variables as standard
               actualize();
               // Spaetestens hier wird die Signalproduktion einbezogen (siehe oben).
            }
         }
      }
   }
}
void sigs::signal_diffuse_QUANTA(const long &i, const signal_molecule &sig_type) {
   // Diffusion nach klassischer Manier:
   short int n;
   // Addiere alle Nachbarn zu deltan
   double deltan = 0.0;
   for (n = 0; n < dim2; n++) {
      if ((knot[i].near_n[n] != -1) && (knot[knot[i].near_n[n]].status != external)) {
         deltan += sigsknot[knot[i].near_n[n]].signal[sig_type];
      } else { deltan += boundary.signal[sig_type]; }
   }
   // Subtrahiere dim2 mal die aktuelle Signaldichte bei i
   if (knot[i].status == external) {
      deltan -= (dim2 * boundary.signal[sig_type]);
   } else {
      deltan -= (dim2 * sigsknot[i].signal[sig_type]);
   }
   // Gewichte den Wert mit der Diffusionskonstante, die bereits dx und dt enthaelt:
   double dn = p_diffuse_signal[sig_type] * deltan;

   // Unterbinde Fluss in diesen Punkt hinein:
   if (dn >= 0) { dn = 0.; }
   // Mache dn positiv
   else { dn *= -1.; }
   // Ganzahliger Teil von dn:
   deltan = int (dn);
   // Rest ist
   dn -= double (deltan);
   /* Das letzte Signalmolekuel wird ausgewuerfelt. Das ist notwendig, da
    * sonst bei kleinen Zeitschritten, wo vielleicht nur einzelne Molekuele
    * hin und wieder wandern, (deterministisch) niemals eine Bewegung stattfaende */
   if (drandom() < dn) { ++deltan; }
   // Jetzt enthaelt deltan die Zahl der zu verteilenden Signalmolekuele
   /* Kleines Problem: Falls die dynamische Rezeptor-Ligand-Bindung aktiviert
    * ist und trotzdem Signal in Quanten produziert und diffundiert wird, kann
    * hier durch die zufaellige Wahl des Restes mehr verbraucht werden als da ist.
    * Das muss korrigiert werden: */
   // if (deltan>sigsknot[i].signal[sig_type]) deltan=sigsknot[i].signal[sig_type];
   // ... nicht mehr noetig, da ich receptor_use=2 mit Quanta-Diffusion verboten habe.

   // Speichere die Veraenderung in _new
   if (knot[i].status != external) { signal_put(i,sig_type,(-1.0) * deltan); }
   if (sigsknot[i].signal_new[sig_type] < 0) {
      cout << "Negative signal-density=" << sigsknot[i].signal_new[sig_type]
           << " of type " << sig_type << " in signal_diffuse_QUANTA!\n";
      exit(1);
   }
   // Verteile die deltan Molekuele auf die Nachbarn:
   // Zunaechst random Reihenfolge der Nachbarn erzeugen:
   nn = nn_permuts.get_rand_set();
   double stmp;
   double scop;
   short int putit;
   long whichone;
   long whichtmp;
   // if (deltan>1) cout<<"punkt="<<i<<","<<knot[i].status<<" deltan="<<deltan<<":\n";
   // Dann verteilen:
   for (int k = 0; k < deltan; k++) {
      stmp = SIGNAL_MAX;
      whichone = -1;
      // Suche den Nachbarn mit der minimalen Menge Signalmolekuele:
      // dabei werden die unveraenderten Werte zugrundegelegt, da die Berechnung
      // auf der Basis der alten Werte erfolgt
      for (n = 0; n < dim2; n++) {
         whichtmp = knot[i].near_n[nn[n]];
         // scop=0.0; // Am Rand wird herausdiffundiert!
         scop = boundary.signal[sig_type];
         if (whichtmp != -1) {
            if (knot[whichtmp].status != external) {
               scop = sigsknot[whichtmp].signal_new[sig_type];
            }
         }
         if (scop < stmp) {
            stmp = scop;
            whichone = whichtmp;
         }
      }
      // whichone enthaelt die Position, stmp den kleinsten Wert. Speichern:
      putit = 0;
      if (whichone != -1) {
         if (knot[whichone].status != external) { putit = 1; }
      }
      if (putit == 1) {
         signal_put(whichone,sig_type,1);
         if (knot[i].status == external) {
            ++signal_total[sig_type];
            --signal_diffout[sig_type];
         }
      } else {
         if (knot[i].status != external) {
            // Aus der totalen Zahl der erzeugten Signale entfernen
            --signal_total[sig_type];
            ++signal_diffout[sig_type];
         }
      }
      /*
       * if (deltan>1) {
       * cout<<"  k="<<k<<": whichtmp="<<whichtmp<<" whichone="<<whichone
       * <<" stmp="<<stmp<<" putit="<<putit;
       * if (whichone!=-1) cout<<" type"<<knot[whichone].status;
       * cout<<"\n";
       * }
       */
   }
}
void sigs::signal_diffuse_EULER(const signal_molecule &sig_type) {
   // New variant in v5.01.1
   for (unsigned int nt = 0; nt < diff_steps[sig_type]; nt++) {
      for (long i = 0; i < pointnumber; i++) {
         if (knot[i].status != external) {
            double deltan = 0.0;
            short degrees_of_freedom = 0;

            // addiere die Signalproduktion in signal_new[]:
            signal_put(i,sig_type,sigsknot[i].signal_tmp[sig_type]);

            // Addiere alle Nachbarn zu deltan
            for (short n = 0; n < dim2; n++) {
               if ((knot[i].near_n[n] != -1)
                   && (knot[knot[i].near_n[n]].status != external)) {
                  deltan += sigsknot[knot[i].near_n[n]].signal[sig_type];
                  ++degrees_of_freedom;
               } else {
                  deltan += boundary.signal[sig_type];
                  ++degrees_of_freedom;
                  signal_diffout[sig_type] += (p_diffuse_signal[sig_type]
                                               * (sigsknot[i].signal[sig_type]
                                                  - boundary.signal[sig_type]));
               }
            }

            // Subtrahiere dim2 mal die aktuelle Signaldichte bei i
            deltan -= (double (degrees_of_freedom) * sigsknot[i].signal[sig_type]);
            // The following variant should hold for no obstacles! But anyway.
            // deltan-=(double(dim2)*sigsknot[i].signal[sig_type]);
            // Gewichte den Wert mit der Diffusionskonstante, die bereits dx und dt enthaelt:
            deltan *= p_diffuse_signal[sig_type];

            // Speichere die Veraenderung in _new
            signal_put(i,sig_type,deltan);
            signal_total[sig_type] += deltan;
            if (sigsknot[i].signal_new[sig_type] < 0) {
               cout << "\nNegative signal-density=" << sigsknot[i].signal_new[sig_type]
                    << " of type " << sig_type << " at point " << i
                    << " in signal_diffuse_EULER!\n";
               exit(1);
            }
         }
      }
      // save the new variables as standard
      actualize();
   }
}
void sigs::signal_diffuse_EULER(const long &i,
                                const signal_molecule &sig_type,
                                space &xyz) {
   double deltan = 0.0;
   short degrees_of_freedom = 0;

   // addiere die Signalproduktion in signal_new[]:
   signal_put(i,sig_type,sigsknot[i].signal_tmp[sig_type]);

   /* ##### ignore_objects ist nun neu zu programmieren:
    * Aktualisiere den grid::status auf dem Signal Gitter von space:: aus.
    * Rufe hier nur noch den Status ab. Beachte, dass das alte Kriterium
    * nur Sinn macht, wenn die Aufloesungen beider Gitter identisch sind.
    * !!! EULER IST IM MOMENT NICHT VERWENDBAR !!! 6.1.2005
    */
   // Falls ignore_objects==0 sind die internen Gitterpunkte eines Objekts
   // gesondert zu behandeln! Pruefe also dann, ob ein interner Punkt vorliegt:
   if ((ignore_objects == 1)
       || (knot[i].status == nothing)
       || (xyz.object_border(i) == 1)) {
      // Addiere alle Nachbarn zu deltan
      for (short n = 0; n < dim2; n++) {
         if ((knot[i].near_n[n] != -1) && (knot[knot[i].near_n[n]].status != external)) {
            if ((ignore_objects == 1)
                || (knot[knot[i].near_n[n]].status == nothing)) {
               // ||
               //   object_border(knot[i].near_n[n])==1)
               deltan += sigsknot[knot[i].near_n[n]].signal[sig_type];
               ++degrees_of_freedom;
            }
         } else {
            deltan += boundary.signal[sig_type];
            ++degrees_of_freedom;
            signal_diffout[sig_type] += (p_diffuse_signal[sig_type]
                                         * (sigsknot[i].signal[sig_type]
                                            - boundary.signal[sig_type]));
         }
      }
      // Subtrahiere dim2 mal die aktuelle Signaldichte bei i
      deltan -= (double (degrees_of_freedom) * sigsknot[i].signal[sig_type]);
      //  deltan-=(double(dim2)*sigsknot[i].signal[sig_type]);
      // Gewichte den Wert mit der Diffusionskonstante, die bereits dx und dt enthaelt:
      deltan *= p_diffuse_signal[sig_type];

      // Speichere die Veraenderung in _new
      signal_put(i,sig_type,deltan);
      signal_total[sig_type] += deltan;
      if (sigsknot[i].signal_new[sig_type] < 0) {
         cout << "\nNegative signal-density=" << sigsknot[i].signal_new[sig_type]
              << " of type " << sig_type << " at point " << i << " in signal_diffuse_EULER!\n";
         exit(1);
      }
   } else {
      // Der Fall don't ignore objects and internal point of an object:
      /* Im Fall von Objektwachstum kann sich auf den internen Punkten eines Objekts
       * eine Restkonzentration von Signal befinden. Das gleiche gilt nach der
       * Bewegung von Objekten. Diese Restkonzentration ist zu eliminieren und
       * an die naechstgelegenen Randpunkte des Objekts abzugeben.
       */
      if (sigsknot[i].signal_new[sig_type] > 0) {
         // merke die aktuelle Signalkonzentration
         deltan = sigsknot[i].signal_new[sig_type];
         // finde naechstgelegenen Randpunkt des gleichen Objekts
         // und addiere die Signale dort
         signal_put(xyz.next_border(i),sig_type,deltan);
         // ##### siehe Kommentar oben zu ignore_object (nur fuer EULER) ...
         // signal_total ist dadurch unveraendert!
         // loesche die Signale vom internen Punkt
         sigsknot[i].signal_new[sig_type] = 0;
      }
   }
}
double sigs::get_rhs_x_ADI(const long &r,
                           const double &alpha,
                           const signal_molecule &sig_type) {
   double rhs;
   if (knot[r].near_n[1] != -1) {
      rhs = sigsknot[knot[r].near_n[1]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs = boundary.signal[sig_type];
   } else { rhs = sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[4] != -1) {
      rhs += sigsknot[knot[r].near_n[4]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   if (knot[r].near_n[2] != -1) {
      rhs += sigsknot[knot[r].near_n[2]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else { rhs += sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[5] != -1) {
      rhs += sigsknot[knot[r].near_n[5]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   rhs *= alpha;
   rhs += (3.0 * sigsknot[r].signal[sig_type]);
   rhs += sigsknot[r].signal_tmp[sig_type];
   // +delta_t*(*ratefield)(i,j,k)); // nur bei Produktion auf dem Gitter
   return rhs;
}
double sigs::get_rhs_y_ADI(const long &r, const double &alpha,
                           const signal_molecule &sig_type) {
   double rhs;
   if (knot[r].near_n[0] != -1) {
      rhs = sigsknot[knot[r].near_n[0]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs = boundary.signal[sig_type];
   } else { rhs = sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[3] != -1) {
      rhs += sigsknot[knot[r].near_n[3]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   if (knot[r].near_n[2] != -1) {
      rhs += sigsknot[knot[r].near_n[2]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else { rhs += sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[5] != -1) {
      rhs += sigsknot[knot[r].near_n[5]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   rhs *= alpha;
   rhs += (3.0 * sigsknot[r].signal[sig_type]);
   rhs += sigsknot[r].signal_tmp[sig_type];
   // +delta_t*(*ratefield)(i,j,k)); // nur bei Produktion auf dem Gitter
   return rhs;
}
double sigs::get_rhs_z_ADI(const long &r, const double &alpha,
                           const signal_molecule &sig_type) {
   double rhs;
   if (knot[r].near_n[0] != -1) {
      rhs = sigsknot[knot[r].near_n[0]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs = boundary.signal[sig_type];
   } else { rhs = sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[3] != -1) {
      rhs += sigsknot[knot[r].near_n[3]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   if (knot[r].near_n[1] != -1) {
      rhs += sigsknot[knot[r].near_n[1]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else { rhs += sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[4] != -1) {
      rhs += sigsknot[knot[r].near_n[4]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   rhs *= alpha;
   rhs += (3.0 * sigsknot[r].signal[sig_type]);
   rhs += sigsknot[r].signal_tmp[sig_type];
   // +delta_t*(*ratefield)(i,j,k)); // nur bei Produktion auf dem Gitter
   return rhs;
}
void sigs::signal_diffuse_ADI(const signal_molecule &sig_type) {
   /* ADI ist fuer Diffusion mit Quellterm gemacht (Heat transfer equation).
    * Hier wird die Produktion durch die Zellen in den Zellklassen gemacht und
    * auf das Gitter als signal_new eingetragen. Entsprechend der Zahl der
    * internen Zeitschritte in der Diffusion wird diese Produktion pro internen
    * Zeitschritt gleichmaessig umgerechnet, in signal_tmp gespeichert
    * und in dieser Routine als Quelle verwendet.
    */

   for (unsigned int loop = 0; loop < diff_steps[sig_type]; loop++) {
      double bet,rhs,bound;
      long r;

      // first third of the timestep: coordinate x is treated implicitly with delta_t/3.0
      // Uebernehme meine alte Nummerierung und ergaenze diese durch zwei temp-Randpunkte.
      for (int j = 0; j < prodimvec[1]; j++) {
         for (int k = 0; k < prodimvec[2]; k++) {
            int i = 0;
            r = Index(i,j,k);

            // repeat the whole procedure for all reaction volumes:
            while (i < prodimvec[0]) {
               // Look for the first non-border point:
               while (i < prodimvec[0] && knot[r].status == external) {
                  ++i;
                  if (i < prodimvec[0]) { r = Index(i,j,k); }
               }
               // (i,j,k) points on the first real point
               // save this one:
               int i0 = i;
               // Note: i0=0 is possible here, implying the boundary outside the lattice!

               // Start matrix inversion if a point was found:
               if (i0 < prodimvec[0]) {
                  // determine the boundary value
                  if (dirichlet[sig_type]) {
                     bound = boundary.signal[sig_type];           // case Dirichlet
                  } else { bound = dx * boundary.signal[sig_type]; }
                  // cerr<<"bound="<<bound<<"\n";
                  // cerr<<"alpha="<<p_diffuse_signal[sig_type]<<"\n";

                  // calculate first real point at Index r=(i,j,k)
                  if (dirichlet[sig_type]) {
                     bet = 3.0 + 2.0 * p_diffuse_signal[sig_type];
                  } else { bet = 3.0 + p_diffuse_signal[sig_type]; }
                  gam[i] = (-p_diffuse_signal[sig_type] / bet);
                  // cerr<<"i="<<i<<", r="<<r<<", beta="<<bet<<", gamma["<<i<<"]="<<gam[i]<<"\n";
                  rhs = get_rhs_x_ADI(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type] * bound) / bet;

                  // go through all non-border points now
                  ++i;
                  r = Index(i,j,k);
                  while (i < prodimvec[0] && knot[r].status != external) {
                     // in the case of dirichlet all points are treated equally
                     // for von Neumann the last point and the point in front of external are
                     // different
                     bet = (3.0 + p_diffuse_signal[sig_type] * (2.0 + gam[i - 1]));
                     gam[i] = (-p_diffuse_signal[sig_type] / bet);
                     // cerr<<"i="<<i<<", r="<<r<<", beta="<<bet<<", gamma["<<i<<"]="<<gam[i]<<"\n";
                     rhs = get_rhs_x_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * sigsknot[knot[r].near_n[0]].signal_new[sig_type]) / bet;
                     ++i;
                     if (i < prodimvec[0]) { r = Index(i,j,k); }
                  }
                  // comes out with i pointing on the first border point
                  //           and r=(i,j,k) or r=(prodimvec[0]-1,j,k)
                  r = Index(i - 1,j,k);
                  // In the last point, the Dirichlet boundary condition has to be added:
                  if (dirichlet[sig_type]) {
                     sigsknot[r].signal_new[sig_type] += p_diffuse_signal[sig_type] * bound
                                                         / bet;
                  } else {
                     // recalculate the last point in the while loop for von Neumann
                     --i;
                     bet = (3.0 + p_diffuse_signal[sig_type] * (1.0 + gam[i - 1]));
                     gam[i] = (-p_diffuse_signal[sig_type] / bet);
                     // cerr<<"i="<<i<<", r="<<r<<", beta="<<bet<<", gamma["<<i<<"]="<<gam[i]<<"\n";
                     rhs = get_rhs_x_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * (sigsknot[knot[r].near_n[0]].signal_new[sig_type]
                              - vonNEUMANNinwards * bound)) / bet;
                     // note that bound enters with negative sign for it was defined with an added -
                     ++i;
                  }

                  // Ruecklauf:
                  for (int ii = (i - 1); ii > i0; ii--) {
                     r = Index(ii,j,k);
                     sigsknot[knot[r].near_n[0]].signal_new[sig_type]
                        -= gam[ii - 1] * sigsknot[r].signal_new[sig_type];
                  }
                  // i points still on the first border point

                  /* use this if considering non-convex or multiple reaction volumes
                   * and restart the whole procedure */
                  if (i < prodimvec[0]) { r = Index(i,j,k); }
               }      // end of if (i0<prodimvec[0]) {
            }     // end of while (i<prodimvec[0]) {
         }    // end of for (k)
      }   // end of for (j)

      // Aktualisiere den Wert von Signal
      for (r = 0; r < pointnumber; r++) {
         if (knot[r].status != external) {
            sigsknot[r].signal[sig_type] = sigsknot[r].signal_new[sig_type];
         }
      }

      // second third of the timestep: coordinate y is treated implicitly with delta_t/3.0
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int k = 0; k < prodimvec[2]; k++) {
            int j = 0;
            r = Index(i,j,k);

            // repeat the whole procedure for all reaction volumes:
            while (j < prodimvec[1]) {
               // Look for the first non-border point:
               while (j < prodimvec[1] && knot[r].status == external) {
                  ++j;
                  if (j < prodimvec[1]) { r = Index(i,j,k); }
               }
               // (i,j,k) points on the first real point
               // save this one:
               int j0 = j;

               // Start matrix inversion if a point has been found:
               if (j0 < prodimvec[1]) {
                  // determine the boundary value
                  if (dirichlet[sig_type]) {
                     bound = boundary.signal[sig_type];           // case Dirichlet
                  } else { bound = dx * boundary.signal[sig_type]; }

                  // calculate first real point at Index r=(i,j,k)
                  if (dirichlet[sig_type]) {
                     bet = 3.0 + 2.0 * p_diffuse_signal[sig_type];
                  } else { bet = 3.0 + p_diffuse_signal[sig_type]; }
                  gam[j] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_y_ADI(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type] * bound) / bet;

                  // go through all non-border points now
                  ++j;
                  r = Index(i,j,k);
                  while (j < prodimvec[1] && knot[r].status != external) {
                     bet = (3.0 + p_diffuse_signal[sig_type] * (2.0 + gam[j - 1]));
                     gam[j] = (-p_diffuse_signal[sig_type] / bet);
                     rhs = get_rhs_y_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * sigsknot[knot[r].near_n[1]].signal_new[sig_type]) / bet;
                     ++j;
                     if (j < prodimvec[1]) { r = Index(i,j,k); }
                  }
                  // comes out with j pointing on the first border point
                  //           and r=(i,j,k) or r=(i,prodimvec[1]-1,k)
                  r = Index(i,j - 1,k);
                  // In the last point, the Dirichlet boundary condition has to be added:
                  if (dirichlet[sig_type]) {
                     sigsknot[r].signal_new[sig_type] += p_diffuse_signal[sig_type] * bound
                                                         / bet;
                  } else {
                     // recalculate the last point in the while loop for von Neumann
                     --j;
                     bet = (3.0 + p_diffuse_signal[sig_type] * (1.0 + gam[j - 1]));
                     gam[j] = (-p_diffuse_signal[sig_type] / bet);
                     rhs = get_rhs_y_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * (sigsknot[knot[r].near_n[1]].signal_new[sig_type]
                              - vonNEUMANNinwards * bound)) / bet;
                     // note that bound enters with negative sign for it was defined with an added -
                     ++j;
                  }

                  // Ruecklauf:
                  for (int jj = (j - 1); jj > j0; jj--) {
                     r = Index(i,jj,k);
                     sigsknot[knot[r].near_n[1]].signal_new[sig_type]
                        -= gam[jj - 1] * sigsknot[r].signal_new[sig_type];
                  }
                  // j points still on the first border point

                  /* use this if considering non-convex or multiple reaction volumes
                   * and restart the whole procedure */
                  if (j < prodimvec[1]) { r = Index(i,j,k); }
               }      // end of if (j0<prodimvec[1]) {
            }     // end of while (j<prodimvec[1]) {
         }
      }
      // Aktualisiere den Wert von Signal
      for (r = 0; r < pointnumber; r++) {
         if (knot[r].status != external) {
            sigsknot[r].signal[sig_type] = sigsknot[r].signal_new[sig_type];
         }
      }

      // last third of the timestep: coordinate z is treated implicitly with delta_t/3.0
      for (int i = 0; i < prodimvec[0]; i++) {
         for (int j = 0; j < prodimvec[1]; j++) {
            int k = 0;
            r = Index(i,j,k);

            // repeat the whole procedure for all reaction volumes:
            while (k < prodimvec[2]) {
               // Look for the first non-border point:
               while (k < prodimvec[2] && knot[r].status == external) {
                  ++k;
                  if (k < prodimvec[2]) { r = Index(i,j,k); }
               }
               // (i,j,k) points on the first real point
               // save this one:
               int k0 = k;

               // Start matrix inversion if a point has been found:
               if (k0 < prodimvec[2]) {
                  // determine the boundary value
                  if (dirichlet[sig_type]) {
                     bound = boundary.signal[sig_type];           // case Dirichlet
                  } else { bound = dx * boundary.signal[sig_type]; }

                  // calculate first real point at Index r=(i,j,k)
                  if (dirichlet[sig_type]) {
                     bet = 3.0 + 2.0 * p_diffuse_signal[sig_type];
                  } else { bet = 3.0 + p_diffuse_signal[sig_type]; }
                  gam[k] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_z_ADI(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type] * bound) / bet;

                  // go through all non-border points now
                  ++k;
                  r = Index(i,j,k);
                  while (k < prodimvec[2] && knot[r].status != external) {
                     bet = (3.0 + p_diffuse_signal[sig_type] * (2.0 + gam[k - 1]));
                     gam[k] = (-p_diffuse_signal[sig_type] / bet);
                     rhs = get_rhs_z_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * sigsknot[knot[r].near_n[2]].signal_new[sig_type]) / bet;
                     ++k;
                     if (k < prodimvec[2]) { r = Index(i,j,k); }
                  }
                  // comes out with k pointing on the first border point
                  //           and r=(i,j,k) or r=(i,j,prodimvec[2]-1)
                  r = Index(i,j,k - 1);
                  // In the last point, the Dirichlet boundary condition has to be added:
                  if (dirichlet[sig_type]) {
                     sigsknot[r].signal_new[sig_type] += p_diffuse_signal[sig_type] * bound
                                                         / bet;
                  } else {
                     // recalculate the last point in the while loop for von Neumann
                     --k;
                     bet = (3.0 + p_diffuse_signal[sig_type] * (1.0 + gam[k - 1]));
                     gam[k] = (-p_diffuse_signal[sig_type] / bet);
                     rhs = get_rhs_z_ADI(r,p_diffuse_signal[sig_type],sig_type);
                     sigsknot[r].signal_new[sig_type]
                        = (rhs + p_diffuse_signal[sig_type]
                           * (sigsknot[knot[r].near_n[2]].signal_new[sig_type]
                              - vonNEUMANNinwards * bound)) / bet;
                     // note that bound enters with negative sign for it was defined with an added -
                     ++k;
                  }

                  // Ruecklauf:
                  for (int kk = (k - 1); kk > k0; kk--) {
                     r = Index(i,j,kk);
                     sigsknot[knot[r].near_n[2]].signal_new[sig_type]
                        -= gam[kk - 1] * sigsknot[r].signal_new[sig_type];
                  }
                  // k points still on the first border point

                  /* use this if considering non-convex or multiple reaction volumes
                   * and restart the whole procedure */
                  if (k < prodimvec[2]) { r = Index(i,j,k); }
               }      // end of if (k0<prodimvec[2]) {
            }     // end of while (k<prodimvec[2]) {
         }
      }
      // Aktualisiere den Wert von Signal
      for (r = 0; r < pointnumber; r++) {
         if (knot[r].status != external) {
            sigsknot[r].signal[sig_type] = sigsknot[r].signal_new[sig_type];
         }
      }
   }
   // Aktualisiere den Wert von Signal
   /*
    * double tmptot=signal_total[sig_type];
    * signal_total[sig_type]=0.0;
    * for (long r=0; r<pointnumber; r++) {
    * if (knot[r].status!=external) {
    *  signal_total[sig_type]+=sigsknot[r].signal_new[sig_type];
    * }
    * }
    * signal_diffout[sig_type]+=(tmptot-signal_total[sig_type]);
    */
}
double sigs::get_rhs_x_CN(const long &r,const double &alpha,
                          const signal_molecule &sig_type) {
   double rhs;
   if (knot[r].near_n[1] != -1) {
      rhs = sigsknot[knot[r].near_n[1]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs = boundary.signal[sig_type];
   } else { rhs = sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[3] != -1) {
      rhs += sigsknot[knot[r].near_n[3]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   rhs *= alpha;
   rhs += (2.0 * sigsknot[r].signal[sig_type]);
   rhs += sigsknot[r].signal_tmp[sig_type];
   // +delta_t*(*ratefield)(i,j,k)); // nur bei Produktion auf dem Gitter
   return rhs;
}
double sigs::get_rhs_y_CN(const long &r,const double &alpha,
                          const signal_molecule &sig_type) {
   double rhs;
   if (knot[r].near_n[0] != -1) {
      rhs = sigsknot[knot[r].near_n[0]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs = boundary.signal[sig_type];
   } else { rhs = sigsknot[r].signal[sig_type] + dx * boundary.signal[sig_type]; }
   rhs -= (2.0 * sigsknot[r].signal[sig_type]);
   if (knot[r].near_n[2] != -1) {
      rhs += sigsknot[knot[r].near_n[2]].signal[sig_type];
   } else if (dirichlet[sig_type]) {
      rhs += boundary.signal[sig_type];
   } else {
      rhs += sigsknot[r].signal[sig_type] - vonNEUMANNinwards * dx
             * boundary.signal[sig_type];
   }
   rhs *= alpha;
   rhs += (2.0 * sigsknot[r].signal[sig_type]);
   rhs += sigsknot[r].signal_tmp[sig_type];
   // +delta_t*(*ratefield)(i,j,k)); // nur bei Produktion auf dem Gitter
   return rhs;
}
void sigs::signal_diffuse_CN(const signal_molecule &sig_type) {
   /* CN steht fuer Cranck-Nickolson und ist analog zu ADI in 3D.
    * Weitere Bemerkungen dort. Behandlung identisch.
    */

   for (unsigned int loop = 0; loop < diff_steps[sig_type]; loop++) {
      double bet,rhs,bound;
      long r;

      // first third of the timestep: coordinate x is treated implicitly with delta_t/2.0
      // Uebernehme meine alte Nummerierung und ergaenze diese durch zwei temp-Randpunkte.
      for (int j = 0; j < prodimvec[1]; j++) {
         int i = 0;
         r = Index(i,j);

         // repeat the whole procedure for all reaction volumes:
         while (i < prodimvec[0]) {
            // Look for the first non-border point:
            while (i < prodimvec[0] && knot[r].status == external) {
               ++i;
               if (i < prodimvec[0]) { r = Index(i,j); }
            }
            // (i,j) points on the first real point
            // save this one:
            int i0 = i;

            // Start matrix inversion if a point has been found:
            if (i0 < prodimvec[0]) {
               // determine the boundary value
               if (dirichlet[sig_type]) {
                  bound = boundary.signal[sig_type];            // case Dirichlet
               } else { bound = dx * boundary.signal[sig_type]; }
               // cerr<<"bound="<<bound<<"\n";
               // cerr<<"alpha="<<p_diffuse_signal[sig_type]<<"\n";

               // calculate first real point at Index r=(i,j)
               if (dirichlet[sig_type]) {
                  bet = 2.0 + 2.0 * p_diffuse_signal[sig_type];
               } else { bet = 2.0 + p_diffuse_signal[sig_type]; }
               gam[i] = (-p_diffuse_signal[sig_type] / bet);
               // cerr<<"i="<<i<<", r="<<r<<", beta="<<bet<<", gamma["<<i<<"]="<<gam[i]<<"\n";
               rhs = get_rhs_x_CN(r,p_diffuse_signal[sig_type],sig_type);
               sigsknot[r].signal_new[sig_type] = (rhs + p_diffuse_signal[sig_type] * bound)
                                                  / bet;

               // go through all non-border points now
               ++i;
               r = Index(i,j);
               while (i < prodimvec[0] && knot[r].status != external) {
                  bet = (2.0 + p_diffuse_signal[sig_type] * (2.0 + gam[i - 1]));
                  gam[i] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_x_CN(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type]
                        * sigsknot[knot[r].near_n[0]].signal_new[sig_type]) / bet;
                  /*
                   * if (sigsknot[r].signal_new[sig_type]<0) {
                   * cout<<"negative signal "<<sig_type
                   * <<". r="<<r<<"=("<<i<<","<<j<<"); rhs="<<rhs
                   * <<"; sig_new="<<sigsknot[r].signal_new[sig_type]<<".\n";
                   * exit(1);
                   * }
                   */
                  ++i;
                  if (i < prodimvec[0]) { r = Index(i,j); }
               }
               // comes out with i pointing on the first border point
               // and r=(i,j) or r=(prodimvec[0]-1,j)
               r = Index(i - 1,j);

               if (dirichlet[sig_type]) {
                  sigsknot[r].signal_new[sig_type] += p_diffuse_signal[sig_type] * bound
                                                      / bet;
               } else {
                  // recalculate the last point in the while loop for von Neumann
                  --i;
                  bet = (2.0 + p_diffuse_signal[sig_type] * (1.0 + gam[i - 1]));
                  gam[i] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_x_CN(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type]
                        * (sigsknot[knot[r].near_n[0]].signal_new[sig_type]
                           - vonNEUMANNinwards * bound)) / bet;
                  // note that bound enters with negative sign for it was defined with an added -
                  ++i;
               }

               // Ruecklauf:
               for (int ii = (i - 1); ii > i0; ii--) {
                  r = Index(ii,j);
                  sigsknot[knot[r].near_n[0]].signal_new[sig_type]
                     -= gam[ii - 1] * sigsknot[r].signal_new[sig_type];
               }
               // i points still on the first border point

               /* use this if considering non-convex or multiple reaction volumes
                * and restart the whole procedure */
               if (i < prodimvec[0]) { r = Index(i,j); }
            }     // end of if (i0<prodimvec[0]) {
         }    // end of while (i<prodimvec[0]) {
      }   // end of for (j)

      // Aktualisiere den Wert von Signal
      for (r = 0; r < pointnumber; r++) {
         if (knot[r].status != external) {
            sigsknot[r].signal[sig_type] = sigsknot[r].signal_new[sig_type];
         }
      }

      // second half of the timestep: coordinate y is treated implicitly with delta_t/2.0
      for (int i = 0; i < prodimvec[0]; i++) {
         int j = 0;
         r = Index(i,j);

         // repeat the whole procedure for all reaction volumes:
         while (j < prodimvec[1]) {
            // Look for the first non-border point:
            while (j < prodimvec[1] && knot[r].status == external) {
               ++j;
               if (j < prodimvec[1]) { r = Index(i,j); }
            }
            // (i,j) points on the first real point
            // save this one:
            int j0 = j;

            // Start matrix inversion if a point has been found:
            if (j0 < prodimvec[1]) {
               // determine the boundary value
               if (dirichlet[sig_type]) {
                  bound = boundary.signal[sig_type];            // case Dirichlet
               } else { bound = dx * boundary.signal[sig_type]; }
               // cerr<<"bound="<<bound<<"\n";
               // cerr<<"alpha="<<p_diffuse_signal[sig_type]<<"\n";

               // calculate first real point at Index r=(i,j)
               if (dirichlet[sig_type]) {
                  bet = 2.0 + 2.0 * p_diffuse_signal[sig_type];
               } else { bet = 2.0 + p_diffuse_signal[sig_type]; }
               gam[j] = (-p_diffuse_signal[sig_type] / bet);
               // cerr<<"j="<<j<<", r="<<r<<", beta="<<bet<<", gamma["<<j<<"]="<<gam[j]<<"\n";
               rhs = get_rhs_y_CN(r,p_diffuse_signal[sig_type],sig_type);
               sigsknot[r].signal_new[sig_type] = (rhs + p_diffuse_signal[sig_type] * bound)
                                                  / bet;

               // go through all non-border points now
               ++j;
               r = Index(i,j);
               while (j < prodimvec[1] && knot[r].status != external) {
                  bet = (2.0 + p_diffuse_signal[sig_type] * (2.0 + gam[j - 1]));
                  gam[j] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_y_CN(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type]
                        * sigsknot[knot[r].near_n[1]].signal_new[sig_type]) / bet;
                  ++j;
                  if (j < prodimvec[1]) { r = Index(i,j); }
               }
               // comes out with j pointing on the first border point
               // and r=(i,j) or r=(i,prodimvec[1]-1)
               r = Index(i,j - 1);

               if (dirichlet[sig_type]) {
                  sigsknot[r].signal_new[sig_type] += p_diffuse_signal[sig_type] * bound
                                                      / bet;
               } else {
                  // recalculate the last point in the while loop for von Neumann
                  --j;
                  bet = (2.0 + p_diffuse_signal[sig_type] * (1.0 + gam[j - 1]));
                  gam[j] = (-p_diffuse_signal[sig_type] / bet);
                  rhs = get_rhs_y_CN(r,p_diffuse_signal[sig_type],sig_type);
                  sigsknot[r].signal_new[sig_type]
                     = (rhs + p_diffuse_signal[sig_type]
                        * (sigsknot[knot[r].near_n[1]].signal_new[sig_type]
                           - vonNEUMANNinwards * bound)) / bet;
                  // note that bound enters with negative sign for it was defined with an added -
                  ++j;
               }

               // Ruecklauf:
               for (int jj = (j - 1); jj > j0; jj--) {
                  r = Index(i,jj);
                  sigsknot[knot[r].near_n[1]].signal_new[sig_type]
                     -= gam[jj - 1] * sigsknot[r].signal_new[sig_type];
               }
               // j points still on the first border point

               /* use this if considering non-convex or multiple reaction volumes
                * and restart the whole procedure */
               if (j < prodimvec[1]) { r = Index(i,j); }
            }     // end of if (j0<prodimvec[1]) {
         }    // end of while (j<prodimvec[1]) {
      }
      // Aktualisiere den Wert von Signal
      for (r = 0; r < pointnumber; r++) {
         if (knot[r].status != external) {
            sigsknot[r].signal[sig_type] = sigsknot[r].signal_new[sig_type];
         }
      }
   }
   /*
    * double tmptot=signal_total[sig_type];
    * signal_total[sig_type]=0.0;
    * for (long r=0; r<pointnumber; r++) {
    * if (knot[r].status!=external) {
    * signal_total[sig_type]+=sigsknot[r].signal[sig_type];
    * }
    * }
    * signal_diffout[sig_type]+=(tmptot-signal_total[sig_type]);
    */
}
void sigs::actualize() {
   // Aktualisiert alle Signal-Variablen:
   for (short a = 0; a < signals; a++) {
      if (signal_use[a] == 1) {
         for (long n = 0; n < pointnumber; n++) {
            sigsknot[n].actualize(signal_molecule(a));
         }
      }
   }
}
void sigs::get_signal_name(const unsigned short &sx, char * name) {
   switch (sx) {
      case sig_differ2CC:
         strcat(name,"si");
         break;

      case sig_proliferate:
         strcat(name,"sipro");
         break;

      case CXCL12:
         strcat(name,"cxcl12");
         break;

      case CXCL13:
         strcat(name,"cxcl13");
         break;

      case antibody:
         strcat(name,"ab");
         break;

      case antigen:
         strcat(name,"ag");
         break;

      case SEMA4D:
         strcat(name,"sema4d");
         break;

      case glucose:
         strcat(name,"glucose");
         break;

      case oxygen:
         strcat(name,"oxygen");
         break;
   }
}
void sigs::write_siglog(double &t) {
   // Signal production
   siglog << t;
   for (short a = 0; a < signals; a++) {
      if (signal_use[a] == 1) {
         signal_total[a] = get_signal_total(signal_molecule(a));
         siglog << "   " << signal_total[a]
                << "   " << signal_diffout[a]
                << "   " << signal_used[a]
                << "   " << signal_produced[a];
         // signal_diffout[a]=0.0;
         signal_produced[a] = 0.0;
      }
   }
   siglog << "\n";
}
void sigs::write_files(suffix tnr, bool forceit) {
   // Find largest value of the signals
   double sigmax[signals];
   for (short a = 0; a < signals; a++) {
      sigmax[a] = 0.;
   }
   for (long i = 0; i < pointnumber; i++) {
      for (short a = 0; a < signals; a++) {
         if ((signal_use[a] == 1) && (sigsknot[i].signal[a] > sigmax[a])) {
            sigmax[a] = sigsknot[i].signal[a];
         }
      }
   }

   // Zeige Ebenen (je nach Dimension unterschiedlich)
   long k[dim];
   for (short d = 0; d < dim; d++) {
      k[d] = 0;
   }
   long d1 = 0;
   long d2 = 1;
   if (dim == 3) { d1 = 1; d2 = 2; k[0] = prodimvec[0] / 2; }

   // Open signal files: ================================
   ofstream sig[signals];
   ofstream sigevery[signals];
   for (unsigned short sx = 0; sx < signals; sx++) {
      if ((signal_use[sx] == 1) && ((forceit == true) || (fix_signal[sx] == false))) {
         char sigdatname[30] = "";
         get_signal_name(sx,sigdatname);
         strcat(sigdatname,tnr);
         strcat(sigdatname,".ppm");
         cout << "Write to " << sigdatname << "\n";

         // Stream oeffnen
         sig[sx].open(sigdatname);
         sig[sx].setf(ios::scientific, ios::floatfield);

         // Kopf:
         sig[sx] << "P3\n";
         sig[sx] << prodimvec[d1] << "\n";
         sig[sx] << prodimvec[d2] << "\n";
         // this is to be tested as it does not follow from the structure before #####
         if (system == 0) { sig[sx] << "255\n"; }

         char sigeverydatname[30] = "";
         get_signal_name(sx,sigeverydatname);
         strcat(sigeverydatname,tnr);
         strcat(sigeverydatname,".dat");
         cout << "Write to " << sigeverydatname << "\n";
         // Stream oeffnen
         sigevery[sx].open(sigeverydatname);
         sigevery[sx].setf(ios::scientific, ios::floatfield);
      }
   }

   // Write to the signal files:
   long i;
   double x,y,z;
   x = 0.;
   y = 0.;
   z = 0.5 * dx;
   for (long i1 = 0; i1 < prodimvec[d2]; i1++) {
      // 1-component in 2D, 2-component in 3D
      for (long i2 = 0; i2 < prodimvec[d1]; i2++) {
         // 0-component in 2D, 1-component in 3D
         k[d1] = i2;
         k[d2] = i1;
         x = (i2 + 0.5) * dx;
         y = (i1 + 0.5) * dx;
         i = Index(k);

         if (knot[i].status == external) {
            for (unsigned short sx = 0; sx < signals; sx++) {
               if ((signal_use[sx] == 1) && ((forceit == true) || (fix_signal[sx] == false))) {
                  sig[sx] << "63 63 63\n";       // dunkelgrau
               }
            }
         } else {
            for (unsigned short sx = 0; sx < signals; sx++) {
               if ((signal_use[sx] == 1) && ((forceit == true) || (fix_signal[sx] == false))) {
                  double tmp = sigsknot[i].signal[sx];
                  // if (sx==CXCL12) {
                  // cout<<"i="<<i<<"=("<<k[0]<<","<<k[1]<<","<<k[2]<<")="<<tmp<<";;; ";
                  // }
                  if (tmp > sigs::SIGNAL_MAX) {
                     cout << "Signal_density " << tmp
                          << " of type " << sx << " too large at point "
                          << i << " in cellman::xfiles()\n";
                     exit(1);
                  }
                  sigevery[sx] << x << " " << y << " " << z << " " << tmp << "\n";
                  if ((sx == CXCL12) || (sx == CXCL13) || (sx == SEMA4D)
                      || (sx == glucose) || (sx == oxygen) || (sx == antibody)) {
                     // rescale with maximum
                     if (sigmax[sx] > 0.) { tmp *= (50. / sigmax[sx]); }
                  }
                  if (tmp > 50.0) { tmp = 50.0; }
                  int tmpi = int ((tmp + 1.0) / 2.0);
                  switch (tmpi) {
                     case 25:
                        sig[sx] << "0 0 0\n";
                        break;

                     case 24:
                        sig[sx] << "0 0 20\n";
                        break;

                     case 23:
                        sig[sx] << "0 0 40\n";
                        break;

                     case 22:
                        sig[sx] << "0 0 60\n";
                        break;

                     case 21:
                        sig[sx] << "0 0 80\n";
                        break;

                     case 20:
                        sig[sx] << "0 0 100\n";
                        break;

                     case 19:
                        sig[sx] << "0 0 120\n";
                        break;

                     case 18:
                        sig[sx] << "0 0 140\n";
                        break;

                     case 17:
                        sig[sx] << "0 0 160\n";
                        break;

                     case 16:
                        sig[sx] << "0 0 180\n";
                        break;

                     case 15:
                        sig[sx] << "0 0 200\n";
                        break;

                     case 14:
                        sig[sx] << "0 0 220\n";
                        break;

                     case 13:
                        sig[sx] << "0 0 255\n";
                        break;

                     case 12:
                        sig[sx] << "20 20 255\n";
                        break;

                     case 11:
                        sig[sx] << "40 40 255\n";
                        break;

                     case 10:
                        sig[sx] << "60 60 255\n";
                        break;

                     case 9:
                        sig[sx] << "80 80 255\n";
                        break;

                     case 8:
                        sig[sx] << "100 100 255\n";
                        break;

                     case 7:
                        sig[sx] << "120 120 255\n";
                        break;

                     case 6:
                        sig[sx] << "140 140 255\n";
                        break;

                     case 5:
                        sig[sx] << "160 160 255\n";
                        break;

                     case 4:
                        sig[sx] << "180 180 255\n";
                        break;

                     case 3:
                        sig[sx] << "200 200 255\n";
                        break;

                     case 2:
                        sig[sx] << "220 220 255\n";
                        break;

                     case 1:
                        sig[sx] << "240 240 255\n";
                        break;

                     case 0:
                        sig[sx] << "255 255 255\n";
                        break;
                  }
                  // if (tmp>0) {
                  // sig<<"0 0 ";
                  // if (s.knot[i].signal[sig_differ2CC]+tmp>6) sig<<"0";
                  // else sig<<40*(7-s.knot[i].signal[sig_differ2CC]-tmp); // blautoene
                  // sig<<"\n";
                  // }
                  // else sig<<"255 255 255\n"; //weis
               }
            }
         }
      }
   }

   // Close signals files:
   for (unsigned short sx = 0; sx < signals; sx++) {
      if ((signal_use[sx] == 1) && ((forceit == true) || (fix_signal[sx] == false))) {
         sig[sx].close();
         sigevery[sx].close();
      }
   }
}
void sigs::write_TEST(double &t) {
   for (short sig_i = 0; sig_i < signals; sig_i++) {
      if (signal_use[sig_i] == 1) {
         char sigfname[30] = "sigtest_";
         get_signal_name(short (sig_i),sigfname);
         strcat(sigfname,TEST_SUF);
         strcat(sigfname,".out");
         cout << "Write signals to " << sigfname << "\n";
         ofstream sigf(sigfname);
         sigf << "! For signal " << sig_i << "\n"
              << "! at time " << t << " hours\n"
              << "! index : " << dim
              << " koords : signal : reference : signal-reference\n";
         for (long ii = 0; ii < pointnumber; ii++) {
            long tmp[dim];
            get_koord(ii,tmp);
            sigf << ii << "    ";
            for (short jj = 0; jj < dim; jj++) {
               sigf << tmp[jj] << "  ";
            }
            double ref = reference(tmp,t,signal_molecule(sig_i));
            sigf << "  "
                 << sigsknot[ii].signal[sig_i] << "   "
                 << ref << "   "
                 << sigsknot[ii].signal[sig_i] - ref
                 << "\n";
         }
         sigf.close();
      }
   }
   addchar(TEST_SUF);
}
void sigs::write_all_signals(double &chemo_max, double &chemo_steep,
                             double &chemo_half) {
   for (short sig_i = 0; sig_i < signals; sig_i++) {
      if (signal_use[sig_i] == 1) {
         char sigfname[30] = "sigs_";
         get_signal_name(sig_i,sigfname);
         strcat(sigfname,".out");
         cout << "Write signals to " << sigfname << "\n";
         ofstream sigf(sigfname);
         // !!! ATTENTION !!!
         // If you change the preamble static const N_PRE_WORDS has to be adapted in signals.h.
         sigf << "! For signal " << sig_i << "\n"
              << "! index : " << dim
              << " koords : conc : gradient : chemo-weight : %of chemo_max : grad-direction\n";
         for (long ii = 0; ii < pointnumber; ii++) {
            long tmp[dim];
            get_koord(ii,tmp);
            sigf << ii << "    ";
            for (short jj = 0; jj < dim; jj++) {
               sigf << tmp[jj] << "  ";
            }
            double grad[dim];
            get_gradient(signal_molecule(sig_i),ii,grad);
            double gradient = get_2norm(grad);
            for (short jj = 0; jj < dim; jj++) {
               if (gradient > 0.) { grad[jj] /= gradient; } else { grad[jj] = 0.; }
            }
            double alpha = chemo_max / (1.0 + exp(chemo_steep * (chemo_half - gradient)));
            sigf << "  "
                 << sigsknot[ii].signal[sig_i] << "  "
                 << gradient << "    "
                 << alpha << "  "
                 << alpha / chemo_max << "    ";
            for (short jj = 0; jj < dim; jj++) {
               sigf << grad[jj] << "  ";
            }
            sigf << "\n";
         }
         sigf.close();
      }
   }
}
