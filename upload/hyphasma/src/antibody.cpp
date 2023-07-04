#include "antibody.h"
#include <math.h>
#include <string.h>
#include <sstream>

int AntibodyDyn::antibodies_resolution = 10;             // Resolution of systemic antibodies in
                                                         // affinity bins (0=no systemic
                                                         // antibodies):
double AntibodyDyn::ag_portion = 1e-8;                      // Threshold Ag-concentration for
                                                            // binding CC (in Mol)"; break;
double AntibodyDyn::antibody_degradation = 30.0;            // Antibody degradation half time
                                                            // (days):
double AntibodyDyn::antibody_production_factor = 1e-17;     // Production rate of systemic
                                                            // antibodies by output population 
                                                            // (mol per hour and cell):
double AntibodyDyn::k_on = 3600 * 1e6;                      // in 1/(Mol hr)
double AntibodyDyn::k_ic_exp_min = 5.5;                     // Exponent of lowest affinity
                                                            // dissociation constant of immune
                                                            // complex (1/Mol):
double AntibodyDyn::k_ic_exp_max = 9.5;                     // Exponent of highest affinity
                                                            // dissociation constant of immune
                                                            // complex (1/Mol):
double AntibodyDyn::inject_ab_affinity = 7.1;               // Inject antibodies: exponent of
                                                            // affinity in l/mol:
long AntibodyDyn::inject_ab_ASindex = -1;                 // Inject antibodies: index in
                                                          // affinityspace (-1 to use exponent of
                                                          // affinity):
double AntibodyDyn::ab_injection_time = 144.0;              // Inject antibodies: time of
                                                            // administration (hours):
double AntibodyDyn::ab_injection_amount = -1.0;             // Inject antibodies: concentration
                                                            // (mol/l; <=0 for none):
double AntibodyDyn::consistency_threshold = 1.e-30;         // This threshold is used in
                                                            // AntibodyDyn::check4consistency(...)

AntibodyBin::AntibodyBin() { }
AntibodyBin::~AntibodyBin() { }

void AntibodyBin::initialise(int &resolution, int &ag,
                             double &kon, double &expmin, double &expmax,
                             double &ag_portion,
                             double &inject_ab_affinity, long &inject_ab_ASindex,
                             AffinitySpace &AS, ofstream &ana) {
   bindim = resolution;
   ag_index = ag;
   antibodies = new double[bindim + 2];
   k_on = new double[bindim + 2];
   k_off = new double[bindim + 2];
   for (int n = 0; n < bindim + 2; n++) {
      antibodies[n] = 0;
      k_on[n] = 0;
      k_off[n] = 0;
   }

   if (bindim > 0) {
      if (ag >= 0) {
         ana << "Install new antibody affinity bins for antigen " << ag_index << "\n";
         ana << "k_offs are derived from the range of dissociation constants:\n";
      }
      for (int t = 0; t <= bindim; t++) {
         k_on[t] = kon;
         // scale k_off with affinities between 5*10^5 and 5*10^10
         // antibodies_k_off[t]=antibodies_k_on[t]/pow(10,5.5+0.5*t); // in 1/hr
         double expsteps = (expmax - expmin) / bindim;
         k_off[t] = k_on[t] / pow(10, expmin + expsteps * t);    // in 1/hr
         if (ag >= 0) {
            ana << "bin=" << t << " --> k_off=" << k_off[t] / 3600.
                << "/s --> K=" << k_on[t] / k_off[t] << "/M\n";
         }
         k_on[t] *= ag_portion;
      }
   }

   // define ab injection affinity bin:
   if (inject_ab_ASindex < 0) {
      // use inject_ab_affinity
      ab_injection_bin = AntibodyDyn::getBinFromRealAffinity(inject_ab_affinity);
   } else {
      // use inject_ab_ASindex
     double injection_affinity = 0;
     if (ag_index >= 0) { injection_affinity = AS.affinity(inject_ab_ASindex,ag_index); }
     ab_injection_bin = AntibodyDyn::getBinFromAbstractAffinity(injection_affinity);
   }

   if (ag >= 0) {
      // open an ofstream for each instantiation
      char datname[50] = "";
      get_filename(datname);
      ofstream ab_file(datname);
      ab_file << "! Time course of homogeneous antibody concentration in GC in M, relative to Ag"
              << ag_index << "\n"
              << "! time . total . bins ... exogenous . ha-frac-without-exo\n";
      ab_file.close();
   }
}
void AntibodyBin::clear_bins() {
   for (int n = 0; n <= bindim + 1; n++) {
      antibodies[n] = 0;
   }
}
void AntibodyBin::set_antigen(int index) {
   ag_index = index;
}
int AntibodyBin::get_antigen() {
   return ag_index;
}
void AntibodyBin::attributeABfromAS(AffinitySpace &AS) {
   /* This routine runs through all Ab producing cell types in AffinitySpace,
    * collects the Ab saved at this position, and attributes it to the
    * corresponding Ab-affinity-bin.
    * Old values on the Ab-bins are overwritten, thus, this routine
    * is made to initialise Ab-bins with Ab as saved in AffinitySpace.
    */
   // get the number of different Ab-producing cells
   int n_producers = AS.get_n_ab_producers();
   // index in AffinitySpace
   long nASindex;
   double Ab_amount;
   clear_bins();
   // go through all types of Ab-producting cells:
   for (int n = 0; n < n_producers; n++) {
      // get the Index of Ab in AffinitySpace associated with Ab-producing cell type <n>
      nASindex = AS.get_AbProducingIndex(n);
      // get the amount of Ab of this type
      Ab_amount = AS.get_AbAmount(nASindex);
      // put this amount to the respective Ab-bins for the new Ab-bin
      production(nASindex,Ab_amount,AS);
   }
}
void AntibodyBin::get_filename(char * datname) {
   strcat(datname, "antibody_ag");
   stringstream itmp;
   itmp << ag_index;
   string stmp = itmp.str();
   strcat(datname, (char*) stmp.c_str());
   strcat(datname, ".out");
}
/*
 * Not needed anymore, as the degradation of all bins is included
 * in AntibodyDyn::produce_antibodies_outside(...).
 * Delete if this holds true for a while.
 */
/*
 * void AntibodyBin::degradation(double degrad_rate) {
 * double d_ab;
 * for (int n = 0; n <= bindim + 1; n++) {
 *  d_ab = antibodies[n] * degrad_rate;
 *  antibodies[n] -= d_ab;
 *  if (antibodies[n] < 0.) antibodies[n] = 0.;
 * }
 * // ... note that the last bin calculates the injected antibodies separately
 * }
 */

void AntibodyBin::production(long n, double &add, AffinitySpace &AS) {
   /* Adds <add> produced/degraded for antibody index <n> on AffinitySpace <AS>
    * add<0 is possible if only degradation is done.
    */
   double aff = AS.affinity(n,ag_index);
   int bin = AntibodyDyn::getBinFromAbstractAffinity(aff);
   antibodies[bin] += add;
}
double AntibodyBin::get_total_antibodies() {
   double t = 0.;
   for (int j = 0; j <= bindim; j++) {
      t += (1.e+15 * antibodies[j]);   // sum in M
   }
   return t;
}
double AntibodyBin::average_ab_affinity() {
   double sum = 0.;
   double average = 0.;
   // bins of different affinity are [0,..,bindim]
   for (int i = 0; i <= bindim; i++) {
      sum += antibodies[i];
      average += antibodies[i] * double (i);
   }
   if (sum > 0.) { average /= sum; } else { average = 0.; }
   /* now average contains the average bin as a double
    * this value has to be transformed into affinity space
    * (see void produce_antibodies_outside(...) below for the definition)
    */
   average /= double (bindim);
   //  cout<<average<<", ";
   return average;
}
void AntibodyBin::do_inject_antibody(double &amount, ofstream &ana) {
   antibodies[ab_injection_bin] += amount;
   antibodies[bindim + 1] += amount;
   cout << "\n"
        << "Inject soluble antibodies: " << amount << "mol/micrometer^3 into affinity bin "
        << ab_injection_bin << " of antigen " << ag_index << ".\n";
   ana << "\n"
       << "Inject soluble antibodies: " << amount << "mol/micrometer^3 into affinity bin "
       << ab_injection_bin << " of antigen " << ag_index << ".\n";
}
void AntibodyBin::show_antibodies(double &time) {
   cout << "Ag=" << ag_index << " t=" << time << "   ";
   for (int n = 0; n <= bindim; n++) {
      cout << antibodies[n] << "   ";
   }
   cout << ", exogenous: " << antibodies[bindim + 1]  // exogenous (or injected antibodies)
        << "\n";
}
void AntibodyBin::write_antibodies(double &time) {
   // open file to append
   char datname[50] = "";
   get_filename(datname);
   ofstream ab_file;
   ab_file.open(datname, ofstream::app);
   // add things
   double total_ab = get_total_antibodies();
   ab_file << time << "   " << total_ab << " ";   // in M
   for (int i = 0; i <= bindim + 1; i++) {
      ab_file << "  " << 1.e+15 * antibodies[i];   // in M
      // the last bin is the exogenous (injected) antibody concentration
   }

   // Now calculate the fraction of high affinity antibodies excluding injected antibodies:
   // ... first define the limit between low and high affinity at 20% // ### could be other value
   int lowafflimit = int (bindim * 0.2);
   // ... add all antibodies above this limit:
   double highaff = 0.;
   for (int i = lowafflimit + 1; i <= bindim; i++) {
      highaff += antibodies[i];
   }
   // remove exogenous antibodies:
   if ((ab_injection_bin >= 0) && (antibodies[bindim + 1] > 0)) {
      // antibody was injected before
      // cout<<"\n\n in mk_cell_sum in if ab_injection_bin>=0 && ... \n\n";
      // reduce the total antibody by the injected one
      total_ab -= 1.e+15 * antibodies[bindim + 1];
      if (ab_injection_bin > lowafflimit) {
         // also reduce the high affinity antibodies by the injected one
         highaff -= antibodies[bindim + 1];
      }
   }
   if (total_ab > 0) {
      ab_file << "   " << 1.e+15 * highaff / total_ab;
   } else {
      ab_file << "   0.0";
   }
   ab_file << "\n";
   // close the file
   ab_file.close();
}
/*******************************************************************************************
*******************************************************************************************/

AntibodyDyn::AntibodyDyn() { }
AntibodyDyn::AntibodyDyn(Parameter &par, AffinitySpace &AS, ofstream &ana) {
   ab_bins.reserve(100);
   ab_bins.clear();
   check4new_antigen(AS,ana);
   if (ab_injection_time >= 0) {
      ask2inject_antibody(par.Value.tmin,par.Value.deltat,ana);
   }
}
void AntibodyDyn::set_statics(const Parameter &par, ofstream &ana) {
   // this is about real antibody production:
   // antibody production rate (unspecific) value from Randall 1992: 2000 antibodies per cell and
   // second
   antibody_production_factor = 0;
   if (par.Value.antibodies_production > 0) {
      antibody_production_factor = par.Value.antibodies_production * par.Value.deltat
                                   * par.Value.N_GC
                                   / (par.Value.V_blood * 1.e+15);
   }
   // unit mol/micrometer^3 per cell and hour (note for V_blood in liter: 1.e+13 micrometer^3, i.e.
   // 10^-2 liter)
   ana << "  Effective Ab production factor = " << antibody_production_factor
       << " mol/micrometer^3\n";

   // antibody degradation
   antibody_degradation = 0;
   // half value is 30 days: prob=log(2.)*timestep(0.002hr)/720.hr
   if (par.Value.antibodies_degradation > 0) {
      antibody_degradation = log(2.) * par.Value.deltat
                             / (24. * par.Value.antibodies_degradation);
   } else {
      // 2.e-06 = half value is 30 days: prob=log(2.)*timestep(0.002hr)/720.hr
      antibody_degradation = 0.;
   }
   ana << "  antibody degradation = " << antibody_degradation
       << " event probability per time step\n";

   // about affinity bins of Abs:
   // resolution of antibody in bins of affinities between 0 and 1
   antibodies_resolution = par.Value.antibodies_resolution;
   ag_portion = par.Value.ag_threshold;
   // note that k_on is now the same for all sets of antibody affinity binning
   k_on = 3600. * par.Value.ic_k_on;  // 3600.*(typical value=1.e+06); // in 1/(M hr)
   k_ic_exp_min = par.Value.k_ic_exp_min;
   k_ic_exp_max = par.Value.k_ic_exp_max;

   // about Ab injection
   ab_injection_time = -1.;
   ab_injection_amount = 0.;
   inject_ab_affinity = par.Value.injected_antibody_affinity;
   inject_ab_ASindex = par.Value.inject_antibody_ASindex;
   if (par.Value.inject_antibody > 0.) {
      ab_injection_time = par.Value.inject_antibody_time;
      ab_injection_amount = par.Value.inject_antibody * 1.e-15;   // mol/l --> mol/micrometer^3
   }
}
AntibodyDyn::~AntibodyDyn() {
   ///Philippe //Michael: why is this needed, vectors are destroyed anyway at the end?
   ab_bins.clear();
}
int AntibodyDyn::getBinFromRealAffinity(double realAffinity) {
   int bin
      = int (antibodies_resolution * (realAffinity - k_ic_exp_min)
             / (k_ic_exp_max - k_ic_exp_min));
   if (bin > antibodies_resolution) { bin = antibodies_resolution; }
   if (bin < 0) { bin = 0; }
   return bin;
}
int AntibodyDyn::getBinFromAbstractAffinity(double abstractAffinity) {
   return int (abstractAffinity * antibodies_resolution);
   /* This defines a map from affinity space to bin space:
    * aff --> bin=int(aff*bindim)
    * inducing the map
    * [0,..,0.1[   --> 0
    * [0.1,..,0.2[ --> 1
    * ...
    * [0.9,..,1.0[ --> 9
    * [1.0]        -->10
    * The corresponding back transformation is then defined by:
    * bin --> aff=double(bin)/double(bindim)
    * This maps the bins to the lower boundary of the respective intervals.
    * The back transformation is used in double& cellman::average_ab_affinity();
    */
}
void AntibodyDyn::check4new_antigen(AffinitySpace &AS, ofstream &ana) {
   /* This routine can be used to scan for the existence of new Ag.
    * This may be at the beginning of a simulation at the time of initialisation
    * or later in the course of the simulation.
    * If a new Ag is detected a new Ab-bin is installed and loaded
    * with Abs pre-existing on the AffinitySpace.
    */
   int n_ags = AS.get_n_Antigen();
   int n_abs = ab_bins.size();
   // If there are more antigens than antibody bin instantiations extend:
   while (n_ags > n_abs) {
      // get the index in AffinitySpace of the next new
      int i_ag = AS.get_Antigen(n_abs);
      AntibodyBin ab_tmp;
      ab_tmp.initialise(antibodies_resolution, i_ag,
                        k_on, k_ic_exp_min, k_ic_exp_max, ag_portion,
                        inject_ab_affinity, inject_ab_ASindex,
                        AS, ana);
      ab_bins.push_back(ab_tmp);
      /* In the case of pre-existing Ab on the AffinitySpace,
       * Ab are attributed to the new set of Ab-bins */
      ab_bins[n_abs].attributeABfromAS(AS);
      ++n_abs;
   }
}
bool AntibodyDyn::check4consistency(AffinitySpace &AS, ofstream &ana) {
   /* Used to compare the content of affinity bins and the amount of
    * antibody stored in AffinitySpace (done for all antigens).
    * Returns true if fine and false if deviation is above threshold.
    * If deviation is above threshold, antibody from AffinitySpace
    * is reattributed to the bins and this event is reported to the logfile.
    */
   cout << "Check for consistency of antibodies on AffinitySpace and in bins ... ";
   bool allfine = true;
   AntibodyBin ab_tmp;
   int i_ag = -1;
   ab_tmp.initialise(antibodies_resolution, i_ag,
                     k_on, k_ic_exp_min, k_ic_exp_max, ag_portion,
                     inject_ab_affinity, inject_ab_ASindex,
                     AS, ana);
   int n_ags = ab_bins.size();
   bool fine;
   for (int a = 0; a < n_ags; a++) {
      fine = true;
      ab_tmp.set_antigen(ab_bins[a].get_antigen());
      ab_tmp.attributeABfromAS(AS);
      // compare ab_tmp to ab_bins[a]
      for (int bin = 0; bin <= antibodies_resolution; bin++) {
         if (ab_tmp.antibodies[bin] - ab_bins[a].antibodies[bin] > consistency_threshold) {
            cerr << "true[" << bin << "]: " << ab_tmp.antibodies[bin]
                 << " bins:" << ab_bins[a].antibodies[bin] << "\n";
            fine = false;
            allfine = false;
         }
      }
      if (not (fine)) {
         for (int bin = 0; bin <= antibodies_resolution; bin++) {
            ab_bins[a].antibodies[bin] = ab_tmp.antibodies[bin];
         }
         ana << "WARNING! "
             << "Reloaded antibodies from AffinitySpace to bins for antigen #"
             << a << " with index " << ab_bins[a].get_antigen() << "\n";
         cout << "WARNING! "
              << "Reloaded antibodies from AffinitySpace to bins for antigen #"
              << a << " with index " << ab_bins[a].get_antigen() << "\n";
      }
   }
   cout << "done.\n";
   return allfine;
}
double AntibodyDyn::get_total_antibodies() {
   double t = 0.;
   // As all elements of ab_bins should contain the same total amount of Abs, just calculate the
   // first:
   if (ab_bins.size() > 0) { t = ab_bins[0].get_total_antibodies(); }
   return t;
}
double AntibodyDyn::average_ab_affinity(int bin_i) {
   /* returns the average affinity of antibodies for antigen bin_i which is also the
    * antigen index on the antigen list in AffinitySpace by construction */
   if (bin_i < 0) { return 1; }
   return ab_bins[bin_i].average_ab_affinity();
}
void AntibodyDyn::ask2inject_antibody(double time, double dt, ofstream &ana) {
   if ((ab_injection_amount >= 0)
       && (time >= ab_injection_time - 0.5 * dt)
       && (time < ab_injection_time + 0.5 * dt)) {
      for (unsigned int a = 0; a < ab_bins.size(); a++) {
         ab_bins[a].do_inject_antibody(ab_injection_amount, ana);
      }
   }
}
void AntibodyDyn::produce_antibodies_outside(AffinitySpace &AS) {
/* This routine is modifying both, the antibody amount in AffinitySpace
 * and its representation in bins within the AntibodyDyn class.
 * Ab is produced according to all output cells
 * that differentiated to the Ab-producing phenotype,
 * which is known from AffinitySpace.
 * Ab and Ab-bins are updated by production and degradation.
 */
   double d_ab = 0.;
   /* Note that the variable as defined in cellman is divided by <ag_threshold>,
    * thus transforming mol to number of ag-portions!?
    */
   // get the number of different Ab-producing cells
   int n_producers = AS.get_n_ab_producers();
   // index in AffinitySpace
   long nASindex;
   // go through all types of Ab-producing cells:
   for (int n = 0; n < n_producers; n++) {
      // get the Index of Ab in AffinitySpace associated with Ab-producing cell type <n>
      nASindex = AS.get_AbProducingIndex(n);
      // calculate the change of Ab-amount for this Ab (production and degradation)
      d_ab = AS.get_AbProducingCells(n) * antibody_production_factor
             - antibody_degradation * AS.get_AbAmount(nASindex);
      // add this amount on AffinitySpace
      AS.put_Ab(nASindex,d_ab);
      // add this amount to the respective Ab-bins for each antigen
      for (unsigned int a = 0; a < ab_bins.size(); a++) {
         ab_bins[a].production(nASindex,d_ab,AS);
      }
   }
}
void AntibodyDyn::show_antibodies(double time) {
   for (unsigned int a = 0; a < ab_bins.size(); a++) {
      ab_bins[a].show_antibodies(time);
   }
}
void AntibodyDyn::write_antibodies(double time) {
   for (unsigned int a = 0; a < ab_bins.size(); a++) {
      ab_bins[a].write_antibodies(time);
   }
}
