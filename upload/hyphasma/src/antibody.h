#ifndef i_ANTIBODYDYN
#define i_ANTIBODYDYN

#include "setparam.h"
#include "affinityspace.h"

/** mmh February 3rd, 2016
 * The following functionality is expected from the antibody class:
 * It allows for antibody production in general (irrespective of any antigen, just derived from PC).
 * --> produce_antibodies_outside() is in the antigen-global class
 * It allows for antibody degradation in general, just every antibody is degraded.
 * --> degradation is part of produce_antibodies_outside()
 * --> ag_portion, antibody_degradation, antibody_production_factor are all in the antigen-global
 * class
 * It allows for antibody injection with a particular affinity for a specified antigen,
 * but while the definition of the antigen determines the type of antibody to be injected,
 * this antibody might also interfere with other antigens. This has to be organised from a
 * meta-level.
 * --> ask2inject_antibody() and do_inject_antibodies() are in the antigen-global class
 * --> ab_injection_xxx variables are in the antigen-global class
 * It allows to return the total antibody and total antibody for each antigen.
 * --> get_total_antibodies() is in the antigen-global class and is the same with respect to any
 * antigen.
 * It allows to return the average_ab_affinity with respect to any specified antigen.
 * --> average_ab_affinity() should be in the antigen-global class and return a vector for all
 * antigens.
 * It allows to generate read-out in out-files. The number of files generated should correspond to
 * the number of antigens. This has to be organised in the global class.
 * --> show_antibodies() is in the antigen-global class.
 * It allows analysis of antibodies for each antigen with discretised affinity classes.
 * Every antigen gets its own binning.
 * --> antibodies[], k_on[], k_off[], ag_index are in an antigen-specific class.
 **/

/** @brief Container class, for an antigen, to store the amount of antibodies produced against it,
 * depending on their affinity */
class AntibodyBin {
  public:
   AntibodyBin();
   ~AntibodyBin();

   /** @brief Initializer function to creates empty bins :
    *
    *  @param resolution     (int    -> bindim)
    *         resolutions of affinity classes for the antibodies     (total number of bins :
    * resolution + 2).
    *  @param ag             (int    -> ag_index)
    *         index of the antigen in the affinity space             (to store it)
    *
    *  parameters to give k_on[] and k_off[] to the affinity bins :
    *  @param kon            (double)  k_on for all antibodies before correction by ag_portion
    *  @param koff           (double)  maximum k_off (for the best bin) in 1/hr
    *  @param expmin         (double)  kon/10^expmin for the best affinity bin
    *  @param expmax         (double)  kon/10^expmax for the worst affinity bin
    *  @param ag_portion     (double)
    *  formula : for each bin t in [|0 .. resolution|]. the bin nr resolution+1 has k_on = k_off =
    *0.
    *       k_on[t]=kon * ag_portion;
    *       k_off[t] = kon / pow(10, expmin + t * (expmax - expmin) / resolution);
    *  formula from real affinity to bin :
    *       bin = int(resolution * (affinity - expmin) / (expmax - expmin));
    *  formula from abstract affinity (between 0 and 1) into bin :
    *       bin = int(aff * resolution);
    */

   void initialise(int &resolution, int &ag,  /* double& ab_degradation,*/
                   double &kon, double &expmin, double &expmax, double &ag_portion,
                   double &injected_ab_affinity, long &inject_ab_ASindex,
                   AffinitySpace &AS, ofstream &ana);

   /** @brief Storage fields : */
   // ### make these private and shift cellFDC::mk_immune_complex into this class
   double * antibodies; // used for antibody affinity bins
   double * k_on;
   double * k_off;
   int ab_injection_bin;  // ### make this private when AntibodyDyn::show_antibodies() was deleted.

   /** @brief Clear all bins */
   void clear_bins();

   /** @brief Attribute antibodies from AffinitySpace from scratch */
   void attributeABfromAS(AffinitySpace &AS);

   /** @brief Attributes/gets an antigen to the Ab-bins */
   void set_antigen(int index);
   int get_antigen();

   /** @brief Apply degradation for a period of dt according to degrad_rate to all the bins */
   // void degradation(double degrad_rate);

   /** @brief Takes the antibody at position n, finds its affinity to this antigen
    *  then deduces in which bin this antibody is, and increase by 'add' its amount */
   void production(long n, double &add, AffinitySpace &AS);

   /** @brief Injects an amount 'amount' of antibodies in the bin 'injection bin'*/
   void do_inject_antibody(double &amount, ofstream &ana);

   /** @brief Output functions */
   double get_total_antibodies();
   double average_ab_affinity();
   void show_antibodies(double &time);
   void get_filename(char * datname);
   void write_antibodies(double &time);

  private:
   /** @brief Index of this antigen in the AffinitySpace */
   int ag_index;
   /** @brief Antibodies_resolution : the number of bins is nibdim+2 */
   int bindim;
};

/** @brief Class containing the AntibodyBins associated with each antigen. Each time a new antigen
 * is created in the Affinity Space, a new AntibodyBin is created for it.
 */
class AntibodyDyn {
  public:
   // Konstruktor
   AntibodyDyn();
   AntibodyDyn(Parameter &par, AffinitySpace &AS, ofstream &ana);
   // Destruktor
   ~AntibodyDyn();

   /** @brief Main storage : 
       ab_bins[i] = AntibodyBin for the antigen nr i in the Affinity Space */
   vector<AntibodyBin> ab_bins;   // ### make this private and particular parts accessible from
                                  // cellFDC::mk_immune_complex

   /** @brief Checks from Affinity space if there are new antigens.
    *  This has to be called from external to generated sets of affinity-bins */
   void check4new_antigen(AffinitySpace &AS, ofstream &ana);

   /** @brief Checks consistency between affinity bins and antibody on AffinitySpace
    * for each antigen */
   bool check4consistency(AffinitySpace &AS, ofstream &ana);

   /** @brief Major function : updates the amount of the antibodies in each bin of each antigen
    * according to their amount of producers from the AffinitySpace */
   void produce_antibodies_outside(AffinitySpace &AS);

   /** @brief checks if antibodies shoulg be injected at this particular time (+/- 0.5 dt)
    *  and then calls do_inject_antibody for each antigen AntibodyBin */
   void ask2inject_antibody(double time, double dt, ofstream &ana);

   /** @brief Functions to show what's happening */
   // show antibody bins for all antigens on the screen
   void show_antibodies(double time);
   double get_total_antibodies();  // returns a value in M
   double average_ab_affinity(int bin_i);
   // writes antibodies to all antigens into files and returns total amount of antibodies
   void write_antibodies(double time);

   ///Philippe : add this functions and put the parameters in static
   static void set_statics(const Parameter &par, ofstream &ana);

   /** @brief gives the bin of an antibody from its affinity in the real scale */
   static int getBinFromRealAffinity(double realAffinity);
   /** @brief gives the bin of an antibody from its affinity in the affinity space (between 0 and
    * 1)*/
   static int getBinFromAbstractAffinity(double abstractAffinity);

   /** static parameters */
   static int antibodies_resolution;
   static double ag_portion;
   static double antibody_degradation;
   static double antibody_production_factor;
   static double k_on,k_ic_exp_min,k_ic_exp_max;

   static double inject_ab_affinity;
   static long inject_ab_ASindex;
   static double ab_injection_time, ab_injection_amount;
   static double consistency_threshold;
};

#endif
