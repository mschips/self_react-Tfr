/**
 * @file kinetics.h
 * @brief class to store information about the GC volume over simulation only at the following
 * time-points : 24, 48, 72, 96, 120, 144, 168, 240, 336, 408, and 504 hours
 *
 * independent of other classes, starts empty, and waits for data being given to it
 * by calling check_GCvolume(const long& vol, const double& t);
 * if t = one of the times points to record (with tolerance [+0 .. +dt]), then it stores the volume
 * if not, it does not store anything.
 *
 * allows to compare with experimental data : liu91 (days 1,2,3,4,6,7,14 21) and hol92 (days 3, 5
 *,10 ,17)
 * and to give the chi2 weights for each time-point
 *
 * @author Michael Meyer-Hermann
 */

#ifndef i_kinetics
#define i_kinetics

#include <iostream>
#include <fstream>
using namespace std;

/**
 * @brief Class for storing information about GC volume over time
 */

class GCkinetics {
  public:
   /**
    * @brief Empty constructor
    * @return GCkinetics (empty) instance
    */
   GCkinetics() { }

   /**
    * @brief Constructor
    * @param par_dt (double) maximum time delay between a roundish time-point and the next
    *simulation
    * time-point
    * @return GCkinetics instance
    *
    * Initial values :
    *    dt = given par_dt			(interval between time-points)
    *    n_vol_max = 20				(allocated size of tables)
    *    n_vol_values= 0				(nb of data-points saved)
    *    GCvolume_t[n_vol_max] = {0, ... 0}      (table of volumes)
    *    t_GCvolume[n_vol_max] = {0, ... 0}	(table of time-points)
    *    max_GCvolume = 0            (max recorded volume)
    */
   GCkinetics(const double &par_dt);

   /**
    * @brief Destructor
    */
   ~GCkinetics() { }

   /**
    * @brief Gives the volume for this time-point (only recorded if it is one of the desired
    * time-points)
    * @param vol Germinal Center volume
    * @param t time
    */
   /// updates max_GCvolume in case vol is bigger than the previous max_GCvolume
   /// and saves vol when time is appropriate:
   void check_GCvolume(const long &vol, const double &t);

   /**
    * @brief Saves the lists of (time, GCvolume, chi2) into "vol_chi2.out" and returns the total
    * sigma-squared from chi2 test for the kinetics
    */
   double get_sigma2_GCvolume();

  private:
   /// Dimension of the allocated volume-array
   static const int n_vol_max = 20;

   /// For calculating the chi2 sigma between simulation and data,
   /// Use division by n (0.) or n-1 (1.)
   static constexpr double use_sigma_n_1 = 0.;

   /// Number of saved volume-values
   int n_vol_values;

   /// Monitores the maximum reached (recorded yet) GC volume
   long max_GCvolume;

   /// Array for the saved volume values
   long GCvolume_t[n_vol_max];

   /// Array with corresponding time points in hours
   /// because they lie between official_time_point and official_point + dt (ex, between 24h and
   // 24h+dt),7
   /// so records the exact value of this time-point
   double t_GCvolume[n_vol_max];

   /// time step width of calling routine (tolerance interval for rounding time-points)
   double dt;

   /// these two variables are computed by get_chi2, and later used by get_sigma2_GCvolume()
   /// sum of chi-squared:
   double chi2_sum;

   /// number of evaluated values
   int n_chi2_values;

   /**
    * @brief Experimental values (days 1,2,3,4,6,7,14 21)
    * @param day time-point of interest
    */
   double liu91(const int &day);

   /**
    * @brief Experimental values (days 3, 5 ,10 ,17)
    * @param day time-point of interest
    */
   double hol92(const int &day);

   /**
    * @brief Adds the chi2 of this time-point to chi2_sum
    * @param a index of the data point in the tables
    *
    * for each data point, checks to which day it corresponds to, computes the chi2 by calling
    * get_chi2(const int& a, const short& liu, const int& day), and adds it to chi2_sum
    */
   double get_chi2(const int &a);

   /**
    * @brief (called by get_chi2()) computes the chi2 of this time-point, from liu's data (liu = 1)
    * or hol (liu != 1)
    * @param a index of the data point in the tables
    * @param liu set to 1 to request data from liu,
    * @param day the equivalent time-point (in days).
    */
   double get_chi2(const int &a, const short &liu, const int &day);

   /// Saves vol at the next free point in GCvolume_t[] and increases n_vol_values
   void save_GCvolume(const long &vol, const double &t);
};

#endif
