/* August 16, 2018 Michael Meyer-Hermann:
 * Class for induction of Cre-recombination by tamoxifen.
 *
 * Initialisation:
 * Define a tamoxifen object and call tamoxifen::set_parameter(Werte& p)
 * (requires an object of class Werte).
 *
 * Usage:
 * Check with tamoxifen::get_tamoxfien_action(time) whether a Cre-recombination shall be induced.
 * If yes, calculate the recombination probability with
 * tamoxifen::get_recombination_probability(time)
 * Then attribute the Cre-dependent property to cell objects with this probability.
 */

#ifndef i_TAMOXIFEN
#define i_TAMOXIFEN

//#include <vector>
//#include <iostream>
//#include <fstream>
//#include <math.h>
#include <setparam.h>
//using namespace std;

class tamoxifen {
 public:
  tamoxifen() {}
  ~tamoxifen() {}
  void set_parameter(Werte& p);
  bool get_tamoxifen_action(double &time);
  double get_recombination_probability(double& time);
  double get_recombination_probability(double& time, long R, long N);
  bool MHC;
 private:
  bool action;
  double t_inject, sim_dt, recombine_fraction;
  double t_decay, t_next, t_stop, recombination_dt, amplitude;
  double get_recombination_amplitude();
};

#endif
