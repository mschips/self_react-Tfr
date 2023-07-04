#include "tamoxifen.h"
#include <math.h>
//#include "random.h"
//#include <string.h>
//#include <cstdlib>

void tamoxifen::set_parameter(Werte& p) {
  // Activity state of tmx: switches on or off the whole tmx class
  action = p.tamoxifen_do;
  // The time point of tmx injection
  t_inject = p.tamoxifen_t_inject;
  // The simulation timestep
  sim_dt = p.deltat / 2.;
  // Fraction of recombined cells in the case of one shot tmx activity (no decay)
  recombine_fraction = p.tamoxifen_recombine_fraction;

  // The half-life to tmx activity
  t_decay = p.tamoxifen_t_decay;
  // A local variable saving the next time of action
  t_next = t_inject;
  // The time point of stop of tmx activity in the case of tmx decay (not a shot)
  t_stop = p.tamoxifen_t_stop;
  // Frequency (time interval) for recombinaion (avoid doing this in every simulation time step)
  recombination_dt = p.tamoxifen_recombination_dt;
  // Target recombined fraction of cells in the case of tmx decay
  amplitude = get_recombination_amplitude();

  // Actions of tamoxifen-induced recombination:
  MHC = p.tmx_MHC;
}

double tamoxifen::get_recombination_amplitude() {
   /* This routine will calculate the probability with exponential decay
    * Relies on the following global variables:
    * uses t_decay: if <0 just return stain_fraction,
    *               use it as half life time in hr otherwise
    * uses recombine_dt as time step for transformation of a rate to a probability
    * uses recombine_fraction as default response or as integrated fraction of stained cells
    */
   if (t_decay < 0) { return recombine_fraction; }
   /* We have:
    * fraction = p_0 \int_0^\infty exp(-t/t_decay) dt = p_0 t_decay
    * For a finite interval the situation is more involved (t_max=t_stop):
    * fraction = p_0 \int_0^{t_max} exp(-t/t_decay) dt
               = p_0 t_decay (1 - exp(-t_max/t_decay) )
    * solved for p_0 this yields:
    * p_0 = (fraction/t_decay) / (1 - exp(-t_max/t_decay) )
    *
    * Note that t_stop is an absolute time parameter.
    * Here it was assumed relative to zero, thus t_stop has to be reset relative to t_inject,
    * because t_inject is the equivalent of zero in this derivation.
    */
   double p_0 = 
     (recombine_fraction / t_decay) / (1.0 - exp(-1.0 * (t_stop - t_inject) / t_decay));
   // p_0 is a rate per hour --> the probability per time step is generated as p_0 recombine_dt:
   p_0 *= recombination_dt; 
   // note that this is a rate, not an inverse half time, so no log(2) needed here.
   return p_0;
}

/* What follows are routines that can be called globally:
 * ======================================================
 */
bool tamoxifen::get_tamoxifen_action(double &time) {
  if (not(action)) { return false; }
  if (t_decay < 0) { 
    if (time >= t_inject - sim_dt && time < t_inject + sim_dt) { return true; }
  } else {
    if (time >= t_next - sim_dt && time < t_next + sim_dt) { 
      // Determine the next time step for tmx action
      t_next += recombination_dt;
      // If it is beyond t_stop, prevent any further tmx action
      if (t_next > t_stop) { t_next = 0; action = false; }
      return true; 
    }
  }
  return false;
}
double tamoxifen::get_recombination_probability(double &time) {
   /* The decay of <amplitude> with an exponential is calculated
    * time is interpreted relative to t_inject (global)
    * the decay constant is t_decay (global)
    */
   if (t_decay < 0) { return recombine_fraction; }
   double value = amplitude * exp((t_inject - time) * log(2) / t_decay);
   return value;
}
double tamoxifen::get_recombination_probability(double& time, long R, long N) {
  double p = get_recombination_probability(time);
  /* The probability returned from get_recombination_probability(double&)
   * assumes that no cell is recombined twice. This implies that the probability
   * corresponds to the fraction of non-recombined cells that have to be
   * recombined in this step. Thus, the actual probability to recombine per cell
   * is higher than the returned value in get_recombination_probability.
   * This is particularly true for high fractions of recombination,
   * when the fraction of non-recombined cells gets small.
   * The probability is rescaled here to account for this under-estimation.
   */
  if (t_decay < 0) { return p; }
  /* Given a total (Cre+) population of N cells, 
   * and R cells being already recombined,
   * and the probability p as derived from get_recombination_probability(double&).
   * Then, the number of cells to be recombined in this time step is pN.
   * These have to be taken from N-R cells, i.e. the not yet recombined cells.
   * Thus, the rescaled probability is p* = pN/(N-R).
   */
  double pstar;
  if (N > R) { pstar = (p * double(N) / double(N-R)); }
  else { pstar = 0.; }
  /* Note that in the limit of R=0 we have p*=p,
   * and in the limit of R=N, p*->\infty.
   * Note also that this is over-estimating the recombination probability
   * when N is growing for other reasons, so for highly dynamic N.
   */
  if(pstar > 1) { pstar = 1.0; }
  return pstar;
}


