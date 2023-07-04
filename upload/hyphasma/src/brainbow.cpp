#include "brainbow.h"
#include "random.h"
#include <string.h>
#include <cstdlib>

brainbow::brainbow() {
   do_clonality = false;
   do_clone_fractions = false;
   do_clone_fractions_rand = false;
   trackfrom = 0;
   trackuntil = -1;
   Ncolours = 10;
   cell_subset_size = 0;
   stain_fraction = 0.4;
   tamoxifen_decay = -1;
   include_late_founders_in_lineages = false;
   include_seclargest_colour = false;
   rb_dim_ini = 1e+05;
   rb_list.reserve(rb_dim_ini);
   merge_identical_origin = false;
   origin_same_colour = false;
   allow_restaining = false;
   dt_Muller = -1;
}
brainbow::brainbow(double from, double until) {
   do_clonality = false;
   do_clone_fractions = false;
   do_clone_fractions_rand = false;
   trackfrom = from;
   trackuntil = until;
   Ncolours = 10;
   cell_subset_size = 0;
   stain_fraction = 0.4;
   tamoxifen_decay = -1;
   include_late_founders_in_lineages = false;
   include_seclargest_colour = false;
   rb_dim_ini = 1e+06;
   rb_list.reserve(rb_dim_ini);
   merge_identical_origin = false;
   origin_same_colour = false;
   allow_restaining = false;
   dt_Muller = -1;
}
brainbow::~brainbow() { }
void brainbow::set_Ncolour(int c) { Ncolours = c; }
void brainbow::set_stain_time(double t) { stain_time = t; }
void brainbow::set_stain_fraction(double f) { stain_fraction = f; }
void brainbow::set_stain_dt(double t) { stain_dt = t; }
void brainbow::set_lineage_time(double t) { lineage_time = t; }
void brainbow::set_tamoxifen_decay(double t) { tamoxifen_decay = t; }
void brainbow::set_include_late_founders_in_lineages(bool t) 
{ include_late_founders_in_lineages = t; }
void brainbow::set_include_seclargest_colour(bool t)
{ include_seclargest_colour = t; }
void brainbow::set_tamoxifen_stop_time(double t) { tamoxifen_stop_time = t; }
void brainbow::set_staining_mode(staining_modes b) {
   how2do_staining = b;
   if (((how2do_staining == gabriel_4colours) && (Ncolours != 4))
       || (((how2do_staining == gabriel_10colours)
            || (how2do_staining == gabriel_10colours_2founder))
           && (Ncolours != 10))) {
      cerr << "ERROR: In brainbow::set_staining_mode(staining_modes):\n"
           << "staining_mode=" << how2do_staining << "\n"
           << "incompatible with Ncolours=" << Ncolours << "\n"
           << "abort program.\n";
      exit(1);
   }
}
void brainbow::set_allow_restaining(bool b) { allow_restaining = b; }
void brainbow::set_origin_same_colour(bool b) { origin_same_colour = b; }
void brainbow::set_merge_identical_origin(bool b) { merge_identical_origin = b; }
void brainbow::set_cell_subset_size(int x) { cell_subset_size = x; }
void brainbow::set_do(bool x, bool y, bool z) {
   do_clonality = x;
   do_clone_fractions = y;
   do_clone_fractions_rand = z;
}
long brainbow::write(double time, bool founder, bool birth, long mother_index, long ss_position) {
   if ((time >= trackfrom) && ((trackuntil == -1) || (time <= trackuntil))) {
      unsigned long event_index = rb_list.size();
      if (event_index >= rb_list.capacity()) { rb_list.reserve(2 * rb_list.capacity()); }
      // cout<<"!!!!  "<<xtime<<"  "<<xfounder<<"  "<<xbirth<<"  "<<xmother_index<<"
      //  "<<xss_position<<"\n";
      brainbow_data x;
      x.time = time;
      x.founder = founder;
      x.lineage = false;   // is defined during analysis only
      x.birth = birth;
      if (founder) { x.mother_index = -1; } else { x.mother_index = mother_index; }
      x.daughter1_index = -1;
      x.daughter2_index = -1;
      x.ss_position = ss_position;
      x.colour = -1;
      rb_list.push_back(x);
      // Now the new event is introduced in the list.
      // In addition, the mother cell now gets a pointer to the daughter cell:
      if (not (founder)) {
         // founder cells don't have a mother cell
         if (rb_list[mother_index].daughter1_index == -1) {
            // write the pointer to this cell to daughter1_index if not yet used
            rb_list[mother_index].daughter1_index = event_index;
         } else {
            // but if used write to daughter2_index
            rb_list[mother_index].daughter2_index = event_index;
         }
         // as every cell has only two daughters, daughter2_index must be free when this point was
         // reached.
      }
      return event_index;
   }
   // else do nothing.
   return -1;
}
void brainbow::save_to_file() {
  /* Saves all events recorded in write(...) so far in the file braibow.out for later analysis.
   * Note that the variable lineage is not saved because it is defined during analysis.
   */
  ofstream fff("brainbow.out");
  // save all data in this file:
  for (unsigned long i = 0; i < rb_list.size(); i++) {
    fff << i << "   "
	<< rb_list[i].time << "\t"
	<< rb_list[i].founder << "\t"
	<< rb_list[i].birth << "\t"
	<< rb_list[i].mother_index << "\t"
	<< rb_list[i].daughter1_index << "\t"
	<< rb_list[i].daughter2_index << "\t"
	<< rb_list[i].ss_position << "\t"
	<< rb_list[i].colour;
    if (i < rb_list.size() - 1) { fff << "\n"; }
  }
  fff.close();
  cout << "\nWrote mutation raw data to brainbow.out.\n";
}
long brainbow::flist_index(long i, vector<long>& v, long max) {
  /** @brief: Returns the index in the flist with v[index] == i.
   ** Exclusively called from save_Muller_diagram and made for flist therein.
   **/
  long ii = 0;
  while ( v[ii] != i && ii < max) { ++ii; }
  if ( ii >= max ) { return -1; } // this means an error and should never happen.
  return ii;
}
void brainbow::set_dt_Muller(double dt) { dt_Muller = dt; }
void brainbow::set_Muller_maxtime(double maxt) { Muller_maxtime = maxt; }
void brainbow::save_Muller_diagram() {
  /** @brief: Muller diagrams are built from two files, attribute and population.
   ** "attribute" contains two columns with an identifier for a cell and its mother cell.
   ** Optionally, a color can be added, which is ignored here.
   ** Also, no mother will be given because Muller.plot generates a new color for
   ** each entry, even if a daughter of a mother. The aim here is to get one color
   ** for each founder clone. Thus, the attribute file will only contain founder cells,
   ** for which "NA" instead of the mother-identifier has to be given.
   ** "population" contains the identifier, the time, and the abundance.
   ** The identifier should be the same as in the attribute-file.
   ** The time is the time of any change of cell numbers, i.e. birth, death, new founder.
   ** Abundance is a relative abundance, so giving simply the number of cells stemming
   ** from a particular founder clone and providing these cell numbers for all founder
   ** clones in regular time steps is sufficient.
   **/
  cout << "Write phylogeny data for Muller-diagram ...\n";
  // find all founders (saved as index of appearance in rb_list):
  vector<long> flist;
  flist.reserve(1000);
  unsigned int Nf = get_all_founders(-1,flist);  // time=-1 infers that all founders are returned
  // Note that these indices are any numbers, not a sequence.
  // Generate a file with all founder cells using the rb_list indices of appearance as names:
  ofstream att("muller_attributes.out");
  for (unsigned long i = 0; i < Nf; i++) { att << flist[i] << "     NA\n"; }
  att.close();
  // Define a vector with associated affinities
  /*
  vector<int> affinities;
  affinities.reserve(Nf+1);
  for ( unsigned long j = 0; j < Nf; j++ ) { 
    affinities.push_back(); 
  }
  */
  // Define a vector with relative abundances:
  vector<double> abundance;
  abundance.reserve(Nf+1);
  for ( unsigned long j = 0; j < Nf; j++ ) { abundance.push_back(0); }
  double dtime = 1.0; // in hours; write abundances every time interval dtime
  double tt = dtime; // This is the first time for writing
  if ( Muller_maxtime < 0 ) {
    // find the latest time point in the rb_list:
    Muller_maxtime = rb_list[rb_list.size() - 1].time;
  }
  ofstream pop("muller_population.out");
  //  ofstream aff("muller_population_aff.out");
  // Go through all events and calculate the relative abundances:
  for (unsigned long i = 0; i < rb_list.size(); i++) {
    if ( tt <= Muller_maxtime ) {
      if ( rb_list[i].time > tt ) {
	// Write all abundancies of existing founder clones to the file
	for ( unsigned long j = 0; j < Nf; j++ ) {
	  if ( abundance[j] > 0 ) {
	    pop << flist[j] << "   " << tt/24. << "   " << abundance[j] << endl;
	    /*
	    // If day 10, write abundancies together with affinities:
	    if ( ( rb_list[i].time > 240. ) && ( rb_list[i].time <= ( 240. + dtime ) ) ) {
	      aff << flist[j] << "   " 
		  << tt/24. << "   " 
		  << abundance[j] << "   "
		  << affinities[j]
		  << endl;
	    }
	    */
	  }
	}
	tt += dtime;
      }
      if ( rb_list[i].founder == 1 ) {
	/* case of a new founder cell:
	 * initialize the abundance as 1 
	 * (abundancies are given as cell numbers stemming from a clone) */
        long iflist = flist_index(i,flist,Nf);
	abundance[iflist] += 1;
        //affinties[iflist] = rb_list[i].affinity;
      } else {
	if ( rb_list[i].birth == 0 ) {
	  /* case of cell death:
	   * reduce the abundance of the founder of this cell by 1.
	   * At first, find the index of the founder cell of the dying cell:
	   * Start from the known mother and recursively follow the mothers: */
	  long ifound = get_founder_cell_index(rb_list[i].mother_index);
	  abundance[flist_index(ifound,flist,Nf)] -= 1;
	} else {
	  /* case of cell birth:
	   * Cell birth happens exclusively by division. The mother cell
	   * is replaced by two daughter cells. Thus, the first appearance
	   * shall be counted as an additional cell, but the second born
	   * daugther should not be counted again. The simplest solution is
	   * to increment the abundance by two portions of 0.5. */
	  long ifound = get_founder_cell_index(rb_list[i].mother_index);
	  abundance[flist_index(ifound,flist,Nf)] += 0.5;
	}
      }
    }
  }
  pop.close();
  //aff.close();
  cout << "\n ... generated muller_{attribute,population,population_aff}.out.\n";
}  
void brainbow::read_from_file() {
  ifstream fff("brainbow.out");
  long i;
  // Get the data until the end of the file:
  while (fff.eof() == 0) {
    brainbow_data x;
    fff >> i
	>> x.time
	>> x.founder
	>> x.birth
	>> x.mother_index
	>> x.daughter1_index
	>> x.daughter2_index
	>> x.ss_position
	>> x.colour;
    x.lineage = false;
    if (fff.fail() == false) {
      rb_list.push_back(x);
    } else { cerr << "ERROR in brainbow::read_from_file()!\n\n"; }
    // cerr<<i<<" ";
  }
  fff.close();
  cout << "\nTotal number of read events = " << rb_list.size() << "\n";
}
// for analysis:

long brainbow::get_founder_cell_index(long &cell_index) {
   // returns the event_index of the founder cell of the cell depicted by the event_index
   // <cell_index>
   long founder_index = cell_index;
   brainbow_data x = rb_list[cell_index];
   while (x.founder == false) {
      founder_index = x.mother_index;
      x = rb_list[founder_index];
   }
   return founder_index;
}
long brainbow::get_lineage_cell_index(long &cell_index) {
  /* returns the event_index of the lineage cell of origin of the cell depicted 
   * by the event_index <cell_index>
   */
  brainbow_data x = rb_list[cell_index];
  long lineage_index = cell_index;
  while (x.lineage == false && x.founder == false) {
    lineage_index = x.mother_index;
    x = rb_list[lineage_index];
  }
  /* If a founder cell was found before the lineage origin was detected, 
   * an error message is returned.
   * This means that founder cells showing up after <lineage_time> 
   * are ignored for the lineage analysis. */
  if (x.lineage == false) { lineage_index = -1; }
  return lineage_index;
}
long brainbow::get_clone_of_founder(long &cell_index) {
  /* returns the position on shape space of the founder cell of the cell depicted 
   * by the event_index <cell_index>
   */
   long founder_index = get_founder_cell_index(cell_index);
   return rb_list[founder_index].ss_position;
}
long brainbow::get_index_for_time(double &time) {
   // returns the smallest event_index (there can be more than one event at a time)
   // corresponding to time <time> or earlier
   if (rb_list.size() == 0) {
      return -1;                    // then there is no event in the list at all
   }
   if (time < 0) {
      return rb_list.size() - 1;         // if time<0 return the last event index
   }
   unsigned long index = 0;
   if (rb_list[index].time > time) {
      return -1;                           // if the first event is already later return -1
   }
   // note that rounding errors require to compare to time + infinitesimal_time 
   // (absolute scale, so not general !!!)
   while (index < rb_list.size() && rb_list[index].time < time + infinitesimal_time) {
      index++;
   }
   // now index points to the first event at a time later than or equal to <time>, so go one back:
   index--;  // it cannot be zero because of the two returns above
   return index;
}
void brainbow::get_cell_indices_at_time(double &time, vector<long> &cellind) {
   /* This returns the event indices of all cells alive and present at time "time".
    * The problem to solve is that there are many cells present at this time
    * but the last event has been before.
    * Only events which are not mother_indices of already identified cells
    * have to be included.
    */
   cellind.clear();  // start from an empty array
   //cout << "  Find cells at time=" << time << " ... ";
   // Get the highest event index for time <time>
   long maxindex = get_index_for_time(time);
   //cout << "maxindex=" << maxindex;
   vector<long> tmp;
   tmp.reserve(maxindex + 1);
   for (long i = 0; i <= maxindex; i++) {
      tmp.push_back(i);
   }
   long tmpi = maxindex;
   long i,j;
   // From here go backwards
   while (tmpi >= 0) {
      if (tmp[tmpi] != -1) {
         i = tmp[tmpi];
         if (rb_list[i].birth) {
            // it has to be a birth event
            // add this cell to the list
            cellind.push_back(i);
         }
         // recursively remove all mothers from the list of possible cells to be included,
         // do this as well when birth==false, which is death.
         j = i;
         while (rb_list[j].mother_index != -1) {
            // if ==-1, the founder cell was reached
            tmp[rb_list[j].mother_index] = -1;
            j = rb_list[j].mother_index;
         }
      }
      --tmpi;
   }
   //cout << ",  found " << cellind.size() << " cells.\n";
   // Note that this returns the event-indices in decreasing order.
}
vector<long> brainbow::get_cell_subset(vector<long> all_cells) {
   // Randomly choses cell_subset_size values from the local copy of cell_indices (all_cells)
   // and saves them in a vector subset, which is returned to the calling routine.
   vector<long> subset;
   subset.reserve(cell_subset_size);
   if (all_cells.size() > 0) {
      for (int i = 0; i < cell_subset_size; ++i) {
         int itmp = irandom(all_cells.size());
         subset.push_back(all_cells[itmp]);
         // remove the value in the local copy to avoid twice the same cell to be included
         if (all_cells.size() > 0) {
            all_cells.erase(all_cells.begin() + itmp);
         } else { i = cell_subset_size; }
      }
   }
   return subset;
}
int brainbow::get_random_colour(int &Ncols) {
   int colour = -1;
   if (how2do_staining == gabriel_4colours) {
      /*
       * colour=rand() % 10; // use 10 colours even though Ncols=4;
       * // now project 0-2 to 0; 3-5 to 1; 6-8 to 2; 9 to 3 to get Gabriels fractions
       * if (colour<3) colour=0; else if (colour<6) colour=1; else if (colour<9) colour=2; else
       * colour=3;
       */
      double pthis = drandom();
      if (pthis < 0.3) { colour = 0; } else if (pthis < 0.6) {
         colour = 1;
      } else if (pthis < 0.9) {
         colour = 2;
      } else { colour = 3; }
   } else if (how2do_staining == gabriel_10colours) {
      // values from Gabriel's email from Nov 5th, 2015 (first value black, others colours):
      double normit = 1.0 - 0.5212;   // 1 - black colour frequency
      double pvec[10] = {
         0.132082352, 0.121671862, 0.105401788, 0.056367556, 0.027128581, 0.016105016,
         0.014065362, 0.003171486, 0.001803856, 0.001002142
      };
      /* For test purpose I used this one and compared it to gabriel_4colours:
       *  double normit=1.0-0.6;
       *  double pvec[10]={ 0.12, 0.12, 0.12, 0.04, 0, 0, 0, 0, 0, 0 };
       */
      for (int i = 0; i < 10; ++i) {
         pvec[i] /= normit;
      }
      for (int i = 1; i < 10; ++i) {
         pvec[i] += pvec[i - 1];
      }
      // for (int i=0; i<10; i++) cout<<pvec[i]<<", "; cout<<"\n";
      if (pvec[9] < 1 - 1.e-08) {
         cout << "WARNING:\n"
	      << "brainbow::get_random_colour(int&,double): total probability of colours < 1\n";
      } else if (pvec[9] > 1 + 1.e-08) {
         cerr << "WARNING:\n"
	      << "brainbow::get_random_colour(int&,double): total probability of colours "
	      << pvec[9] << " > 1\n"
	      << "         abort analysis.\n";
         exit(1);
      }
      double pthis = drandom();
      colour = 0;
      while (pvec[colour] < pthis) {
         ++colour;
      }
      // <colour> contains the colour index to chose
   } else if (how2do_staining == gabriel_10colours_2founder) {
      // values from Gabriel's email from Nov 13th, 2015 (first value black, others colours):
      double normit = 1.0 - 0.1089;   // 1 - black colour frequency
      double pvec[10] = {
         0.172, 0.172, 0.172, 0.0511, 0.08, 0.08,
         0.028, 0.08, 0.028, 0.028
      };
      for (int i = 0; i < 10; ++i) {
         pvec[i] /= normit;
      }
      for (int i = 1; i < 10; ++i) {
         pvec[i] += pvec[i - 1];
      }
      // for (int i=0; i<10; i++) cout<<pvec[i]<<", "; cout<<"\n";
      if (pvec[9] < 1 - 1.e-08) {
         cout << "WARNING:\n" 
	      << "brainbow::get_random_colour(int&,double): total probability of colours < 1\n";
      } else if (pvec[9] > 1 + 1.e-08) {
         cerr
            << "WARNING:\n"
	    << "brainbow::get_random_colour(int&,double): total probability of colours "
            << pvec[9] << " > 1\n"
            << "         abort analysis.\n";
         exit(1);
      }
      double pthis = drandom();
      colour = 0;
      while (pvec[colour] < pthis) {
         ++colour;
      }
      // <colour> contains the colour index to chose
   } else if (how2do_staining == Ncolours_random) {
      colour = rand() % Ncols;   // returns a colour between 0 and Ncolour-1
   } else {
      cerr << "In brainbow::get_random_colour(int&, double):\n"
           << "how2do_staining undefined.\n";
      exit(1);
   }
   return colour;
}
void brainbow::transfer_colour2daughters(long cell_index) {
  // recursively call transfer_colour2daughters for each coloured daughter cell.
  if (rb_list[cell_index].daughter1_index != -1) {
    if (rb_list[rb_list[cell_index].daughter1_index].colour != -1) {
      /* If this happens, a daughter cell exists that acquired a colour from elsewhere.
       * This should never happen unless restaining of cells is allowed. */
      if (not(allow_restaining)) {
	cerr << "COL: mother=" << rb_list[cell_index].colour
	     << ">> daughter=" << rb_list[rb_list[cell_index].daughter1_index].colour
	     << "\n";
	exit(1);
      }
    }
    rb_list[rb_list[cell_index].daughter1_index].colour = rb_list[cell_index].colour;
    transfer_colour2daughters(rb_list[cell_index].daughter1_index);
  }
  if (rb_list[cell_index].daughter2_index != -1) {
    rb_list[rb_list[cell_index].daughter2_index].colour = rb_list[cell_index].colour;
    transfer_colour2daughters(rb_list[cell_index].daughter2_index);
  }
  /* This routine ends without further recursive call when both daughter indices are -1.
   * Both daughter cells are -1 when there is no further progeny in the list or
   * when the cell died.
   * When only one daughter index is !=-1 (which is always the first), this means that
   * the next event is a death or exit event. The colour will also be saved in the
   * death or exit event. */
}
double brainbow::get_stain_amplitude() {
   /* This routine will calculate the probability with exponential decay
    * Relies on the following global variables:
    * uses tamoxifen_decay: if <0 just return stain_fraction,
    *                       use it as half life time in hr otherwise
    * uses stain_dt as time step for transformation of a rate to a probability
    * uses stain_fraction as default response or as integrated fraction of stained cells
    */
   if (tamoxifen_decay < 0) { return stain_fraction; }
   /* We have:
    * stain_fraction = p_0 \int_0^\infty exp(-t/tamoxifen_decay) dt = p_0 tamoxifen_decay
    * For a finite interval the situation is more involved (t_max=tamoxifen_stop_time):
    * stain_fraction = p_0 \int_0^{t_max} exp(-t/tamoxifen_decay) dt
    *                = p_0 tamoxifen_decay (1 - exp(-t_max/tamoxifen_decay) )
    * solved for p_0 this yields:
    * p_0 = (stain_fraction/tamoxifen_decay) / (1 - exp(-t_max/tamoxifen_decay) )
    *
    * Note that tamoxifen_stop_time is an absolute time parameter.
    * Here it was assumed relative to zero, thus tamoxifen_stop_time has to be reset relative to
    * stain_time,
    * because stain_time is the equivalent of zero in this derivation.
    */
   double p_0 = (stain_fraction / tamoxifen_decay)
                / (1.0 - exp(-1.0 * (tamoxifen_stop_time - stain_time) / tamoxifen_decay));

   // p_0 is a rate per hour --> the probability per time step is generated as p_0 stain_dt:
   p_0 *= stain_dt;  // Note that this is a rate, not an inverse half time, so no log(2) needed
                     // here.
   return p_0;
}
double brainbow::get_stain_probability(double &time, double &amplitude) {
   /* The decay of <amplitude> with an exponential is calculated
    * time is interpreted relative to stain_time (global)
    * the decay constant is tamoxifen_decay (global)
    */
   if (tamoxifen_decay < 0) { return stain_fraction; }
   double value = amplitude * exp((stain_time - time) * log(2) / tamoxifen_decay);
   return value;
}
void brainbow::stain(double time, double &stain_prob, int &Ncols, int * Nfc,
                     vector<long> &founder_list,
                     vector<int> &projection) {
   /* Attributes a colour to each cell present at time <time>.
    * The number of colours and their probability of being chosen are hard-coded in
    * get_random_colour(...).
    * <Ncols> is only used if the colours are chosen randomly with same probabilty.
    */
   cout << "Stain fraction " << stain_prob << " at time t=" << (time - stain_time)
        << " hours post tamoxifen ...\n";
   // Get a list of birth_event_indices of cells at time <time>:
   vector<long> stainindices;
   stainindices.reserve(2000);
   /* The following is time consuming and called several times when tamoxifen decay is active.
    * Consider using the previous run as starting point and only make the "delta" to the next
    * time point. This would involve removing events with death events and adding new borns.
    */
   get_cell_indices_at_time(time, stainindices);

   // Save the dimension of the lists
   int Nfounders = founder_list.size();
   // Go through all these cells and attribute a colour:
   for (unsigned long i = 0; i < stainindices.size(); i++) {
      /* If colour was black so far (==-1), staining is done with the given probability.
       * As tamoxifen might recut and change the colour, in the case of active tamoxifen decay,
       * one may set the option allow_restaining=true. If not set, only black cells are stained.
       */
      if ((rb_list[stainindices[i]].colour == -1) || allow_restaining) {
         if (drandom() < stain_prob) {
            if (origin_same_colour) {
               /* every cell is checked for the color attributed to the founder cell in Nfc 
		* and gets this one.
                * Nfc is defined with a dimension corresponding to founder_list or
                * founder_short_list (see below).
                */
               // find the founder cell to event index stainindices[i]
               long founderind = get_founder_cell_index(stainindices[i]);
               // cout<<"staini="<<stainindices[i]<<"; founderind="<<founderind<<", ";
               int founderarraypos = get_origin_arraypos(founderind, Nfounders, founder_list);
               if (merge_identical_origin) {
                  // find the founder cell in the merged short_list, which overwrites
                  // founderarraypos
                  project_origin(founderarraypos,projection);
               }
               // attribute the colour attributed to this founder cell to the cell
               rb_list[stainindices[i]].colour = Nfc[founderarraypos];
               // cout<<"Nfc["<<founderarraypos<<"]="<<Nfc[founderarraypos]<<"; ";
            }     // if (origin_same_colour)
            else {
               rb_list[stainindices[i]].colour = get_random_colour(Ncols);
            }
            // Transfer this to all daughters if staining really took place
            if (rb_list[stainindices[i]].colour != -1) {
               // the coloured cell has to transfer its colour to all of its progeny:
               transfer_colour2daughters(stainindices[i]);
            }
         }    // if (drandom()<stain_prob
      }   // if (rb_list[]...==-1)
   }
   cout << "staining done.\n";
}
void brainbow::stain_all_founder(double stain_prob, int &Ncolours) {
   // find all founders
   vector<long> flist;
   flist.reserve(1000);
   int Nf = get_all_founders(-1,flist);  // time=-1 infers that all founders are returned
   // The probability is now determined on the level of founder cells, thus independent
   // of a tamoxifen decay time (i.e. fixed):
   int attribute_color;
   for (int i = 0; i < Nf; i++) {
      attribute_color = -1;   // stays -1=black if not stained below
      if (drandom() < stain_prob) { attribute_color = get_random_colour(Ncolours); }
      if (attribute_color >= 0) {
         // if color other than black
         rb_list[flist[i]].colour = attribute_color;
         // the coloured cell has to transfer its colour to all of its progeny:
         transfer_colour2daughters(flist[i]);
      }
   }
}
void brainbow::stain(vector<long> &founder_list,
                     vector<int> &projection) {
   /* organises the staining of cells with tamoxifen
    * uses global variables stain_time, stain_dt, Ncolours
    */
   // Load the default staining amplitude (stain_fraction if tamoxifen decay is off).
   // If tamoxifen decay is on, this is the initial probability based on the integrated staining
   // fraction,
   // which is then subject to decay --> stain_probability below.
   double stain_amplitude = get_stain_amplitude();

   if (how2do_staining == gabriel_10colours_2founder) {
      stain_all_founder(stain_amplitude,Ncolours);
   } else {
      /* In the case of origin_same_colour, the colours have to be attributed to all founder cells
       * before staining.
       * These are saved in Nfc[] and can be used later on to attribute colours to the cells.
       */
      // array of colours for all founders (its too long of the founder_short_list is used, but
      // doesn't matter)
      int Nf = founder_list.size();
      int Nfc[Nf];
      // Each founder cell gets a colour with the normal rules (only used if origin_same_colour is
      // true)
      if (origin_same_colour) {
         for (int i = 0; i < Nf; i++) {
            Nfc[i] = get_random_colour(Ncolours);
            // cout<<Nfc[i]<<";";
         }
      }

      // The normal staining procedure is here:
      double cutoff = 1.0e-06;   // this value is empirical and must be < stain_fraction
      double stain_probability;
      double time = stain_time;
      do {
         // cout<<"stain(): time="<<time<<"\n";
         // set stain_probability for time <time>
         stain_probability = get_stain_probability(time,stain_amplitude);
         // cout<<"stain_probability="<<stain_probability<<"\n";
         // stain at time <time>
         stain(time, stain_probability, Ncolours, Nfc, founder_list, projection);
         // set the time for the next staining round
         time += stain_dt;
         // Only continue when tamoxifen_decay is set; stain_amplitude==0 stops the while loop
         if (tamoxifen_decay < 0) { stain_probability = 0; }
      } while (stain_probability > cutoff && time < tamoxifen_stop_time + 1.e-08);
   }
}
int brainbow::get_all_lineages(double time, vector<long> &linlist) {
  /* The lineages at time <time> are fixed once and assumed unchanged.
   * New clones enterring the GC later on are ignored in this analysis.
   * Thus, <lineage_list> is fixed for all analysis time points and is calculated once.
   * All cells, even with same shape space position, are counted as a separate lineage.
   * Therefore, the <linlist> is defined by determination of all cells at <time>.
   * This strategy applies for the cases of with or without tamoxifen decay.
   * Note that a lineage origin cell of a stained cell later on may be unstained
   * when either time<stain_time or when tamoxifen_decay<0.
   */
   // fill <linlist> with all cells existing at time <time>
   get_cell_indices_at_time(time, linlist);
   // attribute lineage=true to all events associated with a lineage origing cell
   for (unsigned int n = 0; n < linlist.size(); n++) {
      rb_list[linlist[n]].lineage = true;
   }
   // return the number of lineages defined
   return linlist.size();
}
int brainbow::get_number_of_different_cell_types(double time, ofstream& diverse) {
  /* At time <time> the number of cells in different shape space positions is determined.
   * All cells are saved in a list.
   * Cells with a shape space position already carried by another cell are removed.
   * The length of the remaining list is returned.
   */
  diverse << time << "   ";
  vector<long> ssposlist;
  // fill <ssposlist> with all cells existing at time <time>
  get_cell_indices_at_time(time, ssposlist);
  // Replace the cell indices in <ssposlist> by their respective shape-space position
  for (unsigned int n = 0; n < ssposlist.size(); n++) {
    ssposlist[n] = rb_list[ssposlist[n]].ss_position;
  }
  long totalcells = ssposlist.size();
  diverse << totalcells << "   ";
  // Remove double entries from the list
  sort( ssposlist.begin(), ssposlist.end() ); 
  ssposlist.erase( unique( ssposlist.begin(), ssposlist.end() ), ssposlist.end() );
  diverse << ssposlist.size() << "   ";
  double fraction = double(ssposlist.size())/double(totalcells);
  diverse << fraction << "   " << 1.-fraction << endl;
  return ssposlist.size();
}
void brainbow::write_number_of_different_cell_types(double deltat, int Nevals) {
  ofstream diverse("Ncelltypes.out");
  diverse << "! Time[h] : number cells : number of types : ratio : 1-ratio\n";
  double tnow = 0;
  for (int i = 0; i < Nevals; ++i) {
    tnow += deltat;
    get_number_of_different_cell_types(tnow, diverse);
  }
  diverse.close();
  cout << "Ncelltypes.out generated.\n";
}
int brainbow::add_late_founders_to_lineages(double fromtime, 
					     vector<long> &foundlist,
					     vector<long> &linlist) {
  /* Add all founders in <foundlist> that entered the GC after time <fromtime> 
   * to <linlist> and returns the updated <linlist> and the increased 
   * number of lineages <totalnumber> to the calling routine.
   */
  for (unsigned int i=0; i<foundlist.size(); i++) {
    if (rb_list[foundlist[i]].time >= fromtime) 
      { linlist.push_back(foundlist[i]); }
  }
  for (unsigned int n = 0; n < linlist.size(); n++) 
    { rb_list[linlist[n]].lineage = true; }
  return linlist.size();
}
int brainbow::get_all_founders(double time, vector<long> &founder_list) {
   // returns the number of founders cells by this time
   // also has to define an array of indices to the founder events
   cout << "Find founder cells at time=" << time << " ... ";
   long timeindex = get_index_for_time(time);
   founder_list.clear();
   for (long i = 0; i <= timeindex; i++) {
      if (rb_list[i].founder) {
         founder_list.push_back(i);    // this one is private in brainbow
      }
   }
   cout << "identified " << founder_list.size() << " founder cells.\n";
   return founder_list.size();
}
int brainbow::merge_origins(int Norigin,
                            vector<long> &origin_list,
                            vector<int> &origin_projection,
                            vector<long> &origin_short_list) {
/** origin_list contains Norigin event_indices (either founders or lineages).
 *  origin_short_list is loaded with those event_indices of founders/lineages
 *                    with different position in affinity space.
 *  origin_projection is loaded with Norigin indices that point to the position in
 *                    orgin_short_list, at which the merged founder/lineage event index is saved.
 **/
   int Nshortlist = 0;
   /* Define a vector with the affinity space positions of all founders/lineages 
    * saved in the short list */
   vector<long> as_position;
   as_position.reserve(Norigin + 1);
   long as_i;
   bool already_there = false;
   for (int i = 0; i < Norigin; i++) {
      // get position in affinity space
      as_i = rb_list[origin_list[i]].ss_position;
      // search whether this position is already present
      already_there = false;
      for (int j = 0; j < Nshortlist; j++) {
         if (as_position[j] == as_i) {
            // save that it was there
            already_there = true;
            // attribute this value to the projection
            origin_projection.push_back(j);
            // don't change origin_short_list
         }
      }
      if (not (already_there)) {
         // add a new origin to the short_list and save the event_index therein
         origin_short_list.push_back(origin_list[i]);
         // add a new affinity space position to the as_position list
         as_position.push_back(as_i);
         // set origin_projection to the current
         origin_projection.push_back(Nshortlist);
         // increment the number of different founders/lineages
         ++Nshortlist;
      }
   }
   return Nshortlist;
}
void brainbow::project_origin(int &listpos, vector<int> &projection) {
   listpos = projection[listpos];
}
int brainbow::get_origin_arraypos(long origin_index, int Norigins, vector<long> &origin_list) {
   // Returns the position of origin_index in the array origin_list of length Norigins
   // This can be the founder_list or the lineage_list
   bool found = false;
   int i = -1;
   while (not (found) && i < Norigins - 1) {
      ++i;
      if (origin_list[i] == origin_index) { found = true; }
   }
   if (not (found)) { return -1; }
   return i;
}
int brainbow::get_dominant_origin(double time, vector<long> &cellind, // time and event indices
                                  int Norigins, vector<long> &origin_list, // list of
                                                                           // founder/lineage and
                                                                           // length
                                  vector<int> &projection, // map from origin_list to its short_list
                                  int largest_colour, bool use_founder, // mode of action
                                  long * freq) {
   // result
   /* Loads freq (of length Norigins) with a distribution of founder/lineage frequency in
    * cells with event indices listed in cellind at time <time>.
    * If largest_colour>=0, the results are restricted to the members of the most frequent colour
    *                       with colour <largest_colour>.
    * Else (i.e. largest_colour==-1) the "real" founder/lineage frequency is saved.
    * The index of the most frequent founder/lineage cell in founder_list/lineage_list is returned.
    * If use_founder==true, Nfounders and founder_list were passed to Norigins and origin_list,
    *                       Nlineages and lineage_list otherwise.
    * if merge_identical_origin==true, <projection> is used to map the resulting origin cell
    *                                  onto the respective short_list (defined outside)
    */
   for (int i = 0; i < Norigins; i++) {
      freq[i] = 0;
   }
   long origin_index;
   for (unsigned long i = 0; i < cellind.size(); i++) {
      if ((largest_colour < 0) || (rb_list[cellind[i]].colour == largest_colour)) {
         // Either do analysis of real dominant founder/lineage or restrict to the dominant colour
         // Get the event index associated with the founder/lineage of cell with even index
         // cellind[i]:
         if (use_founder) {
            origin_index = get_founder_cell_index(cellind[i]);
         } else {
            origin_index = get_lineage_cell_index(cellind[i]);
            // this can return -1, when a founder cell was found at a time after the lineage
            // definition
            // time
         }
         if (origin_index >= 0) {
            // get the position of the origin cell in the origin_list that contains event index
            // origin_index
            int origin_arraypos = get_origin_arraypos(origin_index, Norigins, origin_list);
            if (origin_arraypos >= 0) {
               // eventually project the array position onto the corresponding position in the short
               // list
               if (merge_identical_origin) {
                  project_origin(origin_arraypos,projection);
               }
               // add this to the frequency counter
               ++freq[origin_arraypos];
            }
         }
      }
   }
   // find the most frequent origin cell in origin_frequency:
   int most_frequent_origin = 0;
   for (int i = 1; i < Norigins; i++) {
      if (freq[i] > freq[most_frequent_origin]) { most_frequent_origin = i; }
   }
   /* Note that in the case of lineages not all cells find a lineage because the lineage definition
    * time
    * can be earlier than the end of founder cells enterring the GC. This returns -1 above and
    * prevents
    * the counters to go up. So the sum of all freq[] entries is not equal to the number of cells,
    * i.e. the length of <cellind>. Therefore, in the clone_lineage_fractions.out file, the
    * percentage
    * of lineages does not add up to 100%. The deviation of 100% reflects how many founder cells
    * enterred
    * after the definition time for lineages.
    */
   return most_frequent_origin;
}
unsigned long brainbow::evaluate_colours(long * colour_frequency, int &Ncols,
                                         vector<long> &cellind) {
   // Loads the array <colour_frequency> of dimension Ncols with the number of cells
   // that carry the colour associated with the respective index in the array.
   // Returns the total number of cells
   long Ncells = 0;
   for (int i = 0; i < Ncols; i++) {
      colour_frequency[i] = 0;
   }
   for (unsigned long i = 0; i < cellind.size(); i++) {
      if (rb_list[cellind[i]].colour != -1) {
         ++colour_frequency[rb_list[cellind[i]].colour];
         ++Ncells;
      }
   }
   return Ncells;  // this is the total number of coloured cells in the cell_indices
}
void brainbow::get_fractions(long Ncells, long * frequency, double * fraction, int &Norigins) {
   if (Ncells > 0) {
      for (int i = 0; i < Norigins; i++) {
         fraction[i] = double (frequency[i]) / double (Ncells);
      }
   } else {
      for (int i = 0; i < Norigins; i++) {
         fraction[i] = 0;
      }
   }
}
void brainbow::put_fraction_header(ofstream &xf) {
   xf << "! time  cellnumber  Norigins  percentage_of_this_origin[1 .. Norigins]\n";
}
void brainbow::save_clonality(double &evaluation_time, ofstream &outfile,
                              long unsigned &Nstained_cells, double &largest_colour_fraction,
                              int &dominant_origin_in_colour, long &colour_ss_position,
                              long &Nreal_cells, double &largest_origin_fraction,
                              int &dominant_origin_in_real, long &real_ss_position,
                              bool clonal) {
   // Write everything to the output file
   outfile << evaluation_time << "\t" << evaluation_time - stain_time << "\t"
           << Nstained_cells << "\t"
           << largest_colour_fraction << "\t"
           << dominant_origin_in_colour << "\t"
           << colour_ss_position << "\t"
           << Nreal_cells << "\t"
           << largest_origin_fraction << "\t"
           << dominant_origin_in_real << "\t"
           << real_ss_position << "\n";
   cout << "Largest color occupies " << largest_colour_fraction << " of the stained cells at day "
        << (evaluation_time - stain_time) / 24.0 << " after injection of tamoxifen.\n";
   if (clonal) {
      cout << "Largest founder clone occupies ";
   } else { cout << "Largest lineage occupies "; }
   cout << largest_origin_fraction << " of all cells at day "
        << (evaluation_time - stain_time) / 24.0 << " after injection of tamoxifen.\n";
}
void brainbow::save_fractions(double time, long cell_num, int &Ncolumns, int &maxNcolumns,
                              double * fraction, ofstream &fx) {
   fx << time << "\t" << cell_num << "\t" << Ncolumns << "\t";
   for (int i = 0; i < Ncolumns; ++i) {
      fx << 100 * fraction[i] << "\t";
   }
   for (int i = Ncolumns; i < maxNcolumns; ++i) {
      fx << "0\t";
   }
   fx << "\n";
}
int brainbow::find_largest(long * frequency, int &Nentries) {
   // Returns the colour/clone/lineage_index in frequency with the largest number of cells
   int largest = 0;
   for (int i = 1; i < Nentries; i++) {
      if (frequency[i] > frequency[largest]) { largest = i; }
   }
   return largest;
}
int brainbow::find_seclargest(long * frequency, int ilargest, int &Nentries) {
   // Returns the colour/clone/lineage_index in frequency with the 2nd-largest number of cells
   int seclargest = 0;
   if (seclargest == ilargest) { seclargest = 1; }
   for (int i = seclargest+1; i < Nentries; i++) {
      if (i != ilargest && frequency[i] > frequency[seclargest]) { seclargest = i; }
   }
   return seclargest;
}
void brainbow::make_brainbow(double * evaluation_times, int Nevaluations) {
   cout << "\nBRAINBOW ANALYSIS:\n";

   // open files:
   char clonality_file[50] = "";
   char lineage_clonality_file[50] = "";
   char colour_fractions_file[50] = "";
   char clone_fractions_file[50] = "";
   char lineage_fractions_file[50] = "";
   char clone_fractions_rand_file[50] = "";
   char clone_lineage_fractions_rand_file[50] = "";
   if (do_clonality) {
      strcat(clonality_file,"clonality.out");
      strcat(lineage_clonality_file,"clonality_lineage.out");
      strcat(colour_fractions_file,"clone_colour_fractions.out");
   } else {
      strcat(clonality_file,"clonality_tmp.out");
      strcat(lineage_clonality_file,"clonality_lineage_tmp.out");
      strcat(colour_fractions_file,"clone_colour_fractions_tmp.out");
   }
   ofstream fff(clonality_file);
   fff << "! time  dtime [Ncells_col  %of-largest-colour  dominant_founder_col clone] "
       << "[Ncells_total %of-largest-clone  dominant_founder_real clone]\n";
   ofstream lll(lineage_clonality_file);
   lll << "! time  dtime [Ncells_col  %of-largest-colour  dominant_lineage_col clone] "
       << "[Ncells_total %of-largest-lineage  dominant_lineage_real lineage]\n";
   ofstream colf(colour_fractions_file);
   colf << "! time  stained-cell-number  Ncolours  percentage_of_this_colour[1 .. Ncolours]\n";

   if (do_clone_fractions) {
      strcat(clone_fractions_file,"clone_fractions.out");
      strcat(lineage_fractions_file,"clone_lineage_fractions.out");
   } else {
      strcat(clone_fractions_file,"clone_fractions_tmp.out");
      strcat(lineage_fractions_file,"clone_lineage_fractions_tmp.out");
   }
   ofstream cf(clone_fractions_file);
   put_fraction_header(cf);
   ofstream lf(lineage_fractions_file);
   put_fraction_header(lf);

   if (do_clone_fractions_rand) {
      strcat(clone_fractions_rand_file,"clone_fractions_rand.out");
      strcat(clone_lineage_fractions_rand_file,"clone_lineage_fractions_rand.out");
   } else {
      strcat(clone_fractions_rand_file,"clone_fractions_rand_tmp.out");
      strcat(clone_lineage_fractions_rand_file,"clone_lineage_fractions_rand_tmp.out");
   }
   ofstream cfr(clone_fractions_rand_file);
   put_fraction_header(cfr);
   ofstream lfr(clone_lineage_fractions_rand_file);
   put_fraction_header(lfr);

   // Find out the maximum number of founder cells by calculating the latest evaluation time:
   int Nfounders;
   int Nfounders_short = -1;
   vector<long> founder_list;
   founder_list.reserve(1000);
   vector<long> founder_short_list;
   founder_short_list.reserve(1000);
   vector<int> founder_list_projection;
   founder_list_projection.reserve(1000);
   int maxNfounders = get_all_founders(evaluation_times[Nevaluations - 1],founder_list);
   cout << "Maximum number of founder cells in the considered time frame is " << maxNfounders
        << ".\n";
   if (merge_identical_origin) {
      Nfounders_short = merge_origins(maxNfounders,
                                      founder_list,
                                      founder_list_projection,
                                      founder_short_list);
      cout << "Reduced to " << Nfounders_short << " different founders at max evaluation time "
           << evaluation_times[Nevaluations - 1] << ".\n";
   }
   // ... this is needed in save_fractions(..) to make lines of equal length for R routines.

   // The same for lineages:
   vector<long> lineage_list;
   lineage_list.reserve(5000);
   vector<long> lineage_short_list;
   lineage_short_list.reserve(5000);
   vector<int> lineage_list_projection;
   lineage_list_projection.reserve(5000);
   cout << "Find lineages at time=" << lineage_time << " ...\n";
   int Nlineages = get_all_lineages(lineage_time, lineage_list);  // also sets rb_list[].lineage
   cout << "Found " << Nlineages << " lineages at time " << lineage_time << ".\n";
   if (include_late_founders_in_lineages) {
     cout << "Add founders after time=" << lineage_time << " to the lineages ...\n";
     Nlineages = add_late_founders_to_lineages(lineage_time, founder_list, lineage_list);
     cout << Nlineages << " lineages defined with late founders.\n";
   }
   int Nlineages_short = -1;
   if (merge_identical_origin) {
      Nlineages_short = merge_origins(Nlineages,
                                      lineage_list,
                                      lineage_list_projection,
                                      lineage_short_list);
      cout << "Reduced to " << Nlineages_short
           << " different lineages at lineage definition time " << lineage_time << ".\n";
   }

   stain(founder_list, founder_list_projection);
   long colour_frequency[Ncolours];
   double colour_fraction[Ncolours];

   unsigned long Nstained_cells = 0;
   int i_largest = 0, i_seclargest = 0;
   vector<long> cell_indices;
   cell_indices.reserve(30000);
   cout << "\n";
   // Do evaluations at different times:
   for (int i = 0; i < Nevaluations; i++) {
      double evaluation_time = evaluation_times[i];
      cout << "Evaluation_time=" << evaluation_time << " ... \n";

      // load cell_indices with a list of cell indices in rb_list existing at this time:
      get_cell_indices_at_time(evaluation_time, cell_indices);
      long Ntotal_cells = cell_indices.size();

      if (do_clonality) {
         // find the colour frequencies, save the largest, calculate the fractions:
         Nstained_cells = evaluate_colours(colour_frequency, Ncolours, cell_indices);
         i_largest = find_largest(colour_frequency, Ncolours);
	 i_seclargest = find_seclargest(colour_frequency, i_largest, Ncolours);

         get_fractions(Nstained_cells,colour_frequency,colour_fraction,Ncolours);
         save_fractions(evaluation_time,
                        Nstained_cells,
                        Ncolours,
                        Ncolours,
                        colour_fraction,
                        colf);
      }

      // load the founder_list and return its length
      Nfounders = get_all_founders(evaluation_time,founder_list);
      Nfounders_short = -1;
      if (merge_identical_origin) {
         Nfounders_short = merge_origins(Nfounders,
                                         founder_list,
                                         founder_list_projection,
                                         founder_short_list);
         cout << "Reduced to " << Nfounders_short << " different founders at evaluation time "
              << evaluation_time << ".\n";
      }

      // Generate a distribution of dominant founders in dominant colour and in real:
      long dominant_founder_frequency_in_colour[Nfounders];
      int dominant_founder_in_colour = 0;
      long dominant_lineage_frequency_in_colour[Nlineages];
      int dominant_lineage_in_colour = 0;
      if (do_clonality) {
         dominant_founder_in_colour
            = get_dominant_origin(evaluation_time, cell_indices,
                                  Nfounders, founder_list, founder_list_projection,
                                  i_largest, true,
                                  dominant_founder_frequency_in_colour);
         dominant_lineage_in_colour
            = get_dominant_origin(evaluation_time, cell_indices,
                                  Nlineages, lineage_list, lineage_list_projection,
                                  i_largest, false,
                                  dominant_lineage_frequency_in_colour);
      }

      // Go through all cells and find the dominant founder clone and its fraction:
      long dominant_founder_frequency_in_real[Nfounders];
      int dominant_founder_in_real = 0;
      int clone_largest = 0;
      double clone_fraction[Nfounders];

      long dominant_lineage_frequency_in_real[Nlineages];
      int dominant_lineage_in_real = 0;
      int lineage_largest = 0;
      double lineage_fraction[Nlineages];

      if (do_clonality || do_clone_fractions) {
         dominant_founder_in_real
            = get_dominant_origin(evaluation_time, cell_indices,
                                  Nfounders, founder_list, founder_list_projection,
                                  -1, true,
                                  dominant_founder_frequency_in_real);
         clone_largest = find_largest(dominant_founder_frequency_in_real,Nfounders);
         get_fractions(
            cell_indices.size(),dominant_founder_frequency_in_real,clone_fraction,Nfounders);

         dominant_lineage_in_real
            = get_dominant_origin(evaluation_time, cell_indices,
                                  Nlineages, lineage_list, lineage_list_projection,
                                  -1, false,
                                  dominant_lineage_frequency_in_real);
         lineage_largest = find_largest(dominant_lineage_frequency_in_real,Nlineages);
         get_fractions(
            cell_indices.size(),dominant_lineage_frequency_in_real,lineage_fraction,Nlineages);
      }

      if (do_clone_fractions) {
         save_fractions(evaluation_time,
                        cell_indices.size(), Nfounders, maxNfounders, clone_fraction, cf);
         save_fractions(evaluation_time,
                        cell_indices.size(), Nlineages, Nlineages, lineage_fraction, lf);
      }

      double colfrac = colour_fraction[i_largest];
      if (include_seclargest_colour) { colfrac += colour_fraction[i_seclargest]; }
      if (do_clonality) {
         save_clonality(evaluation_time,fff,
                        Nstained_cells,colfrac,dominant_founder_in_colour,
                        rb_list[founder_list[dominant_founder_in_colour]].ss_position,
                        Ntotal_cells,clone_fraction[clone_largest],dominant_founder_in_real,
                        rb_list[founder_list[dominant_founder_in_real]].ss_position,true);
         save_clonality(evaluation_time,lll,
                        Nstained_cells,colfrac,dominant_lineage_in_colour,
                        rb_list[lineage_list[dominant_lineage_in_colour]].ss_position,
                        Ntotal_cells,lineage_fraction[lineage_largest],dominant_lineage_in_real,
                        rb_list[lineage_list[dominant_lineage_in_real]].ss_position,false);
      }

      if (do_clone_fractions_rand) {
         // Define a subset of the cells by random choice
         // and repeat the clone_fraction analysis with those
         vector<long> cell_subset = get_cell_subset(cell_indices);
         // redefine dominant_founder_frequency_in_real from the scratch (delete old one)
         dominant_founder_in_real
            = get_dominant_origin(evaluation_time, cell_subset,
                                  Nfounders, founder_list, founder_list_projection,
                                  -1, true,
                                  dominant_founder_frequency_in_real);
         clone_largest = find_largest(dominant_founder_frequency_in_real,Nfounders);
         get_fractions(cell_subset_size,
                       dominant_founder_frequency_in_real,
                       clone_fraction,
                       Nfounders);
         save_fractions(evaluation_time,
                        cell_subset_size,
                        Nfounders,
                        maxNfounders,
                        clone_fraction,
                        cfr);
         // redefine dominant_lineage_frequency_in_real from the scratch (delete old one)
         dominant_lineage_in_real
            = get_dominant_origin(evaluation_time, cell_subset,
                                  Nlineages, lineage_list, lineage_list_projection,
                                  -1, false,
                                  dominant_lineage_frequency_in_real);
         lineage_largest = find_largest(dominant_lineage_frequency_in_real,Nlineages);
         get_fractions(cell_subset_size,
                       dominant_lineage_frequency_in_real,
                       lineage_fraction,
                       Nlineages);
         save_fractions(evaluation_time,
                        cell_subset_size,
                        Nlineages,
                        Nlineages,
                        lineage_fraction,
                        lfr);
      }
      cout << "evaluation done.\n\n";
   }
   fff.close();
   lll.close();
   colf.close();
   cf.close();
   cfr.close();
   lf.close();
   lfr.close();
}
void brainbow::analysis() {
   // set analysis parts: do_clonality, do_clone_fractions, do_clone_fractions_rand
   set_do(true,true,true);
   set_Ncolour(10);
   set_stain_time(48);
   set_stain_fraction(0.4);  // (1-0.1089); // (0.4);
   set_stain_dt(2.0);  // hours
   set_lineage_time(stain_time);
   set_tamoxifen_decay(24);  // -1; 24
   set_include_late_founders_in_lineages(false);
   set_include_seclargest_colour(false);
   set_tamoxifen_stop_time(stain_time + 48);
   set_staining_mode(gabriel_10colours);  // (gabriel_10colours_2founder);
   set_origin_same_colour(false);
   set_merge_identical_origin(false);
   set_allow_restaining(false);
   set_cell_subset_size(70);
   double deltat = 24;
   int Nevals = 20;
   double eval_times[Nevals];
   for (int i = 0; i < Nevals; ++i) {
      eval_times[i] = stain_time + (double (i) * deltat);
   }
   // eval_times[0]=3; eval_times[1]=12; eval_times[2]=17;
   // for (int i=0; i<Nevals; ++i) eval_times[i]*=24.;
   make_brainbow(eval_times, Nevals);
}
