#include "cell.h"
#include <math.h>

// ============================================================
// ============================================================
// ============================================================
// ============================================================
// ============================================================

// +++ OPTION ==============================================================================
// sets the time resolution with which the velocity of the cell is analyzed
// double cell::deltat_v=0.166666666667; // the unit is in minutes
double cell::deltat_v = 1.0;  // the unit is in minutes
// =========================================================================================
// Note, the available v_modes that are chosen in the par-files:
// v_modus==1: random choice of v-states out of the N_v, constant \Delta v
// v_modus==2: as ==1 but with weighted probability (weights given within the code!)
// v_modus==3: as ==1 but coupling of v-state changes to changes of polarity
// v_modus==4: velocity states according to adhesion state (not done yet)
// Note, the number of v-states is given by N_v
// end OPTION ==============================================================================

long cell::id_last = 0;
long cell::adhesion_time = 0;
// static chemotaxis properties:
double cell::north_weight = 0.;
double cell::TFR_north_weight = 0.;
double cell::chemo_max = 10.;
double cell::chemo_steep = 1.e+10;
double cell::chemo_half = 2.e-10;
double cell::CXCL12crit = -1.;
double cell::CXCL13crit = -1.;
double cell::CXCL12recrit = -1.;
double cell::CXCL13recrit = -1.;
double cell::use_glucose = 20.e-18;           // 22. in Gernot's continuous model (mol/(cell sec))
double cell::use_oxygen = 95.e-18;            // 34. in Gernot's continuous model
double cell::use_glucose_pro = 20.e-18;       // identical in Gernot's models
double cell::use_oxygen_pro = 95.e-18;        // identical in Gernot's models
double cell::critical_nutrient = 2.5e-08;     // Mol^2
short cell::use_specific_turning_angles = 0;  // use isotropic target polarisation
short cell::lattice_dim = 3;
int cell::ab_resolution = 0;
bool cell::tmx_MHC = false;
bool cell::tmx_MHC_noAgPresentation = false;
bool cell::tmx_MHC_noDivision = false;
bool cell::tmx_MHC_noTfhSignal = false;
//MS
bool cell::exp_stop_apo=false;


cell::cell() {
  // cout<<"In cell default constructor ... ";
  ++id_last;
  id = id_last;
  index = 0;
  brainbow_index = -1;
  status = proliferate;
  trackit = false;
  writethis2track = trackini;
  trackno = -1;
  born_index = 0;
  born_time = 0;
  volume = 1;
  p_move = 0.;
  add_moves = 0.;
  clock = 0;
  age = 0;
  contact_inhibited = false;
  changed_polarity = 0;
  pos_ss = 0;
  p_mutation = 0.;
  v_state = 1.;
  adhesive = 0.;
  pressure = 1.;
  Ki67 = 0;
  BrdU = 0.;
  n_recycling = 0;
  n_mutation = 0;
  n_recandmute = 0;
  n_fdc_encounters = 0;
  n_immobile = 1;
  for (short i = 0; i < 3; i++) {
    polarity[i] = 0.;
    last_position[i] = 0.;
  }
  for (short i = 0; i < signals; i++) {
    responsive2signal[i] = false;
  }
  // cout <<"end of cell default constructor.\n";
}
// ============================================================

cell::cell(const cell &x) {
  // cout<<"In cell::cell(cell&) ...";
  operator =(x);
  // cout<<"end of cell::cell(cell&)\n";
}
cell&cell::operator =(const cell &x) {
  // cout<<"In cell::operator=(cell&) ...";
  index = x.index;
  brainbow_index = x.brainbow_index;
  status = x.status;
  trackit = x.trackit;
  writethis2track = x.writethis2track;
  trackno = x.trackno;
  born_index = x.born_index;
  born_time = x.born_time;
  volume = x.volume;
  clock = x.clock;
  age = x.age;
  contact_inhibited = x.contact_inhibited;
  p_move = x.p_move;
  add_moves = x.add_moves;
  pos_ss = x.pos_ss;
  p_mutation = x.p_mutation;
  v_state = x.v_state;
  adhesive = x.adhesive;
  pressure = x.pressure;
  Ki67 = x.Ki67;
  BrdU = x.BrdU;
  n_recycling = x.n_recycling;
  n_mutation = x.n_mutation;
  n_recandmute = x.n_recandmute;
  n_fdc_encounters = x.n_fdc_encounters;
  n_immobile = x.n_immobile;
  for (short i = 0; i < 3; i++) {
    polarity[i] = x.polarity[i];
    last_position[i] = x.last_position[i];
  }
  for (short i = 0; i < signals; i++) {
    responsive2signal[i] = x.responsive2signal[i];
  }
  // cout<<"end of cell::operator=(cell&)\n";
  return *this;
}
cell::~cell() {
  // cout<<"in ~cell()...\n";
}
// ============================================================

void cell::set_statics(const Parameter &par, ofstream &ana) {
   // adhesion
   adhesion_time = long (par.Value.adhesion_time / (60. * par.Value.deltat) + 0.5);

   // tamoxifen-induced actions
   tmx_MHC = par.Value.tmx_MHC;
   tmx_MHC_noAgPresentation = par.Value.tmx_MHC_noAgPresentation;
   tmx_MHC_noDivision = par.Value.tmx_MHC_noDivision;
   tmx_MHC_noTfhSignal = par.Value.tmx_MHC_noTfhSignal;

   //Bcl2 experiment --> stop apoptosis
   exp_stop_apo = par.Value.exp_stop_apo;

   // lattice dimension
   lattice_dim = par.Value.DimSpace;

   // resolution of affinities for antibodies
   ab_resolution = par.Value.antibodies_resolution;
   cout << "ab_resolution in cell::set_statics = " << ab_resolution << "\n";

   // chemotaxis properties:
   north_weight = par.Value.north_weight;
   TFR_north_weight = par.Value.TFR_north_weight;
   chemo_max = par.Value.chemo_max;
   chemo_steep = par.Value.chemo_steep
                 / (par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15 * par.N_A);
   // in # molecules
   chemo_half = par.Value.chemo_half * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15
                * par.N_A;
   // in # molecules
   CXCL12crit = par.Value.CXCL12crit * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15
                * par.N_A;
   // in # molecules
   CXCL13crit = par.Value.CXCL13crit * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15
                * par.N_A;
   // in # molecules
   CXCL12recrit = par.Value.CXCL12recrit * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15
                  * par.N_A;
   // in # molecules
   CXCL13recrit = par.Value.CXCL13recrit * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15
                  * par.N_A;
   // in # molecules
   ana << "CXCL12crit for desensitisation = " << CXCL12crit << " molecules\n"
       << "CXCL13crit for desensitisation = " << CXCL13crit << " molecules\n";
   ana << "CXCL12recrit for resensitisation = " << CXCL12recrit << " molecules\n"
       << "CXCL13recrit for resensitisation = " << CXCL13recrit << " molecules\n";
   if ((CXCL12crit > 0) && (CXCL12recrit > 0) && (CXCL12recrit > CXCL12crit)) {
      cout << "ERROR: CXCL12 resensitisation threshold concentration is larger than\n"
           << "       CXCL12 sensitisation threshold concentration.\n"
           << "       Please use appropriate values in the parameter file.\n\n\n";
      ana << "ERROR: CXCL12 resensitisation threshold concentration is larger than\n"
          << "       CXCL12 sensitisation threshold concentration.\n"
          << "       Please use appropriate values in the parameter file.\n\n\n";
   }
   if ((CXCL13crit > 0) && (CXCL13recrit > 0) && (CXCL13recrit > CXCL13crit)) {
      cout << "ERROR: CXCL13 resensitisation threshold concentration is larger than\n"
           << "       CXCL13 sensitisation threshold concentration.\n"
           << "       Please use appropriate values in the parameter file.\n\n\n";
      ana << "ERROR: CXCL13 resensitisation threshold concentration is larger than\n"
          << "       CXCL13 sensitisation threshold concentration.\n"
          << "       Please use appropriate values in the parameter file.\n\n\n";
   }

   if (par.Value.use_glucose <= 0.) {
      use_glucose = 0.;
   } else {
      use_glucose = par.Value.use_glucose         // mol/(cell sec)
                    * 3600. * par.Value.deltat    // 3600 hr = sec
                    * par.N_A;                    // /mol
   }
   ana << "Cell glucose consumption per cell (molecules) = " << use_glucose << "\n";
   if (par.Value.use_oxygen <= 0.) {
      use_oxygen = 0.;
   } else {
      use_oxygen = par.Value.use_oxygen          // mol/(cell sec)
                   * 3600. * par.Value.deltat    // 3600 hr = sec
                   * par.N_A;                    // /mol
   }
   ana << "Cell oxygen consumption per cell (molecules) = " << use_oxygen << "\n";
   if (par.Value.use_glucose_pro <= 0.) {
      use_glucose_pro = 0.;
   } else {
      use_glucose_pro = par.Value.use_glucose_pro     // mol/(cell sec)
                        * 3600. * par.Value.deltat    // 3600 hr = sec
                        * par.N_A;                    // /mol
   }
   ana << "Cell glucose consumption per proliferating cell (molecules) = " << use_glucose_pro
       << "\n";
   if (par.Value.use_oxygen_pro <= 0.) {
      use_oxygen_pro = 0.;
   } else {
      use_oxygen_pro = par.Value.use_oxygen_pro      // mol/(cell sec)
                       * 3600. * par.Value.deltat    // 3600 hr = sec
                       * par.N_A;                    // /mol
   }
   ana << "Cell oxygen consumption per proliferating cell (molecules) = " << use_oxygen_pro
       << "\n";
   if (par.Value.critical_nutrient <= 0.) {
      critical_nutrient = 0.;
   } else {
      /* The critical nutrient for necrosis has to be adopted to the resolution
       * of the signal-grid. In 3D it is fully adopted to the signal-grid, while
       * in 2D the third dimension, corresponding to the thickness of the cell-grid,
       * is taken from the cell-grid. */
      cout << "dx=" << par.Value.dx << ", dx_signal=" << par.Value.dx_signal << "\n";
      critical_nutrient = par.Value.critical_nutrient;    // (mol/l)^2
      /* Das "factor" wuerde das kritische Produkt in Molekuelen ausdruecken.
       * Ist aber nicht sinnvoll, da der Vergleich von Konzentrationen Skaleninvariant ist,
       * nicht der von Molekuelzahlen.
       */
      /* The criterion for necrosis is that the concentration-product (glucose
       * and oxygen) at the point on the signal-grid corresponding to the position
       * of the cell is larger than this critical concentration. As this simulation
       * works with numbers of molecules, the nutrient values are given as this
       * within a volume element. Thus the critical_nutrient has to be scaled
       * accordingly, which is ensured by the present factors.
       * -> Changed back to concentrations 10.1.2005
       */
   }
   double factor = par.Value.dx_signal * par.Value.dx_signal * 1.e-15     // dm^2
                   * par.Value.dx_signal * par.Value.dx_signal * 1.e-15   // dm^2
                   * par.N_A * par.N_A;                                   // /mol^2
   if (par.Value.DimSpace == 2) {
      factor *= par.Value.dx * par.Value.dx;
   } else {
      factor *= par.Value.dx_signal * par.Value.dx_signal;
   }

   ana << "Critical nutrient concentration product for necrosis (molecules) = "
       << critical_nutrient << "\n";
   if ((use_oxygen_pro * use_glucose_pro > 0.)
       && (use_oxygen_pro * use_glucose_pro > 0.5 * critical_nutrient * factor)) {
      ana
         <<
         "!!! Attention !!! Nutrient consumption large with respect to 0.5*critical_nutrient.\n";
      cout << "\n!!!!! Attention !!!!!\n"
           << "--> Nutrient consumption " << use_oxygen_pro * use_glucose_pro
           << " large with respect to 0.5*critical_nutrient=" << 0.5 * critical_nutrient
         * factor << ".\n"
           << "Use smaller time steps!!!\n\n";
      exit(1);
   }

   use_specific_turning_angles = par.Value.use_specific_turning_angles;
}
void cell::set_clock() { ++clock; }

void cell::aging(const double &dt) { age += dt; }

void cell::get_adhesion(const double &max) {
   adhesive = max;
   /* ### This is the place to introduce a more sophisticated calculation
    * of the expression of adhesion molecules.
    */
}
// ============================================================

short int cell::do_diffuse(states celltype, const long &li, space &l) {
   short err = 1;
   long j;

   // cout<<"do_diffuse(celltype="<<celltype<<": ";
   // Try to move the cell if it shall according to overall probability:
   if (drandom() < p_move / v_state) {
      // cout<<"p-criterion yes; ";
      double maxprojection = -99;
      long takethisindex = -1;
      // double maxextprojection=-99;
      // long takethisextindex=-1;
      // Preferentially move in direction of polarity:
      for (short k = 0; k < l.dim2; k++) {
         long v_j[l.dim];
         j = l.knot[index].near_n[k];
         if ((j != -1) && (l.cellknot[j].cell == nocell)) {
            // continue only if no cell at knot[j]
            for (short i = 0; i < l.dim; i++) {
               v_j[i] = l.knot[j].x[i] - l.knot[index].x[i];
            }
            // Jetzt enthaelt v_j den Differenzenvektor zwischen dem Nachbarn und index
            // und polarity den Richtungsvektor.
            double scals = l.get_scalarproduct(v_j, polarity);
            if (scals > maxprojection) {
               maxprojection = scals;
               takethisindex = j;
            }
         }
         /* Detect a boundary problem
          * if (j!=-1 && l.knot[j].status==external) {
          * for (short i=0; i<l.dim; i++)
          *  v_j[i]=l.knot[j].x[i]-l.knot[index].x[i];
          * double scals=l.get_scalarproduct(v_j,polarity);
          * if (scals>maxprojection) {
          *  maxextprojection=scals;
          *  takethisextindex=j;
          * }
          * }
          */
      }
      /* This allows to reset the polarity if the cell is captured at the boundary:
       * if (maxextprojection>maxprojection) {
       * l.get_random_direction(polarity);
       * if (trackit) writethis2track=true;
       * }
       */
      // cout<<"maxprojection="<<maxprojection<<"; takethisindex="<<takethisindex<<"; ";
      // takethisindex is positive if a free neighbour place exists
      // maxprojection is positive if a free neighbour place in the direction of polarity exists
      if ((takethisindex >= 0) && (maxprojection >= 0)) {
         // do the movement only if the scalarproduct is positive
         // and an empty neighbour was found:
         l.clear_knot(index);
         // l.cellknot[index].cell=nocell;
         // l.knot[index].listi=-1;
         index = takethisindex;
         l.set_knot(takethisindex, celltype, li);
         // l.knot[takethisindex].listi=li;
         // l.knot[takethisindex].cell=celltype;
         err = 0;
         if (trackit) {
            writethis2track = movement;
         }
      } else {
         // save that the movement was initated but suppressed by lack of space
         contact_inhibited = true;
      }
   }
   // cout<<"do_diffuse-end.\n";
   // No movement analysis is returned as realised in fragmove(..).
   // It is not necessary as any move is a move by one lattice step.
   // So only the information whether a move was performed or not is relevant!
   return err;
}
double cell::get_v_factor(const short &N_v, const double &p_slow_factor) {
   double _v_step = (1. - 1. / p_slow_factor) / (double (N_v) - 1.);
   return 1. / (1. - double (irandom(N_v)) * _v_step);
}
void cell::set_south_polarity(space &l) {
   // directs any predifined polarity towards the south direction
   // by switching the sign of the corresponding compoente of polarity
   if ((l.dim == 2) && (polarity[1] < 0.)) {
      polarity[1] *= -1.;
   }
   if ((l.dim == 3) && (polarity[2] < 0.)) {
      polarity[2] *= -1.;
   }
   // this vector is still normed!
}
void cell::set_north_polarity(space &l, double N_weight) {
   // polarity is already set to contain a new random direction
   double north[l.dim];
   for (unsigned short i = 0; i < l.dim; i++) {
      north[i] = 0.;
   }
   if (l.dim == 2) {
      north[1] = -1.;
   }
   if (l.dim == 3) {
      north[2] = -1.;
   }
   // north[] contains a normed vector pointing to the upper part of the reaction volume
   // now combine this with the north vector with weight
   for (unsigned short i = 0; i < l.dim; i++) {
      polarity[i] = (1.0 - N_weight) * polarity[i] + N_weight * north[i];
//      polarity[i] = (1.0 - north_weight) * polarity[i] + north_weight * north[i];
   }
   // ... and norm it:
   double normit = l.get_2norm(polarity);
   for (unsigned short i = 0; i < l.dim; i++) {
      polarity[i] /= normit;
   }
   if (l.cellknot[index].cell==TFR && TFR_north_weight<0) {
       cerr<< "wrong TFR movement setting\n";
       exit(1);
   }
   /*
    * cout<<"D:(";
    * for (unsigned short i=0; i<l.dim; i++) {
    * cout<<polarity[i]<<",";
    * }
    * cout<<")->"<<l.get_2norm(polarity)<<"/////     ";
    */
}
void cell::set_southTfr_polarity(space &l, double N_weight, const signal_molecule &chemokine, double * pol, sigs &s) {
//    cerr<<N_weight<<endl;
    if (N_weight>0) {
//        set_south_polarity(l);
        double north[l.dim];
        for (unsigned short i = 0; i < l.dim; i++) {
           north[i] = 0.;
        }
        if (l.dim == 2) {
           north[1] = -1.;
        }
        if (l.dim == 3) {
           north[2] = 1.;
        }
        // north[] contains a normed vector pointing to the upper part of the reaction volume
        // now combine this with the north vector with weight
        for (unsigned short i = 0; i < l.dim; i++) {
           polarity[i] = (1.0 - N_weight) * polarity[i] + N_weight * north[i];
     //      polarity[i] = (1.0 - north_weight) * polarity[i] + north_weight * north[i];
        }
        // ... and norm it:
        double normit = l.get_2norm(polarity);
        for (unsigned short i = 0; i < l.dim; i++) {
           polarity[i] /= normit;
        }

    } else {set_Tfrchemotaxis_polarity(chemokine, pol, s);}
}
void cell::set_Tfrchemotaxis_polarity(const signal_molecule &chemokine, double * pol, sigs &s) {
  /* Calculates a normed direction vector resulting from the random direction
   * stored in pol and the chemotaxis direction which is determined here.
   * The resulting direction is returned in pol = cell::polarity = vector to change.
   *
   * This routine is only called when the persistence time is over and reorientation is needed.
   * This philosophy corresponds to the idea that persistence time is a cell internal program,
   * which is fulfilled irrespective of changes in the environment.
   */
    //if chemokine level is below the minimum value we can consider to be in the external region
  if ( ((chemokine == CXCL12) && (CXCL12recrit > 0.)
        && s.undercritical_signal(index, chemokine, CXCL12recrit/2)) ||
       ((chemokine == CXCL13) && (CXCL13recrit > 0.)
        && s.undercritical_signal(index, chemokine, CXCL13recrit)) ) {
    responsive2signal[chemokine] = false;
  } else {
    // Find n_chemotaxis:
    double n_chemotaxis[s.dim];
    s.get_gradient(chemokine, index, n_chemotaxis);
    /* n_chemotaxis contains the UNnormed gradient of the chemokine-field.
     * Thus, only the direction of the gradient not its strength matters for polarity.
     */
    // Get its norm
    double gradient = s.get_2norm(n_chemotaxis);
    // and norm the vector:
    for (unsigned short i = 0; i < s.dim; i++) {
      if (gradient > 0.) {
        n_chemotaxis[i] /= gradient;
        n_chemotaxis[i] *= -1; //no matter CXCL12 or 13 Tfr are in the external area (CD25+)
      } else {
        n_chemotaxis[i] = 0.;
      }
    }

    // Calculate relative weight of chemotaxis and random direction:
    double alpha = chemo_max / (1.0 + exp(chemo_steep * (chemo_half - gradient)));
    // this used a sigmoidal function.
    // +++ OPTION: reduce chemotactic answer of CB to semaphorin
    // if (chemokine==SEMA4D && s.knot[index].cell==CB) alpha*=0.3;
    // +++ OPTION end
    /*
     * cout<<"chemokine "<<chemokine
     * <<" gradient="<<gradient
     * <<"; c_max="<<chemo_max
     * <<"; c_steep="<<chemo_steep
     * <<"; c_half="<<chemo_half
     * <<" --> alpha="<<alpha<<"\n";
     */

    // Combine the random and the chemotaxis vectors with that weight:
    for (unsigned short i = 0; i < s.dim; i++) {
      pol[i] += alpha * n_chemotaxis[i];
    }
    // and norm this resulting vector (caution: re-use alpha!!!):
    alpha = s.get_2norm(pol);
    for (unsigned short i = 0; i < s.dim; i++) {
      pol[i] /= alpha;
    }
    // vector pol is changed in the calling routine and can be used there.
  }
}
void cell::set_ext_polarity(space &l, double N_weight) {
   // polarity is already set to contain a new random direction
   double north[l.dim];
   for (unsigned short i = 0; i < l.dim; i++) {
       north[i]=0.;
//      polarity[i]*=-1;
//      cerr<<"pol"<<polarity[i]<<endl<<i<<" "<<l.knot[index].x[i]<<endl;
   }
   if (l.knot[index].x[2]>32) {
       north[2]=1;
   } else {north[2]=-1;}
//   north[2]=-1;
   if (l.knot[index].x[1]>32) {
       north[1]=1;
   } else {north[1]=-1;}
   if (l.knot[index].x[0]>32) {
       north[0]=1;
   } else {north[0]=-1;}
//   if (l.dim == 2) {
//      north[1] = -1.;
//   }
//   if (l.dim == 3) {
//      north[2] = -1.;
//   }
   // north[] contains a normed vector pointing to the upper part of the reaction volume
   // now combine this with the north vector with weight
//   N_weight=.2;
   for (unsigned short i = 0; i < l.dim; i++) {
      polarity[i] = (1.0 - N_weight) * polarity[i] + N_weight * north[i];
//      polarity[i] = (1.0 - north_weight) * polarity[i] + north_weight * north[i];
   }
   // ... and norm it:
   double normit = l.get_2norm(polarity);
   for (unsigned short i = 0; i < l.dim; i++) {
      polarity[i] /= normit;
   }
   if (l.cellknot[index].cell==TFR && TFR_north_weight<0) {
       cerr<< "wrong TFR movement setting\n";
       exit(1);
   }
   /*
    * cout<<"D:(";
    * for (unsigned short i=0; i<l.dim; i++) {
    * cout<<polarity[i]<<",";
    * }
    * cout<<")->"<<l.get_2norm(polarity)<<"/////     ";
    */
}
void cell::set_chemotaxis_polarity(const signal_molecule &chemokine, double * pol, sigs &s) {
  /* Calculates a normed direction vector resulting from the random direction
   * stored in pol and the chemotaxis direction which is determined here.
   * The resulting direction is returned in pol = cell::polarity = vector to change.
   *
   * This routine is only called when the persistence time is over and reorientation is needed.
   * This philosophy corresponds to the idea that persistence time is a cell internal program,
   * which is fulfilled irrespective of changes in the environment.
   */
  if ( ((chemokine == CXCL12) && (CXCL12crit > 0.)
	&& s.overcritical_signal(index, chemokine, CXCL12crit)) || 
       ((chemokine == CXCL13) && (CXCL13crit > 0.)
	&& s.overcritical_signal(index, chemokine, CXCL13crit)) ) {
    responsive2signal[chemokine] = false;
  } else {
    // Find n_chemotaxis:
    double n_chemotaxis[s.dim];
    s.get_gradient(chemokine, index, n_chemotaxis); 
    /* n_chemotaxis contains the UNnormed gradient of the chemokine-field.
     * Thus, only the direction of the gradient not its strength matters for polarity.
     */
    // Get its norm
    double gradient = s.get_2norm(n_chemotaxis);
    // and norm the vector:
    for (unsigned short i = 0; i < s.dim; i++) {
      if (gradient > 0.) {
	n_chemotaxis[i] /= gradient;
      } else {
	n_chemotaxis[i] = 0.;
      }
      if (chemokine == SEMA4D) {
	n_chemotaxis[i] *= -1.;      // cells run away from semaphorins
      }
    }
    
    // Calculate relative weight of chemotaxis and random direction:
    double alpha = chemo_max / (1.0 + exp(chemo_steep * (chemo_half - gradient)));
    // this used a sigmoidal function.
    // +++ OPTION: reduce chemotactic answer of CB to semaphorin
    // if (chemokine==SEMA4D && s.knot[index].cell==CB) alpha*=0.3;
    // +++ OPTION end
    /*
     * cout<<"chemokine "<<chemokine
     * <<" gradient="<<gradient
     * <<"; c_max="<<chemo_max
     * <<"; c_steep="<<chemo_steep
     * <<"; c_half="<<chemo_half
     * <<" --> alpha="<<alpha<<"\n";
     */
    
    // Combine the random and the chemotaxis vectors with that weight:
    for (unsigned short i = 0; i < s.dim; i++) {
      pol[i] += alpha * n_chemotaxis[i];
    }
    // and norm this resulting vector (caution: re-use alpha!!!):
    alpha = s.get_2norm(pol);
    for (unsigned short i = 0; i < s.dim; i++) {
      pol[i] /= alpha;
    }
    // vector pol is changed in the calling routine and can be used there.
  }
}
void cell::set_polarity_velocity(const double &persistence, space &l, sigs &s) {
   set_polarity_velocity(persistence, 1, 0.0, 1, 1, l, s);
}
void cell::set_polarity_velocity(const double &persistence,
                                 const short &v_modus,
                                 const double &p_switch_v,
                                 const short &N_v,
                                 const double &p_slow_factor,
                                 space &l,
                                 sigs &s) {
   // ----------------------------
   // Polarity and velocity states
   // ----------------------------
   // get a random direction respecting persistence
   if (drandom() < persistence) {
      if (use_specific_turning_angles > 0) {
         changed_polarity = l.get_distribution_direction(polarity);
      } else {
         changed_polarity = l.get_random_direction(polarity);
      }

      // +++ OPTION: chose the cells that are chemotactic active:
      if ((s.signal_use[CXCL13] == 1) && (responsive2signal[CXCL13] == true)
              && l.cellknot[index].cell != TFR) {
          set_chemotaxis_polarity(CXCL13, polarity, s);
      }
      if ((s.signal_use[CXCL12] == 1) && (responsive2signal[CXCL12] == true)
              && l.cellknot[index].cell != TFR) {
          set_chemotaxis_polarity(CXCL12, polarity, s);
      }
      // ##### Switch to novel philosophy here (use responsive2signal[] everywhere!)
      if ((s.signal_use[SEMA4D] == 1)
          && ((l.cellknot[index].cell == CB) || (l.cellknot[index].cell == out))) {
         set_chemotaxis_polarity(SEMA4D, polarity, s);
      }
      if ((l.cellknot[index].cell == TC) && (north_weight > 0)) {
         set_north_polarity(l, north_weight);
      }
      //MSchips
      if (l.cellknot[index].cell == TFR) {
          if (TFR_north_weight==-1 ||TFR_north_weight==-5) { //external area of chemokine-rich
              if ((s.signal_use[CXCL13] == 1) && (responsive2signal[CXCL13] == true)) {
                  if (TFR_north_weight==-5){set_chemotaxis_polarity(CXCL13, polarity, s);}
                  else {set_Tfrchemotaxis_polarity(CXCL13, polarity, s);}
              }
              if ((s.signal_use[CXCL12] == 1) && (responsive2signal[CXCL12] == true)) {
                  set_Tfrchemotaxis_polarity(CXCL12, polarity, s);
              }
          } else if (TFR_north_weight==-2) { //only DZ boundary
//              if (l.knot[index].x[2]<=32){
                  set_southTfr_polarity(l,.2,CXCL12, polarity, s);
//              } else {set_southTfr_polarity(l,TFR_north_weight,CXCL12, polarity, s);}
          } else if (TFR_north_weight>0) { //like TFh
              set_north_polarity(l, TFR_north_weight);
          }
      }

      // +++ OPTION end
      // cout<<"polar-vector=("<<polarity[0]<<","<<polarity[1]<<")\n";

      if (v_modus == 3) {
         // the v-state is changed with orientation change for v_modus==3:
         // cout<<"in reorientation: old vstate="<<v_state;
         v_state = get_v_factor(N_v, p_slow_factor);
         // cout<<"  new vstate="<<v_state<<"\n";
      }
      // remember to write the new polarity to TRACK!
      if (trackit) {
         writethis2track = polarisation;
      }
   }
   // else take previous polarity

   // v_modus==1,2: decoupling of v_state from the change of orientation
   if ((N_v > 1) && (v_modus < 3)
       && ((drandom() < p_switch_v)
           || ((v_modus == 2) && (((v_state == 1.) && (drandom() < p_switch_v))
                                  || ((v_state > 1. + 1.e-08)
                                      && (v_state < p_slow_factor - 1.e-08)
                                      && (drandom()
                                          < p_switch_v / 2.)))))) {
      // cout<<"in random v-state: old vstate="<<v_state;
      v_state = get_v_factor(N_v, p_slow_factor);
      // cout<<"  new vstate="<<v_state<<"\n";
   }
   /*
    * v_modus==4:
    * Define v_state on two criteria:
    * 1) is the current cell velocity below some threshold,
    * 2) is the cell in adhesive contact
    * This is done fragment by fragment, which is intuitive, as adhesion
    * is a fragment property and not a cell property!
    * ==> Nothing to do here (v_state=1.0 & p_slow_factor is given to fragdiffuse())
    * // Still to realize ###
    */

   // +++ OPTION: The following line fixes the orientation of the cell to north!
   // polarity[0]=0.; polarity[1]=-1.; // for running on slit
   // polarity[0]=0.95; polarity[1]=-0.31224989991992; // for crawling on wall
   // polarity[0]=-0.9; polarity[1]=-0.435889894; // angle=25.8 to lattice axis
   // +++ end OPTION.

   // ------------------------------------
   // end of polarity and velocity states.
   // ------------------------------------
}
// ============================================================

short cell::do_mutate(AffinitySpace &shape) {
   // short cell::do_mutate(AffinitySpace& shape, const double& p) {
   // Save old cell position in shape space
   long postmp = pos_ss;
   short int err = 1;
   // probabilistic decision if mutation is done:
   if (drandom() < p_mutation) {
      /* find one of the neighbors (-1 is not excluded):
       * knot[i].pos_ss is old position in shapespace
       * .nn[] are the neighbor coordinates in shapespace
       * a random neighbor is used */
      pos_ss = shape.getMutation(pos_ss);
      // save the number of mutations
      ++n_mutation;
      if (n_recycling > 0) {
         ++n_recandmute;
      }
      // shape space:
      shape.rem_cell(sCB, postmp);
      shape.add_cell(sCB, pos_ss);
      shape.rem_cell(total, postmp);
      shape.add_cell(total, pos_ss);
      // cout<<"previous="<<postmp<<" now="<<pos_ss<<"\n";
      err = 0;
   }
   return err;
}
// ============================================================

long cell::find_mitosis_place(const double &p,
                              bool forceit,
                              const double &dx_max,
                              long * pp,
                              space &l) {
   // Assume that cell li at lattice point index has volume 1 and target volume =1 (==CB_vol)
   short err = 2;
   long int j = 0;
   // probabilistic decision if proliferation is done
   if ((forceit) || (drandom() < p)) {
      // If so, try to find place for new cell
      err = 1;

      // Randomize the sequence of possible neighbors
      l.nn = l.nn_permuts.get_rand_set();

      short int n = 0;
      // repeat until its working but not to often
      while (n < l.dim2 && err == 1) {
         // Chose one of its neighbors
         j = l.knot[index].near_n[int (l.nn[n])];
         // check for space at new position (if ok then err=0)
         if ((j != -1) && (l.cellknot[j].cell == nocell)) {
            err = 0;
         }
         // the calling routine has still to put_CB at position j
         // next try
         ++n;
      }
      // Zaehler aktualisieren
      if (err == 0) {
         ++pp[0];
      }
      // Falls das Speichern nicht moeglich war, probiere bei den Diagonalen Nachbarn:
      if (err == 1) {
         short nmax = 4;
         if (l.dim > 2) {
            nmax = 12;
         }
         dynarray<long> dn(nmax, 0, 1);
         random_sequence(dn, nmax);
         n = 0;
         while (n < nmax && err == 1) {
            // Chose one of its neighbors
            j = l.knot[index].diag_n[int (dn[n])];
            // check for space at new position (if ok then err=0)
            if ((j != -1) && (l.cellknot[j].cell == nocell)) {
               err = 0;
            }
            // calling routine has still to do put_CB(j,newCB,shape);
            // next try
            ++n;
         }
         if (err == 0) {
            ++pp[1];
         }
      }
      if (err == 1) {
         // dann wurde auch dort nichts gefunden: suche im ganzen Gitter
         double tmp, min;
         min = double (l.pointnumber);     // also einfach groesser als jeder Abstand
         for (long k = 0; k < l.pointnumber; k++) {
            if (l.cellknot[k].cell == nocell) {
               tmp = l.Abstand(index, k);
               // Die folgende Auswahl ist nicht isotrop und haengt von der
               // Reihenfolge der geprueften k's ab. ### Aber wohl kleiner Effekt!
               // Pruefe, ob das Ziel nah genug dran ist (Grenze dx_max)
               if ((tmp < dx_max) && (tmp < min)) {
                  min = tmp;
                  j = k;
                  err = 0;
               }
            }
         }
         // j enthaelt den naehesten Index zu i!!
         if (err == 0) {
            ++pp[2];
         }
      }

      if (err == 1) {
         ++pp[3];
      }
   }

   if (err > 0) {
      j = (-1) * err;
   }
   // j=index on the lattice if space was found
   // j=-2 if no try due to small p
   // j=-1 if tried but no space available
   return j;
}
// ============================================================

long cell::get_contact(const states &celltype, space &l) {
   /* Checks if there is a cell of type celltype in contact
    * to self-lattice point index. Returns the index of the cell contact or -1.
    * If more than one contact exists, choose randomly.
    */
   short err = 0;
   long c[l.dim2];
   // Pruefe Zelltyp der Nachbarn
   for (int n = 0; n < l.dim2; n++) {
      long j = l.knot[index].near_n[n];
      if ((j != -1) && (l.cellknot[j].cell == celltype)) {
         c[err] = j;
         ++err;
      }
   }
   if (err == 0) {
      return -1;
   } else if (err > 1) {
      return c[irandom(err)];
   } else {
      return c[0];
   }
}
// ============================================================

short int cell::find_contact(states celltype, const long &i, space &l) {
   /* Checks if there is a cell of type celltype in contact
    * to lattice point i.
    * Note: the cell at point i has to be of type different
    * from celltye */
   // ### May be better within space class!
   short err = 1;
   long j;
   // Pruefe Zelltyp der Nachbarn
   for (int n = 0; n < l.dim2; n++) {
      j = l.knot[i].near_n[n];
      if (j != -1) {
         if (l.cellknot[j].cell == celltype) {
            err = 0;
         }
      }
   }
   return err;
}
// ============================================================

double cell::get_ss_receptor_ligand(const double &K, const double &r0, const double &s) {
   /* Calculates the steady-state fraction of liganded receptors
    * weighted with the total number of receptors r0,
    * using dissociation constant K and signal concentration s. */
   return (r0 * s) / (K + s);
}
double cell::get_receptor_ligand(double &kplus, double &kminus, double &s, double &r0,
                                 double &rold) {
   /* Calculates the steady-state fraction of liganded receptors,
   * using dissociation constant K and signal concentration s. */
   // cout<<"k+="<<kplus<<" k-="<<kminus<<" s="<<s<<" r0="<<r0<<" rold="<<rold<<"\n";
   return kplus * s * (r0 - rold) - kminus * rold;
}
// ============================================================

void cell::signal_secretion(const long &i,
                            const signal_molecule &sig_type,
                            const short &mode,
                            double &p,
                            sigs &l) {
   // cout<<"produce "<<p<<" units of "<<sig_type<<" at index "<<i<<" in mode "<<mode<<".\n";
   if (mode == 0) {
      // continuous
      l.signal_put(i, sig_type, p);
      l.signal_total[sig_type] += p;
      l.signal_produced[sig_type] += p;
   } else {
      // vesicles of molecules
      // Ganzzahliger Teil:
      int ip_mksignal = int (p);
      // Rest
      double dp_mksignal = p - double (ip_mksignal);
      if (drandom() < dp_mksignal) {
         ++ip_mksignal;
      }
      l.signal_put(i, sig_type, double (ip_mksignal));
      // Fuer die Statistik: zaehle alle jemals sekretierten Signalmolekuele
      l.signal_total[sig_type] += double (ip_mksignal);
      l.signal_produced[sig_type] += double (ip_mksignal);
   }
}
void cell::use_nutrient(sigs &s, long &sigindex) {
   // Governs the use of nutrients
   /* #### does not correctly distribute the consumption for multi-subunit-cells!
    * Now the nutrients are taken from the grid point "index", i.e. the soma of the cell
    * in the case of multi-subunit-cells.
    * If nutrients are thought to be distributed on various gridpoints, then an
    * additional routine frag_cell::use_nutrient(sigs&) has to be overlayed.
    */
   double total_nutrient = 0.;
   if ((s.signal_use[glucose] == 1) && (status != necrotic)) {
      double use;
      if (status == proliferate) {
         use = use_glucose_pro;
      } else {
         use = use_glucose;
      }
      total_nutrient = s.sigsknot[sigindex].signal_new[glucose] * s.mk_concentration;
      // cout<<"glucose="<<s.sigsknot[sigindex].signal_new[glucose]<<", use="<<use<<"\n";
      // short all=
      s.signal_get(sigindex, glucose, use);
      // if (all=1) {} // if state=proliferate set state to viable, if viable set to necrotic
      // else {} // set state to proliferate
   }
   if ((s.signal_use[oxygen] == 1) && (status != necrotic)) {
      double use;
      if (status == proliferate) {
         use = use_oxygen_pro;
      } else {
         use = use_oxygen;
      }
      if (s.signal_use[glucose] == 0) {
         total_nutrient = 1.;
      }
      total_nutrient *= s.sigsknot[sigindex].signal_new[oxygen] * s.mk_concentration;
      // short all=
      s.signal_get(sigindex, oxygen, use);
      // if (all=1) {} // if state=proliferate set state to viable, if viable set to necrotic
      // else {} // set state to proliferate
   }
   // cout<<"total="<<total_nutrient<<", critical="<<critical_nutrient<<"\n";
   if (total_nutrient < critical_nutrient) {
      status = necrotic;
   }
}
// ============================================================
// ============================================================
// ====================== FRAG_CELL ===========================
// ============================================================
// ============================================================

double frag_cell::frac_average = 0.;
/* OPTION MOVE_ANALYSIS:
 * long frag_cell::n_try_move=0;
 * long frag_cell::n_move_done=0;
 * long frag_cell::n_move_removed=0;
 * long frag_cell::n_move_forbidden=0;
 * long frag_cell::n_move_self=0;
 * long frag_cell::n_move_forbidden_back=0;
 * long frag_cell::n_move_self_back=0;
 */

frag_cell::frag_cell() {
   // cout<<"in frag_cell default constructor ...";
   volume = 0;
   radius = 0;
   borderpoints = 1;
   // longaxis=0;
   // shortaxis=0;
   elongation = 1.;
   fragments = new long[FRAGMENT_STEP];
   t_immobile = new long[FRAGMENT_STEP];
   for (int j = 0; j < FRAGMENT_STEP; j++) {
      fragments[j] = 0;
      t_immobile[j] = 0;
   }
   for (int i = 0; i < 3; i++) {
      barycenter[i] = 0.0;
   }

   // for deformation properties
   n_moves = 0;
   alpha = 1.;
   alpha_mean10 = 1.;
   alpha_mean = 1.;
   n_directed_moves = 0;
   flag_no_correction = 0.;
   performed2aimed_move = 1.;
   performed2aimed_move_now = 1.;
   // cout<<" end frag_cell default constructor.\n";
}
// ============================================================

frag_cell::frag_cell(const frag_cell &x) : cell(x) {
   /// Philippe Important : cell() is to avoid the warning. I hope it is fine to call this mother
   // constructor in particular
   // cout<<"in frag_cell(const frag_cell&) constructor ...";
   fragments = new long[FRAGMENT_STEP];
   t_immobile = new long[FRAGMENT_STEP];
   operator =(x);
   // cout<<" end of frag_cell(const frag_cell&)\n";
}
frag_cell&frag_cell::operator =(const frag_cell &x) {
   // cout<<"In frag_cell::operator=() ... ";
   cell::operator =(x);
   // volume=x.volume; now in class cell
   radius = x.radius;
   borderpoints = x.borderpoints;
   // longaxis=x.longaxis;
   // shortaxis=x.shortaxis;
   elongation = x.elongation;
   for (int j = 0; j < FRAGMENT_STEP; j++) {
      fragments[j] = x.fragments[j];
      t_immobile[j] = x.t_immobile[j];
   }
   for (int i = 0; i < 3; i++) {
      barycenter[i] = x.barycenter[i];
   }

   // for deformation properties
   n_moves = x.n_moves;
   alpha = x.alpha;
   alpha_mean10 = x.alpha_mean10;
   alpha_mean = x.alpha_mean;
   n_directed_moves = x.n_directed_moves;
   flag_no_correction = x.flag_no_correction;
   performed2aimed_move = x.performed2aimed_move;
   performed2aimed_move_now = x.performed2aimed_move_now;
   // cout<<" end of frag_cell::operator=() \n";
   return *this;
}
// ============================================================

frag_cell::~frag_cell() {
   // cout<<"in frag_cell destructor ...\n";
   // delete[] fragments;
   /* Warum gibt es einen Speicherzugriffsfehler am Programmende, wenn man delete[]
    * hier aufruft. Das gleiche passiert mit dem gleichen deletet
    * in den Destruktoren der abgeleiteten Klassen. ???
    * Ohne diesen Befehl wird der Speicherbedarf erhoeht! Dies habe ich
    * jetzt vermieden, indem die Routine deliberate_memory() nach der
    * Erzeugung von neuen Objekten immer aufgerufen wird.
    */
}
void frag_cell::deliberate_memory() {
   delete[] fragments;
   delete[] t_immobile;
}
// ============================================================

void frag_cell::show_fragments() {
   cout << "fragment-list[]= ";
   for (int a = 0; a < volume; a++) {
      cout << fragments[a] << " ";
   }
   cout << "; length=" << volume << "\n";
}
void frag_cell::show_fragment_immobility() {
   cout << "fragment-immobility[]= ";
   for (int a = 0; a < volume; a++) {
      cout << t_immobile[a] << " ";
   }
   cout << "; length=" << volume << "\n";
}
void frag_cell::show_barycenter(short unsigned &d) {
   cout << "barycenter=(";
   for (short f = 0; f < d; f++) {
      cout << barycenter[f] << ",";
   }
   cout << ")\n";
}
// ============================================================

void frag_cell::set_clock() {
   ++clock;
   for (int a = 0; a < volume; a++) {
      ++t_immobile[a];
   }
}
void frag_cell::reset_clock(const int &i) { t_immobile[i] = 0; }

// ============================================================

int frag_cell::where_fragment(const long &i) {
   // returns the position of the fragment at lattice position i
   // on the fragment list frags
   long n = 0;
   long res = -1;
   while (n < volume && res == -1) {
      if (fragments[n] == i) {
         res = n;
      }
      ++n;
   }
   return res;
}
int frag_cell::where_fragment(const long &i, long * frags, const int &max) {
   // return the position of i on the fragment list frags
   // Falls frags und max die Fragmente und das Volumen zu einer Zelle sind
   // lieber where_fragment(const long&) aufrufen!
   long n = 0;
   long res = -1;
   while (n < max && res == -1) {
      if (frags[n] == i) {
         res = n;
      }
      ++n;
   }
   return res;
}
// ============================================================

void frag_cell::del_fragment(const long &nr, long * frags, int &max) {
   // deletes fragment nr of the fragment list frags* and reduces max by one
   frags[nr] = frags[max - 1];
   --max;
}
void frag_cell::del_fragment(const long &nr) {
   // same as above but acting on the cell fragment list and volume
   fragments[nr] = fragments[volume - 1];
   --volume;
}
// ============================================================

int frag_cell::check_connection(states celltype,
                                const long &i,
                                const long &li,
                                long * unconnected,
                                int &n_unconnected,
                                space &l) {
   /* Achtung: the field "unconnected" and "n_unconnected" are changed here!!!
    * Starts from a list "unconnected" of potentially unconnected fragments.
    * Starts at lattice point "i" and looks at all neighbours of i.
    * Checks if these neighbours belong to the same cell with list index "li".
    * If so and if present on the unconnected list: delete it from there and
    * restart the routine from this new lattice point.
    * If not continue with the next of the previously defined next-neighbors.
    * The result is an unconnected-list of length n_unconnected containing only
    * these fragments with no connection to lattice point i.
    */
   // Go through all neighbors
   long k;
   int fragk;
   for (int n = 0; n < l.dim2; n++) {
      if (n_unconnected > 0) {
         // get index of a neighbor of jcell at lattice point j
         k = l.knot[i].near_n[n];
         // cout<<"  check lattice point k="<<k<<"...\n";
         // Continue if k belongs to cell li only
         if ((k != -1) && (l.cellknot[k].cell == celltype) && (l.cellknot[k].listi == li)) {
            // cout<<"    delete this from fragment list:\n";
            // cout<<"    list= ";
            // for (int a=0; a<n_unconnected; a++) cout<<unconnected[a]<<" ";
            // cout<<"; length="<<n_unconnected<<"\n";
            // Delete this fragment from the unconnected-list
            fragk = where_fragment(k, unconnected, n_unconnected);
            // cout<<"    fragk="<<fragk<<"\n";
            if (fragk != -1) {
               del_fragment(fragk, unconnected, n_unconnected);
               // restart the routine from lattice point k
               // cout<<"    call again check_connection ...\n";
               check_connection(celltype, k, li, unconnected, n_unconnected, l);
            }
         }
      }
   }
   return n_unconnected;
}
int frag_cell::check_connection(states celltype, const long &li, space &l) {
   // As the above check_connection but the original cell remains unchanged!
   // This routine is thought for the check of a final fragment list!
   int tmpvol = volume;
   long * tmpfrags;
   tmpfrags = new long[FRAGMENT_STEP];
   for (int j = 0; j < FRAGMENT_STEP; j++) {
      tmpfrags[j] = fragments[j];
   }
   int fragk = where_fragment(index, tmpfrags, tmpvol);
   if (fragk == -1) {
      cout << "cell does not contain its own center!!!";
      exit(1);
   } else {
      del_fragment(fragk, tmpfrags, tmpvol);
   }
   int unconn = check_connection(celltype, index, li, tmpfrags, tmpvol, l);
   if (unconn != 0) {
      cout << "Unconnected cell (the original) results from proliferation!";
      exit(1);
   }
   delete[] tmpfrags;
   return unconn;
}
// ============================================================

void frag_cell::attribute_neighbors(states celltype,
                                    const long &li,
                                    const long &j,
                                    const long &lj,
                                    frag_cell &targetcell,
                                    space &l) {
   /* Finds those neighbors of the "targetcell" (list index "lj" and type "celltype")
    * at lattice point "j" that belong to the fragment list of the calling cell
    * (list index "li" and same type) and that are nearer to the center of
    * the "targetcell" than to the one of the poolcell.
    * Saves them on the fragment list of the target cell.
    * Deletes them from the fragment list of the calling cell.
    * Recursively looks for neighbors of the new fragments found.
    */
   long k;
   // cout<<"\nAttribute neighbors of lattice point "<<j<<":\n";
   for (int n = 0; n < l.dim2; n++) {
      // get index of a neighbor of the cell at lattice point index
      k = l.knot[j].near_n[n];
      /* Only continue if k is nearer from the center of the cell than from
       * the poolcell, and if k belongs to the poolcell only.
       * Note, that if k belongs to attributed cell, it has already been processed.
       * Otherwise it does not belong to the cell or the poolcell and has not to be processed.
       */
      if ((k != -1) && (l.cellknot[k].cell == celltype)   // is of correct type
          &&
          (l.cellknot[k].listi == li)   // belongs to the calling cell
          &&
          (l.Abstand(k, index) > l.Abstand(k, targetcell.index))) {
         // ### possible anisotropy: comparing to "index" instead of real "barycenter"!
         // Save this fragment on the fragment list of the target cell
         targetcell.fragments[targetcell.volume] = k;
         ++targetcell.volume;
         // Change the map to the cell list
         l.cellknot[k].listi = lj;
         // Take care that the status in the grid class is correct !!! ???
         /* Don't see any problem, as the status of both involved cells is the same
          * and this applies to all subunits of it.
          * --> Nothing to do here!
          */
         // Delete this fragment from the list of the calling cell
         int fragk = where_fragment(k);
         del_fragment(fragk);
         // restart the routine from lattice point k
         attribute_neighbors(celltype, li, k, lj, targetcell, l);
      }
   }
   // cout<<"  ... end of attribute.\n";
}
// ============================================================

long frag_cell::find_nn2cell(space &l, long * k) {
   /*
    * starts from lattice point "k"
    * and find the nearast fragment belonging to the (calling) cell
    */
   double dist = l.dx * l.pointnumber;   // larger than lattice
   // double dist=2.0*l.dx*l.prodim; // larger than lattice
   double newdist;
   long target = -1;
   // long t[l.dim];
   for (int m = 0; m < volume; m++) {
      // l.get_koord(fragments[m],t);
      // newdist=l.get_2norm(k,t);
      newdist = l.get_2norm(k, l.knot[fragments[m]].x);
      if (newdist < dist) {
         dist = newdist;
         target = fragments[m];
      }
   }
   return target;
   // =-1, falls nichts gefunden!
}
// ============================================================

short int frag_cell::contact(states celltype, space &l) {
   /*
    * as "find_contact(states,const long&, space&)" but for
    * cells with more than one fragment. It goes through all
    * fragments.
    */
   short err = 1;
   int f;
   // Pruefe Zelltyp der Nachbarn aller Fragmente
   f = 0;
   while (err > 0 && f < volume) {
      err = find_contact(celltype, fragments[f], l);
      ++f;
   }
   return err;
}
// ============================================================

long frag_cell::get_barycenter(space &l) {
   /*
    * Calculates the barycenter of cell li, saves it (as double value), and
    * determines the lattice point nearest to the barycenter of cell li
    * that belongs to the same cell
    */
   // Define a sum vector s (long) and sd (double)
   long s[l.dim];
   double sd[l.dim];
   int m;
   // initialize the sum vectors:
   for (m = 0; m < l.dim; m++) {
      s[m] = 0;
      sd[m] = 0.;
   }
   // Go through all points of the cell and add
   for (m = 0; m < volume; m++) {
      // add the point-vector of the fragment and add to "s"
      for (short d = 0; d < l.dim; d++) {
         s[d] += l.knot[fragments[m]].x[d];
      }
   }
   // calculate barycenter coordinates as doubles and save them in the cell
   for (short d = 0; d < l.dim; d++) {
      sd[d] = double (s[d]) / double (volume);
      barycenter[d] = sd[d];
   }
   // calculate corresponding integer coordinates
   for (short d = 0; d < l.dim; d++) {
      s[d] = long (sd[d] + 0.5);

      /*
       * // The following avoids possible anisotropies (which did not turn out to be relevant):
       * short knapp=0;
       * if ((double(s[d])-sd[d])<0) { // abgerundet
       * if (double(s[d])-sd[d]+0.5<1.e-08) knapp=-1; //sd knapp drunter
       * }
       * else { // aufgerundet
       * if (sd[d]+0.5-double(s[d])<1.e-08) knapp=1; // sd knapp drueber
       * }
       * if (knapp!=0) {
       * // correct decision according to random decision
       * int limit_decision=irandom(2);
       * if (knapp==-1 && limit_decision==1) ++s[d];
       * if (knapp==1  && limit_decision==0) --s[d];
       * }
       * // end of avoid anisotropy.
       */
   }
   // find nearest neighbour belonging to the cell and save its index on the lattice
   long s_index = find_nn2cell(l, s);
   if (s_index == -1) {
      cout << "No lattice point for cell found in get_barycenter->find_nn2cell.\n";
      exit(1);
   }

   // Older version:
   /*
    * // Restore index-form
    * long s_index=l.Index(s);
    * // Check if the cell at s_index belongs to cell li
    * if (l.knot[s_index].cell!=celltype || l.knot[s_index].listi!=li) {
    * // otherwise find nearest neighbor belonging to the cell
    * // Achtung: Eigentlich koennte man hier find_nn2cell gleich aufrufen.
    * // Dadurch wuerde man die Abfrage mit li sparen, li muesste nicht
    * // uebergeben werden, es muesste Index(s) nicht berechnet werden,
    * // der Abruf von knot.xxx waere beseitigt, und vor allem
    * // die Kontrolle ob Gitterplatz zu CB gehoert wuerde wegfallen.
    * // Dies wuerde die Uebergabe des Zelltyps sparen. CPU-time???
    * // Realisiert am 25.4.2004
    * long s_index=find_nn2cell(l,s);
    * if (s_index==-1) {
    *  cout<<"No lattice point for cell "<<li<<" found in get_barycenter->find_nn2cell.\n";
    *  exit(1);
    * }
    * }
    */
   return s_index;
}
long frag_cell::find_last_fragment(const long &li,
                                   const states &celltype,
                                   space &l,
                                   double * direction) {
   // Looks for the most distant fragments in direction "direction"
   // Start from the barycenter:
   double v_a[l.dim];
   for (short i = 0; i < l.dim; i++) {
      v_a[i] = barycenter[i];
   }
   // Walk through the cell in direction "direction":
   short stop = 0;
   long lastincell = index;
   while (stop == 0) {
      for (short i = 0; i < l.dim; i++) {
         v_a[i] += direction[i];
      }
      long newindex = l.Index(v_a);
      if (newindex == -1) {
         stop = 1;
      } else if (l.self(newindex, li, celltype) == 0) {
         stop = 1;
         /*
          * // Try to find a neighbour of newindex with positive scalarproduct and l.self(..)==1
          * short j=0;
          * while (stop==0 && j<l.dim2) {
          * long tmpindex=l.knot[newindex].near_n[j];
          * if (tmpindex!=-1) {
          *  // check the scalarproduct with "direction"
          *  long tmpvector[l.dim];
          *  l.get_koord(tmpindex,tmpvector);
          *  for (short k=0; k<l.dim; k++) tmpvector[k]-=long(v_a[k]+0.5);
          *  if (l.get_scalarproduct(tmpvector,direction)>0 && l.self(tmpindex,li,celltype)==1) {
          *    stop=0;
          *    newindex=tmpindex;
          *    lastincell=newindex;
          *    for (short k=0; k<l.dim; k++) v_a[k]=tmpvector[k];
          *  }
          * }
          * ++j;
          * }
          * if (stop==1) {
          * // try the next point in direction "direction":
          * for (short i=0; i<l.dim; i++) v_a[i]+=direction[i];
          * newindex=l.Index(v_a);
          * if (l.self(newindex,li,celltype)==1) {
          *  stop=0;
          *  lastincell=newindex;
          * }
          * }
          */
         /* The two variants of refinement of the criterium for finding the last
          * fragment in a specific "direction" turned out to be without major
          * effect on the resulting graph "axis.eps". Thus the last fragment of
          * the cell is now searched with the simplest algorithm in order to save
          * CPU time.
          */
      } else {
         lastincell = newindex;
      }
      /* Note that at this point one may allow some deviation of the optimal path
       * in direction "direction". Then if a point failed to belong to the cell,
       * neighbours of this point may be checked. v_a would have to be adjusted
       * correspondingly as new starting point.
       * ==> A corresponding refinement turned out to be irrelevant (see above)!
       */
   }
   // lastincell is either index if no further fragments exist in direction "direction"
   // or the last point within the cell in direction "direction"
   // if (celltype==out) cout<<"lastincell="<<lastincell<<" and \n";
   return lastincell;
}
double frag_cell::get_axis_length(const long &li,
                                  const states &celltype,
                                  space &l,
                                  double * direction) {
   /*
    * cout<<"frags=(";
    * for (short b=0; b<volume; b++) cout<<fragments[b]<<",";
    * cout<<")\n";
    */

   long lastincell = find_last_fragment(li, celltype, l, direction);
   for (short i = 0; i < l.dim; i++) {
      direction[i] *= -1.;
   }
   long firstincell = find_last_fragment(li, celltype, l, direction);
   // restore "direction" (change relevant in calling routine)
   for (short i = 0; i < l.dim; i++) {
      direction[i] *= -1.;
   }

   // calculate the distance between first- and last-incell and return this
   return l.Abstand(lastincell, firstincell) + 1.;
   /* Note, that "+1." is necessary, as the distance between the two points on the axis
    * are underestimating the length of the cell. In the case of a one fragment cell,
    * the distance would become zero without this. But the extension of the cell is
    * one lattice constant even in this case. For longer distances "+1." signifies
    * that the first as well as the last fragment are included into the returned
    * distance.
    */
}
double frag_cell::get_long2short_axis(const long &li,
                                      const states &celltype,
                                      space &l,
                                      double &long_axis,
                                      double &short_axis) {
   long_axis = get_axis_length(li, celltype, l, polarity);
   // define a set of vector
   double perp[l.dim][l.dim];   // erster Index zaehlt die Vektoren, zweiter die Komponenten
   double lis[l.dim][l.dim];    // linear independent system
   // der erste normierte Vektor ist polarity
   // erstelle ein System linear unabhaengiger Vektoren (unnormiert):
   // get a direction of polarity that is unequal to zero:
   short dimdone = 0;
   for (short i = 0; i < l.dim; i++) {
      if (polarity[i] != 0) {
         dimdone = i;
      }
   }
   // erstelle ein (1,0,0)-System von Vektoren
   for (short n = 0; n < l.dim; n++) {
      for (short i = 0; i < l.dim; i++) {
         lis[n][i] = 0.;
         if (i == n) {
            lis[n][i] = 1.;
         }
      }
   }
   // ersetze den Vektor "n=dimdone" durch "polarity" und setze diesen als ersten Vektor
   for (short i = 0; i < l.dim; i++) {
      lis[dimdone][i] = lis[0][i];
      lis[0][i] = polarity[i];
   }

   // cout<<"pola=("<<polarity[0]<<","<<polarity[1]<<")   ";
   // cout<<"lis0=("<<lis[0][0]<<","<<lis[0][1]<<")   ";
   // cout<<"lis1=("<<lis[1][0]<<","<<lis[1][1]<<")\n";

   // mache mit diesen Vektoren das Schmidtsche Orthogonalisierungsverfahren
   // lis[0][i] is polarity and is normed
   for (short i = 0; i < l.dim; i++) {
      perp[0][i] = lis[0][i];
   }
   // calculate perp[n]:
   for (short n = 1; n < l.dim; n++) {
      double h_n[l.dim];
      for (short i = 0; i < l.dim; i++) {
         h_n[i] = lis[n][i];
      }
      for (short m = 0; m < n; m++) {
         // skalarproduct of lis[n] with perp[m]
         double project = 0.;
         for (short i = 0; i < l.dim; i++) {
            project += lis[n][i] * perp[m][i];
         }
         for (short i = 0; i < l.dim; i++) {
            h_n[i] -= (project) * perp[m][i];
         }
      }
      // save normed h_n in perp[n]:
      double norm_h_n = l.get_2norm(h_n);
      for (short i = 0; i < l.dim; i++) {
         perp[n][i] = h_n[i] / norm_h_n;
      }
   }
   // perp[n][i] contains: for n=0 -> polarity;
   //                      for n>0 and n<l.dim normed vectors perpendicular to polarity

   // cout<<"perp0=("<<perp[0][0]<<","<<perp[0][1]<<")   ";
   // cout<<"perp1=("<<perp[1][0]<<","<<perp[1][1]<<")\n";

   /*
    * // Do a check of the result (deactivated)
    * for (short n=1; n<l.dim; n++) {
    * double prod=0.;
    * for (short i=0; i<l.dim; i++) prod+=polarity[i]*perp[n][i];
    * if (prod>1.e-08 || prod<-1.e-08) {
    *  cout<<"Schmidt's orthogonalisation failed in get_long2short_axis(...)!\n";
    *  exit(1);
    * }
    * if (l.dim==3) {
    *  prod=0.;
    *  for (short i=0; i<l.dim; i++) prod+=perp[1][i]*perp[2][i];
    *  if (prod>1.e-08 || prod<-1.e-08) {
    *    cout<<"Schmidt's orthogonalisation failed in get_long2short_axis(...), 3D!\n";
    *    exit(1);
    *  }
    * }
    * }
    * // end of check
    */

   // cout<<"perp1=("<<perp[1][0]<<","<<perp[1][1]<<")\n";
   short_axis = get_axis_length(li, celltype, l, perp[1]);
   if (l.dim == 3) {
      double short_axis2 = get_axis_length(li, celltype, l, perp[2]);
      if (short_axis2 < short_axis) {
         short_axis = short_axis2;
      }
   }
   // "short_axis" contains the shorter axis perpendicular to "polarity"

   // Comment:
   // But there may be shorter axes other than the one checked! Important?
   // A more sophisticated search for long and short axis did not substantially
   // change the result shown in axis.eps. Thus it may be fine as it is.

   // if the polarity of the cell has just changed, the long axis may be shorter
   // ### check, how often this happens
   /*  if (long_axis<short_axis) return short_axis/long_axis;
    * return long_axis/short_axis;
    */
   return long_axis / short_axis;   // reinstalled in the simple version 2/2007
   /* ### Eventually, better to check more axis in various directions covering the whole cell!?
    * This may become superfluous if in the get_axis_length-routine the direction of walk
    * is defined more flexible.
    */
   /* The whole thing here has to be rethought of for the analysis
    * of shape to velocity correlation. ###
    */
}
double frag_cell::get_radius(const short unsigned &d) {
   // in units of the lattice constant for a hypothetic spherical shape
   if (d == 2) {
      radius = sqrt(double (volume) / 3.141592654);
   } else {
      radius = pow(double (volume) / 4.1888, 1. / 3.);
   }
   return radius;
}
double frag_cell::get_elongation(const long &li, const states &celltype, space &l) {
   long i = find_last_fragment(li, celltype, l, polarity);
   // This calculates the radius of the cell in direction "polarity"!
   elongation = l.Abstand(i, index) + 0.5;
   // if (celltype==out) cout<<"elongation="<<elongation<<" radius="<<radius<<"\n";
   elongation /= radius;
   return elongation;
}
double frag_cell::get_long2short_axis(const long &li, const states &celltype, space &l) {
   double long_axis = 0.;
   double short_axis = 0.;
   return get_long2short_axis(li, celltype, l, long_axis, short_axis);
}
double frag_cell::get_reshaping_force(const double &p_max_tension, const double &K_elongation) {
   double p = 0.;
   if (elongation > 1.) {
      p = p_max_tension * (elongation - 1.) / (elongation + K_elongation - 2.);
   }
   /*
    * cout<<"reshape: p_max="<<p_max_tension
    *  <<", r="<<radius
    *  <<", elongation="<<elongation
    *  <<", p_now="<<p
    *  <<"\n";
    */
   return p;
}
void frag_cell::check_barycenter(space &l) {
   // void frag_cell::check_barycenter(const states& celltype, const long& li, space& l) {
   // check consistency of barycenters (optional control)
   double oldbary[3];
   for (short x = 0; x < 3; x++) {
      oldbary[x] = barycenter[x];
   }
   get_barycenter(l);
   for (short x = 0; x < 3; x++) {
      if (oldbary[x] != barycenter[x]) {
         cout << "Inconsistent barycenter.\n";
         exit(1);
      }
   }
}
short frag_cell::enclosed_object(const long &source,
                                 const long &target,
                                 const long &li,
                                 const states &celltype,
                                 space &l) {
   /*
    * Checks, if the movement of a fragment of the cell defined by (li,celltype)
    * from lattice point source to target leads to the enclosure of a one-fragment
    * object.
    */
   // cerr<<"start enclosed_object(...) ...\n";
   if (target == source) {
      return 0;
   }
   // Go through all neighbours of target
   for (short n = 0; n < l.dim2; n++) {
      // Check only if the neighbour 1) does not belong to the cell itself
      //                         and 2) is not "empty"
      //                         and 3) is not "source"
      long neighbour = l.knot[target].near_n[n];
      // cerr<<"work on neighbour="<<neighbour<<":\n";
      if ((neighbour != -1) && (l.cellknot[neighbour].cell != nocell)
          && (l.self(neighbour, li, celltype) == 0)
          && (neighbour != source)) {
         // so we are dealing with a non-self neighbour:
         /* If "neighbour" has dim2 neighbours belonging to the cell
          * return 1. "target" is included in this check, "source" not.
          */
         short selfnns = 0;
         for (short j = 0; j < l.dim2; j++) {
            long k = l.knot[neighbour].near_n[j];
            if (((l.self(k, li, celltype) == 1) || (k == target)) && (k != source)) {
               ++selfnns;
            }
         }
         // cerr<<"end enclosed_object(...)\n";
         if (selfnns == l.dim2) {
            return 1;
         }
      }
   }
   // cerr<<"end enclosed_object(...)\n";
   return 0;
}
short frag_cell::check_for_ring(const states &celltype, const long &li, space &l) {
   /*
    * Checks if the barycenter belongs to the cell itself.
    * if not, 1 is returned.
    *
    * !!! Not in use !!!
    */
   // calculate integer coordinates of the barycenter:
   long s[l.dim];
   for (short d = 0; d < l.dim; d++) {
      s[d] = long (barycenter[d] + 0.5);
   }
   // get corresponding index, nicht ohnehin bekannt?, nein (siehe naechster Kommentar)
   long i = l.Index(s);
   /* Remark: the result of get_barycenter cannot be used as it may differ
    * from i when the barycenter is not occupied by a cell fragment,
    * which is exactly the interesting case here. */
   if ((l.cellknot[i].cell == celltype) && (l.cellknot[i].listi == li)) {
      return 0;
   } else {
      /* Here one may incorporate a more sophisticatd criterion for
       * the presence of a ring structure. With the present criterion,
       * all movements that lead to a barycenter outside the cell
       * are strictly forbidden. This prevents strongly deformed cells
       * from being built. A test for a really appearing ring structure
       * would enlarge the set of possible cell forms.
       ###
       ##################(but note, that this routine is not in use 25.4.2004)
       */
      return 1;
   }
}
// ============================================================

int frag_cell::get_free_nn(space &l,
                           const long &li,
                           long * n_list,
                           const double &tolerance_min,
                           const double &tolerance_steepness,
                           double * target_point,
                           const long &excluded,
                           const long &include) {
/*
 * Da hier das barycenter der Zellen verwendet wird, kann man diese
 * Routine nicht zu lattice.h tun, wo die Zellen nicht bekannt sind. Damit
 * wird diese Routine eine Zelleigenschaft und sollte in die Zell-Klassen
 * integriert werden.
 */
/* Returns a list "n_list" of points on the lattice that are all neighbours
 *   of points that belong to the cell "li" represented in the list of
 *   fragments "fragments" (of length "volume").
 * "tolerance_min" is a criterion to include the free points on the list.
 *   For tolerance_min=1.e-07 only the nearest free neighbors with respect to
 *   the barycenter are included (non-deforming limit).
 *   For tolerance_min=1 all are included. A double in between denotes
 *   the fraction of farest-nearest that is tolerated for inclusion.
 * "tolerance_steepness" accounts for a adaptive treatment of the range of
 *   neighbours considered as target points. It denotes the deformation
 *   at which the tolerance is set to the average of "1" and "tolerance_min".
 *   For "tolerance_steepness==0.01" the non-deforming limit is taken.
 * It is possible to exclude one lattice point from the fragment list
 *   from considering its neighbors to belong to the list n_list. This is
 *   helpful if this "excluded" point is considered to move somewhere.
 *   In this way the connectivity of the cell is maintained. However, this
 *   is not useful if the cell consists of one fragment only. Then the
 *   excluded point is nevertheless considered. (no exclusion, call with -1!)
 * It is also possible to "include" a specific point as free (avoid by
 *   calling with -1). This is meaningful if a point considered to move
 *   is thought to stay if it is in a reasonable position. This is not
 *   considered for volume==1.
 */
/* 3/2004: "target_point" wird hinzugefuegt. Damit muss die Auswahl nicht
 *   auf der Basis des wirklichen Schwerpunkts geschehen, sondern kann auch
 *   auf der Basis eines anderen Punktes im Raum vorgenommen werden.
 */
// cout<<"In get_free_nn() ... ";

   long j;
   short jj;
   long i, k;
   int used = 0;
   double dist;
   double nearest = l.dx * l.pointnumber;
   // double nearest=2.*l.dx*l.prodim;
   double farest = 0.;
   short ex_holes = 0;
   double center_dist[(l.dim2 - 1) * volume + 2];
   // cout<<"target=("<<target_point[0]<<","<<target_point[1]<<")\n";

   if ((include != -1) && (volume > 1)) {
      // do the same for the point "include" as below for free neighbors
      n_list[used] = include;
      ++used;
      dist = l.get_2norm(target_point, l.knot[include].x);
      center_dist[used - 1] = dist;
      if (dist > farest) {
         farest = dist;
      }
      if (dist < nearest) {
         nearest = dist;
      }
   }
   for (j = 0; j < volume; j++) {
      // get the index of j's fragment
      k = fragments[j];    // is one of the lattice indices

      /*
       * long tmp[l.dim];
       * l.get_koord(k,tmp);
       * cout<<"  consider fragment "<<k<<": ("<<tmp[0]<<","<<tmp[1]<<")\n";
       */

      if ((k != excluded) || (volume == 1)) {
         // Go through its neighbors
         for (jj = 0; jj < l.dim2; jj++) {
            i = l.knot[k].near_n[jj];
            if ((i != -1) && (l.cellknot[i].cell == nocell)) {
               // check if i is already on the result list
               short already_there = 0;
               for (int kkk = 0; kkk < used; kkk++) {
                  if (n_list[kkk] == i) {
                     already_there = 1;
                  }
               }
               // if not ...
               if (already_there == 0) {
                  n_list[used] = i;
                  ++used;
                  // Calculate distance to (target) center
                  dist = l.get_2norm(target_point, l.knot[i].x);

                  /*
                   * l.get_koord(i,tmp);
                   * cout<<"  consider target point "<<i<<": ("<<tmp[0]<<","<<tmp[1]
                   *  <<") dist="<<dist<<"\n";
                   */

                  // dist=Abstand(CB_list[li].index,i);
                  /* Die auskommentierte Variante fuehrte in einigen Konstellation
                   * zu einer Anisotropie in der Bewegung. Dies liegt an einem
                   * konsequenten Fehler bei der Projektion des berechneten Schwerpunkts
                   * auf die Gitterpunkte. */
                  center_dist[used - 1] = dist;
                  if (dist > farest) {
                     farest = dist;
                  }
                  if (dist < nearest) {
                     nearest = dist;
                  }
                  // get the number of self-neighbors (k belongs to the cell!)
                  if (l.get_n_self(i, li, l.cellknot[k].cell, -1) == l.dim2) {
                     ex_holes = 1;
                  }
               }
            }
         }
      }
   }
   // save the total number of neighbours found:
   alpha = double (used);
   // if holes exist select only those neighbors that are holes
   if (ex_holes == 1) {
      j = 0;
      while (j < used) {
         // Delete all places from the n_list that have less than l.dim2 neighbors of the cell:
         if (l.get_n_self(n_list[j], li, l.cellknot[fragments[0]].cell, -1) != l.dim2) {
            --used;
            n_list[j] = n_list[used];
         } else {
            ++j;
         }
      }
      //    cout<<"only holes: ";
   } else {
      // select only those neighbors that are near enough according to the tolerance parameter
      j = 0;
      double range = farest - nearest;
      double deform;    // deformation measure
      if (nearest > 0) {
         deform = range / nearest;     // 0: sphere; infty: strongly deformed
      } else {
         deform = 10000. * range;     // large deformation
      }
      double take_tolerance = (exp(-deform / tolerance_steepness) + tolerance_min)
                              / (1.0 + tolerance_min * exp(-deform / tolerance_steepness));
      range *= take_tolerance;
      range = nearest + range + 1.e-08;
      // range saves the limit up to which distance free neighbors are included
      while (j < used) {
         if (center_dist[j] > range) {
            --used;
            n_list[j] = n_list[used];
            center_dist[j] = center_dist[used];
         } else {
            ++j;
         }
      }
      // this routine includes at least one free point (if existent) on the result list!
      /* ##### This is a problem: If volume==1 a point in the back is also chosen,which
       * is not the case when do_diffuse(...) is used instead! Thus, CB and CC are
       * treated inconsistently.
       */
   }

   // deformation parameter:
   // alpha becomes the fraction of used free neighbours here
   alpha = double (used) / alpha;
   // alpha_mean and alpha_mean10 are changed in fracdiffuse();

   /*
    * cout<<"  take: ";
    * for (int tt=0; tt<used; tt++) cout<<n_list[tt]<<" ";
    * cout<<"\n";
    */

   return used;
}
// ============================================================

double frag_cell::bind_ss_receptor_ligand(const signal_molecule &sig,
                                          double &K,
                                          double &r0,
                                          double &old,
                                          space &l,
                                          sigs &s) {
   double n = 0.0;
   double average = 0.0;
   for (int f = 0; f < volume; f++) {
      if (l.object_border(fragments[f]) == 1) {
         ++n;
         average += s.sigsknot[fragments[f]].signal[sig];
      }
   }
   average /= n;
   // average enthaelt die durchschnittliche Signalkonzentration
   // n enthaelt die Zahl der Randobjektpunkte
   // Berechne den neuen steady state
   double newvalue = get_ss_receptor_ligand(K, r0, average);
   /* Bemerkung: Dies ist eine Naeherung, denn es wird nur der alte
    * Durchschnittswert nach old uebergeben. Dieser wird auf allen Fragmenten
    * als Referenz verwendet und der Verbrauch der Signale lokal entsprechend
    * berechnet. Das ist gut bei schneller Signal-Diffusion und kleinen Zellen.
    ### Bei vielen Fragmenten oder langsamer Signal-Diffusion sollte diese
    ###   Anpassung des Signals auf den lokalen alten Werten von bound_receptor
    ###   in jedem der Fragmente durchgefuehrt werden. Dazu ist ein array
    ###   bound_receptor[] (dimensioniert durch die Zahl der Randfragmente
    ###   eines Objekts) zu definieren und hier zu verwenden.
    */
   // Aktualisiere (siehe Bemerkung in frag_cell::bind_receptor_ligand()
   for (int f = 0; f < volume; f++) {
      if (l.object_border(fragments[f]) == 1) {
         if (s.sigsknot[fragments[f]].signal[sig] - ((newvalue - old) / n) < 0) {
            cout << "\n In frag_cell::bind_ss_receptor_ligand():\n"
                 << "Error: newvalue=" << newvalue << " old=" << old << " frag[" << f
                 << "].sig=" << s.sigsknot[fragments[f]].signal[sig]
                 << " _new=" << s.sigsknot[fragments[f]].signal_new[sig] << "\n";
            exit(1);
         }
         // l.knot[fragments[f]].signal[sig]-=((newvalue-old)/n);
         s.sigsknot[fragments[f]].signal_new[sig] -= ((newvalue - old) / n);
      }
   }
   return newvalue;
}
// ============================================================

double frag_cell::bind_receptor_ligand(const signal_molecule &sig,
                                       double &kplus,
                                       double &kminus,
                                       double &r0,
                                       double &rold,
                                       space &l,
                                       sigs &s) {
   double n = 0.0;
   double average = 0.0;
   for (int f = 0; f < volume; f++) {
      if (l.object_border(fragments[f]) == 1) {
         ++n;
         average += s.sigsknot[fragments[f]].signal[sig];
      }
   }
   average /= n;
   double delta_binding = get_receptor_ligand(kplus, kminus, average, r0, rold);
   // cout<<"in bind: average="<<average<<" delta_binding="<<delta_binding<<" n="<<n<<"\n";
   /* Bemerkung: Dies ist eine Naeherung, denn es wird nur der Durchschnittswert
    * verwendet. Dieser wird auf alle Fragmente gleichmaessig verteilt.
    * Das ist gut bei schneller Signal-Diffusion und kleinen Zellen.
    ### Bei vielen Fragmenten oder langsamer Signal-Diffusion sollte diese
    ###   Anpassung des Signals auf den lokalen alten Werten von bound_receptor
    ###   in jedem der Fragmente durchgefuehrt werden. Dazu ist ein array
    ###   bound_receptor[] (dimensioniert durch die Zahl der Randfragmente
    ###   eines Objekts) zu definieren und hier zu verwenden.
    ###   Aber zuerst muss die raeumliche Ordnung der Fragmente auch nach
    ###   Diffusion des Zellobjekts gewaehrleistet werden. ###
    */
   // Aktualisiere die Zahl der verbrauchten Signale
   s.signal_used[sig] += delta_binding;
   // Reduziere die Signale entsprechend und moeglichst homogen auf Randpunkten
   double rest = delta_binding;
   // cout<<"rest="<<rest<<" sig="<<l.knot[fragments[0]].signal[sig];
   // Redistributing the way the signal is used in the neighboring position (to avoid negative nbs)
   while (rest > 1.0e-08 || rest < -1.0e-08) {
      average = rest / n;
      for (int f = 0; f < volume; f++) {
         if (l.object_border(fragments[f]) == 1) {
            if (s.sigsknot[fragments[f]].signal_new[sig] >= average) {
               // l.knot[fragments[f]].signal[sig]-=average;
               s.sigsknot[fragments[f]].signal_new[sig] -= average;
               rest -= average;
            } else {
               rest -= s.sigsknot[fragments[f]].signal_new[sig];
               // l.knot[fragments[f]].signal[sig]=0.0;
               s.sigsknot[fragments[f]].signal_new[sig] = 0.0;
               --n;
            }
            /* In diese Abschnitt kann man die Veraenderung der Signalmenge
             * auch gleich in .signal[] umsetzen (nicht nur in .signal_new[]).
             * Dann wirkt sich sowohl Zuwachs als auch Abbau direkt auf das Gitter aus.
             * Wenn nicht, wird in lattice.signal_diffuse() der Abbau vor der Diffusion
             * umgesetzt, der Zuwachs auf die dortigen Unterzeitschritte verteilt. */
         }
      }
      if ((n == 0.0) && (rest > 1.0e-08)) {
         cout << "Error in bind_receptor_ligand!\n";
         exit(1);
      }
   }
   // cout<<"  "<<l.knot[fragments[0]].signal[sig]<<"\n";
   return delta_binding;
}
// ============================================================

short frag_cell::find_path(const long &source,
                           const long &target,
                           space &l,
                           const long &li,
                           const states &celltype) {
   double a;
   long b;
   return find_path(source, target, l, li, celltype, a, b);
}
short frag_cell::find_path(const long &source,
                           const long &target,
                           space &l,
                           const long &li,
                           const states &celltype,
                           double &max_dist,
                           long &max_dist_index) {
   /*
    * Looks for a path connecting "source" and "target" within the cell.
    * If found, return 1, else return 0.
    * In addition, measure the distance between the cell and the
    * direct connection line: return this as "max_dist" and remember
    * the most distant fragment index in "max_dist_index".
    */
   // Initialize return variables:
   max_dist = 0.;
   max_dist_index = -2;

   // get source and target vectors:
   long v_s[l.dim];
   l.get_koord(source, v_s);
   long v_t[l.dim];
   l.get_koord(target, v_t);
   // direct connection vector (normed)
   /* v_t_s may be used to call a new routine l.get_distance2line(double* double*)
    * double v_t_s[l.dim];
    * for (short i=0; i<l.dim; i++) v_t_s[i]=double(v_t[i]-v_s[i]);
    * double norm_t_s=l.get_2norm(v_t_s);
    * for (short i=0; i<l.dim; i++) v_t_s[i]/=norm_t_s;
    */
   // define an array to remember if the fragments have already been tried as path
   short triedpath[volume];
   for (int i = 0; i < volume; i++) {
      triedpath[i] = 0;
   }
   // actual point
   long a = source;
   triedpath[where_fragment(a)] = 1;
   long v_a[l.dim];
   for (short i = 0; i < l.dim; i++) {
      v_a[i] = v_s[i];
   }

   // Walk through the cell:
   /*
    * cout<<"frags=(";
    * for (short b=0; b<volume; b++) cout<<fragments[b]<<",";
    * cout<<")\n";
    */
   // cout<<"here: s="<<source<<" t="<<target<<" a="<<a<<"\n";

   while (l.Index(v_a) != target) {
      /*
       * cout<<"here: s="<<source<<" t="<<target<<" a="<<a
       *  <<" v_s=("<<v_s[0]<<","<<v_s[1]<<")"
       *  <<" v_t=("<<v_t[0]<<","<<v_t[1]<<")"
       *  <<" v_a=("<<v_a[0]<<","<<v_a[1]<<")"
       *  <<"\n";
       */
      long take = -1;
      long take1 = -1;
      // direction vector:
      double dir[l.dim];
      for (short i = 0; i < l.dim; i++) {
         dir[i] = double (v_t[i] - v_a[i]);
      }
      double norm_dir = l.get_2norm(dir);
      for (short i = 0; i < l.dim; i++) {
         dir[i] /= norm_dir;
      }
      // find next neighbour of v_a with maximal projection on "dir", save in a
      double tmpdist = l.pointnumber;
      double tmpdist1 = l.pointnumber;
      double project = -1. * l.pointnumber;
      double project1 = -1. * l.pointnumber;
      /*    double tmpdist=10.*l.prodim;
       * double tmpdist1=10.*l.prodim;
       * double project=-10.*l.prodim;
       * double project1=-10.*l.prodim;   */
      short stop = 0;
      for (short i = 0; i < l.dim2; i++) {
         if (stop == 0) {
            if (l.knot[a].near_n[i] == target) {
               // arrived!
               take = l.knot[a].near_n[i];
               tmpdist = 0.;
               stop = 1;
            } else {
               short fragstatus = triedpath[where_fragment(l.knot[a].near_n[i])];
               if ((l.self(l.knot[a].near_n[i], li, celltype) == 1) && (fragstatus < 2)) {
                  long v_n[l.dim];
                  l.get_koord(l.knot[a].near_n[i], v_n);
                  double tmp = l.get_distance2line(v_n, v_s, v_t);
                  for (short j = 0; j < l.dim; j++) {
                     v_n[j] -= v_a[j];
                  }
                  double tmp2 = l.get_scalarproduct(v_n, dir) / l.get_2norm(v_n);
                  if (fragstatus == 0) {
                     if ((tmp2 > project) || ((tmp2 == project) && (irandom(2) == 0))) {
                        tmpdist = tmp;
                        project = tmp2;
                        take = l.knot[a].near_n[i];
                     }
                  } else {
                     // fragstatus==1
                     if ((tmp2 > project1) || ((tmp2 == project1) && (irandom(2) == 0))) {
                        tmpdist1 = tmp;
                        project1 = tmp2;
                        take1 = l.knot[a].near_n[i];
                     }
                  }
               }
            }
         }     // end if (stop==0)
      }        // end for (i)
      if (take == -1) {
         if (take1 == -1) {
            // No additional path found, abort check, forbid move
            // cout<<" unable to find path -> abort!\n";
            return 0;
         } else {
            // case take==-1 && take1!=-1
            take = take1;
            take1 = -1;
            tmpdist = tmpdist1;
            // tmpdist1=10.*l.prodim;
            tmpdist1 = l.pointnumber;
            // i.e. fragments with status 1 are taken only if no fragments with status 0 are left!
         }
      }
      // use "take" one as new starting point
      a = take;
      /* remember that this fragment has already been used,
       * but allow to come back to this one if a neighbour exists
       * that opens a possible new path: */
      int fragind = where_fragment(a);
      if (triedpath[fragind] == 0) {
         short open_paths = 0;
         for (short j = 0; j < l.dim2; j++) {
            if ((l.self(l.knot[a].near_n[j], li,
                        celltype) == 1)
                && (triedpath[where_fragment(l.knot[a].near_n[j])] < 2)) {
               ++open_paths;
            }
         }
         if (open_paths > 1) {
            triedpath[fragind] = 1;
         } else {
            triedpath[fragind] = 2;
         }
      } else {
         triedpath[fragind] = 2;     // try every fragment two times at maximum
      }
      // get the corresponding vector
      l.get_koord(a, v_a);

      if (stop == 0) {
         // remember the distance to the direct connection of v_s and v_t
         /*
          * cout<<"v_a=("<<v_a[0]<<","<<v_a[1]<<") ";
          * cout<<"v_s=("<<v_s[0]<<","<<v_s[1]<<") ";
          * cout<<"v_a_s=("<<v_a_s[0]<<","<<v_a_s[1]<<") ";
          */
         if ((tmpdist > max_dist) || ((tmpdist == max_dist) && (irandom(2) == 0))) {
            max_dist = tmpdist;
            long v_a_s_mean[l.dim];
            for (short j = 0; j < l.dim; j++) {
               v_a_s_mean[j] = long ((3. * double (v_a[j]) + double (v_s[j])) / 4. + 0.5);
            }
            max_dist_index = l.Index(v_a_s_mean);
         }
      }
   }
   // cout<<"max_dist="<<max_dist<<"\n";
   return 1;
}
short frag_cell::narrow_neck(const long &source,
                             long &target,
                             const long &li,
                             const states &celltype,
                             space &l) {
   /*
    * Checks if a neck is went throuh by a fragment displacement from "source" to "target".
    * Criterion: Measure on a straight way over the lattice from "source" to "target"
    *           how far the nearest point of the cell is situated -> "minway".
    *           Compare this distance to the radius of the cell.
    *           if minway=0: no neck, i.e. return 0
    *           if 0<minway/radius<1: return -1 with probability
    *                                 [(radius-minway)/radius]^(dim-1) (note, that pi
    *                                 cancels in this expression for 3D)
    *                                 and index of neck point otherwise.
    *           if minway/radius>1: suppress this move and return the index of neck point
    *           if no path found: return -2
    */

   // remember max_distance to direct connection and corresponding point:
   double max_dist = 0.;
   long max_dist_index = -2;
   short path_exist = find_path(source, target, l, li, celltype, max_dist, max_dist_index);
   if (path_exist == 0) {
      // cout<<" unable to find path -> abort!\n";
      target = source;
      return -4;
   }
   /*
    * Now "max_dist" is the maximal distance to the direct connection line
    * between source and target on the optimal way through the cell.
    * "max_dist_index" is the corresponding index on the lattice.
    */

   // +++ OPTION ================================
   // Critical scale
   double q = 1.4;
   // end OPTION ================================

   // Get radius of the cell in units of the lattice constant on the basis of its "volume":
   // double radius=get_radius(l.dim);

   // cout<<"q="<<q<<" dist="<<max_dist<<" ==> ";
   if ((int (max_dist / 1.4) == 0) || (max_dist / radius < 0.25)) {
      // note this 1.4 is not q!
      // note that the first condition is relevant for smaller cells, the latter for larger ones.
      // cout<<"no problem!\n";
      // target unchanged!
      return 0;
   }

   /*
    * cout<<" d/q*r="<<max_dist/(q*radius)
    *  <<" f="<<pow((q*radius-max_dist)/(q*radius),double(l.dim-1))<<" ==> ";
    */

   if (max_dist / (q * radius) >= 1) {
      // cout<<"huge problem -> abort!\n";
      // target=source;
      // return -3;
      // cout<<"maxdist="<<max_dist<<" radius="<<radius<<" ";
      // cout<<"huge problem -> take another near "<<max_dist_index<<"!\n";
      target = max_dist_index;
      return -2;
   } else {
      // case 0<max_dist/radius<1
      if (drandom() < pow((q * radius - max_dist) / (q * radius), double (l.dim - 1))) {
         // cout<<"little problem -> take another near "<<max_dist_index<<"!\n";
         target = max_dist_index;
         return -1;
         // cout<<"little problem but o.k.!\n";
         // return 0;
      } else {
         // cout<<"maxdist="<<max_dist<<" radius="<<radius<<" ";
         // cout<<"medium problem -> abort!\n";
         target = source;
         return -3;
         // cout<<"medium problem -> take another near "<<max_dist_index<<"!\n";
         // target=max_dist_index;
         // return -2;
      }
   }
}
// ============================================================

double frag_cell::deformable_correction() {
   /* Calculates a reference value (observed in simulations)
    * for the fractions of used free neighbor points. These
    * are fitted to the values seen for K_1/2=1.0 and t_min=0.3.
    * The ratio of the actual fraction (alpha) and this reference
    * value is a suitable approximation to correct for deformation
    * parameters that differ from the above mentioned.
    * The used fit-curves are widely independent of the cells
    * diffusion constant.
    *
    * Used for "use_D_correction==1" only!
    */
   double reference_frac;
   if (volume > 3) {
      reference_frac = (0.44 + exp(-double (volume - 1) / 1.9));
   } else {
      reference_frac = (0.44 + 0.56 * exp(-double (volume - 1)));
   }
   return alpha_mean10 / reference_frac;
}
// ============================================================

double frag_cell::volume_D_factor() {
   // ! if changes are made here, corresponding changes are necessary in cellthis.C
   //   in the set_statics-routines !
   if (volume == 1) {
      return 1.0;
   } else {
      return 0.18 + (2.0 / volume);
   }
}
// ============================================================

double frag_cell::get_distances(char * which, double * ref_point, space &l, double * distance) {
   /*
    * Calculates the distance of all points with index i and "which[i]==1" to "ref_point".
    * The result is saved in "distance[i]" and the average distance is returned.
    * Note: the third potence of the distances is saved and averaged!!!
    */
   double average = 0.;
   int n = 0;
   for (int i = 0; i < volume; i++) {
      if (which[i] == 1) {
         distance[i] = l.get_2norm(ref_point, l.knot[fragments[i]].x);
         // distance[i]*=distance[i];
         // distance[i]*=distance[i];
         average += distance[i];
         ++n;
      }
   }
   if (n == 0) {
      average = -1;
   } else {
      average /= double (n);
   }
   return average;
}
// ============================================================

void frag_cell::set_forbidden_counts(const long &j, double * target, double * old, space &l) {
   /* OPTION: MOVE_ANALYSIS
    * ++n_move_forbidden;
    * double frac_forbid=0.;
    * for (short loop=0; loop<l.dim; ++loop)
    * frac_forbid+=((double(l.knot[j].x[loop])-old[loop])
    *(target[loop]-old[loop]));
    * if (frac_forbid<=0.) ++n_move_forbidden_back;
    */
}
// ============================================================

void frag_cell::fragdiffuse(states celltype,
                            const double &tolerance_min,
                            const double &tolerance_steepness,
                            const long &li,
                            double p,
                            const short &nocutoff,
                            const short &ordered,
                            double * target_point,
                            space &l) {
   /* Note: "celltype" and "li" define the cell under consideration.
    *       "tolerance_min" and "_steepness" define the deformability of the cell.
    *       "p" the probability of move per fragment;
    *          p=-1 infers all fragments to be moved (if possible);
    *       "nocutoff" decides if the movement of fragments is cut off when the barycenter is
    *          moved by p lattice constants (=0 yes, =1 no).
    *       "target_point" is a virtual barycenter to which the fragments movement is related;
    *          if this is the barycenter fragments will basically rearrange,
    *          if it is some other point the whole cell will move towards this point;
    *       "ordered" (0 or 1);
    *          if 0 the border points are moved in random order;
    *          if 1 the border points are ordered according to their distance to the
    *          "target_point", and the move is weighted according to this distance.
    */
   long j;
   if (l.cellknot[index].cell != celltype) {
      cout << "Fehler in fragdiffuse: cell of type " << celltype << " not at index " << index
           << "! ";
      exit(1);
   }

   /* #### Find an algorithm that does not run through all fragments, in particular
    * in case of low probability of fragment movement. This is a long procedure that
    * is called often for reshaping with marginal movement. --> CPU load.
    */

   // cout<<"reshape="<<ordered;

   // save the previous barycenter:
   // ### delete when n_move_xxx_back variables are deleted
   // need it now for cutoff of movement:
   double oldbary[l.dim];
   for (short loop = 0; loop < l.dim; loop++) {
      oldbary[loop] = barycenter[loop];
   }
   // ### isn't that already done in the calling routine fragmove(...)?

   // Declare counter for the number of neighbors of each fragment belonging to the same cell:
   short selfnn;
   // Declare an index var for the target point:
   long newind;
   // Declare and define an array of switches that itemizes those fragments that may be moved:
   char maymoveit[volume];
   // reset the value for borderpoints:
   borderpoints = 0;
   // set this flag to 1 for all border points:
   for (int f = 0; f < volume; f++) {
      // cerr<<"in for (int f=0; ...) { ...\n";
      // Calculate the number of neighbors belonging to the same cell:
      selfnn = l.get_n_self(fragments[f], li, celltype, -1);
      // cerr<<"selfnn="<<selfnn<<"\n";
      // Move only if the fragment is at the border of the cell:
      if (selfnn < l.dim2) {
         // cerr<<"start if (selfnn<l.dim2) { ...\n";
         maymoveit[f] = 1;
         // add the point to borderpoints
         ++borderpoints;

         // ====================================
         // =========== adhesion ===============
         // ====================================
         // check for restriction of movement by adhesion to neighbours
         // Base the restriction of movement on the time of immobility "t_immobile[f]"
         // The adhesive connecion is considered to be established after "adhesion_time"
         /*cout<<"f="<<f<<": t_imm="<<t_immobile[f]<<", ad_time="<<adhesion_time
          * <<", ad="<<adhesive<<"\n";*/
         // if (find_contact(ext,fragments[f],l)==0    // i.e. contact found
         if ((l.no_external(fragments[f]) == 0)    // i.e. contact to external found
             // #### contact to other cells also important !!! Do it soon!!!
             &&
             (t_immobile[f] > adhesion_time)    // adhesive contact assumed to be established
             &&
             (adhesive > 0.)) {
            // the cell is able of adhesion
            /* There are two situations:
             * 1) The fragment is the only with adhesion to some other cell,
             * 2) At least 2 neighbouring fragments of the same cell are adhesive as well.
             * In case 2) the adhesive force is assumed maximal, i.e. no movement.
             * In case 1) a probabilistic decision according to the "adhesive" parameter is done:
             */
            short n_adhesives = 0;
            for (short i = 0; i < l.dim2; i++) {
               long j = l.knot[fragments[f]].near_n[i];
               if (l.self(j, li, celltype) == 1) {
                  // self-neighbour found
                  if (t_immobile[where_fragment(j)] > adhesion_time) {
                     ++n_adhesives;
                  }
               }
            }
            if (n_adhesives > 1) {
               flag_no_correction = 1.;
               maymoveit[f] = 0;
            } else if (drandom() < adhesive) {
               maymoveit[f] = 0;
            }
         }
         // Note, that for contact to other cells, two criteria have to be met:
         // 1) the correct adhesion receptor and ligands are expressed
         // 2) both fragments of the two cells have to be immobile for a while
         // Therefore, introduce receptor and ligand and calculate the product
         // of both immobility times. ###
         // ====================================
         // end of adhesion
         // =================================
         // cerr<<"end if (selfnn< ... {\n";
      } else {
         maymoveit[f] = 0;
      }
   }

   // cerr<<", borderpoints="<<borderpoints<<"\n";
   // long ntryold=n_try_move;

   /* now the fragments are labeled with
    * 0: no move
    * 1: move
    * 2: no move because of adhesion
    * The case 2 is transformed to 0 and the fraction of 2s saved in flag_no_correction:
    */
   /*
    * int n_1s=0,n_2s=0;
    * for (int f=0; f<volume; f++) {
    * if (maymoveit[f]==1) ++n_1s;
    * else if (maymoveit[f]==2) { ++n_1s; ++n_2s; maymoveit[f]=0; }
    * }
    * // fraction of intended moves prevented by adhesion
    * flag_no_correction=double(n_2s)/double(n_1s);
    */

   // =======================================================
   // ======= case of weighted fragments by distance ========
   // =======================================================
   if (ordered == 1) {
      // define distance dependent moving probabilities:
      // cerr<<"start if (ordered==1) { ...\n";
      // define a list of distances
      double distance_list[volume];
      double average_distance = get_distances(maymoveit, target_point, l, distance_list);

      /*
       * cout<<"average_dist(target=("
       *  <<target_point[0]<<","<<target_point[1]<<"))="
       *  <<average_distance<<": ";
       * for (int f=0; f<volume; f++)
       * if (maymoveit[f]==1) {
       *  cout<<"  ("<<l.knot[fragments[f]].x[0]
       *      <<","<<l.knot[fragments[f]].x[1]
       *      <<")->"<<distance_list[f]/average_distance
       *      <<";\n";
       * }
       */

      // the weight is defined as the distance divided by the average distance:
      for (int f = 0; f < volume; f++) {
         if (maymoveit[f] == 1) {
            if ((drandom() < p * distance_list[f] / average_distance) || (p == -1)) {
               /* OPTION: MOVE_ANALYSIS
                * ++n_try_move; // remember the number of tries of fragment movement
                */
            } else {
               maymoveit[f] = 0;
            }
         }
      }
   } else {
      // ie if (ordered==0)
      // cerr<<"start if (ordered==0) { ...\n";
      for (int f = 0; f < volume; f++) {
         /* Note: The following random decision is done for each fragment separately, i.e.
          *       this routine describes a diffusion of fragments within the cell and not
          *       the diffusion of the cell itself!
          *       Alternatively, the probabilistic decision may be done for the whole cell
          *       movement. Afterwards all fragments are moved trying to reach the desired
          *       movement of the cell on the fragment level.
          *       Indeed, the calling routine decides if a whole cell movement is done
          *       (probabilistic decision). If yes, it calls this routine with p=-1,
          *       if no, this routine is called with p>0 and fragment-based diffusion
          *       is performed.
          */
         if (maymoveit[f] == 1) {
            if ((p == -1) || (drandom() < p)) {
               /* OPTION: MOVE_ANALYSIS
                * ++n_try_move; // remember the number of tries of fragment movement
                */
            } else {
               maymoveit[f] = 0;
            }
         }
      }
   }

   // cout<<", old="<<ntryold<<", volume="<<volume<<", move="<<n_try_move-ntryold<<"\n";
   // cout<<", volume="<<volume;

   // Go through all fragments of the cell:
   // cerr<<"vor random2_sequence(...) ...\n";
   long ff[volume];
   random2_sequence(ff, volume);
   /* This is done by randomizing the sequence of fragments, being not clearly
    * necessary as the fragments are not in any specific order and growth of the
    * cell is isotropic, such that the order is not correlated to the position
    * in space. Nevertheless, a slight anisotropic sequence of fragments in the
    * list would break the isotropy of the diffusion movement of the cell and
    * is even enhanced by the deletion of fragments saved at the end of the
    * list at the end of this routine. Therefore, a randomized version of this
    * routine is currently tested. */
   // variable for control of performed move
   short stop_moving = 0;
   int nf = 0;
   // cerr<<"vor while (nf<volume ...) { ...\n";
   while (nf < volume && stop_moving == 0) {
      int f = ff[nf];
      // cerr<<"fragment="<<f<<"... \n";
      // start only if the fragment is itemized:
      if (maymoveit[f] == 1) {
         // reset flag
         //      short flag_forbidden=0;
         // get lattice index of the fragment
         j = fragments[f];
         // cerr<<"Try to move fragment "<<f<<"at index "<<j<<":\n";
         // Make sure that this border point is still a border point
         selfnn = l.get_n_self(j, li, celltype, -1);
         /* Move only if the fragment is still at the border of the cell:
          * Note that has to be checked again, as after the movement of
          * other fragments, fragments initially being border points
          * are not necessarily still. */
         if (selfnn < l.dim2) {
            /* look for free neighbors of the whole cell, i.e. it is allowed that the
             * cell deforms. Delete a fragment on one side and put it to the other is
             * explicitly allowed and corresponds to a shift of cytosol. Note, that
             * the way of switching a fragment over long distances is only correct if
             * the fragments of the same cell are all equivalent to each other. The
             * concept should be revisited if e.g. local adhesion molecule distributions
             * or receptor densities are considered within the same cell! ##
             */
            // Declare list of free neighbors
            long n_list[(l.dim2 - 1) * volume + 2];
            // Store all free neighbors of this cell in n_list:
            int used_nn = get_free_nn(l,
                                      li,
                                      n_list,
                                      tolerance_min,
                                      tolerance_steepness,
                                      target_point,
                                      j,
                                      j);

            // deformation parameter:
            // count the total number of times that this cell has run get_free_nn();
            ++n_moves;
            // change the average fraction using alpha just set in get_free_nn(...)
            alpha_mean10 = (alpha_mean10 * 9.0 + alpha) / 10.;
            alpha_mean = (alpha_mean * double (n_moves - 1) + alpha) / double (n_moves);
            // frac_average (average over all cells) is calculated in cellman.C

            // cerr<<"used_nn="<<used_nn<<"\n";

            /* If no free neighbour was defined and the cell is a one fragment cell
             * the cell-state of contact inhibition is set for use in the calling
             * routine. This allows to define the requirements for a cell-cell exchange.
             */
            if ((volume == 1) && (used_nn == 0)) {
               contact_inhibited = true;
            }

            /* Randomly choose one of those free neighbors in n_list with most self-neighbors
             * and check if the fragment can be moved without separating the cell into
             * two parts:
             */
            if (used_nn > 0) {
               // Choose randomly the new place for fragment f at lattice point j
               newind = n_list[irandom(used_nn)];
               if (newind != j) {
                  /*
                   * show_fragments();
                   * cout<<"n_list[]=";
                   * for (int bb=0; bb<used_nn; ++bb) cout<<n_list[bb]<<",";
                   * cout<<"\n";
                   */

                  /* Before moving, check wether the movement infers a shift of the
                   * cytosol through a narrow neck. Try to avoid this check in clear
                   * cases in order to keep CPU load small.
                   */
                  short neck = narrow_neck(j, newind, li, celltype, l);
                  // cerr<<"narrow_neck --> "<<neck<<"\n";
                  // neck==0 if no narrow neck or move nevertheless, newind unchanged;
                  // neck==-1 move but newind has changed! (neck present, do not move beyond)
                  // neck==-2 move but newind has changed!
                  //          (neck present, do not move beyond, ignore for correction)
                  // neck==-3 newind has changed to j --> do not move, ignore for correction
                  // neck==-4 no path found --> abort movement, ignore for correction
                  if (neck < -1) {
                     // +++ OPTION
                     // suppress any change of the p-diffusion-correction-factor for this move:
                     // flag_no_correction=1.;
                     flag_no_correction += 1. / double (borderpoints);
                     /* This is to be done in order to avoid that the correction factor
                      * counteracts a reduced movement provoked by obstacles.
                      */
                     // end OPTION
                     //	      flag_forbidden=1;
                     /* OPTION: MOVE_ANALYSIS
                      * ++n_move_forbidden;
                      */
                  }
                  if ((neck == -1) || (neck == -2)) {
                     // newind is now the position of an optimal target point
                     // ==> find the nearast to that point:
                     // define vector for neck
                     double v_neck[l.dim];
                     l.get_koord(newind, v_neck);
                     // Declare new list of free neighbors
                     long neck_list[(l.dim2 - 1) * volume + 2];
                     // re-run get_free_nn(...v_neck...) in the non-deforming limit
                     int neck_used
                        = get_free_nn(l, li, neck_list, 1.e-07, 0.01, v_neck, j, j);
                     // take one of the possible target points
                     // cout<<" neck_used="<<neck_used<<" j="<<j;
                     if (neck_used > 0) {
                        newind = neck_list[irandom(neck_used)];
                     } else {
                        newind = j;          // i.e. no movement
                                             //		flag_forbidden=1;
                     }
                     // cout<<" newind="<<newind<<"\n";
                     /* Note, that the target point is not changed here, as the evaluation
                      * of the performed move (in relation to the intended one) has to be
                      * based on the initially fixed target point.
                      */
                  }
               }

               if (newind != j) {
                  // cerr<<"start if (newind!=j) { ...\n";
                  // check a possible decomposition of the cell by moving "j":
                  long decomp = l.decompose(j, li, celltype);
                  // cout<<"j="<<j<<" newind="<<newind<<" decomp="<<decomp<<"\n";
                  if (decomp != j) {
                     /* Get back from l.decompose(...):
                      * "-1": decomposition of the cell --> forbid move!
                      * "j": Either no decomposition of the cell
                      *      (or decomposition but closed ring of width one lattice constant)
                      *      --> seemingly a decomposition, allow the move as it is done!
                      * "other than j": Fragment "j" is part of a tail
                      *                 --> check if newindex is not the end of the tail
                      *                 --> and then move decomp instead!
                      */
                     if (decomp == -1) {
                        /* Possible decompostion. This information is based on the check
                         * of the neighbours of "j": No direct connection of these exists
                         * but via the point "j". This does not exclude a more complicated
                         * path to exist. In particular if "j" is part of a ring structure,
                         * the way to reach one neighbour from another one may be longer
                         * than one step. It is of great interest that such rings if they
                         * have not been prevented from being built are at least destroyed.
                         * To forbid "j" to move in just this situation would invalidate
                         * the cell movement.
                         * So far the problem. Now the solution:
                         */
                        // temporarily delete the fragment:
                        l.cellknot[j].cell = nocell;
                        int left = -1;
                        int right = -1;
                        short foundpath = 1;
                        short loop = 0;
                        while (loop < l.dim2 && foundpath == 1) {
                           if (l.self(l.knot[j].near_n[loop], li, celltype) == 1) {
                              if (left == -1) {
                                 left = l.knot[j].near_n[loop];
                              } else {
                                 right = l.knot[j].near_n[loop];
                                 foundpath = find_path(left, right, l, li, celltype);
                                 // cout<<"path from "<<left<<" to "<<right<<"="<<foundpath<<"; ";
                              }
                           }
                           ++loop;
                        }
                        // cout<<"j="<<j<<" foundpath="<<foundpath<<"\n";
                        // foundpath==1 if all neighbours connected, ==0 otherwise.
                        // restore the fragment on the cell-lattice:
                        l.cellknot[j].cell = celltype;

                        // Suppress the move only if foundpath==0:
                        if (foundpath == 0) {
                           // do not move:
                           //		  flag_forbidden=1; // exclude from "move on itself" counting
                           // cout<<"forbidden to move "<<j<<". decomp="<<decomp<<"\n";
                           /* OPTION: MOVE_ANALYSIS
                            * set_forbidden_counts(j,target_point,oldbary,l);
                            */
                           // Restore newind to j. This prevents the move to be executed.
                           newind = j;
                        }
                        // else move from "j" to "newind" as initially intended
                     }
                     if (decomp != -1) {
                        // case decomp also != j
                        // check that the movement target is not the end of the tail!!!
                        /* Do this by counting the self-neighbours of the target point
                         * and excluding the point decomp of this count. This number
                         * has to be larger than zero!
                         */
                        if (l.get_n_self(newind, li, celltype, decomp) > 0) {
                           // don't try to move fragment f with index j again:
                           maymoveit[f] = 0;
                           // move "decomp" instead of "j" to "newind"
                           j = decomp;
                           /* But the fragment to be changed on the fragment list
                            * has a index different from "f" now!
                            * Look for the fragment-index with value decomp! */
                           f = where_fragment(decomp);
                        } else {
                           // do not move:
                           //		  flag_forbidden=1; // exclude from "move on itself" counting
                           // cout<<"forbidden to move "<<j<<". decomp="<<decomp<<"\n";
                           /* OPTION: MOVE_ANALYSIS
                            * set_forbidden_counts(j,target_point,oldbary,l);
                            */
                           // Restore newind to j. This prevents the move to be executed.
                           newind = j;
                        }
                     }
                  }        // end if (decomp!=j)
                           // else move from j to newind as initially intended
                           // cerr<<"end if (newind!=j) {\n";
               }

               if (newind != j) {
                  // cerr<<"start second if (newind!=j) {...\n";
                  // Check if the movement will enclose some other object
                  // Check one-fragment objects only:
                  if (enclosed_object(j, newind, li, celltype, l) == 1) {
                     // do not move:
                     //	      flag_forbidden=1; // exclude from "move on itself" counting
                     /* OPTION: MOVE_ANALYSIS
                      * set_forbidden_counts(j,target_point,oldbary,l);
                      */
                     // Restore newind to j. This prevents the move to be executed.
                     newind = j;
                  }
                  // cerr<<"end second if (newind!=j) {...\n";
               }

               if (newind != j) {
                  // =====================================
                  // Uff! Do the move from j to newind !!!
                  // =====================================
                  // cerr<<" move from "<<j<<" to -> "<<newind<<":\n";
                  /*
                   * <<"n=("<<n_list[0]<<","<<n_list[1]<<","
                   * <<n_list[2]<<","<<n_list[3]<<","<<n_list[4]<<"...";
                   */
                  /* Save the fragment on the new lattice point and delete from the old*/
                  l.set_knot(newind, celltype, li);
                  // l.knot[newind].cell=celltype;
                  // l.knot[newind].listi=li;
                  l.clear_knot(j);
                  // l.knot[j].cell=empty;
                  // l.knot[j].listi=-1;
                  fragments[f] = newind;
                  reset_clock(f);
                  // Correct for the shifted barycenter
                  index = get_barycenter(l);
                  /* This is done after each moved fragment, which is not without alternative:
                   * Here, the movement of all fragments is thought to be sequential, and each
                   * fragment knows about the state of the cell after previous movements of all
                   * other fragments. Alternatively, all fragments may be moved according to
                   * the fixed initial state of the cell.
                   */
                  // ### Eventually better to calculate relative changes ...

                  // Check how far the cell went:
                  if (nocutoff == 0) {
                     double go = l.get_2norm(barycenter, oldbary);
                     // #####double overshot=get_radius(l.dim)/(2.*volume);
                     double overshot = radius / (2. * volume);
                     // cout<<"nf="<<nf<<"; p="<<p<<"; p-bla="<<p-overshot<<"; go="<<go<<"\n";
                     if ((go > 1. - overshot) || ((go > p - overshot) && (p >= 0.))) {
                        // cout<<"STOP!\n";
                        stop_moving = 1;
                     }
                  }

                  // OPTION: MOVE_ANALYSIS
                  // count the number of successful movements
                  // ++n_move_done;
                  // end OPTION

                  /* Now eliminate further fragments from the maymoveit list if a
                   *  previous border point has been covered by the former movement:
                   *  The old and new index are j and newind!
                   *  Find the neighbour of newind that is most in direction of j: */
                  // define vectors:
                  double vj[l.dim];
                  for (short dd = 0; dd < l.dim; dd++) {
                     vj[dd] = double (l.knot[j].x[dd]) - double (l.knot[newind].x[dd]);
                  }
                  // Get candidates for elimination from movement list:
                  long n_newind[l.dim];
                  short howmany = l.get_nn_directed2(vj, newind, n_newind);
                  /*
                   * cout<<" howmany="<<howmany<<" n_newind=(";
                   * for (short bla=0; bla<howmany; ++bla) cout<<n_newind[bla]<<",";
                   * cout<<")  ";
                   */
                  // "n_newind" contains "howmany" lattice indices
                  short alphaa = 0;
                  short done = 0;
                  // The choice of the fragment that is deactivated
                  // shall happen in a random sequence!!! done:
                  while (howmany > 0 && done == 0) {
                     if (howmany > 1) {
                        alphaa = irandom(howmany);
                     } else {
                        alphaa = 0;
                     }
                     // find the corresponding fragment in the cell
                     int frag2j = where_fragment(n_newind[alphaa]);
                     // if found, remove from move list!
                     if (frag2j != -1) {
                        /* OPTION: MOVE_ANALYSIS
                         * if (maymoveit[frag2j]==1) ++n_move_removed;
                         */
                        maymoveit[frag2j] = 0;
                        done = 1;
                     } else {
                        // remove this one from the list "n_newind"
                        // copy the last to the actual
                        n_newind[alphaa] = n_newind[howmany - 1];
                        // reduce the number of items
                        --howmany;
                     }
                  }
                  /*
                   * for (short bla=0; bla<volume; bla++) cout<<int(maymoveit[bla])<<",";
                   * cout<<"\n";
                   * if (xxxx=='a') cin>>xxxx;
                   */

                  /* Is it possible that this process violates isotropy? We walk
                   * through the fragments always in the same order. Therefore,
                   * always the fragments with higher index in the fragment list
                   * will be deleted from the moving list. However, the
                   * order is not necessarily correlated with the order in space.
                   * This is to be checked!
                   * No problem, as now the fragments are walked through in a random order.
                   */
               }            // end if newind!=j ...
               else {
                  // otherwise the move has been suppressed
                  /* OPTION: MOVE_ANALYSIS
                   * if (flag_forbidden==0) { // case of movement on itself
                   * ++n_move_self;
                   * double frac_self=0.;
                   * for (short loop=0; loop<l.dim; ++loop)
                   *  frac_self+=((double(l.knot[j].x[loop])-oldbary[loop])
                   *(target_point[loop]-oldbary[loop]));
                   * if (frac_self<=0.) ++n_move_self_back;
                   * }
                   */
               }
               // cout<<"end.\n";
            }      // end if used_nn ...
         }         // end if selfnn<l.dim2
         maymoveit[f] = 0;
      }    // end if maymoveit[f]==1
      ++nf;
   }   // end while (...)
       /*
        * cout<<"Check of cell li="<<li<<": ";
        * if (volume>1) {
        *  for (int k=0; k<volume; ++k) {
        *    if (l.get_n_self(fragments[k],li,celltype,-1)==0) {
        *      cout<<"Error: one fragment of cell "<<li<<" has lost contact!\n";
        *      cout<<"Fragment-index = "<<k<<" lattice-index = "<<fragments[k]<<"\n";
        *      exit(1);
        *    }
        *    if (l.knot[fragments[k]].listi!=li) {
        *      cout<<"Error: lattice point of a fragment points to wrong cell!\n";
        *      cout<<"Fragment-index = "<<k<<" lattice-index = "<<fragments[k]
        *          <<" points to "<<l.knot[fragments[k]].listi<<"\n";
        *      exit(1);
        *    }
        *  }
        * }
        * cout<<"done!\n";
        */
       // cout<<"done!\n";
       // fragment_consistency();
}
// ============================================================

// movement of fragments (sum of diffusion and other forces:
double frag_cell::fragmove(states celltype,
                           const long &li,
                           const double &tolerance_min,
                           const double &tolerance_steepness,
                           const double &elongfrac,
                           const double &eta_max,  // max smoothmove
                           const double &p_tension,
                           const double &K_elongation,
                           space &l) {
   // save the previous barycenter:
   double p = p_move;
   double oldbary[l.dim];
   for (short loop = 0; loop < l.dim; loop++) {
      oldbary[loop] = barycenter[loop];
   }

   if (use_D_correction == 1) {
      // Rescale the diffusion probability p for the diffusion of fragments:
      // This is possible by linear scaling as D enters the probability linearly (see cellthis.C).
      p /= volume_D_factor();
      // cout<<"p_1="<<p<<"; alpha_="<<alpha_mean10<<"; v="<<volume
      // <<"; _corr="<<deformable_correction();
      p /= deformable_correction();
      // cout<<"; p_2="<<p<<"\n";
      // it is ensured in cellthis.C that p remains smaller than 1 !!!

      fragdiffuse(celltype, tolerance_min, tolerance_steepness, li, p, 1, 0, barycenter, l);
      // The barycenter is already corrected in fragdiffuse. Nothing to do.
   } else {
      // use_D_correction==0, i.e. use forces and barycenter shift

      // define a target-barycenter
      double wantedbary[l.dim];
      for (short loop = 0; loop < l.dim; loop++) {
         wantedbary[loop] = barycenter[loop];
      }
      // calculate smoothfactor:
      double smoothfactor = eta_max / double (borderpoints);
      if (smoothfactor < 1.) {
         smoothfactor = 1.;
      } else if (smoothfactor > 5 * borderpoints) {
         smoothfactor = double (5 * borderpoints);
      }
      // +++ OPTION
      // Choose the parameter value smoothmove to be interpreted as a constant value:
      // smoothfactor=eta_max;
      // Delete the above line if the parameter value is interpreted as maximum value!
      // end OPTION
      // cout<<"borderpoints="<<borderpoints<<"  smoothfactor="<<smoothfactor<<"\n";

      if (drandom() < smoothfactor * p / v_state) {
         // Shift of the barycenter:
         // #####double elongscale=elongfrac*get_radius(l.dim)*elongation;
         double elongscale = elongfrac * radius * elongation;
         if (elongscale < 1.) {
            elongscale = 1.;
         }
         for (short loop = 0; loop < l.dim; loop++) {
            // Shift the barycenter of the cell by "polarity"
            wantedbary[loop] += (elongscale * polarity[loop]);
         }

         // cerr<<"bary=("<<barycenter[0]<<","<<barycenter[1]<<"); wanted=("
         //  <<wantedbary[0]<<","<<wantedbary[1]<<")\n";

         // move the fragments:
         // version with p-correction factor
         /* Note the deformation parameter are chosen larger than in the non-deforming limit!
          * This is necessary in order to avoid lattice effects that emerge if only the point
          * nearest to the virtual barycenter is considered as possible target. In that case
          * and if the polarity vector opens a small angle to the lattice symmetry axis,
          * the resulting movement becomes parallel to the lattice axis.
          * The values 0.2 and 0.2 are chosen such that the angle is just respected,
          * i.e. as non-deforming as possible. Roughly, the values remain valid for
          * different spatial resolutions.
          */
         // cerr<<"volume="<<volume<<", vor fragdiffuse ...\n";
         fragdiffuse(celltype,
                     tolerance_min,
                     tolerance_steepness,
                     li,
                     (2. - performed2aimed_move_now) / smoothfactor,
                     short (use_dynamic_correction),
                     0,
                     wantedbary,
                     l);
         // cerr<<"p="<<(2.-performed2aimed_move_now)/smoothfactor
         //  <<"; p*volume="<<(2.-performed2aimed_move_now)*volume/smoothfactor<<"\n";

         // Note, that for smoothfactor==1 no accelaration of movement is achieved even for p2a==0!

         // if (use_dynamic_correction==1) {
         // =================================================
         // ============= correction factor =================
         // =================================================
         if (flag_no_correction < 1.) {
            /* Let's compare the performed move compared to the expected one:
             * (wantedbary-oldbary)/smoothfactor is the wanted displacement,
             * barycenter-oldbary is the performed one. Plot the ratio of both:
             * (barycenter-oldbary)/((wantedbary-oldbary)/smoothfactor =
             * (barycenter-oldbary)*smoothfactor
             * as (wantedbary-oldbary) is always 1.
             */
            // +++ OPTION: use the projection instead? No.
            double go = l.get_2norm(barycenter, oldbary);
            /*
             * double barymove[l.dim];
             * for (short a=0; a<l.dim; ++a) barymove[a]=barycenter[a]-oldbary[a];
             * go=l.get_scalarproduct(barymove,polarity);
             */
            // cout<<"go="<<go;
            // add smoothfactor: then go is the distance walked if all would have been tried
            go *= smoothfactor;
            // flag_no_correction is the fraction of moves that have been suppressed but
            // that have not to be treated as intended moves, i.e. not to be included in
            // the correction factor:
            go /= (1. - flag_no_correction);

            // OPTION +++
            // The following is for retarded corrections:
            // The factor "stable" makes the correction smooth
            const double stable = 5.;      // counted in units of one dx displacement
            // The relation of correction with compensation and without is determined here:
            const double compensate = 1.;
            // 0 means no compensation
            // 1 means half compensation half converging to movement 1 dx
            // higher values pronounce the fraction of compensation
            // end OPTION

            // Calculated new correction factor for fragment movement
            /*
             * performed2aimed_move_now=((stable*smoothfactor-1.)*performed2aimed_move_now
             * +go/(2.-performed2aimed_move_now))/(stable*smoothfactor);
             * Note, that this version is metastable! In this way an equilibrium between
             * a "wrong" walking behaviour and the correction factor is found. Instead,
             * the value of "go" has to tend to "1" exactly and always. Wrong walks have
             * also to be equilibrated. Both is ensured by the following variant
             * based on the average of a compensation and a long term target of "1":
             */
            performed2aimed_move_now = ((stable * smoothfactor - 1.) * performed2aimed_move_now
                                        + go * (2. - performed2aimed_move_now + compensate)
                                        / ((compensate + 1.) * (2. - performed2aimed_move_now)))
                                       / (stable * smoothfactor);
            /*
             * cout<<" go/="<<go<<" smooth="<<smoothfactor
             * <<" borderpoints="<<borderpoints
             * <<" flag_no_correction="<<flag_no_correction
             * <<" p2a="<<performed2aimed_move_now<<"\n";
             */
            // For statistics only: Calculate the overall average of this number for ana_file
            ++n_directed_moves;
            performed2aimed_move = ((n_directed_moves - 1) * performed2aimed_move + go)
                                   / n_directed_moves;
            /*// output:
             * cout<<" n="<<n_directed_moves
             * <<" go="<<go
             * <<" p2anow="<<performed2aimed_move_now
             * <<" 2-p2anow="<<2.-performed2aimed_move_now
             * <<" p2a="<<performed2aimed_move
             * <<"\n";
             */
            // end of compare performed to aimed move.
         }
         // =================================================
         // ========= end correction factor =================
         // =================================================
         /*
          * }
          * else {
          * performed2aimed_move=1.;
          * }
          */

         // reset flag_no_correction:
         flag_no_correction = 0.;
         /*
          * //cout<<"alpha="<<alpha<<"; ";
          * cout<<"old=("<<oldbary[0]<<","<<oldbary[1]<<"); "
          *  <<"wanted=("<<wantedbary[0]<<","<<wantedbary[1]
          *  <<"); dist="<<l.get_2norm(oldbary,wantedbary)<<"; "
          *  <<"real=("<<barycenter[0]<<","<<barycenter[1]
          *  <<"); dist="<<l.get_2norm(oldbary,barycenter)<<"\n";
          */
      }

      if (volume > 1) {
         // Do thermic rearrangement of cell fragments only:
         /* 1 One may consider the introduction of a separate fragment diffusion constant, done
          * 2 One may also consider to run this deformation routine
          *   independently of any shift before, done
          * 3 And it seems more appropriate to restrict the fragment movement on movements
          *   to next neighbours only in this case (no cytosol shift) ???
          * Note, that for volume=1 the shift of the barycenter before already
          *       fully accounts for the cell movement. A rearrangement of the fragments
          *       is senseless.
          */

         /* Get new wanted barycenter which is thought to remain unchanged! Therefore,
          * the barycenter is not used itself (which would be dynamical within fragdiffuse)
          * but a fixed point "wantedbary" that is the actual barycenter.
          */
         for (short loop = 0; loop < l.dim; loop++) {
            wantedbary[loop] = barycenter[loop];
         }
         // get actual radius and elongation (basis for reshaping force)
         // #####get_radius(l.dim);
         get_elongation(li, celltype, l);

         /* OPTION: MOVE_ANALYSIS
          * // don't change n_move_done and n_try_move here!
          * long done_old=n_move_done;
          * long try_old=n_try_move;
          * long forbid_old=n_move_forbidden;
          * long forbid_back_old=n_move_forbidden_back;
          * long removed_old=n_move_removed;
          * long self_old=n_move_self;
          * long self_back_old=n_move_self_back;
          */

         // move the fragments:
         fragdiffuse(celltype, 1.e-07, 0.01, li, get_reshaping_force(p_tension,
                                                                     K_elongation), 1, 1, wantedbary,
                     l);

         /* OPTION: MOVE_ANALYSIS
          * // recover statistic variables:
          * n_move_done=done_old;
          * n_try_move=try_old;
          * n_move_forbidden=forbid_old;
          * n_move_forbidden_back=forbid_back_old;
          * n_move_removed=removed_old;
          * n_move_self=self_old;
          * n_move_self_back=self_back_old;
          */
      }

      // Correct for the shifted barycenter at the end of all movements:
      // index=get_barycenter(celltype,li,l);
      // not necessary as done in fragdiffuse(..);
   }

   // Return the distance that the whole cell has moved:
   double moved = l.get_2norm(oldbary, barycenter);
   add_moves += moved;
   if (moved == 0.) {
      ++n_immobile;
   } else if (trackit) {
      writethis2track = movement;
   }
   // if (changed_polarity) cout<<"pol-change="<<changed_polarity<<"\n";
   // if (moved!=0.) cout<<"; moved="<<moved<<"\n";
   return moved;
}
// ============================================================

short frag_cell::do_grow(states celltype,
                         const long &li,
                         const int &target_vol,
                         const double &p_grow,
                         const double &p_shrink,
                         space &l) {
   // cout<<"grow_CB ...";
   // check_listi();
   short err = 1;
   long newind;
   if ((volume < target_vol) && (drandom() < p_grow) && (volume < FRAGMENT_STEP)) {
      err = 0;
      /* look for a free neighbor point on the lattice;
       * we have to go through all fragment points of the actual cell;
       * save all free neighbors of the whole cell in a list */
      long n_list[(l.dim2 - 1) * volume + 2];    // list of free neighbors
      int used = get_free_nn(l, li, n_list, 1.e-07, 0.01, barycenter, -1, -1);

      /* cout<<"In grow: free neighbor list for cell "<<li<<": ";
       * for (int n=0; n<n_list.benutzt(); n++) cout<<n_list[n]<<" ";
       * cout<<"\n";
       */
      // if there has been found at least one neighbor go through the selection:
      if (used > 0) {
         /* Possibility 1: delete those free neighbors that have less than maximum
          * contact to other fragments of the same cell,
          * and choose one of the neighbors randomly */
         newind = n_list[irandom(used)];
         // cout<<"in grow: add new cell at lattice point "<<newind<<"\n";
         /* Possibility 2: choose among all neighbors randomly using a
          * weight function. This is the most general possibility! */
         // The selected neighbor is saved in newind (lattice index)

         /* Save a new fragment on the lattice point found */
         l.set_knot(newind, celltype, li);
         // l.knot[newind].cell=celltype;
         // l.knot[newind].listi=li;
         fragments[volume] = newind;
         ++volume;
         // Correct for the shifted barycenter
         index = get_barycenter(l);
         // ## Eventually better to caculate relative changes ...
         // Correct the changed radius
         get_radius(l.dim);
      } else {
         // if no free neighbor exists break
         err = 1;
      }
   } else if ((volume > target_vol) && (drandom() < p_shrink)) {
      // ###
   }

   // cout<<"after grow of cell "<<li<<" ... ";
   // fragment_consistency();
   if (volume == FRAGMENT_STEP) {
      cout << "Warning: FRAGMENT_STEP=" << FRAGMENT_STEP
           << " has been reached -> change in cells.h!!!\n";
   }

   return err;
}
// ============================================================

short frag_cell::do_mitosis(states celltype,
                            const long &i,
                            const long &li,
                            const long &newli,
                            frag_cell &newCell,
                            space &l) {
   /* Assume that i is a soma point (barycenter of cell)
    * and that the cell is composed of at least two fragments.
    * newCell is a copy of the *this cell on the corresponding cell-list,
    * which will be accordingly adapted.
    * With recursive neighbor control.
    */
   // cout<<"Proliferation ... ";
   short err = 1;
   long j;
   short n;
   int fragk;

   /* define the first fragment of the new cell
    * which is neighbor of the old center at lattice point i
    * and belongs to the same cell*/
   // randomize the neighbors
   l.nn = l.nn_permuts.get_rand_set();
   // go through the neighbors until one belonging to the same cell is found
   n = 0;
   while (n < l.dim2 && err != 0) {
      j = l.knot[i].near_n[int (l.nn[n])];
      if ((j != -1) && (l.cellknot[j].cell == celltype) && (l.cellknot[j].listi == li)) {
         // Initialize the values for the new cell ...
         // Change lattice index in newCell to new center j
         newCell.index = j;
         // initialize fragments (to empty) and add the center only:
         // (so fragment lists handed over in newCB will be destroyed!!!)
         newCell.volume = 1;
         newCell.fragments[0] = j;
         // save the barycenter coordinates of the newCell
         l.get_koord(j, newCell.barycenter);
         // Actualize lattice point
         l.set_knot(j, celltype, newli);
         // l.knot[j].cell=celltype;
         // Write on cell_list und save list-index on the lattice
         // l.knot[j].listi=newli;
         // ... done.
         /*
          * cout<<"initialize newli:\nCell li: frags= ";
          * for (int a=0; a<volume; a++) cout<<fragments[a]<<" ";
          * cout<<"\n";
          * cout<<"Cell newli: frags= ";
          * for (int a=0; a<newCell.volume; a++) cout<<newCell.fragments[a]<<" ";
          * cout<<"\n"; */

         // ... and kill the fragment providing the center of the new cell from the old
         fragk = where_fragment(j);
         del_fragment(fragk);
         /*
          * cout<<"after del_center:\nCell li: frags= ";
          * for (int a=0; a<volume; a++) cout<<fragments[a]<<" ";
          * cout<<"\n";
          * cout<<"Cell newli: frags= ";
          * for (int a=0; a<newCell.volume; a++) cout<<newCell.fragments[a]<<" ";
          * cout<<"\n"; */
         /* now go through all direct neighbors of center j (newCell)
          * and recursively check their neighbors */
         attribute_neighbors(celltype, li, j, newli, newCell, l);
         /* Some fragments of li may be disconnected of the cell now, as
          * the criterion for a switch is the distance to the center only.
          * Therefore fragments nearer to li but connected to newli are
          * still on the li-fragment list.
          * -> Check connection and switch all unconnected fragments to cell newli */
         /*
          * cout<<"after attribute:\nCell li: frags= ";
          * for (int a=0; a<volume; a++) cout<<fragments[a]<<" ";
          * cout<<"\n";
          * cout<<"Cell newli: frags= ";
          * for (int a=0; a<newCell.volume; a++) cout<<newCell.fragments[a]<<" ";
          * cout<<"\n"; */

         long * tmpfrags;
         tmpfrags = new long[FRAGMENT_STEP];
         int tmpvol = volume;
         for (int a = 0; a < FRAGMENT_STEP; a++) {
            tmpfrags[a] = fragments[a];
         }

         // At first delete the center of the cell from the fragment-list
         fragk = where_fragment(i, tmpfrags, tmpvol);
         if (fragk == -1) {
            cout << "cell does not contain its own center!!!";
            exit(1);
         } else {
            del_fragment(fragk, tmpfrags, tmpvol);
         }
         check_connection(celltype, i, li, tmpfrags, tmpvol, l);

         // Unconnected fragments are saved in tmpfrags
         for (int x = 0; x < tmpvol; x++) {
            // save fragment in newli
            newCell.fragments[newCell.volume] = tmpfrags[x];
            ++newCell.volume;
            // Change the map to the cell list
            l.cellknot[tmpfrags[x]].listi = newli;
            // What about the status in the grid class here!
            /* As in the previous similar remark, all initial subunits belong to the same cell
             * and should therefore be in the correct gridpoint-state (i.e. object).
             * But may be I missed something?
             * --> Do an automated grid-state check! Did it 5.1.2005. ! Remove after a while. !
             ###
             */
            if (l.knot[tmpfrags[x]].status != object) {
               cout << "Wrong grid-status " << l.knot[tmpfrags[x]].status << " do_mitosis(..)"
                    << " instead of object=1 \n";
               exit(1);
            }
            // ### Remove up to here.
            // Delete this fragment from the list of cell li
            int fragk = where_fragment(tmpfrags[x]);
            del_fragment(fragk);
         }
         /*
          * cout<<"after saving:\nCell li: frags= ";
          * for (int a=0; a<volume; a++) cout<<fragments[a]<<" ";
          * cout<<"\n";
          * cout<<"Cell newli: frags= ";
          * for (int a=0; a<newCell.volume; a++) cout<<newCell.fragments[a]<<" ";
          * cout<<"\n"; */

         // Calculate new bary-centers
         index = get_barycenter(l);
         newCell.index = newCell.get_barycenter(l);
         // newCell.index=newCell.get_barycenter(celltype,newli,l);

         // Actualise the radii of both cells
         get_radius(l.dim);
         newCell.get_radius(l.dim);

         err = 0;
         delete[] tmpfrags;
      }
      // next try
      ++n;
   }
   // cout<<" done.\n";
   return err;
}
// ============================================================
// ============================================================
// ============================================================
// ===========================================================
// ======= Vergleichsoperatoren ==============================
// ============================================================

char operator ==(const cell &a, const cell &b) {
   char back = 1;
   if ((a.index != b.index) || (a.pos_ss != b.pos_ss)) {
      back = 0;
   }
   return back;
}
// ============================================================

char operator !=(const cell &a, const cell &b) { return !((a == b) == 1); }

// ============================================================
