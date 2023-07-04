#include "cellthis.h"
#include <string.h>
#include <math.h>

// ============================================================
// ============================================================
// ============================================================
// ============================================================
// ============================================================
// ============================================================

immunoglobulin_class::immunoglobulin_class() { Ig_class = IgM; }

immunoglobulin_class::~immunoglobulin_class() { }

short immunoglobulin_class::do_switch_classes = 0;
double immunoglobulin_class::switch_matrix[immunoglobulin_class::matrix_dimension];

void immunoglobulin_class::load_matrix(const Parameter &par, ofstream &ana) {
   do_switch_classes = par.Value.do_switch_classes;
   for (int i = 0; i < matrix_dimension; i++) {
      switch_matrix[i] = par.Value.switch_matrix[i];
   }
   if (do_switch_classes == 0) {
      ana << "Class switching is off.\n";
      cout << "Class switching is off.\n";
   } else {
      ana << "Class switching is on. ";
      cout << "Class switching is on --> " << do_switch_classes << "\n";
      if (do_switch_classes == 1) {
         ana << "BCs switch upon selection by TFH.\n";
      } else if (do_switch_classes == 2) {
         ana << "BCs switch upon division.\n";
      }
      ana << "Switch probabilities:\n";
      for (short i = 0; i < nIg_classes; i++) {
         for (short j = 0; j < nIg_classes; j++) {
            Ig_classes ic = Ig_classes(i);
            Ig_classes jc = Ig_classes(j);
            ana << switch_matrix[get_matrix_index(ic,jc)] << " ";
            //	cout<<switch_matrix[get_matrix_index(ic,jc)]<<" ";
         }
         ana << "\n";
         //      cout<<"\n";
      }
   }
}
int immunoglobulin_class::get_matrix_index(const Ig_classes &i, const Ig_classes &j) {
   // The switch_matrix is 2D with dimensions nIg_classes x nIg_classes.
   // Each line contains the switch probability to any of the target classes.
   // The switch probability from class i to class j is in line i and column j.
   // The numbering is line by line
   return nIg_classes * i + j;
}
void immunoglobulin_class::class_switch() {
   // switch with a given probability derived from the switch_matrix
   /* switching matrix:
    * from   IgM   IgG   ....  IgA
    * IgM   s_MM  s_MG        s_MA
    * IgG   s_GM  s_GG        s_GA
    * ...
    * IgA   s_AM  s_AG        s_AA
    * where s_ii is the probability for not switching,
    * s_ij with j\neq i is the probability of switching from i to j,
    * \sum_j s_ij = 1.
    */
   // cerr<<"CS: "<<Ig_class<<" -> ";
   bool done = false;
   Ig_classes j = IgM;
   double rval = drandom();  // get value between [0,..,1].
   double prob = switch_matrix[get_matrix_index(Ig_class,j)];
   //cout<<"Switch from "<<Ig_class<<" with rval="<<rval<<" to ";
   while (not (done) && j < nIg_classes) {
     //cout<<"(j="<<j<<",p="<<prob<<");";
      if (rval < prob) { done = true; } else {
         j = Ig_classes(int (j) + 1);
         prob += switch_matrix[get_matrix_index(Ig_class,j)];
      }
   }
   // save the new class
   //if (done) { Ig_class = j; }
   Ig_class = j;
   //cout<<" --> result class ="<<Ig_class<<", with j="<<j<<endl;
   // note that if no switch is chosen, then i==Ig_class anyway
   // such that no if condition is needed for this attribution.
   // cerr<<Ig_class<<". ";
}
void immunoglobulin_class::set_class(const immunoglobulin_class &c) { Ig_class = c.Ig_class; }

// ============================================================
// ============================================================
// ============================================================
// ============================================================
// ============================================================
// ============================================================

double cellCB::p_pro = 0.;
double cellCB::delta_p_pro = 0.;
double cellCB::limit_volume = 1.0;
double cellCB::p_mut = 0.;
double cellCB::p_mut_after_selection = 0.;
double cellCB::p_mut_after_dec_selection = 0.;
double cellCB::p_mut_affinity_exponent = 0.;
bool cellCB::p_mut_affinity_dependent = false;
double cellCB::p_difu = 0.;
double cellCB::p_difu_width = -1.;
double cellCB::p_tension = 0.;
double cellCB::p_grow = 0.;
double cellCB::p_shrink = 0.;
double cellCB::p_dif = 0.;
double cellCB::p_dif_target = 0.;
double cellCB::start_differentiate = 72.;
double cellCB::max_pro_distance = 0.;
double cellCB::diffusion_tolerance_min = 0.01;
double cellCB::diffusion_tolerance_steepness = 0.01;
double cellCB::elongation = 1.0;
double cellCB::K_elongation = 2.0;
double cellCB::smoothmove = 1.0;
double cellCB::persistence = 0.0;
int cellCB::target_volume = 1;
short cellCB::receptor_use = 0;
double cellCB::receptors = 0.0;
double cellCB::receptor_activation = 0.0;
double cellCB::receptor_binding = 0.0;
double cellCB::receptor_dissociation = 0.0;
double cellCB::p_CXCR4down = 0.;
short cellCB::v_modi = 1;
short cellCB::n_v_states = 1;
double cellCB::v_slow_factor = 1.0;
double cellCB::p_switch_v = 0.0;
double cellCB::max_adhesion = 0.0;
double cellCB::average_seeder_affinity = 0.0;
bool cellCB::SMOOTH_DIFFERENTIATION = false;
double cellCB::smooth_differentiation_time = 3.0;
double cellCB::ag_preloaded = 0;
bool cellCB::transmit_CCdelay2cellcycle = false;
double cellCB::total_n_of_divisions = 0;
double cellCB::total_n_of_DEC205_divisions = 0;
double cellCB::dtphase[cb_statenumber] = { 0 };
int cellCB::dtphase_frequency[phase_number][dtphase_resolution] = {
   { 0 }
};
double cellCB::fraction_of_phase_width = 0.0;
bool cellCB::ag_loaded_CB_diff2output = false;
bool cellCB::ag_loaded_CB_stop_mutation = false;
double cellCB::asymmetric_polarity_index = 1.0;
double cellCB::smooth_PI = 0.0;
double cellCB::IgE_factor_cellcycle = 1.0;
double cellCB::IgE_factor_divisions = 1.0;
long cellCB::attributed_n_of_divisions[max_n_of_divisions + 1] = { };
long cellCB::cummulative_attributed_n_of_divisions[max_n_of_divisions + 1] = { };
long cellCB::attributed_mutation_prob[mutation_bins] = { };
long cellCB::cummulative_attributed_mutation_prob[mutation_bins] = { };
//MSchips
long cellCB::selfBC_attributed_n_of_divisions[max_n_of_divisions + 1] = { };
long cellCB::attributed_mutation_prob_SELF[mutation_bins] = { };
long cellCB::cummulative_attributed_mutation_prob_SELF[mutation_bins] = { };
///MODIFY THIS
short cellCB::TFR_CC_interaction_mode;

cellCB::cellCB() {
   // cout<<"in cellCB default constructor ...";
   state = cb_normal;
   receptor_ligand = 0.0;
   time_of_cycle_state_switch = 0.;
   cycle_state_time = 0.;
   // IgX.Ig_class=IgM;
   retained_ag = 0.;
   collected_ag_portions.clear();
   iamhighag = false;
   n_divisions2do = 0;
   DEC205 = false;
   DEC205_ova = false;
   MHCdeficient = false;
   diff2output = false;
   set_p_move();
   get_adhesion(max_adhesion);

   //MSchips
   selfMutation = 0;
   redeemed = 0;
   inheritance = 0;
   // cout<<"end of cellCB default constructor.\n";
}
// ============================================================

cellCB::cellCB(const cellCB &x) : frag_cell(x) {
   /// Philippe Important : I put the constructor but it seems dangerous because will first do the
   // frag_cell copy!!
   // cout<<"in cellCB::cellCB(cellCB&) ...";
   operator =(x);
   // cout<<"end of cellCB::cellCB(cellCB&)\n";
}
cellCB::cellCB(const cellCC &x) {
   // cout<<"in cellCB::cellCB(cellCC&) ...";
   operator =(x);
   // cout<<"end of cellCB::cellCB(cellCC&)\n";
}
// ============================================================

cellCB&cellCB::operator =(const cellCC &x) {
   // cout<<"in cellCB::=(cellCC&) ... ";
   cell::operator =(x);
   volume = 1;
   make_CB_new();

   DEC205 = x.DEC205;
   DEC205_ova = x.DEC205_ova;
   MHCdeficient = x.MHCdeficient;
   IgX.set_class(x.IgX);
   //MSchips
   selfMutation = x.selfMutation;
//   redeemed = 0; // redemption is forgotten once cell has gone through a GC-cycle
   redeemed = x.redeemed; // redemption is passed CB-->CC-->CB
   ccl3KO = x.ccl3KO;
   inheritance = 1;
   nTFRcontacts = x.nTFRcontacts;

   double ndivtmp = total_n_of_divisions;
   if (x.DND >= 0) {
      ndivtmp = x.DND;
      /* This includes the case of DEC205_ova, where DND was 
       * already attributed to the CC at the time of selection.
       */
   } else {
     // The case of DEC205_ova needs only to be treated if DND is off.
     if (DEC205_ova) { ndivtmp = total_n_of_DEC205_divisions; }
   }
   if (IgX.Ig_class == IgE) { ndivtmp *= IgE_factor_divisions; }
   n_divisions2do = int (ndivtmp);
   ndivtmp -= (double (n_divisions2do));
   if (drandom(1.) < ndivtmp) { ++n_divisions2do; }
   if (tmx_MHC && tmx_MHC_noDivision && MHCdeficient) { n_divisions2do = 0; }

   // show and statistics
   // cout<<"n_divisions2do="<<n_divisions2do<<" --> ";
   int ndivtmpi = n_divisions2do;
   if (ndivtmpi > max_n_of_divisions) { ndivtmpi = max_n_of_divisions; }
   attributed_n_of_divisions[ndivtmpi]++;
   cummulative_attributed_n_of_divisions[ndivtmpi]++;

   //MSchips
   if (selfMutation) selfBC_attributed_n_of_divisions[ndivtmpi]++;

   
   receptor_ligand = 0.0;  // ### 0.0 is the only way as long as CC are without receptors!
                           // ### May be useful to write a receptor-class.

   retained_ag = double (x.nFDCcontacts);  // transmit the number of collected ag to the CB
   collected_ag_portions = x.collected_ag_portions;
   diff2output = false;  // standard no forced differentiation to output
   fragments[0] = x.index;

   set_p_move();
   t_immobile[0] = 0;
   // for (short i=0; i<3; i++) barycenter[i]=0.0;
   /* !!! Achtung: die aufrufende Funktion muss die Koordinaten noch berechnen!
    * Der Befehl fuer hier ist
    * for (short d=0; d<l.dim; d++) barycenter[d]=double(l.knot[x.index].x[d]);
    * (x[d] auf dem lattice ist der einzige Ort wo die Koordinate einer single
    * subunit Zelle gespeichert ist). Da CC immer ein single subunit-Objekt ist,
    * genuegt es einerseits diese Ortskoordinate zu kennen, es ist aber andererseits
    * ein lattice-object zu uebergeben, um diese zugaenglich zu machen.
    */

   // cout<<"end of cellCB::=(cellCC&)\n";
   return *this;
}
// ============================================================

cellCB&cellCB::operator =(const cellCB &x) {
   // cout<<"in cellCB::=(cellCB&) ... ";
   frag_cell::operator =(x);
   state = x.state;
   // put the new cell in state G1 or cb_normal in the calling routine by evoking make_CB_new();
   cycle_state_time = x.cycle_state_time;
   time_of_cycle_state_switch = x.time_of_cycle_state_switch;
   // eventually reset the cell cycle clock to zero in the calling routine
   // ++++++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++++++++++++++++++++++
   n_divisions2do = x.n_divisions2do;  // inherit the number of remaining divisions to be done!
   // ++++++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++++++++++++++++++++++++
   IgX.set_class(x.IgX);
   receptor_ligand = x.receptor_ligand;
   iamhighag = x.iamhighag;
   retained_ag = x.retained_ag;  // has to be handled in the calling routine
   collected_ag_portions = x.collected_ag_portions;
   DEC205 = x.DEC205;
   DEC205_ova = x.DEC205_ova;
   MHCdeficient = x.MHCdeficient;
   set_p_move();
   diff2output = x.diff2output;
   //MSchips
   selfMutation = x.selfMutation;
   redeemed = x.redeemed;
   inheritance = x.inheritance;
   nTFRcontacts = x.nTFRcontacts;
   ccl3KO = x.ccl3KO;
   // cout<<"end of cellCB::=(cellCB&)\n";
   return *this;
}
void cellCB::make_CB_new() {
   /* This routine may be called after the operator= and adapts the copy of a CB
    * to the state associated with a daugther cell after division.
    * It is automatically called when a CB is made from a CC.
    */
   state = cb_normal;
   if (fixed_number_of_divisions()) {
      // cout<<"set new CB .. ";
      state = cb_G1;
      time_of_cycle_state_switch = set_cycle_state_duration(state);
   }
   cycle_state_time = 0.;
}
void cellCB::set_remaining_divisions() {
   if (fixed_number_of_divisions()) {
      n_divisions2do--;
      if (n_divisions2do <= 0) {
         state = cb_stop_dividing;
      }
   }
}
centroblasts cellCB::set2differentiation() {
   if (fixed_number_of_divisions()) { return cb_stop_dividing; } else { return cb_differentiate; }
}
bool cellCB::fixed_number_of_divisions() {
   if (total_n_of_divisions > 0) { return true; } else { return false; }
}
double cellCB::total_cell_cycle_duration() {
   double ccd = 0.;
   ccd += dtphase[cb_G1];
   ccd += dtphase[cb_S];
   ccd += dtphase[cb_G2];
   ccd += dtphase[cb_M];
   return ccd;
}
// ============================================================

cellCB::~cellCB() {
   // cout<<"in ~cellCB()...\n";
}
// ============================================================

void cellCB::set_statics(const Parameter &par, space &l, ofstream &ana) {
   // get target volume of CBs in number of lattice points
   if (l.dim
       == 2) {
      target_volume = int (3.1415 * pow(par.Value.CB_radius,2.) / pow(l.dx,2.) + 0.5);
   } else { target_volume = int (4.1888 * pow(par.Value.CB_radius,3.) / pow(l.dx,3.) + 0.5); }
   if (target_volume == 0) { target_volume = 1; }
   ana << "  CB target volume = N = " << target_volume << "\n";
   // Maximaler Abstand fuer die CB-Proliferation in Gittereinheiten
   max_pro_distance = par.Value.dx_CB / par.Value.dx;

   ana << "Calculate action probabilities ... \n";
   ana << "  CB somatically hypermutate from t=" << par.Value.Start_Mutation << " with "
       << "probability " << par.Value.mutation << "\n";
   ana << "  CB differentiation to CC (from t=" << par.Value.Start_Differentiation << ") = "
       << par.Value.tolight * par.Value.deltat << "\n";
   // for centroblasts
   p_pro = par.Value.proliferate * par.Value.deltat;
   // +++ OPTION: diminish proliferation rate with time
   delta_p_pro = 0.;
   // delta_p_pro=p_pro/((par.Value.tmax-par.Value.tmin)/par.Value.deltat);
   // end OPTION
   ana << "  CB proliferation = " << p_pro << "\n";

   if (par.Value.CB_fixed_times_of_divisions > 0.) {
      /* The following is from the time before resolution of cell cycle phases:
       * It allows to prolong the duration of cycling of CB.
       * A similar thing could be used when combining CB differentiation by rate
       * with prolonged duration of division phase in respose to DEC205-OVA.
       * This is not yet programmed (see get_new_state to do so). ### */
      /*
       * cell_cycle_delay=long((log(2.)*par.Value.CB_fixed_times_of_divisions/
       *         (par.Value.proliferate*par.Value.deltat))+0.5);
       * ana<<"  CB differentiation is delayed by "<<par.Value.CB_fixed_times_of_divisions
       * <<" cell cycle times,\n"
       * <<" corresponding to "<<cell_cycle_delay<<" time steps of "
       * <<par.Value.deltat<<" hours.\n";
       * DEC205_cell_cycle_delay=long(par.Value.DEC205_p_factor*double(cell_cycle_delay));
       */

      dtphase[cb_G0] = par.Value.CB_dt_G0;
      dtphase[cb_G1] = par.Value.CB_dt_G1;
      dtphase[cb_G2] = par.Value.CB_dt_G2;
      dtphase[cb_S] = par.Value.CB_dt_S;
      dtphase[cb_M] = par.Value.CB_dt_M;
      fraction_of_phase_width = par.Value.CB_dtphase_width;
      total_n_of_divisions = par.Value.CB_fixed_times_of_divisions;
      total_n_of_DEC205_divisions = total_n_of_divisions * par.Value.DEC205_p_factor;
      ana << "  Number of divisions fixed to " << total_n_of_divisions << ".\n"
          << "  Proliferation rate is ignored.\n";
   } else { 
     ana << "  No delay of CB differentiation for division, thus, rate based division.\n"; 
   }

   IgE_factor_cellcycle = par.Value.IgE_factor_cellcycle;
   if (immunoglobulin_class::do_switch_classes && (IgE_factor_cellcycle != 1.0)) {
      ana << "IgE BC have a cell cycle duration corrected by the factor "
          << IgE_factor_cellcycle << ".\n";
   }
   IgE_factor_divisions = par.Value.IgE_factor_divisions;
   if (immunoglobulin_class::do_switch_classes && (IgE_factor_divisions != 1.0)) {
      ana << "IgE BC have a number of divisions per round corrected by the factor "
          << IgE_factor_divisions << ".\n";
   }

   transmit_CCdelay2cellcycle = par.Value.transmit_CC_delay_to_CB_cycle;
   ag_preloaded = par.Value.BC_ag_preloaded;
   ag_loaded_CB_diff2output = par.Value.ag_loaded_CB_diff2output;
   ag_loaded_CB_stop_mutation = par.Value.ag_loaded_CB_stop_mutation;
   asymmetric_polarity_index = par.Value.asymmetric_polarity_index;
   smooth_PI = par.Value.smooth_PI;

   limit_volume = par.Value.CB_maxvolume4differ2CC;
   ana << "  CB limit volume 4 differentiation = " << limit_volume << " x " << target_volume
       << "\n";
   if (fixed_number_of_divisions()) {
      if (target_volume != 1) {
         cout
         <<
         "ERROR: Growth rate in the case of fixed cell cycle phase duration not defined.\n"
         << "       See cellCB-constructor in cellthis.C.\n"
         << "ABORT!\n\n";
         exit(1);
         // ##### Think of what is the right abort criterion here.
      }
   } else {
      p_grow = par.Value.grow * par.Value.deltat;
      ana << "  CB growth = " << p_grow << "\n";
      p_shrink = par.Value.shrink * par.Value.deltat;
      ana << "  CB shrink = " << p_shrink << "\n";
      // get the corrected proliferation rate including mitosis only (no growth)
      // p_proliferate=p_proliferate/(1.-growth_fraction);
      if (target_volume > 1) {
         // ignore this in the case of one cell = one lattice point
         p_pro = p_pro * p_grow / (p_grow - target_volume * p_pro);
      }
      ana << "   Probability of mitosis after completed growth = " << p_pro << "\n";
      // calculate growth-phase fraction
      // ana<<par.Value.proliferate<<"  "<<par.Value.grow<<"   "<<l.dx<<"   "<<CB_vol<<"\n";

      // double growth_fraction=0.5*CB_vol*par.Value.proliferate/par.Value.grow;
      double growth_fraction = target_volume * par.Value.proliferate / par.Value.grow;
      if (target_volume == 1) { growth_fraction = 0.; }
      ana << "   Fraction of growth phase in the total cell cycle "
      // <<"p*V_CB/(2 p_growth) = "
          << growth_fraction << "\n";
      if ((growth_fraction > 0.5) || (p_pro < 0.)) {
         cout << "Proliferation rate too large or CB growth rate to small!\n";
         exit(1);
      }
   }

   ///// PHILIPPE 13/04/2016
   cout << "The mutation probability is taken from :";
   if (par.Value.use_sequence_space) {
      p_mut = par.Value.sequence_mut_per_base * par.Value.size_sequences;
      cout << "Sequence Space ";
   } else if (par.Value.use_arup_space) {
      p_mut = par.Value.arup_mutation * par.Value.arup_length_sequences;
      cout << "Arup Space ";
   } else {
      p_mut = par.Value.mutation;
      cout << "Shape Space ";
   }
   cout << " and = " << p_mut << " (per total sequence per division)\n";

   // Mutation //#####MMH: This overwrites what was done above!!! ERROR!? Philippe?
   p_mut = par.Value.mutation;
   p_mut_after_selection = par.Value.mutation_after_tc;
   if (p_mut_after_selection < 0) { p_mut_after_selection = p_mut; }
   // ++++++dec++++++++++
   p_mut_after_dec_selection = par.Value.mutation_after_dec_tc;
   // keep the value smaller 0 of p_mut_after_dec_selection
   // ++++++dec++++++++++
   p_mut_affinity_exponent = par.Value.mutation_affinity_exponent;
   if (p_mut_affinity_exponent > 1.e-08) {
      p_mut_affinity_dependent = true;
   } else {
      p_mut_affinity_dependent = false;
      p_mut_affinity_exponent = 0.0;
   }

   // Differentiation to CC
   start_differentiate = par.Value.Start_Differentiation;
   p_dif_target = par.Value.tolight * par.Value.deltat;
   if (fixed_number_of_divisions()) { p_dif = p_dif_target; }

   SMOOTH_DIFFERENTIATION = par.Value.smooth_differentiation;
   smooth_differentiation_time = par.Value.smooth_differentiation_time;
   if (SMOOTH_DIFFERENTIATION && (fixed_number_of_divisions() == false)) {
      ana << "Smooth onset of CB differentiation to CC (width="
          << smooth_differentiation_time << " hours).\n";
   } else { ana << "Abrupt onset of CB differentiation to CC.\n"; }

   // Adhesion:
   // =========
   max_adhesion = par.Value.CB_max_adhesion;

   // Motility:
   // =========
   // Toleranzparameter fuer die Abweichung von der Kugelform
   diffusion_tolerance_min = par.Value.distance_tolerance;
   ana << "  CB diffusion_tolerance_min  = t_{\rm min} = " << diffusion_tolerance_min << "\n";
   /* get the exponential factor in the formula for the tolerance
    * for which the half-tolerance-deformation given in the par-file is respected: */
   diffusion_tolerance_steepness = -par.Value.half_tolerance_deformation
                                   / (log((1.0 - diffusion_tolerance_min)
                                          / (2.0 - diffusion_tolerance_min
                                             * (1.0 + diffusion_tolerance_min))));
   ana << "  CB diffusion_tolerance_steepness = K_{1/2} = " << diffusion_tolerance_steepness
       << "\n";
   // elongation=par.Value.CB_elongation*par.Value.CB_radius/l.dx;
   // minimal elongation is 1 lattice constant!
   // if (elongation<1.) elongation=1.;
   elongation = par.Value.CB_elongation;
   ana << "  CB elongation = epsilon = " << elongation << "\n";
   K_elongation = par.Value.CB_K_elongation;
   ana << "  CB elongation for half reshaping force = K_epsilon = " << K_elongation << "\n";
   /*
    * if (par.Value.CB_radius/l.dx<elongation) {
    * ana<<"!!!Warning!!! Elongation leads to barycenter shifts that are larger than r_CB!\n";
    * cout<<"!!!Warning!!! Elongation leads to barycenter shifts that are larger than r_CB!\n";
    * }
    */
   smoothmove = par.Value.CB_smoothmove;
   ana << "  CB smoothness of active walk = eta_max = " << smoothmove << "\n";
   if (par.Value.CB_persistence > 60. * par.Value.deltat) {
      persistence = 60. * par.Value.deltat / par.Value.CB_persistence;
   } else {
      persistence = 1.;   // i.e. change polarity in every time step!
   }
   ana << "  CB polarity persistence = " << par.Value.CB_persistence
       << " i.e. p_{change-polarity} = " << persistence << "\n";

   if (par.Value.CB_D_cytosol >= 0.) {
      p_tension = double (l.dim2) * par.Value.CB_D_cytosol * par.Value.deltat
                  / (par.Value.dx * par.Value.dx);
   } else {
      p_tension = 60. * par.Value.v_CB_cytosol * par.Value.deltat / par.Value.dx;
   }
   ana << "  CB surface tension (tau): Fragment movement probability = " << p_tension << "\n";
   if (p_tension > 0.5) {
      cout << "Centroblast-Cytosol-Fragment-Movement-Probability (p="
           << p_tension
           << ") is too large for dx and dt !!!\n";
      exit(1);
      /* Hier wird die Wahrscheinlichkeit berechnet, mit der ein Centroblast
       * in irgendeine Richtung bewegt. Daher keine Normierung mit 1/(2*l.dim)! */
   }

   if (par.Value.D_CB >= 0.) {
      p_difu = double (l.dim2) * par.Value.D_CB * par.Value.deltat
               / (par.Value.dx * par.Value.dx);
   } else {
      p_difu = 60. * par.Value.v_CB * par.Value.deltat / par.Value.dx;
   }
   p_difu_width = par.Value.v_CB_width * p_difu;
   double dfactor = 0.7;  // approximate value for use_D_correction==0
   if ((use_D_correction == 1) && (target_volume > 1)) { 
     dfactor = (0.18 + (2.0 / double (target_volume))); 
   }
   // Calculate smoothfactor for a spherical cell with target volume:
   double smoothfactor;
   if (l.dim == 2) {
      smoothfactor = smoothmove * par.Value.dx / (6.28 * par.Value.CB_radius);
   } else {
      smoothfactor = smoothmove * par.Value.dx * par.Value.dx
                     / (12.56 * par.Value.CB_radius * par.Value.CB_radius);
   }
   if (p_difu * smoothfactor / (dfactor) > 1.0) {
      cout << "Centroblast-Diffusion (p="
           << p_difu * smoothfactor / (dfactor)
           << ") is too large for dx and dt !!!\n";
      exit(1);
      /* Hier wird die Wahrscheinlichkeit berechnet, mit der ein Centroblast
       * in irgendeine Richtung diffundieren kann. Daher keine Normierung
       * mit 1/(2*l.dim)! */
   }
   ana << "  CB movement probability = "
       << p_difu
       << ", with additional factors = "
       << p_difu * smoothfactor / dfactor
       << "\n";
   double tmpdt = par.Value.CB_persistence;
   if (tmpdt < 60. * par.Value.deltat) {
      tmpdt = 60. * par.Value.deltat;
      ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
   }
   if (par.Value.D_CB < 0.) {
      ana << "  Expected effective diffusion constant for CB = "
          << 60. * par.Value.v_CB * par.Value.v_CB * tmpdt / double (l.dim2)
          << " microns^2/hr = "
          << par.Value.v_CB * par.Value.v_CB * tmpdt / double (l.dim2)
          << " microns^2/min\n";
   } else {
      ana << "  Expected mean cell velocity for CB = "
          << sqrt(double (l.dim2) * par.Value.D_CB / (60. * tmpdt))
          << " microns/min\n";
   }
   v_modi = par.Value.CB_v_modi;
   n_v_states = par.Value.CB_n_v_states;
   if (n_v_states > 1) {
      v_slow_factor = par.Value.v_CB_factor;
      p_switch_v = 60. * par.Value.deltat / par.Value.v_CB_switch_deltat;
      if (p_switch_v < 0.) { p_switch_v = 0.; }
   }

   p_CXCR4down = par.Value.CXCR4down * par.Value.deltat;
   if (p_CXCR4down < 0) { p_CXCR4down = 0.; }
   ana << "  CB CXCR4 down regulation probability = " << p_CXCR4down << "\n";

   // end of motility.
   // ================

   receptor_use = par.Value.CBreceptor_use;
   if (receptor_use > 0) {
      ana << "Initialize receptors ... ";
      receptor_activation = par.Value.CBreceptor_activation;
      receptor_binding = par.Value.CBreceptor_binding;
      receptor_dissociation = par.Value.CBreceptor_dissociation;
      receptor_activation *= par.Value.CBreceptor_total;
      if (receptor_use == 2) {
         receptor_binding *= par.Value.deltat;
         receptor_dissociation *= par.Value.deltat;
         if ((receptor_binding > 1.0) || (receptor_dissociation > 1.0)) {
            cout << "Receptor dynamics: deltat*k+=" << receptor_binding
                 << " deltat*k-=" << receptor_dissociation << " too large!\n";
            exit(1);
         }
      }
      /* Note: if receptor_use=1 then receptor_binding unimportant
       *                          and receptor_dissociation=K
       *              receptor_activation is fraction
       *   if receptor_use=2 then receptor_binding=k_+*deltat
       *                          and receptor_dissociation=k_-*deltat
       *              receptor_activation is total number
       */
      // fixed total number of receptors
      receptors = par.Value.CBreceptor_total;
      ana << "done.\n";
   }


   //MSchips
   ///MODIFY THIS
   TFR_CC_interaction_mode = par.Value.TFR_mode;


   // p_mutate wird durch Aufruf von set_pars von der Zeit abhaengig gesetzt
   // p_differ2CC wird durch Aufruf von set_pars von der Zeit abhaengig gesetzt
}
void cellCB::set_statics(const double &time, const Parameter &par) {
   if (fixed_number_of_divisions()) { p_dif = p_dif_target; } //0.02
   else if (SMOOTH_DIFFERENTIATION) {
      set_differentiation(time);
   } else if (time >= par.Value.Start_Differentiation - 1.0e-09) {
      p_dif = p_dif_target;
      cout << "Differentiation to centrocytes (t=" << time << " hr) (discrete), ";
      // Random generator initialization after start of differentiation
      for (long n = 0; n < par.Value.late_ini_random; n++) {
         irandom(1);
      }
   } else { p_dif = 0.; }
}
void cellCB::set_differentiation(const double &time) {
   /* This routine has to be called at the beginning of a time step only once
    * because it sets a static variable. Note that the starting time
    * of differentiation is the half value of the sigmoidal.
    */
   if (fixed_number_of_divisions()) { p_dif = p_dif_target; }//0.02
   else if (SMOOTH_DIFFERENTIATION) {
      /// Philippe : this is redundant with CellCB::setstatics
      // Use a sigmoidal to switch differentiation on:
      // time is in hours
      p_dif = p_dif_target
              / (1.0 + exp((start_differentiate - time) / smooth_differentiation_time));
   }
}
// ============================================================

void cellCB::ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape) {
   /* Diese ini-Routine geht von einer bereits vordefinierten Zelle
    * aus, die in dem Objekt gespeichert ist. Es werden an dem Objekt
    * die Aenderungen vorgenommen, die es zu einem neuen Zellobjekt
    * machen. Ausserdem werden lattice und shapespace aktualisiert.
    * Achtung: state wird nicht gesetzt!
    */
   // Change lattice index in newCB to i
   index = i;
   born_index = i;
   born_time = t;
   //MS
   ccl3KO=0;
   // initialize fragments (to empty) and add the CB center only:
   // (so fragment lists handed over in newCB will be destroyed!!!)
   volume = 1;
   fragments[0] = i;
   l.get_koord(i,barycenter);
   for (short a = 0; a < l.dim; a++) {
      last_position[a] = barycenter[a];
   }
   // Actualize lattice point
   l.set_knot(i,CB,li);
   // l.knot[i].cell=CB;
   // Write on cell_list und save list-index on the lattice
   // l.knot[i].listi=li;
   // Actualize shape space statistics
   shape.add_cell(sCB,pos_ss);
   shape.add_cell(total,pos_ss);
   // set mutation to standard mutation probability
   p_mutation = p_mut;
   // get initial polarity vector
   l.get_random_direction(polarity);
   // get initial radius
   get_radius(l.dim);
   // cout<<"polar-vector=("<<polarity[0]<<","<<polarity[1]<<")\n";
}
// ============================================================

void cellCB::preload_with_ag() {
   if (ag_preloaded > 0) {
     double width = 0.2 * double(ag_preloaded);
     retained_ag = -1;
      while (retained_ag < 0) {
         retained_ag = get_positive_sample_from_normal(ag_preloaded,width);
      }
      // or use a fixed number:
      // retained_ag=100.;
      /* Attribute this to the vector collected_ag_portions:
       *  This is not differentiated for multiple antigens.
       *  Here, the whole preloaded antigen is attributed to the first antigen.
       *  There is no clear rule how to distribute it better, so any assumption
       *  may be tested here.
       */
      if (collected_ag_portions.size() == 0) { 
	collected_ag_portions.push_back(0); 
      }
      collected_ag_portions[0] += int (retained_ag + 0.5);
   }
}
// =======================================================================
void cellCB::attribute_DEC205(double fracofpos) {
  /* the default value in the CB constructor is DEC205 == false;
     here the value is set true with a probability (parameter) corresponding
     to the fraction of DEC205 positive BCs in the experimental setup
  */
  if (drandom(1.) < fracofpos) { DEC205 = true; }
}
// =======================================================================

double cellCB::set_cycle_state_duration(centroblasts &s) {
  /* This follows
   * https://controls.engin.umich.edu/wiki/index.php/SPC:_random_sampling_from_a_stationary_Gaussian_process
   */
   // The error function is calculated based on ln and sqrt
   if (dtphase[s] <= 0.) { return 0.0; }
   double arg;
   double thisphase = dtphase[s];
   double width = fraction_of_phase_width * thisphase;
   if (IgX.Ig_class == IgE) { thisphase *= IgE_factor_cellcycle; }
   double phaselength = 3.0 * thisphase;
   while (phaselength <= 0. || phaselength >= 2. * thisphase) {
      // This samples from a normal distribution with mean thisphase and width width
      arg = (2. * drandom() - 1.);
      phaselength = thisphase + sqrt(2.) * width * inverse_erf(arg);
   }
   // cout<<"("<<thisphase<<","<<shift<<"); ";
   // double sign=drandom();
   // if (sign<0.5) shift*=-1.0;
   int dt_index = int (phaselength / (2.0 * thisphase / double (dtphase_resolution - 1)) + 0.5);
   if ((dt_index < 0) || (dt_index > dtphase_resolution - 1))
   { cout << "Error in cellCB::set_cycle_state_duration(...) \n"; exit(1); }
   ++dtphase_frequency[s - cb_G1][dt_index];
   // if (s==4) cout<<".......................... phase "<<s<<"  mean="<<dtphase[s]<<"
   //  length="<<phaselength<<";\n";
   return phaselength;
}
centroblasts cellCB::progress_cycle_phase() {
  if (state == cb_G1) { 
    return cb_S; 
  } else if (state == cb_S) {
    return cb_G2;
  } else if (state == cb_G2) {
    return cb_M;
  } else if (state == cb_M) {
    return cb_divide;
  } else if (state == cb_G0) {
    return cb_S;
  }
  return cb_G0;
}
bool cellCB::shiftCCdelay2CBcycle() {
  return (fixed_number_of_divisions() && transmit_CCdelay2cellcycle);
}
centroblasts cellCB::get_virtual_cell_cycle_phase(double waited) {
  centroblasts phase = cb_G1;
  while (phase != cb_divide && waited >= dtphase[phase]) {
    // note that the waiting time is compared to the mean values of the cycle durations
    // which are saved in dtphase[] (and static)
    waited -= dtphase[phase];
    phase = centroblasts(int(phase) + 1);
    if (phase == cb_G0) phase = cb_S;
  }
  if (phase == cb_divide) { phase = cb_G0; }
  return phase;
}
void cellCB::transmit_CCdelay2cycle(double waited_time) {
  if ( shiftCCdelay2CBcycle() ) {
      // cout<<"waited="<<waited_time<<"-> ";
      // find out the position of the cell cycle
      // also do for DEC205_ova positive
      /* states are enum centroblasts{cb_normal,cb_differentiate,
       *    cb_G1,cb_G0,cb_S,cb_G2,cb_M,cb_divide,cb_stop_dividing,
       *    cb_statenumber}; */
      // move forward cell cycle state if waited_time is longer
      // than the duration of the respective phase
      while (state != cb_M && waited_time >= dtphase[state]) {
         waited_time -= dtphase[state];
         state = progress_cycle_phase();
      }
      // the cycle_state_time (so the clock) is set to the remaining of the <waited_time>:
      cycle_state_time = waited_time;
      // the duration of the starting state has to be determined
      time_of_cycle_state_switch = set_cycle_state_duration(state);
      if (cycle_state_time > time_of_cycle_state_switch) {
         cycle_state_time = time_of_cycle_state_switch;
      }
      // cout<<" (state="<<state<<", cycle-time="<<cycle_state_time<<");  ";
   }
}
void cellCB::get_new_state(const long &i, double &dt, space &l, sigs &s) {
   // in the case of differentiation by rate and still in state cb_normal
   if (state == cb_normal) {
      /* This might be used for the case rate differentation
       * plus DEC205 induced prolongation of division period.
       * Needs to be programmed everywhere (not just activating). */
      /*
       * if (cell_cycle_delay==0 ||
       * (DEC205_ova==false && cycling_time>cell_cycle_delay) ||
       * // +++++DEC+++++
       * (DEC205_ova==true && cycling_time>DEC205_cell_cycle_delay)
       * // +++++DEC+++++
       * ) {
       */
      // Ohne Signalstoff wird alles auf Differenzierung gestellt
      if (s.signal_use[sig_differ2CC] == 0) { state = cb_differentiate; }
      /// Philippe: should break here, no ?

      // Der Fall CB im Normalzustand und Differenzierung durch Signalquanta:
      if ((receptor_use == 0)
          && (s.signal_use[sig_differ2CC] == 1)
          && (s.sigsknot[i].signal[sig_differ2CC] >= 1.0)) {
         // cout<<"switch to differentiate\n";
         // Falls Signalmolekuel da, Differenzierung ermoeglichen
         state = cb_differentiate;
         // responsive2signal[CXCL13]=true;
         --s.sigsknot[i].signal[sig_differ2CC];
         --s.sigsknot[i].signal_new[sig_differ2CC];
         ++s.signal_used[sig_differ2CC];

         /* // Fuer fragments: ###
          * // ### Spaeter nicht nur im Zentrum Signal verbrauchen!!!
          * // Vorschlag:
          * // Eliminate one quantum of signal molecules for each CB-fragment
          * for (int ff=0; ff<CB_list[li].volume; ff++)
          * if (l.knot[CB_list[li].fragments[ff]].signal[sig_differ2CC]>0) {
          * // Das Signalmolekuel wird verbraucht:
          * l.knot[CB_list[li].fragments[ff]].signal[sig_differ2CC]--;
          * l.knot[CB_list[li].fragments[ff]].signal_new[sig_differ2CC]--;
          * // Hier ist es wichtig, das Molekuel vollstaendig zu eliminieren, da sonst
          * // das gleiche Molekuel noch diffundieren koennte!!!
          * // ### Vorsicht: Bei Aktivierung ist zu beachten, dass dies die
          * // Ueberpruefung der verbrauchten Signale im Vergleich zu den
          * // Differenzierungsbefehlen modifiziert
          * // --> Folgeaenderung in mk_cell_sum()
          * }
          */
      }
      /// Philippe : suggestion: redefine the options 0,1,2 by enum names ?
      // Der Fall CB im Normalzustand und Ligandierung der Rezeptoren durch Signale:
      if ((receptor_use == 1) && (s.signal_use[sig_differ2CC] == 1)) {
         // Berechne den Anteil der ligandierten Rezeptoren im steady state
         double br = receptor_ligand;
         receptor_ligand = bind_ss_receptor_ligand(sig_differ2CC,
                                                   receptor_dissociation,receptors,
                                                   br,l,s);
         // change mode if above-threshold receptors are liganded
         if (receptor_ligand > receptor_activation) {
            state = cb_differentiate;
            receptor_ligand = 0.0;
         }
      }

      if ((receptor_use == 2) && (s.signal_use[sig_differ2CC] == 1)) {
         /// Philippe : should it be done for multiple fragments as well?
         // Berechne den Anteil der ligandierten Rezeptoren
         double br = receptor_ligand;
         receptor_ligand += bind_receptor_ligand(sig_differ2CC,
                                                 receptor_binding,receptor_dissociation,
                                                 receptors,br,l,s);
         // change mode if above-threshold receptors are liganded
         if (receptor_ligand > receptor_activation) {
            state = cb_differentiate;
            receptor_ligand = 0.0;     // ligand is really used here!
         }
      }

      /// Philippe : don't get that
      // Der Fall CB im Normalzustand und Differenzierung durch FDC-Kontakt
      if (s.signal_use[sig_differ2CC] == -1) {
         // Dies ist die Variante: Differenzierung nach Kontakt zu FDC
         if (contact(FDC,l) == 0) { state = cb_differentiate; }
      }
   }

   // in the case of differentiation by rate and already in state cb_differentiate
   // --> nothing to do

   // in the case of fixed number of divisions:
   if (fixed_number_of_divisions()) {
      // cout<<"in get_new_state(...) --> case fixed_number_of_divisions ... state="<<state<<"\n";
      //  if (state>=cb_G1 && <cb_M) {
      cycle_state_time += dt;
      if (cycle_state_time >= time_of_cycle_state_switch) {
         // cout<<"time="<<cycle_state_time<<", switch at="<<time_of_cycle_state_switch<<",
         // state="<<state<<"; ";
         if (state < cb_divide) {
            // ### here a checkpoint for G0 might be introduced
            state = progress_cycle_phase();
            if (state != cb_divide) {
               time_of_cycle_state_switch = set_cycle_state_duration(state);
            }
            cycle_state_time = 0.;
            // cout<<" new switch="<<time_of_cycle_state_switch<<"; \n";
         }
      }
   }
}
void cellCB::set_mutation_after_TC(AffinitySpace &shape) {
   // cout<<"mutation-aff-dependent="<<p_mut_affinity_dependent<<";
   // mut-freq-after-TC="<<p_mut_after_selection<<"\n";
   if ((p_mut_after_dec_selection >= 0.0) && DEC205_ova) {
      // if p_mut_after_dec_selection gets a specific value and the cell is dec205-ova-activated
      p_mutation = p_mut_after_dec_selection;
   } else {
      if (p_mut_affinity_dependent) {
	// ++++++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++
	/* This version is the default version, which is based on a Hill-function.
	 * With increasing affinity, the mutation probability goes down. */
 	 double rhoaff = shape.best_affinity_norm(pos_ss);
         p_mutation = p_mut + (p_mut_after_selection - p_mut) * pow(rhoaff,p_mut_affinity_exponent);
	 /* In addition, one might activate a reduction of the mutation rate for low affinity cells.
	  * Thus, the mutation probability goes down for large and for low affinity cells. */
         //p_mutation = p_mut * exp(-1.0 * pow((log(rhoaff) + 2.303)/1.0, 2));
	 // ### switch between those with a parameter in the future. 
	 // ++++++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++
         // cout<<"p="<<p_mutation<<";   ";

         if ((TFR_CC_interaction_mode==30||TFR_CC_interaction_mode==31||TFR_CC_interaction_mode==34) && nTFRcontacts>0) {
             if (TFR_CC_interaction_mode==30 || selfMutation) {
                 p_mutation=0.5;
//                    pred=1;
             }
         }

      } else {
         p_mutation = p_mut_after_selection;
         // cout<<"affinity="<<shape.best_affinity_norm(pos_ss)<<"; p_mutation="<<p_mutation<<"\n";
      }
   }

   // save this for the statistics
   int mutbin = int (p_mutation * double (mutation_bins));
   if (mutbin >= mutation_bins) { mutbin = mutation_bins - 1; }
   attributed_mutation_prob[mutbin]++;
   cummulative_attributed_mutation_prob[mutbin]++;
   //MS
   attributed_mutation_prob_SELF[mutbin]++;
   cummulative_attributed_mutation_prob_SELF[mutbin]++;
}
void cellCB::set_CXCR4expression() {
   if (responsive2signal[CXCL12] && (p_CXCR4down > 0) && (drandom() < p_CXCR4down)) {
      responsive2signal[CXCL12] = false;
   }
}
void cellCB::resensitise4CXCL12(sigs &s) {
   if ((s.signal_use[CXCL12] == 1)  // CXCL12 has to be used at all
       && not (responsive2signal[CXCL12])  // otherwise its sensitive anyway
       && (CXCL12recrit > 0)  // re-sensitisation has to be on
       && s.undercritical_signal(index,CXCL12,CXCL12recrit)
       // CXCL12 has to be below the threshold value for re-sensitisation
       ) {
      responsive2signal[CXCL12] = true;
      // cout<<"resensitised CB with index "<<index<<" for CXCL12.\n";
   }
}
void cellCB::adapt_specific_ag(double factor) {
   int sum = 0;
   // Reduce collected antigens by the factor <factor>
   for (unsigned int f = 0; f < collected_ag_portions.size(); f++) {
     collected_ag_portions[f] = int (factor * double(collected_ag_portions[f]) + 0.5);
     sum += collected_ag_portions[f];
   }
   // Now double check whether total collected antigen is still correct
   int diff = sum - int (retained_ag + 0.5);
   if (diff != 0) { retained_ag += double (diff); }
   /* This is the easiest solution to make the representations by retained_ag
    *  and collected_ag_portions consistent again -- if a deviation occured.
    *  However, this might induce a kind of stochastic loss or growth of antigen
    *  in the retained_ag representation. As this will be equilibrated by
    *  large numbers, no effect of this approximation is to be expected.
    */
}
// ============================================================

short cellCB::ask_differentiate() {
   short int err = 1;
//   double rr = drandom();
   // if CB is in the right state
   // then probabilistic decision if differentiation to CC is done:
   if (((state == cb_differentiate) || (state == cb_stop_dividing))
       && (double (volume) / double (target_volume) <= limit_volume)
       && (drandom() < p_dif)) {
//           && ( (selfMutation && rr<100*p_dif)||(!selfMutation && rr<p_dif) ) ) {
//       && selfMutation ) {
      // && volume>0.9*target_volume)
//       if(selfMutation)cerr<<"rr: "<<rr<<" 100*pdif:"<<100*p_dif<<" pdfi: "<<p_dif<<endl;
      err = 0;
      // cout<<".";
   } /*else if ((TFR_CC_interaction_mode==-5 && selfMutation)
              || (TFR_CC_interaction_mode==0 && selfMutation && nTFRcontacts==0)) {err=0;}*/
   if (err == 0) {
      // i.e. if CB stops dividing and starts differentiating
      if (ag_loaded_CB_diff2output && (retained_ag > 0.) && iamhighag) { diff2output = true; }

      //MSchips
      if ( (TFR_CC_interaction_mode==2 && selfMutation && nTFRcontacts>0)
           /*|| (TFR_CC_interaction_mode==3 && nTFRcontacts==0)*/ ) {
          diff2output = false;
      } /*else if ((TFR_CC_interaction_mode==-5 && selfMutation)
                 || (TFR_CC_interaction_mode==0 && selfMutation && nTFRcontacts==0)) {
          diff2output = true;
      }*/
      // cout<<"ag="<<retained_ag<<" induces diff2output="<<diff2output<<"\n";
      // if (DEC205_ova) cout<<"retaind_ag="<<retained_ag<<"; diff2output="<<diff2output<<"\n";
   }
   // if (err==0) cout<<"vol="<<volume<<" target="<<target_volume<<" limit="<<limit_volume<<"\n";
   return err;
}
long cellCB::ask_mitosis(long * pp, space &l) {
   // This routine also makes the probabilistic decision according to the proliferation rate
   // whether division is initiated or not. The probability is p_pro.
   // For IgE BC the probability is adapted with IgE_factor_cellcycle.
   // Shorter cell cycle (factor<1) leads to higher probability p_pro:
   double p_now = p_pro;
   // This allows to test a hypothetic disadvantage of IgE BCs in their division potential:
   if (IgX.Ig_class == IgE) { p_now /= IgE_factor_cellcycle; }
   // Now use this value for either one point objects or multi-fragment objects:
   if ((volume == 1) && (target_volume == 1)) {
      return find_mitosis_place(p_now,state == cb_divide,max_pro_distance,pp,l);
   } else {
     /* case of more fragment object:
      * no need to search for new space, as the volume is conserved during division.
      * proliferation is allowed if the total volume is near the target-volume of a cell.
      * then probabilistic decision whether proliferation will be done: 
      */
      if ((volume > 0.9 * target_volume)
          && ((state == cb_divide) || (drandom() < p_now))) {
         return 0;
      }
      // 
      return 1;
   }
}
void cellCB::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
  if (p_difu_width > 0) {
    p_move = get_positive_sample_from_normal(p_difu, p_difu_width);
  }
}
double cellCB::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   // call polarisation
   // cerr<<"vor set_polarity_velocity(...) ...\n";
   set_polarity_velocity(persistence,v_modi,p_switch_v,n_v_states,v_slow_factor,l,s);
   if (writethis2track == polarisation) {
      td.Write_movement(trackno,CB,time,barycenter,polarity,writethis2track);
      writethis2track = trackini;
   }
   // cerr<<"vor fragmove(...) ...\n";
   return fragmove(CB,li,diffusion_tolerance_min,diffusion_tolerance_steepness,
                   elongation,smoothmove,p_tension,K_elongation,l);
   // Note that p_move is not needed in the parameter list, as it is a property of class cell.
}
void cellCB::show_number_of_divisions(double &time, ofstream &ndivtime) {
   // get average and sd since last writing
   double ndiv_av = 0, ndiv_sd = 0;
   int ndiv_nn = 0;
   for (int i = 0; i <= max_n_of_divisions; i++) {
      ndiv_av += double (attributed_n_of_divisions[i] * i);
      ndiv_nn += attributed_n_of_divisions[i];
   }
   if (ndiv_nn > 0) { ndiv_av /= double (ndiv_nn); } else { ndiv_av = 0; }
   for (int i = 0; i <= max_n_of_divisions; i++) {
      for (int j = 0; j < attributed_n_of_divisions[i]; j++) {
         ndiv_sd += (double (i) - ndiv_av) * (double (i) - ndiv_av);
      }
   }
   if (ndiv_nn > 1) {
      ndiv_sd /= double (ndiv_nn - 1);
      ndiv_sd = sqrt(ndiv_sd);
   } else { ndiv_sd = 0; }

   // write to file
   ndivtime << time << "   " << ndiv_av << "   " << ndiv_sd << "\n";
   // clear the array for the next time interval
   for (int i = 0; i <= max_n_of_divisions; i++) {
      attributed_n_of_divisions[i] = 0;
   }
}
//MSchips
void cellCB::selfBC_number_of_divisions(double &time, ofstream &ndivtime) {
   // get average and sd since last writing
   double ndiv_av = 0, ndiv_sd = 0;
   int ndiv_nn = 0;
   for (int i = 0; i <= max_n_of_divisions; i++) {
      ndiv_av += double (selfBC_attributed_n_of_divisions[i] * i);
      ndiv_nn += selfBC_attributed_n_of_divisions[i];
   }
   if (ndiv_nn > 0) { ndiv_av /= double (ndiv_nn); } else { ndiv_av = 0; }
   for (int i = 0; i <= max_n_of_divisions; i++) {
      for (int j = 0; j < selfBC_attributed_n_of_divisions[i]; j++) {
         ndiv_sd += (double (i) - ndiv_av) * (double (i) - ndiv_av);
      }
   }
   if (ndiv_nn > 1) {
      ndiv_sd /= double (ndiv_nn - 1);
      ndiv_sd = sqrt(ndiv_sd);
   } else { ndiv_sd = 0; }

   // write to file
   ndivtime << time << "   " << ndiv_av << "   " << ndiv_sd << "\n";
   // clear the array for the next time interval
   for (int i = 0; i <= max_n_of_divisions; i++) {
      selfBC_attributed_n_of_divisions[i] = 0;
   }
}
//MS
void cellCB::show_mutation_prob_SELF(double &time,
                                double &muta_avS,
                                double &muta_sdS,
                                ofstream &Smutation_time) {
   // get average and sd since last writing
   double mut_av = 0, mut_sd = 0;
   int mut_nn = 0;
   for (int i = 0; i < mutation_bins; i++) {
      mut_av += double (attributed_mutation_prob_SELF[i]) * double (double (i) + 0.5)
                / double (mutation_bins);
      mut_nn += attributed_mutation_prob_SELF[i];
   }
   if (mut_nn > 0) { mut_av /= double (mut_nn); } else { mut_av = 0; }
   for (int i = 0; i < mutation_bins; i++) {
      double muttmp = 0;
      for (int j = 0; j < attributed_mutation_prob_SELF[i]; j++) {
         muttmp = (double (i) + 0.5) / double (mutation_bins);
         mut_sd += (muttmp - mut_av) * (muttmp - mut_av);
      }
   }
   if (mut_nn > 1) {
      mut_sd /= double (mut_nn - 1);
      mut_sd = sqrt(mut_sd);
   } else { mut_sd = 0; }

   // write to file
   Smutation_time << time << "   "
                 << mut_av << "  " << mut_sd << "     "  // version based on selected BCs since last
                                                         // writing
                 << muta_avS << "  " << muta_sdS << "\n";  // version based on all current CBs (got it
                                                         // from cellman)
   // clear the array for the next time interval
   for (int i = 0; i < mutation_bins; i++) {
      attributed_mutation_prob_SELF[i] = 0;
   }
}


void cellCB::show_cummulative_number_of_divisions() {
   ofstream ndivhisto("ndivhisto.out");
   ndivhisto << "! Histogram for the number of divisions attributed to selected BCs:\n"
             << "! number of divisions . number of occurrences in BCs\n";
   for (int i = 0; i <= max_n_of_divisions; i++) {
      ndivhisto << i << "   " << cummulative_attributed_n_of_divisions[i] << "\n";
   }
   ndivhisto.close();
}
void cellCB::show_mutation_prob(double &time,
                                double &muta_av,
                                double &muta_sd,
                                ofstream &mutation_time) {
   // get average and sd since last writing
   double mut_av = 0, mut_sd = 0;
   int mut_nn = 0;
   for (int i = 0; i < mutation_bins; i++) {
      mut_av += double (attributed_mutation_prob[i]) * double (double (i) + 0.5)
                / double (mutation_bins);
      mut_nn += attributed_mutation_prob[i];
   }
   if (mut_nn > 0) { mut_av /= double (mut_nn); } else { mut_av = 0; }
   for (int i = 0; i < mutation_bins; i++) {
      double muttmp = 0;
      for (int j = 0; j < attributed_mutation_prob[i]; j++) {
         muttmp = (double (i) + 0.5) / double (mutation_bins);
         mut_sd += (muttmp - mut_av) * (muttmp - mut_av);
      }
   }
   if (mut_nn > 1) {
      mut_sd /= double (mut_nn - 1);
      mut_sd = sqrt(mut_sd);
   } else { mut_sd = 0; }

   // write to file
   mutation_time << time << "   "
                 << mut_av << "  " << mut_sd << "     "  // version based on selected BCs since last
                                                         // writing
                 << muta_av << "  " << muta_sd << "\n";  // version based on all current CBs (got it
                                                         // from cellman)
   // clear the array for the next time interval
   for (int i = 0; i < mutation_bins; i++) {
      attributed_mutation_prob[i] = 0;
   }
}
void cellCB::show_cummulative_mutation_prob() {
   ofstream mutset("mutation_set.out");
   mutset << "! Histogram for the probability of mutation attributed to selected BCs:\n"
          << "! bin . mutation probability . number of occurrences in BCs\n";
   for (int i = 0; i < mutation_bins; i++) {
      mutset << i << "   "
             << (double (i) + 0.5) / double (mutation_bins) << "   "
             << cummulative_attributed_mutation_prob[i] << "\n";
   }
   mutset.close();
}
void cellCB::show_cell_cycle_phase_duration() {
   // Writes the histogramms saved in dtphase_frequency[phase][dtphase_resolution]
   // filled according to
   // dt_index=int(phaselength/(2.0*dtphase[phase]/(dtphase_resolution-1.0))+0.5);
   // ++dtphase_frequency[s-cb_G1][dt_index];
   // the mean phase duration is saved in dtphase[phase]
   // and the phases are cb_G1, cb_G0, cb_S, cb_G2, cb_M (in this order)
   ofstream ccp;
   double duration;
   for (short p = cb_G1; p <= cb_M; p++) {
      if (p == cb_G1) { ccp.open("cellcycle_duration_g1.out"); } else if (p == cb_G0) {
         ccp.open("cellcycle_duration_g0.out");
      } else if (p == cb_S) {
         ccp.open("cellcycle_duration_s.out");
      } else if (p == cb_G2) {
         ccp.open("cellcycle_duration_g2.out");
      } else if (p == cb_M) {
         ccp.open("cellcycle_duration_m.out");
      }
      for (int i = 0; i < dtphase_resolution; i++) {
         duration = 2.0 * double (i) * dtphase[p] / double (dtphase_resolution - 1);
         ccp << duration << "   " << dtphase_frequency[p - cb_G1][i] << endl;
      }
      ccp.close();
   }
}
// ============================================================
// ============================================================
// ==================== cellCC ================================
// ============================================================
// ============================================================
// ============================================================

cellCC::cellCC() {
   state = unselected;
   came_from_FDCselected = false;
   affinity = 0.;
   tc_clock = 0.;
   FDCselected_clock = 0.;
   ICOSL = 1;
   hadcontact2Tfh = false;
   mTORC1 = 0;
   FoxO = FoxO_ini;
   FoxO_upregulation = FoxO_production;
   tc_interaction_time = 0.;
   tc_search_duration = 0.;
   individual_dif_delay = 0.;
   bound_ag_index = -1;
   selected_clock = 0.;
   nFDCcontacts = 0;
   BCRsignal = 0.;
   cMyc = 0.;
   nTCcontacts = 0;
   TimeOfLastTbinding = -1;
   collected_ag_portions.clear();
   tc_signal = 0.;
   TFHsignalAUC = 0.;
   SSTactive = false;
   tc_index = -1;
   last_tc_index = -1;
   last_tc_id = -1;
   selected4output = false;
   DEC205 = false;
   DEC205_ova = false;
   MHCdeficient = false;
   fdc_clock = 0.;
   selectable = 1;
   mobile = 1;
   set_p_move();
   //  Ig_class=IgM;
   BCRexpression = 1.0;
   CXCR5failure = 0;
   DND = -1;
   fatetracker_n = 0;

   //MSchips
   //this is 1 until selfBC interacts with TFR
   tfhHelpReduction = 1;
   nTFRcontacts = 0;
   tfr_index = -1;
   last_tfr_index = -1;
   last_tfr_id = -1;
   tfr_tmp_index = -1;
   ccInd_TFRbound = -1;
   ccInd_TFHbound = -1;
   tfr_interaction_time = 0;
   alterednFDCcontacts = 0;
   CD138 = 0;
   tfr_clock = 0;
   tfr_signal = 0;
   came_from_selected=0;
   time_before_tfr=0;
}
cellCC::cellCC(const cellCC &x) : cell(x) {
   operator =(x);
}
cellCC&cellCC::operator =(const cellCC &x) {
   cell::operator =(x);
   state = x.state;
   came_from_FDCselected = x.came_from_FDCselected;
   affinity = x.affinity;
   selectable = x.selectable;
   mobile = x.mobile;
   set_p_move();
   tc_clock = x.tc_clock;
   FDCselected_clock = x.FDCselected_clock;
   ICOSL = x.ICOSL;
   hadcontact2Tfh = x.hadcontact2Tfh;
   mTORC1 = x.mTORC1;
   FoxO = x.FoxO;
   FoxO_upregulation = x.FoxO_upregulation;
   tc_interaction_time = x.tc_interaction_time;
   tc_search_duration = x.tc_search_duration;
   individual_dif_delay = x.individual_dif_delay;
   bound_ag_index = x.bound_ag_index;
   selected_clock = x.selected_clock;
   nFDCcontacts = x.nFDCcontacts;
   BCRsignal = x.BCRsignal;
   cMyc = x.cMyc;
   nTCcontacts = x.nTCcontacts;
   TimeOfLastTbinding = x.TimeOfLastTbinding;
   collected_ag_portions = x.collected_ag_portions;
   tc_signal = x.tc_signal;
   TFHsignalAUC = x.TFHsignalAUC;
   SSTactive = x.SSTactive;
   tc_index = x.tc_index;
   last_tc_index = -1; 
   last_tc_id = -1;
   selected4output = x.selected4output;
   DEC205 = x.DEC205;
   DEC205_ova = x.DEC205_ova;
   MHCdeficient = x.MHCdeficient;
   fdc_clock = x.fdc_clock;
   IgX.set_class(x.IgX);
   BCRexpression = x.BCRexpression;
   CXCR5failure = x.CXCR5failure;
   DND = x.DND;
   fatetracker = x.fatetracker;
   fatetracker_n = x.fatetracker_n;
   //MSchips
   selfMutation = x.selfMutation;
   redeemed = x.redeemed;
   ccl3KO = x.ccl3KO;
   tfhHelpReduction = x.tfhHelpReduction;
   nTFRcontacts = x.nTFRcontacts;
   tfr_index = x.tfr_index;
   last_tfr_index = -1;
   last_tfr_id = -1;
   tfr_tmp_index = x.tfr_tmp_index;
   ccInd_TFRbound = x.ccInd_TFRbound;
   ccInd_TFHbound = x.ccInd_TFHbound;
   tfr_interaction_time = x.tfr_interaction_time;
   alterednFDCcontacts = x.alterednFDCcontacts;
   CD138 = x.CD138;
   tfr_clock = x.tfr_clock;
   tfr_signal = x.tfr_signal;
   came_from_selected = x.came_from_selected;
   time_before_tfr = x.time_before_tfr;
   return *this;
}
cellCC&cellCC::operator =(const cellCB &x) {
   cell::operator =(x);
   state = unselected;
   came_from_FDCselected = false;
   affinity = 0.;
   selectable = 1;
   mobile = 1;
   set_p_move();
   individual_dif_delay = 0.;
   bound_ag_index = -1;
   tc_clock = 0.;
   FDCselected_clock = 0.;
   ICOSL = 1.;
   hadcontact2Tfh = false;
   mTORC1 = 0.;
   FoxO = FoxO_ini;
   FoxO_upregulation = FoxO_production;
   tc_interaction_time = 0.;
   tc_search_duration = 0.;
   selected_clock = 0.;
   nFDCcontacts = int (x.retained_ag + 0.5);
   BCRsignal = 0.;
   cMyc = 0.;
   nTCcontacts = 0;
   TimeOfLastTbinding = -1;
   collected_ag_portions = x.collected_ag_portions;
   tc_signal = 0.;
   TFHsignalAUC = 0.;
   SSTactive = false;
   tc_index = -1;
   last_tc_index = -1;
   last_tc_id = -1;
   selected4output = false;
   DEC205 = x.DEC205;
   DEC205_ova = x.DEC205_ova;
   MHCdeficient = x.MHCdeficient;
   fdc_clock = 0.;
   IgX.set_class(x.IgX);
   BCRexpression = 1.0;
   if (IgX.Ig_class == IgE) { BCRexpression = IgE_BCRlevel; }
   CXCR5failure = 0;
   DND = -1;
   fatetracker_n = 0;
   //MSchips
   selfMutation = x.selfMutation;
   redeemed = x.redeemed;
   ccl3KO = x.ccl3KO;
   nTFRcontacts = 0;
   tfr_index = -1;
   last_tfr_index = -1;
   last_tfr_id = -1;
   tfr_tmp_index=-1;
   ccInd_TFRbound = -1;
   ccInd_TFHbound = -1;
   tfr_interaction_time = 0;
   CD138 = 0;
   tfr_clock = 0;
   tfr_signal = 0;
   came_from_selected=0;
   time_before_tfr=0;
   return *this;
}
cellCC::~cellCC() {
   // cout<<"in ~cellCC() ...\n";
}
// ============================================================

double cellCC::p_apo = 0.;
double cellCC::p_apo4FDCselected = 0.;
double cellCC::p_mph = 0.;
double cellCC::p_FDCdetach = 0.;
double cellCC::TCselect_prob = 1.;
double cellCC::p_FDCsignalling = 1.;
double cellCC::p_dif = 0.;
double cellCC::dif_delay = -1.;
double cellCC::dif_delay_DEC = -1.;
double cellCC::p_difu = 0.;
double cellCC::p_difu_width = -1.;
double cellCC::p_CXCR5down = 0.;
double cellCC::p_dif2out = 0.;
double cellCC::p_dif2out_DEC = 0.;
double cellCC::p_dif2out_target = 0.;
double cellCC::p_dif2out_DEC_target = 0.;
double cellCC::p_final_differentiation = -1.;
double cellCC::start_differentiate = 120.;
long cellCC::test_delay = -1;
long cellCC::ICAM_delay = -1;
double cellCC::persistence = 0.;
short cellCC::v_modi = 1;
short cellCC::n_v_states = 1;
double cellCC::v_slow_factor = 1.0;
double cellCC::p_switch_v = 0.0;
short cellCC::mode_of_setting_tc_time = 0;
double cellCC::tc_time = 1.5;
double cellCC::tc_time_width = 0.5;
double cellCC::BTtime_K = 1.0;
double cellCC::BTtime_min = 0.0;
double cellCC::BTtime_max = 2.0;
double cellCC::BTtime_n = 2.0;
short cellCC::TFHsignal_delivery_mode = 0;
double cellCC::TFHsignal_delivery_min = 0;
double cellCC::TFHsignal_delivery_max = 2;
double cellCC::TFHsignal_delivery_KpMHC = -1;
double cellCC::TFHsignal_delivery_n = 2;
double cellCC::TFHsignal_delivery_factor[TFHsignal_delivery_max_pMHC] = { 0.0 };
bool cellCC::TFHadapt2pMHCseen = 0;
double cellCC::TFHsignal_decay = -1;
double cellCC::cMyc_decay_factor = -1;
bool cellCC::cMyc_FoxO_break = false;
double cellCC::cMyc_FoxO_break_K = 0.5;
double cellCC::cMyc_FoxO_break_n = 1;
double cellCC::rescue_signal_threshold = 1.0;
double cellCC::AgThreshold4Selection = -1.0;
short cellCC::tc_search_duration_mode = 0;
short cellCC::tc_search_time_determination_mode = 0;
bool cellCC::FCRB = false;
double cellCC::tc_search_duration_per_FDCcontact = 0.75;
double cellCC::tc_search_duration_fixed = 3.0;
double cellCC::tc_dec205ova_binding_time = 0.0;
double cellCC::max_tc_signal_portion = -1.0; // defined in constructor
double cellCC::SST_tc_signal = 0.;
unsigned short cellCC::CC_FDC_selection = 1;
short cellCC::TC_CC_selection = 0;
bool cellCC::negativeTCselection = true;
bool cellCC::BCstaysonTCbyTCtime = false;
bool cellCC::force_selection_at_TFHthreshold = false;
int cellCC::fdc_encounters[fdc_max_encounters];
bool cellCC::SMOOTH_DIFFERENTIATION = false;
double cellCC::smooth_differentiation_time = 3.0;
bool cellCC::collectFDCsignals = false;
bool cellCC::multipleTFHcontacts = 0;
int cellCC::AgThreshold4Diff = 0;
double cellCC::collectFDCperiod = 0.;
bool cellCC::inhibitFDCapoptosis = false;
bool cellCC::simultaneousTfhFDC = false;
double cellCC::prob2bindFDCfirst = 0.5;
double cellCC::prob2kill_noFDCcontactBCs = 1.0; 
int cellCC::reset_antigen_after_collection = -1;
double cellCC::ignore_affinity = -1.0;
bool cellCC::ag_loaded_CC_directly2TFH = false;
bool cellCC::ag_deleted_in_fresh_CC = true;
double cellCC::IgE_BCRlevel = 0.3;
double cellCC::IgE_prob_CXCR5down = 0.;
short cellCC::apoptotic_motility_mode = 1;
double cellCC::p_apo_randomwalk = 0.0;
short cellCC::mode_of_apoptosis = 0;
bool cellCC::pMHC_dependent_division = false;
bool cellCC::SIND = false;
bool cellCC::use_cMyc4DND = false;
double cellCC::DND_P_standard = 2.0;
double cellCC::DND_P_min = 1.0;
double cellCC::DND_P_max = 6.0;
double cellCC::pMHC_dependent_K = 8.0;
double cellCC::cMyc_dependent_K = 1.5;
double cellCC::TFHsignal_dependent_K = 1.5; 
double cellCC::TFHgradient_dependent_K = 1.5; 
double cellCC::DND_nHill = 1.0;
double cellCC::pMHC_dependent_pMHC_of_2divisions = 2.0;
short cellCC::mode_of_DND_of_Tfh = 0;
double cellCC::TFHsignal_of_P0divisions = 0.75;
double cellCC::TFHgradient_of_P0divisions = 0.75;
double cellCC::pMHC_of_DEC205_ova = 0;
bool cellCC::ICOSL_dependent_Tfh_signals = false;
bool cellCC::ICOSL_memory = false;
short cellCC::ICOSL_upregulation_mode = 0;
double cellCC::ICOSL_upregulation_time = 0.;
short cellCC::FoxO_mode = 0;
short cellCC::dT_FoxO_start = 0;
short cellCC::dT_FoxO_reg = 3;
bool cellCC::stopFoxOonTFH = false;
double cellCC::FoxO_ini = 1.;
double cellCC::FoxO_production = 1.0;
double cellCC::FoxOup_min = 0.5;
double cellCC::FoxOup_max = 12.0;
double cellCC::FoxOup_K = 2.0;
double cellCC::FoxOup_n = 2.0;
double cellCC::KFoxO = 5.0;
double cellCC::nFoxO = 1.0;
short cellCC::mTOR_mode = 0;
double cellCC::mTORC1_production = 0.;
double cellCC::mTOR_Tfh_production = 0.;
double cellCC::mTOR_BCR_K = 20.;
double cellCC::mTOR_BCR_n = 1.;
long cellCC::cummulative_ag_collection_all[max_n_of_ag_portions + 1] = { 0 };
long cellCC::cummulative_ag_collection_selected[max_n_of_ag_portions + 1] = { 0 };
short cellCC::present_specific_ag2TC = 0;
short cellCC::outputfiles = 1;
double cellCC::write_trackfate = 0;
bool cellCC::no_rebinding_to_moved_TCs = false;

//MSchips
short cellCC::TFR_CC_interaction_mode = -1;
double cellCC::TfhHelp_reduction_factor = 0;
double cellCC::time_firstTFRinteraction = 0;
short cellCC::TFR_bind_onlyself = 0;
double cellCC::tfr_time = 0;
bool cellCC::exp_CCL3KO = 0;
long cellCC::num_tfr_contacts_self[ntfr + 1] = {0};
long cellCC::num_tfr_contacts_nonself[ntfr + 1] = {0};
//
histo_monitor cellCC::monitor_pureFDCsearch;
histo_monitor cellCC::monitor_nTCcontacts_selected;
histo_monitor cellCC::monitor_nTCcontacts_deleted;
histo_monitor cellCC::monitor_BsearchTtime_selected;
histo_monitor cellCC::monitor_BsearchTtime_deleted;
histo_monitor cellCC::monitor_BinteractTtime;
histo_monitor cellCC::monitor_TimeBetweenTbindings;
histo_monitor cellCC::monitor_TfhSignalAtSelection;
histo_monitor cellCC::monitor_TfhSignalAtSelection_selected;
histo_monitor cellCC::monitor_TfhSignalAtSelection_deleted;
histo_monitor cellCC::monitor_pMHCatFateDecision_selected;
histo_monitor cellCC::monitor_pMHCatFateDecision_deleted;
histo_monitor cellCC::monitor_TfhSignalSpeed_selected;
histo_monitor cellCC::monitor_TfhSignalSpeed_deleted;
histo_monitor cellCC::monitor_TFHsignalAUC;
histo_monitor cellCC::monitor_TFHsignalAUC_selected;
histo_monitor cellCC::monitor_TFHsignalAUC_deleted;
histo_monitor cellCC::monitor_cMycAtSelection;
histo_monitor cellCC::monitor_cMycAtSelection_selected;
histo_monitor cellCC::monitor_cMycAtSelection_deleted;
histo_monitor cellCC::monitor_mTORatSelection;
histo_monitor cellCC::monitor_mTORatSelection_selected;
histo_monitor cellCC::monitor_mTORatSelection_deleted;
histo_monitor cellCC::monitor_TFHintensity;
//MSchips
histo_monitor cellCC::monitor_nTFRcontacts_CCselected;
histo_monitor cellCC::monitor_nTFRcontacts_CCdeleted;
histo_monitor cellCC::monitor_nTFRcontacts_SelfCCselected;
histo_monitor cellCC::monitor_nTFRcontacts_SelfCCdeleted;
histo_monitor cellCC::monitor_TFHsig_postTFR_CCselected;
histo_monitor cellCC::monitor_TFHsig_postTFR_SelfCCselected;
histo_monitor cellCC::monitor_BinteractTFRtime_self;
histo_monitor cellCC::monitor_BinteractTFRtime_nonself;
histo_monitor cellCC::monitor_TFRsig_select;
histo_monitor cellCC::monitor_TFRsig_apopt;
histo_monitor cellCC::monitor_nMut_nonSelfCD138;
histo_monitor cellCC::monitor_nMut_SelfCD138;
histo_monitor cellCC::monitor_nMut_RedeemedCD138;
histo_monitor cellCC::monitor_time_first_TFRcontact_Selfselected;
histo_monitor cellCC::monitor_time_first_TFRcontact_NonSelfselected;
//
int cellCC::fatetracker_ndt = 8; // number of timesteps between two write events

void cellCC::set_statics(const Parameter &par, space &l, ofstream &ana) {
  //
  // For centrocytes
  mode_of_apoptosis = par.Value.mode_of_apoptosis;
  if (mode_of_apoptosis == 0) {
    p_apo = par.Value.apoptosis * par.Value.deltat;
    ana << "  CC apoptosis = " << p_apo << "\n";
    if (par.Value.apoptosis4FDCselected > 0) {
      p_apo4FDCselected = par.Value.apoptosis4FDCselected * par.Value.deltat;
      ana << "  CC apoptosis after FDC selection = " << p_apo4FDCselected << "\n";
    } else { ana << "  no CC apoptosis after FDC selection.\n"; }
  } else { // if mode_of_apoptosis == 1
    // Just put the life-times in both variables:
    p_apo = par.Value.apoptosis;
    ana << "  CC life time = " << p_apo << " h\n";
    if (par.Value.apoptosis4FDCselected > 0) {
      p_apo4FDCselected = par.Value.apoptosis4FDCselected;
      ana << "  CC life time after FDC selection = " << p_apo4FDCselected << " h\n";
    }
  }
  p_FDCdetach = par.Value.selection * par.Value.deltat;
  ana << "  CC selection = " << p_FDCdetach << "\n";
  p_FDCsignalling = par.Value.FDCsignalling;
  ana << "  Prob of FDC signalling to CC in contact = " << p_FDCsignalling << "\n";
  IgE_BCRlevel = par.Value.IgE_BCRlevel;
  if (IgE_BCRlevel
      < 1.0) { ana << "  BCR expression of IgE BC reduced to " << IgE_BCRlevel << "\n"; }
  IgE_prob_CXCR5down = par.Value.CC_IgE_prob_CXCR5down;
  ana << "  Prob of failure of IgE-BC to upregulate CXCR5 upon CB2CC differentiation "
      << IgE_prob_CXCR5down << "\n";
  
  ag_loaded_CC_directly2TFH = par.Value.ag_loaded_CC_directly2TFH;
  CC_FDC_selection = par.Value.CC_FDC_selection;
  collectFDCsignals = par.Value.collectFDCsignals;
  multipleTFHcontacts = par.Value.multipleTFHcontacts;
  present_specific_ag2TC = par.Value.present_specific_ag2TC;
  if (present_specific_ag2TC==2) {
      cerr<<present_specific_ag2TC<<endl;exit(1);
  }
  reset_antigen_after_collection = par.Value.reset_antigen_after_collection;
  ignore_affinity = par.Value.ignore_affinity;
  ag_deleted_in_fresh_CC = par.Value.ag_deleted_in_fresh_CC;
  collectFDCperiod = par.Value.collectFDCperiod;
  AgThreshold4Diff = par.Value.AgThreshold4Diff;
  prob2kill_noFDCcontactBCs = par.Value.prob2kill_noFDCcontactBCs;
  inhibitFDCapoptosis = false;
  if (prob2kill_noFDCcontactBCs < 1 && prob2kill_noFDCcontactBCs >= 0) {
    inhibitFDCapoptosis = true;
  }

  // cout<<"in cellthis: collect-dt="<<collectFDCperiod<<"\n\n\n";
  AgThreshold4Selection = par.Value.AgThreshold4Selection;
  tc_search_duration_mode = par.Value.tc_search_duration_mode;
  tc_search_time_determination_mode = par.Value.tc_search_time_determination_mode;
  if (tc_search_duration_mode == 5) { FCRB = true; }
  no_rebinding_to_moved_TCs = par.Value.no_rebinding_to_moved_TCs;
  tc_search_duration_per_FDCcontact = par.Value.tc_search_duration_per_FDCcontact;
  tc_search_duration_fixed = par.Value.tc_search_duration_fixed;
  if (tc_search_duration_mode == 0) { tc_search_duration_fixed = 0.; }
  TC_CC_selection = par.Value.TC_CC_selection;
  negativeTCselection = par.Value.negativeTCselection;
  BCstaysonTCbyTCtime = par.Value.BCstaysonTCbyTCtime;
  force_selection_at_TFHthreshold = par.Value.force_selection_at_TFHthreshold;
  tc_time = par.Value.TC_time;
  tc_time_width = par.Value.TC_time_width;
  BTtime_K = par.Value.BTtime_K;
  BTtime_min = par.Value.BTtime_min;
  BTtime_max = par.Value.BTtime_max;
  BTtime_n = par.Value.BTtime_n;
  //
  TFHsignal_delivery_mode = par.Value.TFHsignal_delivery_mode;
  TFHsignal_delivery_min = par.Value.TFHsignal_delivery_min;
  TFHsignal_delivery_max = par.Value.TFHsignal_delivery_max;
  TFHsignal_delivery_KpMHC = par.Value.TFHsignal_delivery_KpMHC;
  TFHsignal_delivery_n = par.Value.TFHsignal_delivery_n;
  if (TFHsignal_delivery_mode == 1 && TFHsignal_delivery_KpMHC < 0) {
    /* Determine the K-value from the condition that 
     * X0 is the amount of pMHC at which Hill(X0) = 1,
     * yielding   K = X0 [\frac{Pmax - Pmin}{1 - Pmin} - 1]^(1/n)
     */
    TFHsignal_delivery_KpMHC 
      = get_Hill_K(par.Value.TFHsignal_delivery_pMHCof1,
		   TFHsignal_delivery_min, TFHsignal_delivery_max,
		   TFHsignal_delivery_n, 1.0);
    /*      = par.Value.TFHsignal_delivery_pMHCof1
     *      * pow( ( (TFHsignal_delivery_max - 1.0) / (1.0 - TFHsignal_delivery_min) ), 
     *      (1.0 / TFHsignal_delivery_n) );
     */
    ana << "  ... calculated KpMHC for TFH signal delivery to BCs to " 
	<< TFHsignal_delivery_KpMHC << "\n";
  }
  // Load an array with integer pMHC values:
  for (int i = 0; i < TFHsignal_delivery_max_pMHC; i++) {
    TFHsignal_delivery_factor[i] = get_TFHsignal_delivery_factor(double(i), true);
  }
  TFHadapt2pMHCseen = par.Value.TFHadapt2pMHCseen;
  // Default value of TFHsignal_decay is -1.
  if (par.Value.TFHsignal_decay > 0) {
    TFHsignal_decay = get_decay_factor(par.Value.TFHsignal_decay, par.Value.deltat);
    ana << "Decay of tc_signal with factor " << TFHsignal_decay << " per time step\n"
	<< "corresponding to a half life of " << par.Value.TFHsignal_decay << " hours.\n";
  } else {
    ana << "No decay of tc_signal.\n";
  }
  //
  // Same decay for cMyc
  if (par.Value.cMyc_halflife > 0) {
    cMyc_decay_factor = get_decay_factor(par.Value.cMyc_halflife, par.Value.deltat);
    ana << "cMyc decays with factor " << cMyc_decay_factor
	<< " per time step, i.e. with a half life of "
	<< par.Value.cMyc_halflife << " hours.\n";
  } else {
    ana << "No decay of cMyc.\n";
  }
  cMyc_FoxO_break = par.Value.cMyc_FoxO_break;
  cMyc_FoxO_break_K = par.Value.cMyc_FoxO_break_K;
  cMyc_FoxO_break_n = par.Value.cMyc_FoxO_break_n;
  //
  mode_of_setting_tc_time = par.Value.mode_of_setting_TC_time;
  tc_dec205ova_binding_time = tc_time;
  if (par.Value.TC_dec205ova_time > 1.e-08) {
    tc_dec205ova_binding_time = par.Value.TC_dec205ova_time;
  }
  /* The max portion of TC signals transferred to CCs is just the time step.
     Up to v20160927 simply the duration of CC-TC interaction was used, thus,
     the signal portion was the time step. This can still be used if one
     thinks of a rate of an intracellular signal, which is then rescaled to 
     1 per timestep, which would then also lead to the following definition.
     Note that any suboptimal TC property might reduce this max portion. 
     Signal portions are handled in cellCC::add_tc_signal(...);
     #####
     However, the rescue_signal_threshold might be set differently
     depending on the kind of signal used.
  */
  rescue_signal_threshold = par.Value.TC_rescue_time;
  if (FCRB) { rescue_signal_threshold = par.Value.cMyc_selection_threshold; }
  max_tc_signal_portion = par.Value.deltat;
  SST_tc_signal = par.Value.SST_tc_signal;
  if (TC_CC_selection == 0) {
    ana << "Use selection of CC in interaction with FDC only!\n";
  } else { 
    ana << "Use selection of CC in interaction with FDC and TC!\n";
    if (negativeTCselection) {
      ana << "Activate negative selection of CC by TC.\n";
    } else { ana << "No negative selection of CC by TC.\n"; }
    if (BCstaysonTCbyTCtime) {
      ana << "BC stay in contact to TC when the received above threshold signals.\n";
    } else {
      ana << "BC release contact to TC right when the threshold signalling is reached.\n";
    }
  }
  if (CC_FDC_selection == 1) {
    ana << "CC have to see FDC for selection!\n";
  } else { ana << "CC may be selected without interaction with FDC!\n"; }
  TCselect_prob = par.Value.TCell;
  if (TCselect_prob < 1) {
    ana << "TFH positively selected CC survive with probability " << TCselect_prob << "!\n";
  } else { ana << "All TFH positively selected CC survive.\n"; }
  simultaneousTfhFDC = par.Value.simultaneousTfhFDC;
  if (simultaneousTfhFDC) {
    ana << "BCs continue collecting antigen from FDC during search for Tfh signals.\n";
  } else {
    ana << "BCs stop collecting antigen from FDC during search for Tfh signals.\n";
  }
  prob2bindFDCfirst = 1. - par.Value.prob2bindTFHfirst;
  
  if (((1. / par.Value.ccdiff) > 0.5) && (par.Value.ccdiff_delay > 0.)) {
    cout << "WARNING: Delay of CC differentiation and rate of differentiation are both on.\n"
	 << "         Delay of " << par.Value.ccdiff_delay << " and rate of 1/" << 1.
      / par.Value.ccdiff << " hours.\n"
	 << "         If delay is on, differentiation by rate shoud be fast.\n\n";
    ana << "WARNING: Delay of CC differentiation and rate of differentiation are both on.\n"
	<< "         Delay of " << par.Value.ccdiff_delay << " and rate of 1/" << 1.
      / par.Value.ccdiff << " hours.\n"
	<< "         If delay is on, differentiation by rate shoud be fast.\n\n";
  }
  p_dif = par.Value.ccdiff * par.Value.deltat;
  ana << "  CC differentiation = " << p_dif << "\n";
  p_final_differentiation = par.Value.final_differentiation_rate * par.Value.deltat;
  ana << "  CC final differentiation to output = " << p_final_differentiation << "\n";
  dif_delay = par.Value.ccdiff_delay;
  dif_delay_DEC = par.Value.ccdiff_delay_DEC;
  //  p_dif_DEC=par.Value.ccdiff_DEC*par.Value.deltat;
  //  if (p_dif_DEC<0.) p_dif_DEC=p_dif;
  if (dif_delay < 0.) {
    ana << "  No delay of CC2CB differentation upon selection\n";
  } else {
    ana << "  CC differentiation delayed after selection by = " << dif_delay << "\n";
  }
  if (dif_delay_DEC < 0.) {
    ana << "  No delay of CC2CB differentation upon with DEC205-ova bound\n";
  } else {
    ana << "  CC differentiation delayed after selection of CC with DEC205-ova bound by = "
	<< dif_delay_DEC << "\n";
  }
  if (par.Value.CC_test_delay < 0) {
    test_delay = -1;
  } else { test_delay = long (par.Value.CC_test_delay / par.Value.deltat + 0.5); }
  ana << "  Gap between CC affinity tests = " << test_delay << "\n";
  if (par.Value.CC_ICAM_delay < 0) {
    ICAM_delay = -1;
  } else { ICAM_delay = long (par.Value.CC_ICAM_delay / par.Value.deltat + 0.5); }
  ana << "  CC immobile after affinity tests for = " << ICAM_delay << "\n";
  start_differentiate = par.Value.StartOutput;
  p_dif2out_target = par.Value.output;
  SMOOTH_DIFFERENTIATION = par.Value.smooth_dif2out;
  smooth_differentiation_time = par.Value.smooth_dif2out_time;
  if (SMOOTH_DIFFERENTIATION) {
    ana << "Smooth onset of CC differentiation to output (width = "
	<< smooth_differentiation_time << " hours).\n";
  } else { ana << "Abrupt onset of CC differentiation to output.\n"; }
  p_dif2out_DEC_target = par.Value.output_DEC;
  
  if (par.Value.CC_persistence > 60. * par.Value.deltat) {
    persistence = 60. * par.Value.deltat / par.Value.CC_persistence;
  } else {
    persistence = 1.;   // i.e. change polarity in every time step!
  }
  ana << "  CC polarity persistence = " << par.Value.CC_persistence
      << " i.e. p_{change-polarity} = " << persistence << "\n";
  
  if (par.Value.D_CC >= 0.) {
    p_difu = double (l.dim2) * par.Value.D_CC * par.Value.deltat
      / (par.Value.dx * par.Value.dx);
  } else {
    p_difu = 60. * par.Value.v_CC * par.Value.deltat / par.Value.dx;
  }
  if (p_difu > 0.5) {
    cout << "Centrocyte-Motility = " << p_difu << " to large for dx and dt !!!\n";
    exit(1);
    // Siehe Kommentar bei Centroblasten!
  }
  p_difu_width = par.Value.v_CC_width * p_difu;
  ana << "  CC movement probability = " << p_difu << "\n";
  
  p_CXCR5down = par.Value.CXCR5down * par.Value.deltat;
  if (p_CXCR5down < 0) { p_CXCR5down = 0.; }
  ana << "  CC CXCR5 down regulation probability = " << p_CXCR5down << "\n";
  apoptotic_motility_mode = par.Value.CC_apoptotic_motility_mode;
  ana << "  apoptotic CC move in mode " << apoptotic_motility_mode
      << " (0=no; 1=CXCL13; 2=random; 3=CXCL12)\n";
  if ((apoptotic_motility_mode == 3) && (par.Value.p_apo_randomwalk > 0.)) {
    ana << "  CXCL12 sensitivity of apoptotic CC is replaced by random walk with half life "
	<< log(2.) / par.Value.p_apo_randomwalk << " hours.\n";
  }
  p_apo_randomwalk = par.Value.p_apo_randomwalk * par.Value.deltat;
  
  double tmpdt = par.Value.CC_persistence;
  if (tmpdt < 60. * par.Value.deltat) {
    tmpdt = 60. * par.Value.deltat;
    ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
  }
  if (par.Value.D_CC < 0.) {
    ana << "  Expected effective diffusion constant for CC = "
	<< 60. * par.Value.v_CC * par.Value.v_CC * tmpdt / double (l.dim2)
	<< " microns^2/hr = "
	<< par.Value.v_CC * par.Value.v_CC * tmpdt / double (l.dim2)
	<< " microns^2/min\n";
  } else {
    ana << "  Expected mean cell velocity for CC = "
	<< sqrt(double (l.dim2) * par.Value.D_CC / (60. * tmpdt))
	<< " microns/min\n";
  }
  v_modi = par.Value.CC_v_modi;
  n_v_states = par.Value.CC_n_v_states;
  if (n_v_states > 1) {
    v_slow_factor = par.Value.v_CC_factor;
    p_switch_v = 60. * par.Value.deltat / par.Value.v_CC_switch_deltat;
    if (p_switch_v < 0.) { p_switch_v = 0.; }
  }
  
  p_mph = par.Value.macrophage * par.Value.deltat;
  ana << "  CC macrophage deletion = " << p_mph << "\n";
  
  for (int a = 0; a < fdc_max_encounters; a++) {
    fdc_encounters[a] = 0;
  }
  
  // set variables associated with tc induced number of divisions
  pMHC_dependent_division = par.Value.pMHC_dependent_division;
  SIND = par.Value.signal_dependent_number_of_divisions;
  use_cMyc4DND = par.Value.use_cMyc4DND;
  mode_of_DND_of_Tfh = par.Value.mode_of_DND_of_Tfh;
  if (TC_CC_selection == 0) {
    if (pMHC_dependent_division || SIND) {
      ana << "WARNING: pMHC/signal_dependent_division is switched off: "
	  << "CC selection independent of TFH was chosen.\n";
      cerr << "WARNING: pMHC/signal_dependent_division is switched off: "
	   << "CC selection independent of TFH was chosen.\n";
      pMHC_dependent_division = false;
      SIND = false;
    }
  }
  if (pMHC_dependent_division || SIND || use_cMyc4DND) {
    if (SIND) { ana << "Tfh-signal-dependent number of divisions is set on:\n"; }
    else if (use_cMyc4DND) { ana << "cMyc-dependent number of divisions is set on:\n"; }
    else { ana << "pMHC-dependent number of divisions is set on:\n"; }
    DND_P_standard = par.Value.pMHC_dependent_P_standard;
    DND_P_min = par.Value.pMHC_dependent_P_min;
    DND_P_max = par.Value.pMHC_dependent_P_max;
    DND_nHill = par.Value.pMHC_dependent_nHill;
    pMHC_dependent_K = par.Value.pMHC_dependent_K;
    if (pMHC_dependent_K < 0) {
      // Determine the K-value from the condition that Hill(pMHC_0) = P0 (= P_standard)
      pMHC_dependent_K = get_Hill_K(par.Value.pMHC_dependent_pMHC_of_2divisions,
				    DND_P_min, DND_P_max, DND_nHill, DND_P_standard);
      ana << "  ... calculated K_pMHC of the Hill-function to " << pMHC_dependent_K << "\n";
    } else {
      ana << "  ... ignored P_0 and A_0 and used K=" << pMHC_dependent_K
	  << " from parameter file.\n";
    }
    cMyc_dependent_K = par.Value.cMyc_dependent_K;
    if (cMyc_dependent_K < 0) {
      /* Determine the K-value from the condition Hill(cMyc_0) = P0 (= P_standard) .
       * Note that all variables in the Hill-function starting with pMHC
       * are based on division numbers and are thus valid for pMHC as for TFHsignal.
       */
      cMyc_dependent_K = get_Hill_K(par.Value.cMyc_of_P0divisions, 
				    DND_P_min, DND_P_max, DND_nHill, DND_P_standard);
      ana << "  ... calculated K_cMyc of the Hill-function to " << cMyc_dependent_K << "\n";
    } else {
      ana << "  ... ignored P_0 and A_0 and used K=" << cMyc_dependent_K
	  << " from parameter file.\n";
    }
    TFHsignal_dependent_K = par.Value.TFHsignal_dependent_K;
    if (TFHsignal_dependent_K < 0) {
      // Determine the K-value from the condition analogue to cMyc above.
      TFHsignal_dependent_K = get_Hill_K(par.Value.TFHsignal_of_P0divisions,
					 DND_P_min, DND_P_max, DND_nHill, DND_P_standard);
      ana << "  ... calculated K_TFHsignal of the Hill-function to " 
	  << TFHsignal_dependent_K << "\n";
    } else {
      ana << "  ... ignored P_0 and A_0 and used K=" << TFHsignal_dependent_K
	  << " from parameter file.\n";
    }
    TFHgradient_dependent_K = par.Value.TFHgradient_dependent_K;
    if (TFHgradient_dependent_K < 0) {
      // Determine the K-value from the condition analogue to cMyc above.
      TFHgradient_dependent_K = get_Hill_K(par.Value.TFHgradient_of_P0divisions,
					   DND_P_min, DND_P_max, DND_nHill, DND_P_standard);
      ana << "  ... calculated K_signal_gradient of the Hill-function to " 
	  << TFHgradient_dependent_K << "\n";
    } else {
      ana << "  ... ignored P_0 and A_0 and used K=" << TFHgradient_dependent_K
	  << " from parameter file.\n";
    }
    /* All pMHC_dependent parameters are now set or calculated.
     * For DEC205_ova, the amount of pMHC attributed to DEC205 positive BCs has to be set.
     * It is defined as DEC205_ova_pMHC_factor times the pMHC corresponding to 
     * meandiv = DND_P_standard divisions.
     * This is calculated from the condition Hill(pMHC_0) = 2 = meandiv, yielding
     * pMHC_0 = K {(meandiv-P_min)/(P_max-meandiv)}^(1/n)
     * which is the global function get_Hill_K(...) with min and max values exchanged.
     */
    // ++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    double DEC205_ova_pMHC_factor = 5.0; // factor read-off from Fig S5E in Victora Cell 2010
    double meandiv = DND_P_standard; // A value of 2 is given in Gitlin 2014
    // ++++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++	 
    if (meandiv <= DND_P_min || meandiv >= DND_P_max) {
      cerr << "ERROR: meandiv = " << meandiv << " in cellCC::set_statics(...) "
	   << "is required to be in the interval ]" << DND_P_min << ","
	   << DND_P_max << "[. Exit simulation.\n";
      exit(1);
    }
    /* The effect of anti-DEC205-OVA is assumed to increase pMHC presentation,
     * thus, only an increased pMHC is calculated, nothing to do for TFH-signals.
     * However, in the case of Tfh-signal-depdendent division, the parameter
     * pMHC_dependent_K is calculated from an otherwise unused parameter.
     * Eventually, <pMHC-dependent division number Hill: pMHC for 2 divisions>
     * has to be set to a different value, in particular when simultaneous
     * search for Tfh and FDC is set on.
     */
    pMHC_of_DEC205_ova = get_Hill_K(DEC205_ova_pMHC_factor * pMHC_dependent_K,
				    DND_P_max, DND_P_min, DND_nHill, meandiv);
    ana << "pMHC amount for DEC205_ova positive cells is set to " 
	<< pMHC_of_DEC205_ova << "\n";
  } else { 
    ana << "pMHC-dependent number of divisions is set off.\n"; 
  }
  //
  //
  // Histograms based on the class histo_monitor:
  double max_value;
  /* Initialize the counters for the frequency of periods of pure antigen search on FDC */
  max_value = 4. * par.Value.collectFDCperiod;
//  max_value = 100. * par.Value.collectFDCperiod;
  
  monitor_pureFDCsearch.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
				   "pureFDCsearchDuration.out",
				   "pureFDCsearchDuration");
  /* Initialise the counters for the frequency of number of CC-TC-contacts:
   * Monitors the number directly, thus, max_value = dimension of array.
   */
  max_value = 100; 
  monitor_nTCcontacts_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					  "nTCcontacts_selected.out", 
					  "nTCcontacts_selected");
  monitor_nTCcontacts_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					  "nTCcontacts_deleted.out", 
					  "nTCcontacts_deleted");
  /* Initialise the counters for the frequency of B search Tfh times:
   * Eventually introduce a factor on the RHS if using DEC205-OVA experiments
   * The factor 1.5 in mode 2 is just for security reasons to avoid error messages on the 
   * screen. This is in particular important since the model was extended to make the search
   * duration depend on FoxO levels, which can generate longer durations.
   */
  max_value = 72.0; // used for tc_search_duration_mode == 0
  if (tc_search_duration_mode == 1) {
    max_value = 2.0 * tc_search_duration_fixed;
  } else if (tc_search_duration_mode == 2) {
    max_value = 1.5 * pMHC_of_DEC205_ova * tc_search_duration_per_FDCcontact;
  } else if (tc_search_duration_mode == 3) {
    max_value = 15.0 * tc_search_duration_fixed;
  }
  monitor_BsearchTtime_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					   "BsearchTtime_selected.out", 
					   "BsearchTtime_selected");
  monitor_BsearchTtime_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					  "BsearchTtime_deleted.out", 
					  "BsearchTtime_deleted");
  /* Initialise the counters for the frequency of B-Tfh interaction times:
   */
  max_value = 3.0; // used for tc_search_duration_mode == 0
  monitor_BinteractTtime.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
				    "BinteractTtime.out", "BinteractTtime");
  /* Initialise the counters for the times between two T-B-bindings
   */
  max_value = 20.0;
  monitor_TimeBetweenTbindings.initialize(max_value, 400, par.Value.tmin, par.Value.tmax,
	 "TimeBetweenTbindings.out", "TimeBetweenTbindings");
  /* Initialise the counters for the frequency of integrated Tfh signal levels:
   * Eventually introduce a factor on the RHS if using DEC205-OVA experiments
   * The factor 1.5 in mode 2 below is just for security reasons to avoid error messages 
   * on the screen. This is in particular important since the model was extended to make 
   * the search duration depend on FoxO levels, which can generate longer durations.
   */
  max_value = 10.0; // used for tc_search_duration_mode == 0
  if (tc_search_duration_mode == 1) {
    max_value = 2.0 * tc_search_duration_fixed;
  } else if (tc_search_duration_mode == 2) {
    max_value = 1.5 * pMHC_of_DEC205_ova * tc_search_duration_per_FDCcontact;
  } else if (tc_search_duration_mode == 3) {
    max_value = 15.0 * tc_search_duration_fixed;
  }
  if (force_selection_at_TFHthreshold) {
    max_value = 5. * rescue_signal_threshold;
  }
  monitor_TfhSignalAtSelection.initialize(max_value, 1000, par.Value.tmin, par.Value.tmax,
					  "TfhSignalAtSelection.out", 
					  "TfhSignalAtSelection");
  monitor_TfhSignalAtSelection_selected.initialize(max_value, 1000, 
                                                   par.Value.tmin, par.Value.tmax,
                                                   "TfhSignalAtSelection_selected.out", 
                                                   "TfhSignalAtSelection_selected");
  max_value = 1.5 * rescue_signal_threshold;
  monitor_TfhSignalAtSelection_deleted.initialize(max_value, 100, 
                                                  par.Value.tmin, par.Value.tmax,
                                                  "TfhSignalAtSelection_deleted.out", 
                                                  "TfhSignalAtSelection_deleted");
  max_value = 500;
  monitor_pMHCatFateDecision_selected.initialize(max_value, max_value, 
                                                 par.Value.tmin, par.Value.tmax,
                                                 "pMHCatFateDecision_selected.out", 
                                                 "pMHCatFateDecision_selected");
  monitor_pMHCatFateDecision_deleted.initialize(max_value, max_value, 
                                                par.Value.tmin, par.Value.tmax,
                                                "pMHCatFateDecision_deleted.out", 
                                                "pMHCatFateDecision_deleted");
  /* Initialise the counters for the frequency of B Tfh-signal-speed:
   * Here the range of possible values has to be estimated:
   * The time the BC spends to search and interact with Tfh is always longer
   * than the received signal. Thus, the ratio signal/FDCselected_clock < 1.
   * Use this for now; alternatively,
   * the K-value for attribution of DND would be a good orientation for the range.
   * Eventually introduce a factor on the RHS if using DEC205-OVA experiments.
   */
  max_value = 1.0; // used for tc_search_duration_mode == 0
  if (TFHsignal_delivery_mode == 1) { max_value *= TFHsignal_delivery_max; }
  monitor_TfhSignalSpeed_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					     "TfhSignalSpeed_selected.out", 
					     "TfhSignalSpeed_selected");
  monitor_TfhSignalSpeed_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
					    "TfhSignalSpeed_deleted.out", 
					    "TfhSignalSpeed_deleted");
  max_value = par.Value.TFHsignal_delivery_max;
  monitor_TFHintensity.initialize(max_value, 200, 0, 0, 
				  "TfhSignalIntensity.out", "Tfh signal intensity");
  /* Initialise the counters for the frequency of AUC of Tfh signal at the time of selection:
   * Eventually introduce a factor on the RHS if using DEC205-OVA experiments
   */
  max_value = 40.0; // used for tc_search_duration_mode == 0
  if (TFHsignal_delivery_mode == 1) { max_value *= par.Value.TFHsignal_delivery_max; }
  monitor_TFHsignalAUC.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				  "TFHsignalAUC.out", "THFsignalAUC");
  monitor_TFHsignalAUC_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
					   "TFHsignalAUC_selected.out", "THFsignalAUCinsel");
  monitor_TFHsignalAUC_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
					  "TFHsignalAUC_deleted.out", "THFsignalAUCindel");
  // cMyc and mTOR
  max_value = 5.;
  monitor_cMycAtSelection.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "cMycAtSelection.out", "cMycAtSelection");
  monitor_cMycAtSelection_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "cMycAtSelection_selected.out", "cMycAtSelectionINsel");
  monitor_cMycAtSelection_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "cMycAtSelection_deleted.out", "cMycAtSelectionINdel");
  monitor_mTORatSelection.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "mTORatSelection.out", "mTORatSelection");
  monitor_mTORatSelection_selected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "mTORatSelection_selected.out", "mTORatSelectionINsel");
  monitor_mTORatSelection_deleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax, 
				    "mTORatSelection_deleted.out", "mTORatSelectionINdel");

  //
  // ICOSL:
  ICOSL_dependent_Tfh_signals = par.Value.ICOSL_dependent_Tfh_signals;
  ICOSL_memory = par.Value.ICOSL_memory;
  ICOSL_upregulation_mode = par.Value.ICOSL_upregulation_mode;
  ICOSL_upregulation_time = par.Value.ICOSL_upregulation_time;
  if (ICOSL_dependent_Tfh_signals) {
    ana << "Tfh signals to BCs are modulated with ICOSL expression on BCs.\n"
	<< "ICOSL on BCs is upregulated in mode = " << ICOSL_upregulation_mode << ".\n"
	<< "ICOSL upregulation time is " << ICOSL_upregulation_time << "hours.\n";
  }
  // FoxO-mTORC1
  FoxO_mode = par.Value.FoxO_mode;
  dT_FoxO_start = par.Value.dT_FoxO_start;
  dT_FoxO_reg = par.Value.dT_FoxO_reg;
  stopFoxOonTFH = par.Value.stopFoxOonTFH;
  FoxO_ini = double(par.Value.FoxO_ini);
  FoxO_production = par.Value.deltat / par.Value.dT_FoxO;
  FoxOup_min = par.Value.deltat / par.Value.dT_FoxO_max;
  FoxOup_max = par.Value.deltat / par.Value.dT_FoxO_min;
  // for the Hill function regulating FoxO growth (mode 1,2) or reduction (mode 4)
  FoxOup_K = par.Value.dT_FoxO_K;
  FoxOup_n = par.Value.dT_FoxO_n;
  // and for the Hill function regulating Tfh-dependent FoxO-deflections:
  KFoxO = par.Value.KFoxO;
  nFoxO = par.Value.nFoxO;
  mTOR_mode = par.Value.mTOR_mode;
  mTORC1_production = par.Value.deltat / par.Value.dT_mTORC1;
  mTOR_Tfh_production = par.Value.deltat / par.Value.mTOR_dTfh;
  mTOR_BCR_K = par.Value.mTOR_BCR_K;
  mTOR_BCR_n = par.Value.mTOR_BCR_n;

  //MSchips
  //TFR -- modify this !!
  TFR_CC_interaction_mode = par.Value.TFR_mode;
  time_firstTFRinteraction = par.Value.time_firstTFRinteraction;
  TfhHelp_reduction_factor = par.Value.TfhHelp_reduction_factor;
  TFR_bind_onlyself = par.Value.TFR_bind_onlyself;
  tfr_time = par.Value.TFR_time;
  exp_CCL3KO = par.Value.exp_CCL3KO;

  max_value = 100;
  monitor_nTFRcontacts_SelfCCselected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "nTFRcontacts_SelfCCselected.out",
                                          "nTFRcontacts_SelfCCselected");
  monitor_nTFRcontacts_SelfCCdeleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "nTFRcontacts_SelfCCdeleted.out",
                                          "nTFRcontacts_SelfCCdeleted");
  monitor_nTFRcontacts_CCselected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "nTFRcontacts_CCselected.out",
                                          "nTFRcontacts_CCselected");
  monitor_nTFRcontacts_CCdeleted.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "nTFRcontacts_CCdeleted.out",
                                          "nTFRcontacts_CCdeleted");
  monitor_TFHsig_postTFR_CCselected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "TFHsig_postTFR_CCselected.out",
                                          "TFHsig_postTFR_CCselected");
  monitor_TFHsig_postTFR_SelfCCselected.initialize(max_value, 100, par.Value.tmin, par.Value.tmax,
                                          "TFHsig_postTFR_SelfCCselected.out",
                                          "TFHsig_postTFR_SelfCCselected");
  monitor_BinteractTFRtime_self.initialize(3, 1000, par.Value.tmin, par.Value.tmax,
                                    "BinteractTFRtime_self.out", "BinteractTFRtime_self");
  monitor_BinteractTFRtime_nonself.initialize(3, 1000, par.Value.tmin, par.Value.tmax,
                                    "BinteractTFRtime_nonself.out", "BinteractTFRtime_nonself");
  monitor_TFRsig_apopt.initialize(3, 1000, par.Value.tmin, par.Value.tmax,
                                    "TFRsig_apopt.out", "TFRsig_apopt");
  monitor_TFRsig_select.initialize(3, 1000, par.Value.tmin, par.Value.tmax,
                                    "TFRsig_select.out", "TFRsig_select");

  short maxV=20, stepV=20;
  monitor_nMut_nonSelfCD138.initialize(maxV,stepV,par.Value.tmin,par.Value.tmax,
                                       "nMut_nonselfCD138.out","nMut_nonselfCD138");
  monitor_nMut_SelfCD138.initialize(maxV,stepV,par.Value.tmin,par.Value.tmax,
                                       "nMut_selfCD138.out","nMut_selfCD138");
  monitor_nMut_RedeemedCD138.initialize(maxV,stepV,par.Value.tmin,par.Value.tmax,
                                       "nMut_redeemedCD138.out","nMut_redeemedCD138");
  double maxxV=20, stepxV=100;
  monitor_time_first_TFRcontact_NonSelfselected.initialize(maxxV, stepxV,par.Value.tmin,par.Value.tmax,
                                                           "time_first_TFRcontact_NonSelfselected.out","time_first_TFRcontact_NonSelfselected");
  monitor_time_first_TFRcontact_Selfselected.initialize(maxxV, stepxV,par.Value.tmin,par.Value.tmax,
                                                        "time_first_TFRcontact_Selfselected.out","time_first_TFRcontact_Selfselected");


  // Fate tracking frequency:
  write_trackfate = par.Value.write_trackfate; // is given in minutes, thus the factor 60:
  fatetracker_ndt = int( write_trackfate / ( par.Value.deltat * 60 ) );
  // Output settings
  outputfiles=par.Value.outputfiles;
  if (outputfiles == 0) {
    ofstream apolog;
    apolog.open("apolog.out");
    apolog << "! time   state   nFDCcontacts  tc_clock  tc_signal   "
	   << "best_affinity  FDCselected_clock\n";
    apolog.close();
  }
}
void cellCC::trackfate_initialize() {
  fatetracker.mk_new_fate_track();
}
void cellCC::trackfate_show() {
  fatetracker.show_track();
}
//MSchips here we might want to add ntfrcontacts
///(TODO!)
void cellCC::trackfate(double time, bool forceit) {
  if ((fatetracker_n == fatetracker_ndt && state < selected) || forceit) {
    //  if ((fatetracker_n == fatetracker_ndt) || forceit) {
    fatetracker.write_fate(short(state), pos_ss, affinity, time, 
			   nFDCcontacts, get_pMHC_presentation(),
			   FDCselected_clock, tc_search_duration, tc_interaction_time,
			   FoxO, FoxO_upregulation,
			   mTORC1, ICOSL, nTCcontacts, tc_signal, TFHsignalAUC,
			   DND, BCRsignal, cMyc);
    fatetracker_n = 0;
  } 
  ++fatetracker_n;
}
void cellCC::writetrackfate() {
  if (write_trackfate > 0) { fatetracker.write2file(); }
  fatetracker.clear_fate_track();
}
// ==============================================================================
double cellCC::calc_TfhSignalSpeed() {
  /** @brief: calculates the ratio of Tfh-signals to the time in the state FDCselected.
   ** As the B cell only spends a part of his time bound to Tfh and receiving signals,
   ** this ratio is < 1. However, pMHC-dependent Tfh signalling intensity might bring 
   ** this up to larger values. 
   ** This ratio is considered as a measure for the speed of signal upregulation
   ** in the sense of \Delta signal per time. This interpretation allows to use the
   ** units signals per hour. However, the signal is correlated with the time of 
   ** B-Tfh-interactions, such that the signal comes in hours as well.
   **/
  double v = 0.;
  if (FDCselected_clock > 0) {
    v = tc_signal / FDCselected_clock;
  }
  return v;
}
double cellCC::get_decay_factor(double T, double dt) {
  /* @brief: Calculates the factor of reducing the tc_signal per time step:
   * dx/dt = - g x, 
   * with g=ln(2)/T the decay rate and T the half life given as parameter.
   * This is solved by x = x_0 exp(-gt) .
   * Find a factor of reduction per time step F, such that x(t+dt) = F x(t):
   * F = exp(-g dt) = exp(- dt ln(2)/T), which is used here.
   */
  double F = exp(-1.0 * dt * log(2.) / T);
  return F;
}
// ============================================================================
void cellCC::set_statics(const double &time, const Parameter &par) {
   if (SMOOTH_DIFFERENTIATION) {
      set_differentiation(time);
   } else if (time >= par.Value.StartOutput - 1.0e-09) {
      p_dif2out = p_dif2out_target;
      p_dif2out_DEC = p_dif2out_DEC_target;
      cout << "Output-Production (t=" << time << " hr) (discrete).";
   } else {
      p_dif2out = 0.;
      p_dif2out_DEC = 0.;
   }
}
void cellCC::set_differentiation(const double &time) {
   /* This routine has to be called at the beginning of a time step only once
    * because it sets a static variable. Note that the starting time
    * of differentiation is taken as half value of the sigmoidal.
    */
   if (SMOOTH_DIFFERENTIATION) {
      // Use a sigmoidal to switch differentiation on:
      // time is in hours
      double weight = (1.0 + exp((start_differentiate - time) / smooth_differentiation_time));
      p_dif2out = p_dif2out_target / weight;
      p_dif2out_DEC = p_dif2out_DEC_target / weight;
   }
}
double cellCC::set_selected_CC_delay() {
   if ((not (DEC205_ova) && (dif_delay <= 0.)) || (DEC205_ova && (dif_delay_DEC <= 0.))) {
      return 0.0;
   }
   double shift = 0.;
   double sigmoi;
   double delaytmp = dif_delay;
   if (DEC205_ova) { delaytmp = dif_delay_DEC; }
   if (!CD138) {delaytmp /= 2;}
   // +++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++
   // set the fraction here to determine the width of variation
   double width = delaytmp * 0.1;
   // +++++++++++++++++++ end OPTION ++++++++++++++++++++++++++
   double delaylength = 3.0 * delaytmp;
   while (delaylength <= 0. || delaylength >= 2. * delaytmp) {
      sigmoi = drandom();
      if ((sigmoi == 1.) || (sigmoi == 0.)) { shift = 0; } else {
         shift = width * log(
            (1. - sigmoi) / sigmoi);
      }
      // cout<<"("<<log((1.-gaussf)/gaussf)<<","<<shift<<") ";
      delaylength = (delaytmp + shift);
   }
   //  cout<<".......................... length="<<delaylength<<";\n";
//   cerr<<"TFR "<<TFR_CC_interaction_mode<<endl;

   return delaylength;
}
double cellCC::get_pMHC_presentation() {
  /* pMHC presentation is assumed to be proportional to nFDCcontacts,
   * not to BCRsignal. Thus, retained antigen form the previous round
   * of selection would lead to pMHC presentation in the next round.
   */
   double pMHCp = double(nFDCcontacts); // this is OK if (present_specific_ag2TC == 0)
   if (present_specific_ag2TC > 0) {
       cerr << "if here there is a mistake\n"<<present_specific_ag2TC<<endl;exit(1);
      pMHCp = double (get_max_collected_ag(false));
   }
   /* Increase Ag-presentation for B cells treated with anti-DEC205-OVA.
      In order to make this compatible with present_specific_ag2TC > 0,
      the pMHCp defined above has to be incremented by the pMHC based
      on the antigen provided via DEC205. Changed in v181101 on 181209.
      Used pMHCp = double(pMHC_of_DEC205_ova + nFDCcontacts); before.
    */
   if (DEC205_ova) { pMHCp += double(pMHC_of_DEC205_ova); }
   /* MHC deficient cells do not present anything:
    * Note that MHCdeficient can only be true in DEC205+ BCs by definition in 
    * cellman::recombine_cells(time). */
   if (tmx_MHC && tmx_MHC_noAgPresentation && MHCdeficient) { pMHCp = 0; }
   return pMHCp;
}
double cellCC::get_tc_search_duration() {
  double searchdt = tc_search_duration_fixed;
  /* Note that tc_search_duration_fixed was loaded with 0 for mode==0, and
   * with tc_search_duration_fixed from parameter file for mode==1 and 3.
   * So mode in [0,1,3] are done.
   * In mode 4 tc_search_duration is not used, but FoxO and tc_signal.
   * In mode 5 (FCRB model), the tc_search_duration is determined according to
   *   tc_search_time_determination_mode:
   *   If 0 the time is also is also fixed (done).
   *   If 1 the time is pMHC dependent.
   */
  if ( ( tc_search_duration_mode == 5 && tc_search_time_determination_mode == 1 ) ||
       tc_search_duration_mode == 2 ) {
    /* derive the duration of search for TC help from the amount of pMHC, which
     * is saved in nFDCcontacts (or in pMHC_of_DEC205_ova if DEC205_ova==true).
     * Assumes a linear increase with each portion of collected antigen.
     * Note that BCRsignal is not used here, because the long LZ duration
     * in BCR-independent provision of antigen could not be explained then.
     * So the idea of this choice is that it is antigen processing rather than
     * BCR signalling that prolongs the LZ phase.
     */
     searchdt = get_pMHC_presentation() * tc_search_duration_per_FDCcontact;
  }
  return searchdt;
}
void cellCC::attribute_tc_search_duration() {
  tc_search_duration = get_tc_search_duration();
}
double cellCC::get_tc_interaction_time() {
  double tt = tc_time;
  if (DEC205_ova) { 
     /* In case of DEC205_ova == true, the cell has max possible pMHC presentation
      * and the number of contacts is zero. Here, affinity dependent interaction
      * times are not useful and <tc_dec205ova_binding_time> is always used.
      */
     tt = tc_dec205ova_binding_time; 
  }
  if (mode_of_setting_tc_time == 1) { // gaussian interaction time
     /* The static variable <tc_time> is derived from parameter file and
      * interpreted as the mean interaction time. The static <tc_time_width> 
      * is taken from the parameter file as well.
      */
     double ttwidth = tc_time_width;
     if (DEC205_ova) { // adapt the width to be the same relative width:
	ttwidth = (tc_time_width / tc_time) * tc_dec205ova_binding_time;
     }
     tt = get_positive_sample_from_normal(tt, ttwidth);
  } else if (mode_of_setting_tc_time == 2) { // pMHC-dependent
     /* Here, the number of successful contacts with FDCs (leading to antigen uptake)
      * are used as a measure for affinity and determines the interaction time.
      * The number of collected pMHC portions has to be scaled with a typical measure,
      * which is the parameter BTtime_K. This is also used when the DND-mechanism is
      * switched off. The static variable <tc_time> now corresponds to the interaction
      * time adopted in the case of <BTtime_K> successful contacts with FDCs.
      * This time is made smooth on the last integer with a uniform distribution.
      * In the case of DEC205_ova==true, get_pMHC_presentation() returns the
      * number of FDC contacts <pMHC_of_DEC205_ova> adapted to this case.
      */
     //cout << tt << "(" << get_pMHC_presentation() << ") --> ";
     tt *= ( get_pMHC_presentation() - 0.5 + drandom() ) / BTtime_K;
     //cout << tt << endl;
  } else if (mode_of_setting_tc_time == 3) { // pMHC-dependent Hill function
    double pMHClevel = get_pMHC_presentation() - 0.5 + drandom();
    tt = get_Hill(pMHClevel, BTtime_min, BTtime_max, BTtime_K, BTtime_n);
  }
  //cout<<"tt="<<tt<<"; ";
  return tt;
}
double cellCC::get_tfr_interaction_time() {
    double tt = tfr_time;
    double Ttime=tt;
    if (tt<0) { //interaction as it happens for tfh but half time
        Ttime=get_tc_interaction_time();
        Ttime/=2;
//        if (tt==0) {Ttime=tt;}
    } else if (tt==0) {
        tt = tc_time;
        double ttwidth = tc_time_width;
        tt = get_positive_sample_from_normal(tt, ttwidth);
        Ttime=tt/2;
    }
    if (TFR_bind_onlyself==2 && !selfMutation && !came_from_selected) {Ttime=0.01;}
    return Ttime;
}
// ============================================================
double cellCC::get_ICOSL_expression() {
   /* Returns the level of ICOSL expression of the BC on the PM
    * dependent on the mode of action.
    */
   if (ICOSL_upregulation_mode == 1) {
      if (ICOSL_memory == 1 && n_recycling > 0) {
	 /* When ICOSL_memory is on, the cell that has upregulated ICOSL in the
	  * previous round of selection now uses it. Thus, ICOSL is up right away.
	  */
	 return 1.0;
      }
      /* This uses a deterministic up-regulation with time
       * starting from the time of enterring the state <FDCselected>,
       * which is saved in <FDCselected_clock>.
       * A Hill function is used with parameters set here
       * and with a time constant taken from the parameter file.
       */
      double hillcoefficient = 2.0;
      return get_Hill(FDCselected_clock, ICOSL_upregulation_time, hillcoefficient);
   }
   return 1.0; // case of ICOSL_upregulation_mode == 0
}
      
// ============================================================

void cellCC::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
  if (p_difu_width > 0) {
    p_move = get_positive_sample_from_normal(p_difu, p_difu_width);
  }
}
short cellCC::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   set_polarity_velocity(persistence,v_modi,p_switch_v,n_v_states,v_slow_factor,l,s);
   if (writethis2track == polarisation) {
      double tmpi[l.dim];
      l.get_koord(index,tmpi);
      td.Write_movement(trackno,CC,time,tmpi,polarity,writethis2track);
      writethis2track = trackini;
   }
   short err = do_diffuse(CC,li,l) == 1;
   // do_diffuse==0 if movement happened, i.e. err==1 if no movement happened
   if ((err == 0) && (test_delay == -1)) {
      selectable = 1;
      // cout<<"clock="<<clock<<"; ";
   }
   return err;
}
// ============================================================

void cellCC::set_selectable() {
   set_clock();
   // check for the end of motility suppression after each FDC test:
   if ((mobile == 0) && (ICAM_delay > 0) && (clock > ICAM_delay)) { mobile = 1; }
   // if CC is not selectable but inhibition time is over, make it selectable:
   if ((selectable == 0) && (test_delay > -1) && (clock > test_delay)) {
      selectable = 1;
      // cout<<"clock="<<clock<<", ";
   }
}
// ============================================================

void cellCC::set_CXCR5expression() {
   if (responsive2signal[CXCL13] && (p_CXCR5down > 0) && (drandom() < p_CXCR5down)) {
      responsive2signal[CXCL13] = false;
   }
}
void cellCC::resensitise4CXCL13(sigs &s) {
   if ((s.signal_use[CXCL13] == 1)  // CXCL13 has to be used at all
       && not (responsive2signal[CXCL13])  // otherwise its sensitive anyway
       && (CXCL13recrit > 0)  // re-sensitisation has to be on
       && s.undercritical_signal(index,CXCL13,CXCL13recrit)
       // CXCL13 has to be below the threshold value for re-sensitisation
       ) {
      responsive2signal[CXCL13] = true;
      // cout<<"resensitised CC with index "<<index<<" for CXCL13.\n";
   }
}
void cellCC::set_CXCR4expression() {
   if (responsive2signal[CXCL12] && (cellCB::p_CXCR4down > 0)
       && (drandom() < cellCB::p_CXCR4down)) {
      responsive2signal[CXCL12] = false;
   }
}
void cellCC::resensitise4CXCL12(sigs &s) {
   if ((s.signal_use[CXCL12] == 1)  // CXCL12 has to be used at all
       && not (responsive2signal[CXCL12])  // otherwise its sensitive anyway
       && (CXCL12recrit > 0)  // re-sensitisation has to be on
       && s.undercritical_signal(index,CXCL12,CXCL12recrit)
       // CXCL12 has to be below the threshold value for re-sensitisation
       ) {
      responsive2signal[CXCL12] = true;
      // cout<<"resensitised CC with index "<<index<<" for CXCL12.\n";
   }
}
bool cellCC::set_apoptotic_motility(sigs &s) {
   // returns true if apoptotic cells move at all
   if (apoptotic_motility_mode == 0) { return false; }
   // further the mode of motility is set depending on static short apoptotic_motility_mode
   else if (apoptotic_motility_mode == 1) {
      // continue to be sensitive for CXCL13 ... (but not if CXCR5failure before)
      if (CXCR5failure == 0) {
         set_CXCR5expression();
         resensitise4CXCL13(s);
      }
   } else if (apoptotic_motility_mode == 3) {
      // continue to be sensitive for CXCL12, but down-regulate with rate
      if (drandom() < p_apo_randomwalk) { responsive2signal[CXCL12] = false; }
      if (responsive2signal[CXCL12]) {
         set_CXCR4expression();
         resensitise4CXCL12(s);
      }
      /* Note that in resensitise4CXCL12(s), responsive2signal[CXCL12] might become true
       * again. This is not wanted for apoptotic cells which are supposed to down-regulate
       * everything. So once, responsive2signal[CXCL12] is false, it will remain false.
       */
   }
   // for random walk (apoptotic_motility_mode=2) responsiveness is false but still walk:
   return true;
}
// ============================================================

void cellCC::go2TCselection(AffinitySpace &shape) {
  /* This routine is used to define tc_search_duration instead of progress_selection_state()
   * when antigen search is prevented, ag is preloaded, or anti-DEC205-OVA exps are done.
   */
   if (state != unselected) {
      cout << "Wrong CC-type calling go2TCselection(AffinitySpace&)!\n";
      exit(1);
   } else {
      shape.add_cell(sCCFDCselected,pos_ss);
      shape.rem_cell(sCCunselected,pos_ss);
      state = FDCselected;
      came_from_FDCselected = true;
      FDCselected_clock = 0.;
      tc_search_duration = get_tc_search_duration();
      if (dT_FoxO_reg == FDCselected) { 
	double pMHClevel = get_pMHC_presentation();
	FoxO_upregulation = get_FoxO_rate(pMHClevel); 
      }
   }
}
// ============================================================
void cellCC::delete_antigen() {
  /* Resets nFDCcontacts and collected_ag_portions to zero.
   * However, BCRsignal is kept as it is.
   */
  nFDCcontacts = 0;
  collected_ag_portions.clear();
}
void cellCC::process_antigen() {
  if (ag_deleted_in_fresh_CC) { delete_antigen(); }
}
// ============================================================

long cellCC::contact2FDC(space &l) {
   /* Searches for a neighbour of the object with type FDC and 
    * returns its index on the space-lattice.
    * This cannot be done with cell::find_contact(FDC...) because
    * FDCs might be transparent and would only be identified by also 
    * checking for the lattice property lattice.knot[x].FDClisti.
    *
    * Note that contact2FDC is called in different context such that the test
    * whether the cell is selectable has to stay outside this routine.
    */
   // cout<<"in selection ... ";
   long found_i = -1;
   long j;
   // Pruefe Zelltyp der Nachbarn
   for (int n = 0; n < l.dim2; n++) {
      j = l.knot[index].near_n[n];
      if (j != -1) {
         if ((l.cellknot[j].cell == FDC) || (l.cellknot[j].FDClisti > -1)) { 
	   found_i = j; 
	 }
      }
   }
   return found_i;
}
// ============================================================

short cellCC::bind_antigen(cellFDC &fdc,
                           int frag_index,
                           int &ag_index,
                           AffinitySpace &shape,
                           double &threshold) {
   /* antigen is the shape space position of the antigen on FDC,
    * FDCposition ist der Index des Fragments auf dem benachbarten FDC
    */
   // determine the interaction partner Ag on FDC fragment:
   bool suppress_next_interaction = false;
   if (ag_index < 0) {
      suppress_next_interaction = true;
      // drandom();
      /* This addition call of random number was required for comparison of the results
       * with results in hyphasma15.12.1. Used until 16.03.2 and then deleted.
       * New reference simulations generated in 16.03.2. 
       */
   } else {
      // Only continue if ag is actually there.
      // Determine the position in AffinitySpace of this Ag:
      long agASpos = shape.get_Antigen(ag_index);
      // get the binding probability
      double bindprobability
         = shape.affinity(pos_ss, agASpos, threshold) * BCRexpression;
      // If set, just ignore bindingprobabilities:
      if (ignore_affinity >= 0) { bindprobability = ignore_affinity; }
      // Randomly chose whether the interaction takes place
      if (drandom() < bindprobability) {
         short success = fdc.consume_ag(frag_index,ag_index);
         if (success == 1) {
            // cout<<"affinity o.k.";
            bound_ag_index = ag_index; // save the type of Ag that was bound
            state = contact;
            n_fdc_encounters = 0;
            // shape space:
            shape.add_cell(sCCcontact,pos_ss);
	    if (came_from_FDCselected) { shape.rem_cell(sCCFDCselected, pos_ss); }
	    else { shape.rem_cell(sCCunselected,pos_ss); }
            return 0;     // induces cell immobility in the calling routine!
         } else {
            suppress_next_interaction = true;
         }
      } else {
         suppress_next_interaction = true;
      }
   }
   if (suppress_next_interaction) {
      // Verhindere weitere Tests bis sich die Zelle bewegt hat:
      selectable = 0;
      if (ICAM_delay > 0) { mobile = 0; }
      ++n_fdc_encounters;
      /* Diese Massnahme macht die Affinitaetsabhaengige Auswahl der Zellen vom
       * der jeweiligen Feinheit der Zeitschritte unabhaengig: jede Zelle kann
       * am FDC nur einmal probieren, ob sie passt, danach muss sie sich zunaechst
       * bewegen. Die Haeufigkeit des Testens wird daher primaer durch die
       * Diffusionskonstante der Centrocyten bestimmt. Jetzt neu: internal clock: */
      // Reset internal clock to suppress further affinity tests
      clock = 0;
      // cout<<"bind-->selectable=0";
      return 1;   // allows the cell to move afterwards
   }

   // note : this line should never be called
   return 0;  // induces cell immobility in the calling routine!
}
// ============================================================

void cellCC::progress_selection_state(AffinitySpace &shape) {
  /* This is only called in the state unselected and is associated
   * with single or repeated interactions of CC with FDCs.
   * Called from cellCC::FDCdetach(...) which is called in state contact(to FDC).
   */
  if (TC_CC_selection == 0) {
    // If no selection by TCs is needed go right away to state selected.
    // In case of simultaneousTfhFDC==true, this happens only once and the BC never comes back.
    state = selected;   // that's it
    cerr<<"toadaptfrt"<<endl; exit(1);
    antigen_collection_statistics(nFDCcontacts, true);
    if (immunoglobulin_class::do_switch_classes == 1) { IgX.class_switch(); }
    selected_clock = 0.;
    individual_dif_delay = set_selected_CC_delay();
    shape.add_cell(sCCselected,pos_ss);
  } else {
    // If TC help is needed go to state FDCselected first.
    state = FDCselected;
    // FDCselected_clock is set only at first passage:
    if (not(came_from_FDCselected)) { FDCselected_clock = 0.; }
    // Remember that you have been in state FDCselected before:
    came_from_FDCselected = true;
    tc_search_duration = get_tc_search_duration();
    if (dT_FoxO_reg == FDCselected) { 
      double pMHClevel = get_pMHC_presentation();
      FoxO_upregulation = get_FoxO_rate(pMHClevel); 
    }
    shape.add_cell(sCCFDCselected,pos_ss);
  }
}
void cellCC::return2unselected(AffinitySpace &shape) {
  if (came_from_FDCselected) {
    // parameters relying on pMHC might have changed and have to be adapted:
    progress_selection_state(shape);
  } else {
    state = unselected;
    shape.add_cell(sCCunselected, pos_ss);
  }
  clock = 0;
  selectable = 0;
}
void cellCC::add_collected_ag() {
   // in case that the vector didn't have entries at agi by now, extend the vector
   unsigned int uagi = bound_ag_index;
   if (uagi >= collected_ag_portions.size()) {
      collected_ag_portions.resize(bound_ag_index + 1); // adds 0 up to agi
   }
   if (bound_ag_index == -1) {
      // delete this security question after a test phase
      cerr << "ERROR cellCC::add_collected_ag(): bound_ag_index == -1. Abort.\n";
      exit(1);
   }
   // add one Ag portions at agi
   ++collected_ag_portions[bound_ag_index];
   // reset bound_ag_index:
   bound_ag_index = -1;
}
short cellCC::FDCdetach(AffinitySpace &shape) {
   // Assumes that cell==CC, is called in the state contact to FDC in cellman::calc_CC(...)
   short int err = 1;
   // Probability of selection:
   if (drandom() < p_FDCdetach) {
      // Was FDC-signalling sufficient ?
      if (drandom() < p_FDCsignalling) {
	/* The bound antigen has index bound_ag_index; add one portion in collected_ag_portions;
         * do this even if collectFDCsignals==false to memorize for the interaction with TFH */
         add_collected_ag();
         // Evolve the CC state
         if (collectFDCsignals) {
            // If one arrives here, an additional Ag was collected via BCR from FDC; count it
            ++nFDCcontacts;
             ++alterednFDCcontacts;
	    BCRsignal += 1.;
	    // Return to state un(FDC)selected to continue Ag-collection(Tfh-search)
	    return2unselected(shape);
         } else {
            progress_selection_state(shape);
         }
      } else {
         // No FDC-rescue signals -> go back to unselected (or FDCselected)
         return2unselected(shape);
         /// Philippe : should do bind_ag_index = -1, no ?
      }
      // shape space: remove from contact pool anyway
      shape.rem_cell(sCCcontact,pos_ss);
      err = 0;
   }
   return err;
}
// ============================================================
void cellCC::antigen_collection_statistics(int Ncontacts, bool count_selected) {
  // save for ag_portions-statistic:
  int agportionstmp = Ncontacts;
  if (agportionstmp > max_n_of_ag_portions) { agportionstmp = max_n_of_ag_portions; }
  cummulative_ag_collection_all[agportionstmp]++;
  if (count_selected) {
    cummulative_ag_collection_selected[agportionstmp]++;
  }
}

// ============================================================

short cellCC::stop_collecting_FDCsignals(AffinitySpace &shape, const double& time, double &dt) {
  /* Finishes the state unselected and induces the next state.
   * Assumes CC in state unselected; called from cellman::calc_CC(...)
   * Never called from state FDCselected even when simultaneousTfhFDC==true.
   */
  if (collectFDCsignals) {
    selected_clock += dt;
    //if ( (AgThreshold4Diff > 0 && nFDCcontacts >= AgThreshold4Diff) ||
    /* MMH removed this line above january 9th, 2019.
     * It was used in bcintime37a, which was still in testing phase, only.
     * So no harm to other parameter files.
     * This line was not a requirement of a minimum amount of collected antigen for
     * further state progression but a possibility of further state progression.
     * It was meant to kill the cells if they don't meet this requirement.
     * This is now in the condition below deciding on apoptosis.
     */
    if (selected_clock > collectFDCperiod) {
      // Monitor this event
      monitor_pureFDCsearch.add(selected_clock);
      /* Either kill or select:
       * inhibitFDCapoptosis is true for prob2kill_noFDCcontatBCs in [0,1[, i.e.
       * inhibitFDCapoptosis is false for prob2kill...==1 or <0, i.e. kill always.
       */
      if ( ((nFDCcontacts == 0) || (AgThreshold4Diff > 0 && nFDCcontacts < AgThreshold4Diff)) && 
	   (inhibitFDCapoptosis == false || drandom()<prob2kill_noFDCcontactBCs)
           ) {
	/* Case of no FDC contact up to end of collection period and no retained antigen
	 * from previous rounds. BCRsignal is not used here to allow for passing this
	 * step even if no new BCR signal was induced but antigen is kept from previous rounds.
	 */
          if (!exp_stop_apo) {
              make_apoptosis(time, shape, 0);
          } else {return 0;} // Bcl2 experiment on-- cells continue collecting ag --> can't switch to FDCselected
      } else {
	// Case of at least {one, more than AgThreshold4Diff} contacts to FDCs
	progress_selection_state(shape);
	shape.rem_cell(sCCunselected,pos_ss); 
	/* Philippe: suggestion: put that inside the progress_selection_state 
	 *           so the add and rem are together.
	 * MMH: is not as simple, as progress_selection_state without shape.rem_cell
	 *      is called as well. Need a deeper thought to do it right.
	 */
	if (reset_antigen_after_collection > 0) {
	  /* This in silico experiment allows to fix the amount of antigen after
	   * the phase of antigen collection, thus, imposing a particular amount
	   * of collected antigen (saved in nFDCcontacts). The amount of removed
	   * antigen is made consistent with this imposed value.
	   * Note that in the case of simultaneousTfhFDC, the BC in state FDCselected
	   * continues to collect antigen and this resetting is not called again,
	   * such that the imposed value only applies to the first phase of antigen collection.
	   */
            //MS control
            if (TFR_CC_interaction_mode==7||TFR_CC_interaction_mode==70) {
                cerr<<"trogocytosis is on and could cause problems in this experiment";exit(1);
            }
	  double factor = double (reset_antigen_after_collection) / double (nFDCcontacts);
	  nFDCcontacts = reset_antigen_after_collection;
	  /* ================================= OPTION =======================================
	   * As this is meant to reset the amount of collected antigen for pMHC presentation
	   * to Tfh, BCR signal is not set equal to nFDCcontacts. However, depending on the
	   * idea of this in silico experiment, one might synchronize BCRsignal with 
	   * nFDCcontacts by activating the next line.
	   */
	  // BCRsignal = nFDCcontacts;
	  // ============================= end OPTION =======================================
	  for (unsigned int f = 0; f < collected_ag_portions.size(); f++) {
	    collected_ag_portions[f] = int (factor * double(collected_ag_portions[f]) + 0.5); 
	  }
	  /* ### eventually better to use the function adapt_specific_ag(factor)
	   *  and to put collected_ag_portions as well as cellCB::adapt_specific_ag(double)
	   *  into the cell class.
	   */
	}
      }
      // either if apoptotic or if progression of state:
      return 1;
    }
  }
  // either not in mode <collectFDCsignals> or still in collection period:
  return 0;
}
// ============================================================
int cellCC::get_max_collected_ag(bool returnindex) {
  int maxcol = 0, maxindex = 0;
  // cout << "collected=(";
  for (unsigned int f = 0; f < collected_ag_portions.size(); f++) {
    // cout<<collected_ag_portions[f]<<",";
    if (collected_ag_portions[f] > maxcol) {
      maxcol = collected_ag_portions[f];
      maxindex = f;
    }
  }
  // cout<<")\n";
  if (returnindex) { return maxindex; }
  return maxcol;
}
bool cellCC::same_TC_as_before(long new_tc_lattice_index, long new_tc_id) {
  if (no_rebinding_to_moved_TCs) {
    if (last_tc_id == new_tc_id) { return true; }
  } else {
    if (last_tc_index == new_tc_lattice_index) return true;
  }
  return false;
}
bool cellCC::try2findTFH() {
  /* Determines whether it is allowed to bind to TFH:
   * Returns TRUE unless, depending on tc_search_duration_mode,
   * the criteria of interacting with TFH are not fulfilled.
   * This is relevant for the models with multiple BC-TFH interactions,
   * where the period of search for TFH is either fixed, determined
   * by the number of antigen collection events, or by the strength
   * of BC-TFH interactions, or by Sustained-Signalling-Threshold (SST).
   */
  bool try2find = true;
  if (state == FDCselected) { 
    /* if tc_search_duration_mode == 0, it is always allowed to bind TFH.
     * The idea here is a competition between apoptosis and binding TFH.
     * No additional exceptions.
     */
    if ( (tc_search_duration_mode == 1 || tc_search_duration_mode == 2) &&
	 FDCselected_clock > tc_search_duration ) 
      { 
	/* For fixed pMHC-dependent tc_search_duration, the tc_search_duration
	 * is the only determinant of allowing for interactions with Tfh.
	 * It is always allowed to bind Tfh by the end of tc_search_duration.
	 */
	try2find = false; 
      }
    if ( tc_search_duration_mode == 3 && 
	 ( get_pMHC_presentation() == 0 || FoxO > 1.) ) 
      { 
	/* Note that this returns true when tc_search_duration_mode == 3,
	 * provided some antigen was collected and processed.
	 * Thus, the first interaction with Tfh requires some processed antigen.
	 * When FoxO > 1, the BC is stopped from interacting with Tfh.
	 */
	try2find = false; 
      }
    if ( tc_search_duration_mode == 4 &&
	 ( FoxO > 1. || ( SSTactive && tc_signal < SST_tc_signal ) ) ) 
      { 
	/* When Sustained-Signalling-Threshold (SST) is used for differentiation,
	 * FoxO is also active as upper limit of tc_search_duration.
	 * BCs are stopped from interacting with further Tfh when
	 * either FoxO passes the value 1 or when the BC passed the SST once
	 * (so SSTactive==true) and the tc_signal dropped again below the SST.
	 */
	try2find = false; 
      }
    if ( tc_search_duration_mode == 5 &&
	 ( FDCselected_clock > tc_search_duration || mTORC1 >= 1.) ) {
      /* In the FCRB model cells are finishing Tfh search either when the
       * fixed time period of Tfh search is over or when mTORC1 reached
       * its threshold value of 1. */
      try2find = false;
    }
  }
  if ((TFR_CC_interaction_mode==3)
          && selfMutation && nTFRcontacts>0) {try2find=false;} //for apoptosis model no further interaction with tfh if self & tfrbound
  return try2find;
}
double cellCC::get_FoxO() {
  return FoxO;
}
double cellCC::get_prob2bindFDCfirst() {
  /* Here one may introduce a more complicated rule why BCs would
   * prefer to bind a Tfh versus a FDC if both are available.
   * For example one might imagine that a BC with high BCR signalling level
   * would prefer to bind the Tfh, while a BC with low BCR level would go for FDCs.
   */
  return prob2bindFDCfirst;
}
void cellCC::inhibitFoxO(double& pMHClevel) {
  /* This is only active if FoxO_mode == 0.
   * Deflects FoxO depending on the pMHC level.
   * Requires tc_search_duration_mode == 3.
   * This is only called for collectFDCsignals==true.
   * Should the call of this routine be shifted to the BC state TCcontact,
   * one needs to take care that collectFDCsignals==true.
   */
  // Calculate the downwards deflection of FoxO induced by TCR-pMHC signals:
  // At first a Hill-function in dependence of pMHC is calculated
  double pfac = pow(pMHClevel, nFoxO);
  double deflect = pfac / ( pfac + pow(KFoxO, nFoxO));
  // Secondly, the deflections is inhibited by the current mTORC1 level
  FoxO -= deflect * (1. - mTORC1);
  // As the starting FoxO can be < 1 from previous B-T interactions, FoxO < 0 can happen:
  if (FoxO < 0) { FoxO = 0.; }
  /* From now on, FoxO and mTORC1 signals have to evolve.
   * The flag hadcontact2Tfh is set in bind_TC(..) */
}
double cellCC::get_FoxO_rate(double& pMHClevel) {
  double foxrate = FoxO_production; // works for FoxO_mode==3,4
  // Put the formula here in dependence on the modes 1 or 2 
  if (FoxO_mode == 1) { // reduce FoxO growth rate with higher pMHC
    // min + (max - min) (1-H) = max + (min - max) H:
    foxrate = get_Hill(pMHClevel, FoxOup_max, FoxOup_min, FoxOup_K, FoxOup_n);
  } else if (FoxO_mode == 2) { // increase FoxO growth rate with higher pMHC
    foxrate = get_Hill(pMHClevel, FoxOup_min, FoxOup_max, FoxOup_K, FoxOup_n);
  }
  return foxrate;
}
void cellCC::bind_TC(const double& time, cellTC &tcell, space &l, AffinitySpace &shape) {
  if (state==selected) { // security check; not needed, never showed up sinc 2016nov14
    cerr<<"wrong state "<<state<<" of cell_ss_pos="<<pos_ss<<" in bind_TC\n"; 
  }
  ccInd_TFHbound = index;

  state = TCcontact;
  hadcontact2Tfh = true;
  ++nTCcontacts;
  if (TimeOfLastTbinding > 0) // meaning there was a contact before
    { monitor_TimeBetweenTbindings.add(time - TimeOfLastTbinding); }
  // cout<<"bind2TC!\n";
  shape.add_cell(sCCTCcontact,pos_ss);
  shape.rem_cell(sCCFDCselected,pos_ss);
  tc_clock = 0.;
  tc_interaction_time = get_tc_interaction_time();
  /* Setting tc_signal to zero here is not necessary because it is set
     to zero at the time of CB2CC differentation. With only a single CC-TFH interaction
     it is always zero. With multiple CC-TFH interactions it must not be set to zero.
     if (tc_signal > 0) { cout<<tc_signal<<"! "; exit(1); }
     tc_signal = 0.;
  */
  /* The current index on the space lattice is memorized. When the BC detaches again,
   * this very same TC index is excluded from binding. However, when the TC moves
   * it can also bind the same BC again.
   * ##### The latter property might be avoided by introducing an absolute
   *       TC identifier, or generally an identifier in the cell-class
   */
  // save TC-position on the lattice (not that TC remain immobile during contacts)
  tc_index = tcell.index;
  // save TC-identifier to avoid rebinding to the same TC repeatedly.
  last_tc_index = tcell.index;
  last_tc_id = tcell.id;
  // cout<<nFDCcontacts<<", ";
  // Change the state of the TC if necessary and add index to its contact list
  /* ### This is the place to introduce a possible option with which
   *  not the whole amount of collected antigen is transferred to the TFH interaction
   *  but only the part which is specific for one particular antigen with index <ag_index>,
   *  i.e. collected_ag_portions[ag_index] */
  if (collectFDCsignals) {
    /* note: there are two functions make_tc_cc_link, 
     * and in one case they use the TCR sequence or not (MMH: what does "or not" mean here?).
     * #####
     * The following can be simplified:
     * use get_pMHC_presentation() for present_specific_ag2TC in [0,1],
     * remove DEC205_ova as argument of make_tc_cc_link (done in get_pMHC_presentation),
     * check below for the case of not collecting FDC signals
     * #####end
     */
    if (present_specific_ag2TC == 2) {
      // if present_specific_ag2TC==2
      /* ### At this point one may get more specific and determine the ag_index
       *  to which collected_ag_portions[ag_index] is returned from a specifificty
       *  of the <tcell>. Note that the cellTC-class already has a property <pos_ss>,
       *  which is randomly loaded from the available antigens.
       *  For now, the highest amount of collected Ag is provided. 
       *  ISSUE: DEC205_ova positive cells might interact with Tfh of wrong specificity,
       *         which leads to unclear behaviour!
       */
      cerr << "\nWARNING: cellCC::bind_TC(...) with present_specific_ag2TC = "
	   << present_specific_ag2TC << " not yet programmed.\n"
	   << "present_specific_ag2TC == 1 is set now. Continue ...\n\n";
      exit(1);
      present_specific_ag2TC = 1;
    }
    // get the amount of pMHC presented on the surface:
    double pMHClevel = get_pMHC_presentation();
    /* transfer total number of collected Ag to TC or 
     * the amount of collected Ag of the Ag-type with highest amount */
    tcell.make_tc_cc_link(index,pMHClevel);
    /* Now make a downwards deflection of FoxO or adaption of FoxO_upregulation.
     * Only do if the corresponding Tfh search duration mode is set (save CPU time otherwise)
     * This is repeated with different strength for each new contact with Tfh.
     * ### Note that the deflection does not depend on the polarisation 
     *     of the Tfh towards this BC in this version. If this was wanted,
     *     one would have to shift this call to the state TCcontact and
     *     introduce an extra flag to prevent multiple FoxO deflections
     *     from the same B-Tfh contact.
     */
    if (tc_search_duration_mode == 3 || tc_search_duration_mode == 4) {
      if (FoxO_mode == 0) { inhibitFoxO(pMHClevel); }
      else { // if FoxO_mode in [1..3]
	if (dT_FoxO_reg == TCcontact) {
	  FoxO_upregulation = get_FoxO_rate(pMHClevel);
	}
      }
    }
    /* If tc_search_duration_mode == 5 nothing to do here. FoxO is reduced in
     * progress_signalling() in this setting. */
  } else {
    /* If collectFDCsignals==false, the amount of Ag is not relevant, just the type.
     *  The first successful BC-FDC interaction leads to Ag-uptake and progressing selection.
     *  As above, we have different choices controled by <present_specific_ag2TC>:
     *  0: Use affinity of the BCR to the antigen associated with the TC-type
     *  1: Use the BCR-affinity to the ag-type encountered at the time of interaction with FDC
     *  2: no meaning yet, projected on 0
     */
    if (present_specific_ag2TC == 1) {
      // In this version, BCR and antigen index are passed to the TC selection procedure
      tcell.make_tc_cc_link(index,pos_ss,get_max_collected_ag(true),shape,DEC205_ova);
    } else {
      // In this version, the affinity of pos_ss to the TC.pos_ss is used for selection
      tcell.make_tc_cc_link(index,pos_ss,-1,shape,DEC205_ova);
      /* Here, different TC-subpopulations may be distinguished and it gets important
       *  whether you find the right TFH.*/
    }
  }
}
void cellCC::decay_TFHsignal() {
  if (TFHsignal_decay > 0) {
    if (tc_signal > 0) {
      /* Not sure what is faster here:
       * if (tc_signal > 0) which is a really sure variant
       * if (state >= FDCselected) which might get wrong when the CC returns to unselected
       *                           which never happens in the current setting
       * No if at all, which makes a frequent superfluous multiplication with zero.
       * ### TEST speed.
       */
      /* TFHsignal_decay is the factor per time step that corresponds to
       * half life of Werte::TFHsignal_decay */
      tc_signal *= TFHsignal_decay;
    }
  }
}
void cellCC::decay_cMyc() {
  if (cMyc_decay_factor > 0) { cMyc *= cMyc_decay_factor; }
}
void cellCC::addTFHsignalAUC(double& dt) {
  TFHsignalAUC += tc_signal * dt; 
}
void cellCC::progress_signalling() {
  /* FoxO and mTOR signals are progressed here depending on the FoxO_mode.
   * In FoxO_mode == [1,2], the cells come out of the DZ with low FoxO and
   * FoxO is upregulated with a rate determined elsewhere.
   * In FoxO_mode == 0, the cells might either come with low or high FoxO
   * out of the DZ. The FoxO signal is deflected by interactions with Tfh
   * in both cases elsewhere. Regrowth is regulated here.
   *
   * In Foxo_mode == 4, the cells come from the DZ with FoxO=1 and FoxO is
   * reduced in dependence on BCRsignal. 
   * This uses the K and n value dT_FoxO_{K,n} (FoxOup_{K,n}) from the parameter file. 
   * FoxO_upregulation is used as the maximum rate of down-regulation of FoxO 
   * (i.e. rate with highest possible BCRsignal).
   *
   * Note that mTOR_mode > 0 requires FoxO_mode > 3 (forced in setparam.cpp)
   */
  if (FoxO_mode == 4) { // FCRB model
    double FoxOrelease = get_Hill(BCRsignal, FoxOup_K, FoxOup_n);
    FoxO -= FoxO_upregulation * FoxOrelease;
    if (FoxO < 0) { FoxO = 0.; }
    if (mTOR_mode > 0) { // i.e. 1 or 2
      mTOR_BCRcontrol = 1.; // good for mTOR_mode == 1
      if (mTOR_mode == 2) { mTOR_BCRcontrol = get_Hill(BCRsignal, mTOR_BCR_K, mTOR_BCR_n); }
      mTORC1 += mTORC1_production * mTOR_BCRcontrol;
    }
  } else if (FoxO_mode > 0) { // and FoxO_mode <= 3 and mTOR_mode == 0
    // Stop increasing FoxO when already selected or apoptotic:
    if (short(state) >= dT_FoxO_start && state < selected) {
      // Don't increase FoxO when in contact to Tfh when in FoxO_mode==3:
      // if (not(FoxO_mode == 3 && state == TCcontact)) {
      if (not(state == TCcontact && stopFoxOonTFH)) {
	FoxO += FoxO_upregulation;
      }
    }
  } else { // if (FoxO_mode == 0) and (mTOR_mode == 0)
    if (hadcontact2Tfh) { 
      /* Note that hadconact2Tfh is only used in FoxO_mode==0, however I still keep 
       * the FoxO_mode == 0 condition here in case hadcontact2Tfh would be extended to other
       * FoxO_mode-values, which is the case meanwhile. */
      // ... only signal if a first contact has happened
      // ... and tc_search_duration_mode == 3 
      // ... or FoxO_mode == [1,2], i.e. FoxO grows irrespective of 
      mTORC1 += mTORC1_production;
      // mTORC1 is limited to its max value 1 to prevent negative values in inhibitFoxO(double&)
      if (mTORC1 > 1) { mTORC1 = 1.; }
      if (state < selected) { 
	// stop increasing FoxO when already selected or apoptotic
	FoxO += FoxO_upregulation; 
      }
    }
  }
  // FoxO > 1 is tolerated, as this induces differentiation to the DZ phenotype anyway.
}
// ============================================================

void cellCC::make_apoptosis(const double& time, AffinitySpace &shape, bool fromTfh) {
  //   if (outputfiles == 0) {
  ofstream apolog;
  apolog.open("apolog.out", ofstream::app);
  apolog << time << "  " 
	 << state << "  "
	 << nFDCcontacts << "  " 
	 << tc_clock << "  "
	 << tc_signal << "  "
	 << shape.best_affinity_norm(pos_ss) << "  "
	 << FDCselected_clock
	 << "\n";
  apolog.close();
  //}
  shape.add_cell(sCCapoptosis,pos_ss);
  shape.add_cell(sallapoptosis,pos_ss);
  shape.rem_cell(cells(sCCunselected + state),pos_ss);
  shape.rem_cell(total,pos_ss);
  //MS
//  if (selfMutation) {shape.rem_cell(sCCself,pos_ss);}
  antigen_collection_statistics(nFDCcontacts, false);
  if (fromTfh) {
    monitor_nTCcontacts_deleted.add(nTCcontacts);
    monitor_BsearchTtime_deleted.add(FDCselected_clock);
    monitor_TfhSignalAtSelection.add(tc_signal);
    monitor_TfhSignalAtSelection_deleted.add(tc_signal);
    monitor_pMHCatFateDecision_deleted.add(get_pMHC_presentation());
    monitor_TfhSignalSpeed_deleted.add(calc_TfhSignalSpeed());
    monitor_TFHsignalAUC.add(TFHsignalAUC);
    monitor_TFHsignalAUC_deleted.add(TFHsignalAUC);
    monitor_cMycAtSelection.add(cMyc);
    monitor_cMycAtSelection_deleted.add(cMyc);
    monitor_mTORatSelection.add(mTORC1);
    monitor_mTORatSelection_deleted.add(mTORC1);
    //MSchips
    if (selfMutation) {
        monitor_nTFRcontacts_SelfCCdeleted.add(nTFRcontacts);
        if (nTFRcontacts>0) monitor_time_first_TFRcontact_Selfselected.rm(FDCselected_clock);
    } else {
        monitor_nTFRcontacts_CCdeleted.add(nTFRcontacts);
        if (nTFRcontacts>0) monitor_time_first_TFRcontact_NonSelfselected.rm(FDCselected_clock);
    }
    monitor_TFRsig_apopt.add(tfr_signal);
  }

  // Note that TFHsignalAUC and integratedTFHsignal do not distinguish deleted and selected
  state = apoptosis;
  if (apoptotic_motility_mode == 1) {
    responsive2signal[CXCL12] = false;
    responsive2signal[CXCL13] = true;
  } else if (apoptotic_motility_mode == 3) {
    responsive2signal[CXCL12] = true;
    responsive2signal[CXCL13] = false;
  } else {
    // just do random walk ...
    responsive2signal[CXCL12] = false;
    responsive2signal[CXCL13] = false;
  }
}
// ============================================================
void cellCC::make_apoptosis2(const double& time, AffinitySpace &shape, bool fromTfh) {
    cerr<<"dontuse!!!!!\n";exit(1);
  //   if (outputfiles == 0) {
  ofstream apolog;
  apolog.open("apolog.out", ofstream::app);
  apolog << time << "  "
         << state << "  "
         << nFDCcontacts << "  "
         << tc_clock << "  "
         << tc_signal << "  "
         << shape.best_affinity_norm(pos_ss) << "  "
         << FDCselected_clock
         << "\n";
  apolog.close();
  //}
  shape.add_cell(sCCapoptosis,pos_ss);
  shape.add_cell(sallapoptosis,pos_ss);
  shape.rem_cell(sCCTFRcontact,pos_ss);
  shape.rem_cell(total,pos_ss);
  antigen_collection_statistics(nFDCcontacts, false);
  if (fromTfh) {
    monitor_nTCcontacts_deleted.add(nTCcontacts);
    monitor_BsearchTtime_deleted.add(FDCselected_clock);
    monitor_TfhSignalAtSelection.add(tc_signal);
    monitor_TfhSignalAtSelection_deleted.add(tc_signal);
    monitor_pMHCatFateDecision_deleted.add(get_pMHC_presentation());
    monitor_TfhSignalSpeed_deleted.add(calc_TfhSignalSpeed());
    monitor_TFHsignalAUC.add(TFHsignalAUC);
    monitor_TFHsignalAUC_deleted.add(TFHsignalAUC);
    monitor_cMycAtSelection.add(cMyc);
    monitor_cMycAtSelection_deleted.add(cMyc);
    monitor_mTORatSelection.add(mTORC1);
    monitor_mTORatSelection_deleted.add(mTORC1);
  }
  //MSchips
  if (selfMutation) {monitor_nTFRcontacts_SelfCCdeleted.add(nTFRcontacts);}
  else {monitor_nTFRcontacts_CCdeleted.add(nTFRcontacts);}

  // Note that TFHsignalAUC and integratedTFHsignal do not distinguish deleted and selected
  state = apoptosis;
  if (apoptotic_motility_mode == 1) {
    responsive2signal[CXCL12] = false;
    responsive2signal[CXCL13] = true;
  } else if (apoptotic_motility_mode == 3) {
    responsive2signal[CXCL12] = true;
    responsive2signal[CXCL13] = false;
  } else {
    // just do random walk ...
    responsive2signal[CXCL12] = false;
    responsive2signal[CXCL13] = false;
  }
}
// ============================================================

double cellCC::get_pMHC_dependent_division() {
   /* For the momement, this only makes sense together with collectFDCsignals==true.
    * Eventually think of introducing an affinity-dependent measure of DND
    * which would be applied here to get a dynamics number of divisions.
    * For the moment, DND is automatically switched off in setparam.cpp when
    * collectFDCsignals is false.
    * Note that a separate treatment of the DEC205ova case is not necessary,
    * as the amount of pMHC used is adapted in get_pMHC_presentation().
    */
   double pMHC = get_pMHC_presentation();
   /* ### see comment in cellCC::bind_TC and distinguish present_specific_ag2TC = 1,2 later.
    * which extends to the treatment of DEC205_ova when multiple Tfh specificities and
    * antigens are used.
    */
   return get_Hill(pMHC, DND_P_min, DND_P_max, pMHC_dependent_K, DND_nHill);
}
// ============================================================

double cellCC::get_cMyc_dependent_DND() {
  /* This calculates the DND based on the current level of cMyc (instead of Tfh signal or pMHC).
   * It uses a Hill-function very much as it was used for DND and Tfh signals.
   * The function uses the same min, max, and the same Hill-coefficient
   * as set for pMHC-dependent, but an extra K specific to cMyc-dependent DND.
   */
   return get_Hill(cMyc, DND_P_min, DND_P_max, cMyc_dependent_K, DND_nHill);
}
// ============================================================

double cellCC::get_signal_induced_number_of_divisions() {
   /* Calculates the number of divisions based on the amount or AUC of collected
    * signals from Tfh, saved in <tc_signal> or <TFHsignalAUC>, respectively.
    * A Hill function with parameters from the parameter file is used.
    * Note that a separate treatment of the DEC205ova case is not necessary,
    * as the BCs which got antigen via DEC205 already have a longer period
    * of search for Tfh, which should be translated into a higher signal level
    * when this routine is called.
    */
  double ndiv = 0;
  if (mode_of_DND_of_Tfh == 0) { // i.e. DND is calculated from absolute signal levels
    ndiv = get_Hill(tc_signal, DND_P_min, DND_P_max, TFHsignal_dependent_K, DND_nHill);
  } else if (mode_of_DND_of_Tfh == 2) { // use AUC
    ndiv = get_Hill(TFHsignalAUC, DND_P_min, DND_P_max, TFHsignal_dependent_K, DND_nHill);
  } else { // mode_of_DND_of_Tfh == 1, i.e. DND is calculated from the gradient of signal
    double sig_speed = tc_signal / FDCselected_clock;
    ndiv = get_Hill(sig_speed, DND_P_min, DND_P_max, TFHgradient_dependent_K, DND_nHill);
  }
  return ndiv;
}
// ============================================================

short cellCC::get_tc_selected(AffinitySpace &shape) {
  /** @brief: Progresses the B cell state to <selected> provided conditions are met.
   ** The used condition depends on the settings in the parameter file.
   ** In case of selection, also DND is calculated.
   **/
  // cerr<<"..............................("<<tc_signal<<","<<rescue_signal_threshold<<")\n";
  if ( (AgThreshold4Selection > 0 && nFDCcontacts > AgThreshold4Selection) ||
       (tc_search_duration_mode <= 3 && tc_signal >= rescue_signal_threshold) ||
       (tc_search_duration_mode == 4 && TFHsignalAUC > rescue_signal_threshold) ||
       (tc_search_duration_mode == 5 && mTORC1 >= 1. && cMyc > rescue_signal_threshold) ||
       ((TFR_CC_interaction_mode==300||TFR_CC_interaction_mode==500) && (tc_signal >= rescue_signal_threshold && tfr_signal<tc_signal) ) ) {
    //cerr<<"Select it: dt="<<tc_signal<<"\n";
    if (immunoglobulin_class::do_switch_classes == 1) { IgX.class_switch(); }
    if ((TFR_CC_interaction_mode==21 && (nTFRcontacts==0 || !selfMutation))
            || TFR_CC_interaction_mode!=21) {
        if (drandom()<p_dif2out_target/* 0.1*/) {
            CD138=1;
            if (selfMutation) {monitor_nMut_SelfCD138.add(n_mutation);}
            else if (redeemed) {monitor_nMut_RedeemedCD138.add(n_mutation);}
            else {monitor_nMut_nonSelfCD138.add(n_mutation);}
        }
    }
    individual_dif_delay = set_selected_CC_delay();
    //cerr<<"DEC="<<DEC205_ova<<", tc_signal="<<tc_signal<<", pMHC="<<get_pMHC_presentation();
    if (pMHC_dependent_division) {
      DND = get_pMHC_dependent_division();
    } else if (SIND) { // Signal Induced Number of Divisions
      DND = get_signal_induced_number_of_divisions();
    } else if (FCRB) { // pathway model involving cMyc, mTORC1, and FoxO1
      DND = get_cMyc_dependent_DND();
    }
    if ((TFR_CC_interaction_mode==16 || TFR_CC_interaction_mode==20 || TFR_CC_interaction_mode==34) && nTFRcontacts>0) {
        //these are the models where pre- & post-selection effect of contact with TFR is different
        //pre-selection constact affects only self-reactive cells which won't reach this state
        //post-selection will affect both types differently (see dif2outorcb fct)
        //therefore ntffrcontacts needs to be reset
        nTFRcontacts=0;
    }
    //cerr<<", NDIV="<<DND<<"\n";
    shape.add_cell(sCCselected,pos_ss);
    if (state == TCcontact) {
      shape.rem_cell(sCCTCcontact,pos_ss);
    } else { //if (state == FDCselected)
      shape.rem_cell(sCCFDCselected,pos_ss);
    }	
    // For the statistics:
    antigen_collection_statistics(nFDCcontacts, true);
    monitor_nTCcontacts_selected.add(nTCcontacts);
    monitor_BsearchTtime_selected.add(FDCselected_clock);
    monitor_TfhSignalAtSelection.add(tc_signal);
    monitor_TfhSignalAtSelection_selected.add(tc_signal);
    monitor_pMHCatFateDecision_selected.add(get_pMHC_presentation());
    monitor_TfhSignalSpeed_selected.add(calc_TfhSignalSpeed());
    monitor_TFHsignalAUC.add(TFHsignalAUC);
    monitor_TFHsignalAUC_selected.add(TFHsignalAUC);
    monitor_cMycAtSelection.add(cMyc);
    monitor_cMycAtSelection_selected.add(cMyc);
    monitor_mTORatSelection.add(mTORC1);
    monitor_mTORatSelection_selected.add(mTORC1);
    //MSchips
    if (selfMutation) {
        monitor_nTFRcontacts_SelfCCselected.add(nTFRcontacts);
        if (nTFRcontacts>0) monitor_TFHsig_postTFR_SelfCCselected.add(tc_signal);
    } else {
        monitor_nTFRcontacts_CCselected.add(nTFRcontacts);
        if (nTFRcontacts>0) monitor_TFHsig_postTFR_CCselected.add(tc_signal);
    }
    monitor_TFRsig_select.add(tfr_signal);
    // change states ...
    state = selected;
    if ((TFR_CC_interaction_mode==13 || TFR_CC_interaction_mode==17)&&nTFRcontacts>0) {cerr<<"shouldnothaveinteracted\n";exit(1);}
    selected_clock = 0.;
    //cout<<nTCcontacts<<", ";
    return 1;
  }
  return 0;
}
// ============================================================
bool cellCC::do_selection() {
  /** @brief: Called by try_selection(...) only (see next routine).
   ** Returns true if a decision for selection or apoptosis can be made.
   ** Depending on the mode set in tc_search_duration_mode, conditions are different.
   **/
  if (AgThreshold4Selection > 0 && nFDCcontacts > AgThreshold4Selection) { return true; }
  if (tc_search_duration_mode == 0) { return false; }
  if (tc_search_duration_mode < 3 && FDCselected_clock > tc_search_duration) { return true; }
  if (tc_search_duration_mode == 3 && FoxO > 1.0) { return true; }
  if (tc_search_duration_mode == 4 && 
      (FoxO > 1.0 || ( SSTactive && tc_signal < SST_tc_signal) ) ) { return true; }
  if (tc_search_duration_mode == 5 &&
      (FDCselected_clock > tc_search_duration || mTORC1 >= 1.)) { return true; }
  return false;
}
short cellCC::try_selection(const double &time, const double &dt, AffinitySpace &shape) {
  // Is called from calc_CC in the state FDCselected.
  FDCselected_clock += dt;
  if (do_selection()) {
     /* We come here when the tc_search_duration_mode is either determined by a fixed
      * duration or by a duration determined by the number of successful contacts with FDCs,
      * or by dynamic attribution of Tfh search periods in dependence on FoxO signals. 
      * In this case, selection at the time when the signal threshold is passed
      * during an interaction with Tfh is not helpful. Selection has to be done
      * at a time when the BC is not bound to a Tfh. This is done here.
      */
     /* If the integrated signal is below the survival threshold do apoptosis
      * else, if the integrated signal is above the threshold calculate the number of divisions.
      * check whether the selection signal threshold was met.
      * gotit==1 if the selection threshold was met, ==0 otherwise. 
      */
     short gotit = get_tc_selected(shape);
     // if (gotit == 0) { cerr<<"no success in get_tc_selected ...................\n"; }
     /* The parameter TCell (=TCselect_prob) in setparam controls if all successful BC are
      * surviving (default=1) */
     if ( ((gotit == 1) && (TCselect_prob < 1) && (drandom() > TCselect_prob)) ) {
       make_apoptosis(time, shape, 1); 
       cerr<<"here11";exit(1);
       // Note that when this happens, the cell was already counted in nTCcontacts_selected:
       monitor_nTCcontacts_selected.rm(nTCcontacts);
       monitor_BsearchTtime_selected.rm(FDCselected_clock);
       // was added in get_tc_selected and in make_apoptosis:
       monitor_TfhSignalAtSelection.rm(tc_signal); 
       // was wrongly added in get_tc_selected:
       monitor_TfhSignalAtSelection_selected.rm(tc_signal); 
       monitor_pMHCatFateDecision_selected.rm(get_pMHC_presentation()); 
       monitor_TfhSignalSpeed_selected.rm(calc_TfhSignalSpeed());
       monitor_TFHsignalAUC.rm(TFHsignalAUC);
       monitor_TFHsignalAUC_selected.rm(TFHsignalAUC);
       monitor_cMycAtSelection.rm(cMyc);
       monitor_cMycAtSelection_selected.rm(cMyc);
       monitor_mTORatSelection.rm(mTORC1);
       monitor_mTORatSelection_selected.rm(mTORC1);
       //MSchips
       if (selfMutation) {
           monitor_nTFRcontacts_SelfCCselected.rm(nTFRcontacts);
           if (nTFRcontacts>0) monitor_TFHsig_postTFR_SelfCCselected.rm(tc_signal);
       } else {
           monitor_nTFRcontacts_CCselected.rm(nTFRcontacts);
           if (nTFRcontacts>0) monitor_TFHsig_postTFR_CCselected.rm(tc_signal);
       }
       monitor_TFRsig_select.rm(tfr_signal);
     }
     return gotit; // also returns 1 when apoptosis was done!
   }
   return 0;
}
// ============================================================
double cellCC::get_TFHsignal_delivery_factor(double pMHC, bool calculate) {
  /* Returns the value of the Hill-function for pMHC-dependent Tfh signalling strength. 
     Note that the usage of the array is not compatible with double pMHC values.
     Should those be introduced at some point, the calculation has to be imposed
     by setting the flag <calculate==true>.
   */
  if (not(calculate) && pMHC < TFHsignal_delivery_max-pMHC) { 
    return TFHsignal_delivery_factor[int(pMHC)]; 
  }
  return get_Hill(pMHC, TFHsignal_delivery_min, TFHsignal_delivery_max,
		  TFHsignal_delivery_KpMHC, TFHsignal_delivery_n);
}
double cellCC::get_TFHsignal_delivery_factor(double pMHC, cellTC& TC, bool calculate) {
  /* Same as above but allowing for an adapted K-value provided by TC.
   * Note that when the K-value is supposed to be adapted by the TFH memory,
   * the array cannot be used (to save time) and the calculation has to be
   * performed even when called with calculate==false.
   */
  if (not(calculate) && not(TFHadapt2pMHCseen) && pMHC < TFHsignal_delivery_max-pMHC) { 
    return TFHsignal_delivery_factor[int(pMHC)]; 
  }
  double K_value = TFHsignal_delivery_KpMHC;
  if ( TFHadapt2pMHCseen ) { K_value = TC.K_signal_intensity_pMHC; }
  return get_Hill(pMHC, TFHsignal_delivery_min, TFHsignal_delivery_max,
		  K_value, TFHsignal_delivery_n);
}
void cellCC::add_tc_signal(double& intensity) {
  /* Add signals to the CC. The addition depends on the kind of signals used.
     Up to v20160927 this is always time (duration of interaction).
     Then ICOSL was used to modulate the signal strength.
     Multiplication with ICOSL is 1 normally, unless ICOSL_dependent_Tfh_signals == true.
     Since v20180130, the signal strength can also be modified irrespective of ICOSL.
  */
  double portion = max_tc_signal_portion;
  if (ICOSL_dependent_Tfh_signals) {
    ICOSL = get_ICOSL_expression(); 
    // returns a value in [0,1] based on FDCselected_clock
    portion *= ICOSL;
    //   cerr<<"<"<<ICOSL<<"> ";
  }

  if (TFHsignal_delivery_mode == 1) { portion *= intensity; }
  if (tmx_MHC && tmx_MHC_noTfhSignal && MHCdeficient) { portion = 0; }
  tc_signal += portion;
  // Set SSTactive if the cell gets running (used in tc_search_duration_mode == 4):
  if (tc_signal > SST_tc_signal) { SSTactive = true; }
  /* One may add alternatives here.
   * Note that there are parameters prepared in setparam.h which are TC properties
   * which might be used to modulate TC signalling.
   */
}
double cellCC::get_cMyc_FoxO_break() {
  return 1. - get_Hill(FoxO, cMyc_FoxO_break_K, cMyc_FoxO_break_n);
}
void cellCC::add_cMyc_signal(double& intensity) {
  /* Add cMyc signals to the CC.
   * cMyc signals might be modulated by ICOSL or by pMHC-dependent signalling strength.
   * In addition, cMyc is modulated by the Foxo level.
   */
  double portion = max_tc_signal_portion;
  if (ICOSL_dependent_Tfh_signals) {
    ICOSL = get_ICOSL_expression(); 
    // returns a value in [0,1] based on FDCselected_clock
    portion *= ICOSL;
    //   cerr<<"<"<<ICOSL<<"> ";
  }
  if (TFHsignal_delivery_mode == 1) { portion *= intensity; }
  if (tmx_MHC && tmx_MHC_noTfhSignal && MHCdeficient) { portion = 0; }
  // Calculate the Foxo break of cMyc signalling according to Luo et al 2018
  if (cMyc_FoxO_break) { portion *= get_cMyc_FoxO_break(); }
  cMyc += portion;
  //ms
  ///check to remove at some point --never appeared
  if ((TFR_CC_interaction_mode==3) && selfMutation && nTFRcontacts>0) {
      if (cMyc>0) {cerr<<cMyc<<"WRONG MYC LEVEL!"<<endl; exit(1);}
  }
}
void cellCC::add_mTOR_signal(double& intensity) {
  if (mTOR_mode > 0) { // i.e. 1 or 2
    /* Calculate \alpha_R H_R(BCRsignal) T(pMHC)
     * mTOR_BCRcontrol was already calculated in progress_signalling() called in every time step
     */
    double portion = mTOR_Tfh_production; // this is \alpha_R T(pMHC)
    if (mTOR_mode == 2) { 
      portion *= mTOR_BCRcontrol; // this is H_R(BCRsignal)
    }
    /* Note that T(pMHC) is included in mTOR_Tfh_production because this
     * is repeated for every time step as long the Tfh is polarised towards the BC.
     * One may think of adding a real pMHC-dependent Tfh signalling level:
     * This has to be set in the parameter file in the <Mode of Tfh signal transmission>. */
    if (TFHsignal_delivery_mode == 1) { portion *= intensity; }
    // Optionally, ICOSL-regulation of Tfh signals might be used:
    if (ICOSL_dependent_Tfh_signals) {
      ICOSL = get_ICOSL_expression(); 
      // returns a value in [0,1] based on FDCselected_clock
      portion *= ICOSL;
      //   cerr<<"<"<<ICOSL<<"> ";
    }
    // In case of MHC-deficiency, no signal is transmitted:
    if (tmx_MHC && tmx_MHC_noTfhSignal && MHCdeficient) { portion = 0; }
    mTORC1 +=  portion;
  }
}
// ============================================================
short cellCC::got_tc_signals(const double &time, const double &dt, cellTC &tcell, 
			     space &l, AffinitySpace &shape) {
  /** @brief: The name of the routine has to be read as a question.
   ** It is called from calc_CC in cellman.cpp for B cells in the state <TCcontact>.
   ** Note, that selection from the state <FDCselected> is controlled in <try_selection(...)>.
   ** At first, Tfh signals are added and clocks are progressed.
   ** If the duration of the Tfh-B-interaction was reached, the signal level
   ** is checked for positive or negative selection (depending on the settings),
   ** and the B cell is deliberated from the Tfh.
   ** Positive selection is then induced in <get_TC_selected(...)> and negative
   ** selection is induced in <make_apoptosis(...)>.
   ** It returns 1, when the B cell shall stay bound to the Tfh, 0 otherwise.
   **/
   short staybound = 1;
   // Try to get selected by TC attention:
   tc_clock += dt;
   FDCselected_clock += dt; // also progress this one in case of multiple TC signal integration
   // cout<<"tc_clock="<<tc_clock<<"; dt="<<dt<<"\n";
   // Checks if polarity is set on this cell:
   long ix[l.dim];
   short howmany = l.get_nn_directed2( tcell.polarity, tc_index, ix );
   /* double check if it's working:
    * cout<<"pol=("<<tcell.polarity[0]<<","<<tcell.polarity[1]<<") ";
    * double tmp[l.dim];
    * l.get_diff_vector(tc_index,index,tmp);
    * cout<<"dif=("<<tmp[0]<<","<<tmp[1]<<") ";
    */
   if ((howmany == 1) && (ix[0] == index) && (tcell.able2signal) ) {
     // only add signal if TC polarisation is pointing to a single cell
     double intensity = 1.;
     if (TFHsignal_delivery_mode == 1) {
       intensity = get_TFHsignal_delivery_factor( get_pMHC_presentation(), tcell, false );
     }
     if (TFR_CC_interaction_mode==33||TFR_CC_interaction_mode==35) {
         FoxO=0.0;
         if (TFR_CC_interaction_mode==35) {
             add_mTOR_signal(tfr_interaction_time);
         }
     }
     monitor_TFHintensity.add(intensity);
     add_tc_signal(intensity);
     add_cMyc_signal(intensity);
     add_mTOR_signal(intensity);
     // cout<<" get_signal!";
   }
   // cout<<"\n";
   bool try2select = false;
   if (BCstaysonTCbyTCtime) {
     if ( tc_clock > tc_interaction_time ) { 
       try2select = true; 
     }
   } else {
     if (not (DEC205_ova) || (tc_clock > tc_interaction_time)) { 
       // If DEC205_ova==false, try2select=true is always set.
       try2select = true; 
     }
   }
   if (try2select) {
     /* check whether the selection signal threshold was met.
      * gotit==1 if the selection threshold was met, ==0 otherwise. */
     short gotit = 0;
     /* If tc_search_duration_mode>0, this means that the time of search for Tfh is determined
      * otherwise. Thus, every started BT interaction is finished and the decision about selection
      * is taken back in the state FDCselected in try_selection(...).
      */
     if (tc_search_duration_mode == 0 || force_selection_at_TFHthreshold) {
       gotit = get_tc_selected(shape);
       /* The parameter TCell (=TCselect_prob) in setparam controls if all successful BC are
	* surviving (default=1) */
       if ((gotit == 1) && (TCselect_prob < 1) && (drandom() > TCselect_prob)) {
	 /* TCselect_prob is now interpreted as negative selection, which is a feature
	  * dealt with a few lines below. In case of single T-B-interactions this is fine,
	  * because the BC which is not positively selected will be killed anyway.
	  * But in case of multiple T-B-interactions, TCselect_prob should only reduce
	  * the probability of positive selection and must not kill the BC.
	  * In settings with multiple T-B-interactions TCselect_prob doesn't make sense and
	  * Werte::TCell=1 -> cellCC::TCselect_prob is forced during reading parameters.
          */
	 make_apoptosis(time, shape, 1); 
         cerr<<"comehere2"<<endl;exit(1);
	 /* Note that when this happens, the cell was already counted in nTCcontacts_selected.
	  * It doesn't matter for tc_search_duration_mode == 0.
	  * But when force_selection_at_TFHthreshold is true, tc_search_duration_mode can
	  * be different from 0. Still correct. 
	  * See do_selection for further explanations (same applies here).
	  */
	 monitor_nTCcontacts_selected.rm(nTCcontacts);
	 monitor_BsearchTtime_selected.rm(FDCselected_clock);
	 monitor_TfhSignalAtSelection.rm(tc_signal);
	 monitor_TfhSignalAtSelection_selected.rm(tc_signal);
	 monitor_pMHCatFateDecision_selected.rm(get_pMHC_presentation());
	 monitor_TfhSignalSpeed_selected.rm(calc_TfhSignalSpeed());
	 monitor_TFHsignalAUC.rm(TFHsignalAUC);
	 monitor_TFHsignalAUC_selected.rm(TFHsignalAUC);
	 monitor_cMycAtSelection.rm(cMyc);
	 monitor_cMycAtSelection_selected.rm(cMyc);
	 monitor_mTORatSelection.rm(mTORC1);
         monitor_mTORatSelection_selected.rm(mTORC1);
         //MSchips
         if (selfMutation) {
             monitor_nTFRcontacts_SelfCCselected.rm(nTFRcontacts);
             if (nTFRcontacts>0) monitor_TFHsig_postTFR_SelfCCselected.rm(tc_signal);
         } else {
             monitor_nTFRcontacts_CCselected.rm(nTFRcontacts);
             if (nTFRcontacts>0) monitor_TFHsig_postTFR_CCselected.rm(tc_signal);
         }
         monitor_TFRsig_select.rm(tfr_signal);
         staybound = 0;
       }
     }
     /* Check whether the cell was too long on the TC, then go for apoptosis if
      * negative selection is on.
      * Note that negativeTCselection is always false when tc_search_duration>0
      * or multipleTFHcontacts==true.
      */
     if (negativeTCselection && (gotit == 0) && (tc_clock > tc_interaction_time)) {
       make_apoptosis(time, shape, 1);
       gotit = 1;
     }
     /* Another reason to get off is that the interaction time is over
      * and not sufficient signal yet and (multiple interactions with Tfh on
      * or tc_search_duration_mode>0)
      */
     if (gotit == 0 
	 && ( multipleTFHcontacts || tc_search_duration_mode > 0) // works for 1..5
	 && tc_clock > tc_interaction_time) {
       state = FDCselected;
       // do not reset FDCselected_clock here!
       shape.add_cell(sCCFDCselected,pos_ss);
       shape.rem_cell(sCCTCcontact,pos_ss);
       gotit = 1;
     }
      /* gotit == 1 if CC passed selection threshold, got apoptosis signal,
       * or reached the end of the interaction time with this Tfh in order
       * to get more signals from another ones.
       * Deliberate the cell from the Tfh */
     if (gotit == 1) { 
       monitor_BinteractTtime.add(tc_clock);
       TimeOfLastTbinding = time;
       tcell.liberateCC(index);
       tcell.adapt_K_signal_intensity(int(get_pMHC_presentation()));
       staybound = 0;
     }
     /* The behaviour in the case gotit==0 (didn't get sufficient signal so far) 
      * and negativeTCselection==false is that the BC stays in contact to the TFH. 
      * In the next time step, it might get further signals from the TFH.
      * Every cell will then get sufficient signal to continue with 1 division at some point.
      * What happens if in addition BCstaysonTCbyTCtime==true? The same.
      * It can only be prevented by setting multipleTFHcontacts==true.
      if (gotit==1) {
      cout<<nFDCcontacts<<", "<<tc_clock<<", "<<tc_signal<<", "<<tc_time
      <<", "<<gotit<<"\n";
      }
      */
   }
   return staybound;
}
// ============================================================
//ms
short cellCC::tfh_pointing(const double &time, const double &dt, cellTC &tcell,
                           space &l, AffinitySpace &shape) {
    // this fct returns 1 if this cc is the one receiving the tfh signal, 0 otherwise
    short point = 0;
    long ix[l.dim];
    short howmany = l.get_nn_directed2( tcell.polarity, tc_index, ix );
    if ((howmany == 1) && (ix[0] == index) && (tcell.able2signal) ) {
        point=1;
    }
    return point;
}
// ============================================================
bool cellCC::canbindTFR(const double &time) {
    /** @brief: this functions returns true by default
     * conditions are set for stopping CC-Tfr interaction **/
    bool try2find = true;
    if (exp_CCL3KO && ccl3KO) {try2find=false;} //highest in hierarchy

    //BCs can bind TFR given they still have time
    if (state == FDCselected) {
        if (time<time_firstTFRinteraction // time-dep Tfr action
                || FDCselected_clock > tc_search_duration
                || TFR_CC_interaction_mode==13 //SemiGate
                || TFR_CC_interaction_mode==19 /*SemiGate38*/) {
            //In the gateing models only TC-selected (state=selected) CCs can bind Tfrs
            try2find = false;
        }
    }
    if (TFR_CC_interaction_mode<0
            || (TFR_bind_onlyself==1 && !selfMutation)) {
        try2find = false;
    }
    if ((TFR_CC_interaction_mode==3)
            && selfMutation && nTFRcontacts>0) {
        //In the Apoptosis model s--BCs do not rebind
        try2find=false;
    }
    return try2find;
}
//for the time being tc_clock and tc_interaction_time are kept as for tfh -->CHECK THIS!!!
void cellCC::bind_TFR(const double& time, cellTFR &tcell, space &l, AffinitySpace &shape) {
//  if (state!=FDCselected) { // security check
//    cerr<<"in bind_tfr wrong state "<<state<<" num cont "<<nTFRcontacts<<" of cell_ss_pos="<<pos_ss<<"\n";
//    exit(1);
//  }
    if (state==FDCselected) { shape.rem_cell(sCCFDCselected,pos_ss);}
    else if (state==TCcontact) { shape.rem_cell(sCCTCcontact,pos_ss); cerr<<"not yet tested!"<<endl; exit(1);}
    else if (state==selected) { shape.rem_cell(sCCselected,pos_ss); came_from_selected=1;}


    if (TFR_CC_interaction_mode==13|| TFR_CC_interaction_mode==17||TFR_CC_interaction_mode==19) {
        //security check
        if (state!=selected) {cerr<<state<<" wrong state for gate models"<<endl;exit(1);}
    }

  state = TFRcontact;
  shape.add_cell(sCCTFRcontact,pos_ss);

  ++nTFRcontacts;

  int tmpcont=nTFRcontacts;
  if (tmpcont>ntfr) {tmpcont=ntfr;}
  if (selfMutation) {
      if (nTFRcontacts==0) monitor_time_first_TFRcontact_Selfselected.add(FDCselected_clock);
      num_tfr_contacts_self[tmpcont]++;
  } else {
      if (nTFRcontacts==0) monitor_time_first_TFRcontact_NonSelfselected.add(FDCselected_clock);
      num_tfr_contacts_nonself[tmpcont]++;
  }

  ccInd_TFRbound=index;

  tfr_clock = 0.;
  tfr_interaction_time = get_tfr_interaction_time();

  tfr_index = tcell.index;

  last_tfr_index = tcell.index;
  last_tfr_id = tcell.id;
  tcell.make_tfr_cc_link(index);
  tfr_tmp_index=-1;

}
//ms
bool cellCC::same_TFR_as_before(long new_tfr_lattice_index, long new_tfr_id) {
    ///no_rebinding parameter is the same used for TC -- no need for another one
    ///1=yes; 0=it can once the Tcell moved
  if (no_rebinding_to_moved_TCs) {
    if (last_tfr_id == new_tfr_id) { return true; }
  } else {
    if (last_tfr_index == new_tfr_lattice_index) return true;
  }
  return false;
}
short cellCC::got_tfr_interaction(const double &time, const double &dt, cellTFR &tcell,
                             space &l, AffinitySpace &shape) {
    short staybound=1;

    tfr_clock+=dt;
    FDCselected_clock+=dt;
    tfr_signal+=dt;

    if (index!=ccInd_TFRbound) {
        cerr<<"in gottfr fct index now: "<<index
           <<" index at binding "<<ccInd_TFRbound<<endl;
        cerr<<"forced out"<<endl;
        exit(1);
    }
    bool try2unbind = false;
    if (BCstaysonTCbyTCtime) {
        if (tfr_clock > tfr_interaction_time || tcell.to_delete) { try2unbind = true; }
    }
    if (try2unbind) {
        if (TFR_CC_interaction_mode!=-1) {
            if (selfMutation) {
                if (TFR_CC_interaction_mode==3) {
                    tfhHelpReduction = TfhHelp_reduction_factor;
                    tc_signal=0;
                    cMyc=0; //In DisseD and MiXed cMyc is used for apoptosis vs selection
                }
            }
            if (came_from_selected) {
                state = selected;
                shape.add_cell(sCCselected,pos_ss);
            } else {
                state = FDCselected;
                shape.add_cell(sCCFDCselected,pos_ss);
            }
            shape.rem_cell(sCCTFRcontact,pos_ss);
            staybound = 0;
        }
        if (staybound==0) {
            if (selfMutation) monitor_BinteractTFRtime_self.add(tfr_clock);
            else monitor_BinteractTFRtime_nonself.add(tfr_clock);
            tcell.liberateCC(index);
        }
    }
    return staybound;
}
// ==================================================================================
short cellCC::dif2OUTorCB(double &dt) {
  /* returns 0 if no differentiation
   * returns 1 if differentiation to OUT-cell
   * returns 2 if differentiation to CB
   * returns 3 if differentiation to OUT-cell but with a special delay
   */
  
  /*
   * double p_differentiate=p_dif;
   * if (DEC205_ova) p_differentiate=p_dif_DEC;
   * //  cout<<"(ova="<<DEC205_ova<<"p="<<p_differentiate<<"); ";
   * if (drandom()<p_differentiate) {
   * // In this old version above, the delay of differentiation to CB upon DEC205-ova binding was
   * done as a rate
   * // now it is done with a real delay and then doing normal differentiation with a rate as
   * without DEC:
   */

    //MSchips
    bool force_diff=false;
    bool removecd138=false;

        if ((TFR_CC_interaction_mode==13) && nTFRcontacts>0) { //SemiGate
            individual_dif_delay=0;
            if (CD138) {
                force_diff=true;
                if (selfMutation) {CD138=0; removecd138=true;}
            }
        }
        if ((TFR_CC_interaction_mode==19) && nTFRcontacts>0) { //SemiGate.38
            individual_dif_delay=0;
            if (!CD138) {cerr<<"WRONG TFR ITN MOD19\n";exit(1);}
            if (CD138) {
                force_diff=true;
            }
        }
        if (removecd138) {
            if (selfMutation) {monitor_nMut_SelfCD138.rm(n_mutation);}
            else if (redeemed) {monitor_nMut_RedeemedCD138.rm(n_mutation);}
            else {monitor_nMut_nonSelfCD138.rm(n_mutation);}
        }

  selected_clock += dt;
  if (DEC205_ova) {
    //    if (selected_clock>dif_delay_DEC) {
    if (selected_clock > individual_dif_delay) {
      if (drandom() < p_dif) {
	// Entweder Differentiation zu out oder recycling zu CB
	if (drandom() < p_dif2out_DEC) {
	  selected4output = true;
	  if (p_final_differentiation > 0) { return 3; } else { return 1; }
	} else { return 2; }
      }
    }
  } else {
    // if (DEC205_ova==false)
    // if (selected_clock>dif_delay) {
    if (selected_clock > individual_dif_delay) {
        if (drandom() < p_dif || force_diff) {
            if (TFR_CC_interaction_mode==13 ||
                    (TFR_CC_interaction_mode==19 && CD138)) {
                if (selfMutation) {monitor_nTFRcontacts_SelfCCselected.add(nTFRcontacts);}
                else {monitor_nTFRcontacts_CCselected.add(nTFRcontacts);}
            }
            // Either differentiate to OUT or recycle to CB:
            if (CD138) {
                selected4output = true;
                if (p_final_differentiation > 0) { return 3; } //delayed diff2OUT
                else { return 1; } //immediate diff2OUT
            } else { return 2; } //diff2CB
        }
    }
  }
  return 0;
}
// ============================================================

bool cellCC::final_differentiation() {
   if (selected4output && (p_final_differentiation > 0)) {
      // second condition is redundant!?
      if (drandom() < p_final_differentiation) {
         return true;
      }
   }
   return false;
}
  
// ============================================================
short cellCC::apoptose(const double &time, AffinitySpace &shape) {
   // cout<<"\napoptose";
   short int err = 1;

   /* Try apoptosis if state==unselected. If it's not, it is FDCselected.
    * Then try only when tc_search_duration_mode==0, i.e. do apoptosis during search for Tfh,
    * or when the protected tc_search_duration is over in mode==[1,2],
    * or when the protected tc_search_duration is over in mode==3 without contact to Tfh,
    * or when tc_search_duration_mode==3 FoxO>1 is reached,
    * or when tc_search_duration_mode==4 and FoxO>1 and tc_signal never reached SST_tc_signal,
    * or when tc_search_duration_mode==5 and mTORC1>=1 and cMyc below threshold (tested before),
    * 
    */
   if (( state != FDCselected ) || // apoptosis rate of unselected BCs
       ( tc_search_duration_mode == 0 ) || // apoptosis rate as FDCselected: compete for survival
       ( tc_search_duration_mode < 3 && FDCselected_clock > tc_search_duration ) ||
       // apoptosis rate when the Tfh search period is over
       ( tc_search_duration_mode == 5 &&
	 (FDCselected_clock > tc_search_duration || mTORC1 >= 1.0) ) ||
       ( tc_search_duration_mode == 3 && 
	( ( FoxO_mode == 0 && FDCselected_clock > tc_search_duration && not(hadcontact2Tfh) ) || 
	  FoxO > 1.0 ) ) ||
       // apoptosis rate when FoxO>1; also when no Tfh found for fixed Tfh search period
       /* Is this latter condition needed? 
	* If FoxO_mode=1,2 no, because FoxO starts at low values and grows anyway. 
	* If FoxO_mode==0, FoxO==1 holds until the first Tfh is found, thus, the cells may 
	* never encounter a Tfh and never fulfill FoxO>1. It makes sense to limit this. */
       ( tc_search_duration_mode == 4 && 
	 ( FoxO > 1.0 || ( SSTactive && tc_signal < SST_tc_signal ) ) )
       ) {
     if (mode_of_apoptosis == 1) {
       if ( (state != FDCselected && age > p_apo) ||
	    (state == FDCselected && FDCselected_clock > p_apo4FDCselected) )
	 { err = 0; }
     } else {     
       double p = p_apo;
       if (state == FDCselected) { 
	 /* If tc_search_duration_mode == 0 (a single T-B-interaction),
	  * p_apo4FDCselected can define a limited time during which the BC has the chance
	  * to get Tfh help. If it's set to -1 (=\infty) cells may never die, when they
	  * never pass to the TCcontact state.
	  * If tc_search_duration_mode > 0 && < 3, this time is anyway defined
	  * by tc_search_duration and cells arriving here have a FDCselected_clock
	  * larger than this search period, i.e. the search period is over. 
	  * If tc_search_duration_mode == 3, the time is irrelevant and the criterion
	  * FoxO>1 is used.
	  * In both latter cases, p_apo4FDCselected is used, 
	  * which now denotes the time the unselected BC takes to die 
	  * after the Tfh-search-period or after FoxO>1 was reached.
	  */
	 p = p_apo4FDCselected; 
       }
       if ( (p > 0) && (drandom() < p) ) {
	 /* put p>0 in order to keep the number of drandom calls the same after introduction
	  * of p_apo4FDCselected */
	 err = 0;
       }
     }
   }
   //MS --> Bcl2 exp --> stop apoptosis
   if (err == 0 && !exp_stop_apo) {
     make_apoptosis(time, shape, state==FDCselected);
     // Von der Liste wird die CC erst geloescht, wenn von Makrophagen gefressen
   } else {err=1;}
   // cout<<"errapop="<<err;
   return err;
}
// ============================================================

short cellCC::macrophagocyte(AffinitySpace &shape) {
   if (drandom() < p_mph) {
      /// Philippe : don't get it : n_FDC_encounters is the number of failed interactions ???
      if (n_fdc_encounters < fdc_max_encounters) {
         ++fdc_encounters[n_fdc_encounters];
      } else { ++fdc_encounters[fdc_max_encounters - 1]; }
      shape.rem_cell(sCCapoptosis,pos_ss);
      shape.rem_cell(sCC,pos_ss);
      // total wurde bereits erhoeht als apoptose stattfand
      return 1;
   } else { return 0; }
}
// ============================================================
void cellCC::update_histograms() {
  monitor_pureFDCsearch.update_day();
  monitor_nTCcontacts_selected.update_day();
  monitor_nTCcontacts_deleted.update_day();
  monitor_BsearchTtime_selected.update_day();
  monitor_BsearchTtime_deleted.update_day();
  monitor_BinteractTtime.update_day();
  monitor_TimeBetweenTbindings.update_day();
  monitor_TfhSignalAtSelection.update_day();
  monitor_TfhSignalAtSelection_selected.update_day();
  monitor_TfhSignalAtSelection_deleted.update_day();
  monitor_pMHCatFateDecision_selected.update_day();
  monitor_pMHCatFateDecision_deleted.update_day();
  monitor_TfhSignalSpeed_selected.update_day();
  monitor_TfhSignalSpeed_deleted.update_day();
  monitor_TFHsignalAUC.update_day();
  monitor_TFHsignalAUC_selected.update_day();
  monitor_TFHsignalAUC_deleted.update_day();
  monitor_cMycAtSelection.update_day();
  monitor_cMycAtSelection_selected.update_day();
  monitor_cMycAtSelection_deleted.update_day();
  monitor_mTORatSelection.update_day();
  monitor_mTORatSelection_selected.update_day();
  monitor_mTORatSelection_deleted.update_day();
  //MS
  monitor_nTFRcontacts_SelfCCselected.update_day();
  monitor_nTFRcontacts_SelfCCdeleted.update_day();
  monitor_nTFRcontacts_CCselected.update_day();
  monitor_nTFRcontacts_CCdeleted.update_day();
  monitor_TFHsig_postTFR_CCselected.update_day();
  monitor_TFHsig_postTFR_SelfCCselected.update_day();
  monitor_BinteractTFRtime_self.update_day();
  monitor_BinteractTFRtime_nonself.update_day();
  monitor_TFRsig_select.update_day();
  monitor_TFRsig_apopt.update_day();
  monitor_nMut_nonSelfCD138.update_day();
  monitor_nMut_RedeemedCD138.update_day();
  monitor_nMut_SelfCD138.update_day();
  monitor_time_first_TFRcontact_NonSelfselected.update_day();
  monitor_time_first_TFRcontact_Selfselected.update_day();
}
void cellCC::write_histograms() {
  monitor_pureFDCsearch.write2file();
  monitor_nTCcontacts_selected.write2file();
  monitor_nTCcontacts_deleted.write2file();
  monitor_BsearchTtime_selected.write2file();
  monitor_BsearchTtime_deleted.write2file();
  monitor_BinteractTtime.write2file();
  monitor_TimeBetweenTbindings.write2file();
  monitor_TfhSignalAtSelection.write2file();
  monitor_TfhSignalAtSelection_selected.write2file();
  monitor_TfhSignalAtSelection_deleted.write2file();
  monitor_pMHCatFateDecision_selected.write2file();
  monitor_pMHCatFateDecision_deleted.write2file();
  monitor_TfhSignalSpeed_selected.write2file();
  monitor_TfhSignalSpeed_deleted.write2file();
  monitor_TFHsignalAUC.write2file();
  monitor_TFHsignalAUC_selected.write2file();
  monitor_TFHsignalAUC_deleted.write2file();
  monitor_cMycAtSelection.write2file();
  monitor_cMycAtSelection_selected.write2file();
  monitor_cMycAtSelection_deleted.write2file();
  monitor_mTORatSelection.write2file();
  monitor_mTORatSelection_selected.write2file();
  monitor_mTORatSelection_deleted.write2file();
  //MS
  monitor_nTFRcontacts_SelfCCselected.write2file();
  monitor_nTFRcontacts_SelfCCdeleted.write2file();
  monitor_nTFRcontacts_CCselected.write2file();
  monitor_nTFRcontacts_CCdeleted.write2file();
  monitor_TFHsig_postTFR_CCselected.write2file();
  monitor_TFHsig_postTFR_SelfCCselected.write2file();
  monitor_BinteractTFRtime_self.write2file();
  monitor_BinteractTFRtime_nonself.write2file();
  monitor_TFRsig_select.write2file();
  monitor_TFRsig_apopt.write2file();
  monitor_nMut_nonSelfCD138.write2file();
  monitor_nMut_RedeemedCD138.write2file();
  monitor_nMut_SelfCD138.write2file();
  monitor_time_first_TFRcontact_NonSelfselected.write2file();
  monitor_time_first_TFRcontact_Selfselected.write2file();
}
void cellCC::show_cummulative_ag_collection() {
   ofstream aghisto("ag_endhisto.out");
   aghisto << "! Histogram of number of BCs having collected a number of antigen portions:\n"
           << "! number of antigen portions . occurrence in TC selected BCs ."
	   << " occurrence in BCs after ag-collection\n";
   for (int i = 0; i <= max_n_of_ag_portions; i++) {
      aghisto << i
              << "   " << cummulative_ag_collection_selected[i]
              << "   " << cummulative_ag_collection_all[i]
              << "\n";
   }
   aghisto.close();
}
void cellCC::show_freq_of_TfrCont(double &time, ofstream &freq_nonself, ofstream &freq_self) {
    double fre_av_self=0, fre_sd_self=0,
            fre_av_nonself=0, fre_sd_nonself=0;
    int num_self=0, num_nonself=0;
    for (int i=0; i<=ntfr; i++) {
        fre_av_nonself += double(num_tfr_contacts_nonself[i] * i);
        num_nonself += num_tfr_contacts_nonself[i];
        fre_av_self += double(num_tfr_contacts_self[i] * i);
        num_self += num_tfr_contacts_self[i];
    }
    if (num_nonself>0) {fre_av_nonself/=double(num_nonself);}
    else {num_nonself=0;}
    if (num_self>0) {fre_av_self/=double(num_self);}
    else {num_self=0;}

    for (int i = 0; i <= ntfr; i++) {
       for (int j = 0; j < num_tfr_contacts_nonself[i]; j++) {
          fre_sd_nonself += (double (i) - fre_av_nonself) * (double (i) - fre_av_nonself);
       }
       for (int j = 0; j < num_tfr_contacts_self[i]; j++) {
          fre_sd_self += (double (i) - fre_av_self) * (double (i) - fre_av_self);
       }
    }
    if (num_nonself > 1) {
       fre_sd_nonself /= double (num_nonself - 1);
       fre_sd_nonself = sqrt(fre_sd_nonself);
    } else { fre_sd_nonself = 0; }
    if (num_self > 1) {
       fre_sd_self /= double (num_self - 1);
       fre_sd_self = sqrt(fre_sd_self);
    } else { fre_sd_self = 0; }

    //nonself file
    freq_nonself << time << " " << fre_av_nonself << " " << fre_sd_nonself << endl;
    //self file
    freq_self << time << " " << fre_av_self << " " << fre_sd_self << endl;

    for (int i = 0; i <= ntfr; i++) {
       num_tfr_contacts_nonself[i] = 0;
       num_tfr_contacts_self[i] = 0;
    }
}
// ============================================================
// ============================================================
// ============================================================
// ======================== TC ================================
// ============================================================
// ============================================================
// ============================================================

double cellTC::p_difu = 0.;
double cellTC::p_difu_width = -1.;
double cellTC::persistence = 0.0;
bool cellTC::do_division = 0;
double cellTC::proliferation = -1; // rate to be defined in set-statics
double cellTC::meancycle = 10.0; // hours
double cellTC::cyclewidth = 0.5; // fraction
int cellTC::Ndivisions = 1; // number
double cellTC::max_distance2divide = 0; // micron
double cellTC::dTrepolarise = 0;
double cellTC::dTsignal = 0;
double cellTC::pMHCsignal_adaptation_weight = 1;
double cellTC::K_signal_intensity_pMHC_default = 0;

cellTC::cellTC() {
   // cout<<"in cellTC default constructor ...";
   state = TCnormal;
   cell_cycle_clock = 0;
   cell_cycle_duration = 0;  // has to be defined later
   set_p_move();
   n_CC_nn = 0;
   CC_nn = new long[6];
   CC_affinity = new double[6];
   ini_polarisation_state();
   set_changed_nn(0);
   // cout<<"end of cellTC default constructor.\n";
}
// ============================================================

//cellTC::cellTC(const cellTC &x) : frag_cell(x) { // MMH: activate this, it works
cellTC::cellTC(const cellTC &x) : cell(x) {
   // cout<<"in cellTC::cellTC(cellTC&) ...";
   operator =(x);
   // cout<<"end of cellTC::cellTC(cellTC&)\n";
}
// ============================================================

cellTC&cellTC::operator =(const cellTC &x) {
   /* This generates a real copy of the TC. Don't use this for TC division! */
   // cout<<"in cellTC::=(cellTC&) ... ";
   cell::operator =(x);
   state = x.state;
   pos_ss = x.pos_ss;
   cell_cycle_clock = x.cell_cycle_clock;
   cell_cycle_duration = x.cell_cycle_duration;
   set_p_move();
   n_CC_nn = x.n_CC_nn;
   for (short a = 0; a < n_CC_nn; a++) {
      CC_nn[a] = x.CC_nn[a];
      CC_affinity[a] = x.CC_affinity[a];
   }
   able2repolarise = x.able2repolarise;
   able2signal = x.able2signal;
   CurrentSignalTarget = x.CurrentSignalTarget;
   K_signal_intensity_pMHC = x.K_signal_intensity_pMHC;
   TsinceBCpolarisation = x.TsinceBCpolarisation;
   changed_nn = x.changed_nn;
   // cout<<"end of cellTC::=(cellTC&)\n";
   return *this;
}
// ============================================================

cellTC::~cellTC() {
   // cout<<"in ~cellTC()...\n";
   delete[] CC_nn;
   delete[] CC_affinity;
}
// ============================================================

void cellTC::set_statics(const Parameter &par, space &l, ofstream &ana) {
   if (par.Value.TC_persistence > 60. * par.Value.deltat) {
      persistence = 60. * par.Value.deltat / par.Value.TC_persistence;
   } else {
      persistence = 1.;   // i.e. change polarity in every time step!
   }
   ana << "  TC polarity persistence = " << par.Value.TC_persistence
       << " i.e. p_{change-polarity} = " << persistence << "\n";

   p_difu = 60. * par.Value.v_TC * par.Value.deltat / par.Value.dx;
   if (p_difu > 0.5) {
      cout << "TC-motility (p="
           << p_difu
           << ") is too large for dx and dt !!!\n";
      exit(1);
   }
   p_difu_width = par.Value.v_TC_width * p_difu;
   ana << "  TC movement probability = "
       << p_difu
       << "\n";
   //
   double tmpdt = par.Value.TC_persistence;
   if (tmpdt < 60. * par.Value.deltat) {
      tmpdt = 60. * par.Value.deltat;
      ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
   }
   ana << "  Expected effective diffusion constant for TC = "
       << 60. * par.Value.v_TC * par.Value.v_TC * tmpdt / double (l.dim2)
       << " microns^2/hr = "
       << par.Value.v_TC * par.Value.v_TC * tmpdt / double (l.dim2)
       << " microns^2/min\n";
   //
   dTrepolarise = par.Value.dTrepolarise / 60.; // transform minutes in par-file to hours
   dTsignal = par.Value.dTsignal / 60.; // transform minutes in par-file to hours
   pMHCsignal_adaptation_weight = 1. / par.Value.TFHpMHCadaptation_smoothness;
   K_signal_intensity_pMHC_default = cellCC::TFHsignal_delivery_KpMHC;
   cout << "K_signal_intensity_pMHC_default = " << K_signal_intensity_pMHC_default << "\n";
   if (par.Value.TFHsignal_delivery_mode > 0 && K_signal_intensity_pMHC_default < 0) {
     cerr << "K_signal_intensity_pMHC_default is negative in cellTC::set_statics(...).\n"
	  << "Exit.\n";
     exit(1);
   }
   ana << "  TC can repolarise to a different BC after " << par.Value.dTrepolarise 
       << " minutes of polarisation.\n"
       << "  TC start providing signal after " << par.Value.dTsignal 
       << " minutes of polarisation.\n";
   //
   do_division = par.Value.do_TC_division;
   proliferation = par.Value.TC_doubling * par.Value.deltat;  // rate times timestep
   meancycle = par.Value.TC_meancycle;
   cyclewidth = par.Value.TC_cyclewidth;
   Ndivisions = par.Value.TC_Ndivisions;
   max_distance2divide = par.Value.dx_TC;
}
// ============================================================

void cellTC::ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape) {
   /* Diese ini-Routine geht von einer bereits vordefinierten Zelle
    * aus, die in dem Objekt gespeichert ist. Es werden an dem Objekt
    * die Aenderungen vorgenommen, die es zu einem neuen Zellobjekt
    * machen. Ausserdem werden lattice und shapespace aktualisiert.
    * Achtung: state wird nicht gesetzt!
    */
   // Change lattice index in newTC to i
   index = i;
   born_index = i;
   born_time = t;
   l.get_koord(i,last_position);
   // Actualize lattice point
   l.set_knot(i,TC,li);
   // l.cellknot[i].cell=TC;
   // Write on cell_list und save list-index on the lattice
   // l.knot[i].listi=li;
   // define pos_ss randomly from any of the antigens:
   pos_ss = shape.get_Antigen();
   // Actualize shape space statistics
   // shape.add_cell(sTC,pos_ss);
   // shape.add_cell(total,pos_ss);
   // get initial polarity vector
   l.get_random_direction(polarity);
   // Set default value for the signalling intensity:
   K_signal_intensity_pMHC = K_signal_intensity_pMHC_default;
   // cout<<"polar-vector=("<<polarity[0]<<","<<polarity[1]<<")\n";
}
void cellTC::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
  if (p_difu_width > 0) {
    p_move = get_positive_sample_from_normal(p_difu, p_difu_width);
  }
}
short cellTC::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   set_polarity_velocity(persistence,l,s);
   if (writethis2track == polarisation) {
      double tmpi[l.dim];
      l.get_koord(index,tmpi);
      td.Write_movement(trackno,TC,time,tmpi,polarity,writethis2track);
      writethis2track = trackini;
   }
   short err = do_diffuse(TC,li,l);
   return err;
}

// ============================================================

void cellTC::set_CXCR5expression() {
   if (responsive2signal[CXCL13] && (cellCC::p_CXCR5down > 0)
           && (drandom() < cellCC::p_CXCR5down)) {
      responsive2signal[CXCL13] = false;
   }
}
void cellTC::resensitise4CXCL13(sigs &s) {
   if ((s.signal_use[CXCL13] == 1)  // CXCL13 has to be used at all
       && not (responsive2signal[CXCL13])  // otherwise its sensitive anyway
       && (CXCL13recrit > 0)  // re-sensitisation has to be on
       && s.undercritical_signal(index,CXCL13,CXCL13recrit)
       // CXCL13 has to be below the threshold value for re-sensitisation
       ) {
      responsive2signal[CXCL13] = true;
      // cout<<"resensitised CC with index "<<index<<" for CXCL13.\n";
   }
}

void cellTC::ini_polarisation_state() {
  able2repolarise = true;
  able2signal = false;
  CurrentSignalTarget = -1;
  TsinceBCpolarisation = 0;
}
void cellTC::make_tc_cc_link(const long &index, const double &pMHCpresentation) {
  // This is the version in case of cellCC::collectFDCsignals==true
  /* Switch state to "in contact":
   * if first BC bound -> allow for MTOC polarisation: */
  if (state != TC_CCcontact) { ini_polarisation_state(); }
  state = TC_CCcontact;
  // save the index of the CC on the lattice
  CC_nn[n_CC_nn] = index;
  CC_affinity[n_CC_nn] = pMHCpresentation;
  // increase the number of used list members
  ++n_CC_nn;
  set_changed_nn(1);
}
void cellTC::make_tc_cc_link(const long &index,
			     const long &CCpos,
			     int ag_index,
			     AffinitySpace &shape,
			     const bool &highag) {
  // This is the version in case of cellCC::collectFDCsignals==false
  /* Switch state to "in contact":
   * if first BC bound -> allow for MTOC polarisation: */
  if (state != TC_CCcontact) { ini_polarisation_state(); }
  state = TC_CCcontact;
  // save the index of the CC on the lattice
  CC_nn[n_CC_nn] = index;
  if (highag) {
    CC_affinity[n_CC_nn] = 1.0;
  } else {
    if (ag_index == -1) {
      CC_affinity[n_CC_nn] = shape.affinity(CCpos, pos_ss);
    } else {
      CC_affinity[n_CC_nn] = shape.affinity(CCpos, shape.get_Antigen(ag_index));
    }
  }
  // increase the number of used list members
  ++n_CC_nn;
  set_changed_nn(1);
}
double cellTC::adapt_K_signal_intensity(int last_pMHC) {
  /* Adapts the K-value for the signalling intensity used in cellCC-class.
   * Takes the old K-value and adapts it to last_pMHC, which is the value
   * seen on the last detached B cell.
   * Adaption is dampened with pMHCsignal_adaptation_weight, which is the inverse 
   * number of B cell contacts until adaptation, to avoid to nervous changes.
   */ 
  /*
  cerr << "(index=" << index << ", id=" << id << ", K=" << K_signal_intensity_pMHC
       << ", pMHCnew=" << last_pMHC << ", Knew=";
  */
  K_signal_intensity_pMHC = ( (1. - pMHCsignal_adaptation_weight) * K_signal_intensity_pMHC + 
			      pMHCsignal_adaptation_weight * double(last_pMHC) );
  //cerr << K_signal_intensity_pMHC << "); ";
  return K_signal_intensity_pMHC;
}

short cellTC::get_n_boundCC() {
//    cerr <<"n bound CC " <<n_CC_nn_TFR<<"\n";
    for (short jj=0;jj<n_CC_nn;jj++) {
//        cerr<<jj<<"  ind  "<<CC_nn_TFR[jj]<<"\n";
    }
    return n_CC_nn;
}

void cellTC::liberateCC(const long &index) {
   short where = -1;
   for (short w = 0; w < n_CC_nn; w++) {
      if (CC_nn[w] == index) {
         where = w;
      }
   }
   if (where < n_CC_nn - 1) {
      CC_nn[where] = CC_nn[n_CC_nn - 1];
      CC_affinity[where] = CC_affinity[n_CC_nn - 1];
   }
   if (where==-1) {cerr<<"wrong cc index"<<endl; exit(1);}
   --n_CC_nn;
   // If this was the last contact the TC had to CC, return to normal state.
   if (n_CC_nn == 0) { state = TCnormal; }
   // If the signal target BCs is deliberated, a new target will have to be defined.
   if (index == CurrentSignalTarget) { ini_polarisation_state(); }
   set_changed_nn(1);
}
void cellTC::set_changed_nn(short x) {
   changed_nn = x;
}
void cellTC::evolve_polarisation_state(double& dt) {
  // ... concerning the polarisation state for signalling to BCs.
  TsinceBCpolarisation += dt;
  if (TsinceBCpolarisation > dTrepolarise) { able2repolarise = true; }
  if (TsinceBCpolarisation > dTsignal) { able2signal = true; }
}
void cellTC::set_polarity(space &l) {
  /* Sets the polarity not for movement but for secretion of signals to specific CC
   * Find the bound CC with highest affinity:
   */
  if (changed_nn == 1 && able2repolarise) {
    /* Any BC arriving here has bound at least one portion of antigen!
     * Not true anymore since 2018: Depending on the setting,
     * no-Ag-BCs can have a small chance to bind TC.
     */
    double haff = 0.;
    short a = 0;
    for (short i = 0; i < n_CC_nn; i++) {
      // cout<<"aff="<<CC_affinity[i]<<", ";
      if (CC_affinity[i] > haff) {
	a = i;
	haff = CC_affinity[i];
      } 
      /* How to treat equal affinities of different BCs? 
       * Do random choice? Not needed:
       * As is means: The BC with equal affinity, which bounds first, keeps polarity.
       */
    }
    /*if (CurrentSignalTarget == CC_nn[a]) { cerr << CurrentSignalTarget << "! "; }
      else { if (CurrentSignalTarget != -1) { cerr << CurrentSignalTarget << "? "; } } */
    // Reset polarity if either first bound BC or changed target BC:
    if (CurrentSignalTarget == -1 || CurrentSignalTarget != CC_nn[a]) {
      // Set polarity vector such that it points to this neighbour, i.e. CC_nn[a]-index
      l.get_diff_vector(index,CC_nn[a],polarity);
      able2repolarise = false; // able2signal == false anyway
      TsinceBCpolarisation = 0;
      CurrentSignalTarget = CC_nn[a];
    }
    set_changed_nn(0);
  }
  // else keep old polarity.
}
double cellTC::set_cell_cycle_duration() {
   double arg;
   double cyclelength = 3.0 * meancycle;
   while (cyclelength <= 0. || cyclelength >= 2. * meancycle) {
      // This samples from a normal distribution with mean meancycle and width cyclewidth
      arg = (2. * drandom() - 1.);
      cyclelength = meancycle + sqrt(2.) * cyclewidth * meancycle * inverse_erf(arg);
   }
   return cyclelength;
}
void cellTC::reset_cycle_times() {
   cell_cycle_duration = 0;
   cell_cycle_clock = 0;
}
void cellTC::ask_enter_cell_cycle() {
   /* Takes a decision whether this cell shall enter the cell cycle.
    * This routine is called by TC in state==TCnormal exclusively.
    * If the cell cycle is enterred, the cycle_duration is set,
    * the clock is initialised and the state is set to state==TCdivide.
    */
   // More mechanistic models of TC division might be used here.

   // For the moment it is just homeostatic division with a constant rate:
   if (drandom() < proliferation) {
      cell_cycle_duration = set_cell_cycle_duration();
      cout<<cell_cycle_duration<<"::";
      cell_cycle_clock = 0;
      state = TCdivide;
   }
}
bool cellTC::progress_cell_cycle(double &dt) {
   /* Adds a time step to the cycle clock.
    * If the end of the cell cycle is reached,
    * true is returned.
    */
   bool endofcycle = false;
   // Add time step deltat to the internal clock:
   cell_cycle_clock += dt;
   if (cell_cycle_clock > cell_cycle_duration) { endofcycle = true; }
   return endofcycle;
}
long cellTC::ask_mitosis(long * pp, space &l) {
   /* Because the same routine is used to determine a target position of a new cell,
    * <proliferation> is transported through, but ignored in the <find_mitosis_place(...)>
    * because the state of the TC is always <TCdivide> when this routine is called.
    * The place for the target cell is returned or -1 if none was found.
    */
   //
   //  if (volume==1 && target_volume==1) {
   return find_mitosis_place(proliferation,state == TCdivide,max_distance2divide,pp,l);
   /* state==TCdivide is used to force division irrespective of probabilities.
    * find_mitosis_place allows fpr two modes of division:
    * (i) Division based on a division rate. Here, a TC in the state TCnormal would call this
    *     routine and would then induce division with a probability as given by <proliferation>.
    * (ii) Division based on a cell cycle. Here, the decision for division is taken elsewhere.
    *      When the TC is mature to divide (to do mitosis), the state is set to TCdivide and
    *  only then this routine is called, which then imposes division.
    * Here, only (ii) was programmed and (i) is just included for symmetry with cellCB class.
    */
   /* pp is used for statistics of BC divisions, now being disturbed by TC division events.
    * Think of how this is needed and eventually suppress TC division counts in pp.
    * Actually, no because this is meant to count how often inhibition by lack of space is
    *occurring.
    * ### check whether there is any use of <pp> and <prolog> specific to BCs.
    * It is not used now anyway, pp is updated but only written for outputfiles==0. Not urgent.
    * See also cellman::divide_TC()
    */
   // }
   // Note that this object is used as single node.
   // Upon extension the following might be activated:
   /*
    * else { // case of more fragment object
    * // probabilistic decision if proliferation is done
    * if ( volume>0.9*target_volume
    * && ((state==tc_divide) || (drandom() < proliferation)) )
    *  return 0;
    * // Proliferation is allowed if the total volume is near the target-volume of a cell!
    * return 1;
    * }
    */
}


//MSchips
// ============================================================
// ============================================================
// ============================================================
// ===================== TFRcells =============================
// ============================================================
// ============================================================
// ============================================================

double cellTFR::p_difu = 0.;
double cellTFR::p_difu_width = -1.;
double cellTFR::persistence = 0.0;
bool cellTFR::do_division = 0;
double cellTFR::proliferation = -1; // rate to be defined in set-statics
double cellTFR::meancycle = 10.0; // hours
double cellTFR::cyclewidth = 0.5; // fraction
int cellTFR::Ndivisions = 1; // number
double cellTFR::max_distance2divide = 0; // micron
double cellTFR::dTrepolarise = 0;
double cellTFR::dTsignal = 0;
double cellTFR::pMHCsignal_adaptation_weight = 1;
double cellTFR::K_signal_intensity_pMHC_default = 0;


cellTFR::cellTFR() {
//    cerr<<"in cellTFR default constructor ...";
   state = TFRnormal;
   cell_cycle_clock = 0;
   cell_cycle_duration = 0;  // has to be defined later
   set_p_move();
   n_CC_nn_TFR = 0;
   CC_nn_TFR = std::array<long, 6> { -1, -1, -1, -1, -1, -1 };
   ini_polarisation_state();
   tmp_blocked=0;
   to_delete=0;
   last_cc_ind=-1;
   cc_ind=-1;
//   set_changed_nn(0);
//    cerr<<"end of cellTFR default constructor.\n";
}
// ============================================================

cellTFR::cellTFR(const cellTFR &x) : cell(x) {
   // cout<<"in cellTFR::cellTFR(cellTFR&) ...";
   operator =(x);
   // cout<<"end of cellTFR::cellTFR(cellTFR&)\n";
}
// ============================================================

cellTFR&cellTFR::operator =(const cellTFR &x) {
   // cout<<"in cellTFR::=(cellTFR&) ... ";
   cell::operator =(x);
   state = x.state;
   pos_ss = x.pos_ss;
   cell_cycle_clock = x.cell_cycle_clock;
   cell_cycle_duration = x.cell_cycle_duration;
   set_p_move();
   n_CC_nn_TFR = x.n_CC_nn_TFR;
   for (short a = 0; a < n_CC_nn_TFR; a++) {
      CC_nn_TFR[a] = x.CC_nn_TFR[a];
   }
   able2repolarise = x.able2repolarise;
   able2signal = x.able2signal;
   CurrentSignalTarget = x.CurrentSignalTarget;
   K_signal_intensity_pMHC = x.K_signal_intensity_pMHC;
   TsinceBCpolarisation = x.TsinceBCpolarisation;
   tmp_blocked=x.tmp_blocked;
   to_delete=x.to_delete;
   last_cc_ind=x.last_cc_ind;
   cc_ind=x.cc_ind;
//   changed_nn = x.changed_nn;
   // cout<<"end of cellTFR::=(cellTFR&)\n";
   return *this;
}
// ============================================================

cellTFR::~cellTFR() {
//    cout<<"in ~cellTFR()...\n";
   // <sb>
   //delete[] CC_nn_TFR;
   //delete[] CC_affinity;
   // </sb>
}
// ============================================================

void cellTFR::set_statics(const Parameter &par, space &l, ofstream &ana) {
   if (par.Value.TC_persistence > 60. * par.Value.deltat) {
      persistence = 60. * par.Value.deltat / par.Value.TC_persistence;
   } else {
      persistence = 1.;   // i.e. change polarity in every time step!
   }
   ana << "  TFR polarity persistence = " << par.Value.TC_persistence
       << " i.e. p_{change-polarity} = " << persistence << "\n";

   p_difu = 60. * par.Value.v_TC * par.Value.deltat / par.Value.dx;
   if (p_difu > 0.5) {
      cout << "TFR-motility (p="
           << p_difu
           << ") is too large for dx and dt !!!\n";
      exit(1);
   }
   p_difu_width = par.Value.v_TC_width * p_difu;
   ana << "  TC movement probability = "
       << p_difu
       << "\n";
   //
   double tmpdt = par.Value.TC_persistence;
   if (tmpdt < 60. * par.Value.deltat) {
      tmpdt = 60. * par.Value.deltat;
      ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
   }
   ana << "  Expected effective diffusion constant for TFR = "
       << 60. * par.Value.v_TC * par.Value.v_TC * tmpdt / double (l.dim2)
       << " microns^2/hr = "
       << par.Value.v_TC * par.Value.v_TC * tmpdt / double (l.dim2)
       << " microns^2/min\n";
   //
   dTrepolarise = par.Value.dTrepolarise / 60.; // transform minutes in par-file to hours
   dTsignal = par.Value.dTsignal / 60.; // transform minutes in par-file to hours
   pMHCsignal_adaptation_weight = 1. / par.Value.TFHpMHCadaptation_smoothness;
   K_signal_intensity_pMHC_default = cellCC::TFHsignal_delivery_KpMHC;
   cout << "K_signal_intensity_pMHC_default = " << K_signal_intensity_pMHC_default << "\n";
   if (par.Value.TFHsignal_delivery_mode > 0 && K_signal_intensity_pMHC_default < 0) {
     cerr << "K_signal_intensity_pMHC_default is negative in cellTFR::set_statics(...).\n"
          << "Exit.\n";
     exit(1);
   }
   ana << "  TFR can repolarise to a different BC after " << par.Value.dTrepolarise
       << " minutes of polarisation.\n"
       << "  TFR start providing signal after " << par.Value.dTsignal
       << " minutes of polarisation.\n";
   //
   do_division = par.Value.do_TC_division;
   proliferation = par.Value.TC_doubling * par.Value.deltat;  // rate times timestep
   meancycle = par.Value.TC_meancycle;
   cyclewidth = par.Value.TC_cyclewidth;
   Ndivisions = par.Value.TC_Ndivisions;
   max_distance2divide = par.Value.dx_TC;
}
// ============================================================

void cellTFR::ini(const long &i, const long &li, const double &t, space &l, AffinitySpace &shape) {
   /* Diese ini-Routine geht von einer bereits vordefinierten Zelle
    * aus, die in dem Objekt gespeichert ist. Es werden an dem Objekt
    * die Aenderungen vorgenommen, die es zu einem neuen Zellobjekt
    * machen. Ausserdem werden lattice und shapespace aktualisiert.
    * Achtung: state wird nicht gesetzt!
    */
   // Change lattice index in newTFR to i
   index = i;
   born_index = i;
   born_time = t;
   tmp_blocked=0;
   to_delete=0;
   l.get_koord(i,last_position);
   // Actualize lattice point
   l.set_knot(i,TFR,li);

   // define pos_ss randomly from any of the antigens:
   pos_ss = shape.get_Antigen();

   // get initial polarity vector
   l.get_random_direction(polarity);
   // Set default value for the signalling intensity:
   K_signal_intensity_pMHC = K_signal_intensity_pMHC_default;
   // cout<<"polar-vector=("<<polarity[0]<<","<<polarity[1]<<")\n";
}
void cellTFR::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
  if (p_difu_width > 0) {
    p_move = get_positive_sample_from_normal(p_difu, p_difu_width);
  }
}
short cellTFR::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   set_polarity_velocity(persistence,l,s);
   if (writethis2track == polarisation) {
      double tmpi[l.dim];
      l.get_koord(index,tmpi);
      ///Warning: TRACKing needs to be checked (TODO!)
      td.Write_movement(trackno,TFR,time,tmpi,polarity,writethis2track);
      writethis2track = trackini;
   }
   short err = do_diffuse(TFR,li,l);
   return err;
}
// ============================================================
void cellTFR::set_CXCR5expression() {
   if (responsive2signal[CXCL13] && (cellCC::p_CXCR5down > 0)
           && (drandom() < cellCC::p_CXCR5down)) {
      responsive2signal[CXCL13] = false;
   }
}
void cellTFR::resensitise4CXCL13(sigs &s) {
if ((s.signal_use[CXCL13] == 1)  // CXCL13 has to be used at all
        && not (responsive2signal[CXCL13])  // otherwise its sensitive anyway
        && (CXCL13recrit > 0)  // re-sensitisation has to be on
        && s.undercritical_signal(index,CXCL13,CXCL13recrit)
        // CXCL13 has to be below the threshold value for re-sensitisation
        ) {
    responsive2signal[CXCL13] = true;
    // cout<<"resensitised CC with index "<<index<<" for CXCL13.\n";
}
}
void cellTFR::set_CXCR4expression() {
   if (responsive2signal[CXCL12] && (cellCB::p_CXCR4down > 0)
       && (drandom() < cellCB::p_CXCR4down)) {
      responsive2signal[CXCL12] = false;
   }
}
void cellTFR::resensitise4CXCL12(sigs &s) {
if ((s.signal_use[CXCL12] == 1)  // CXCL12 has to be used at all
        && not (responsive2signal[CXCL12])  // otherwise its sensitive anyway
        && (CXCL12recrit > 0)  // re-sensitisation has to be on
        && s.undercritical_signal(index,CXCL12,CXCL12recrit)
        // CXCL12 has to be below the threshold value for re-sensitisation
        ) {
    responsive2signal[CXCL12] = true;
    // cout<<"resensitised CC with index "<<index<<" for CXCL12.\n";
}
}
// ============================================================

void cellTFR::ini_polarisation_state() {
  able2repolarise = true;
  able2signal = false;
  CurrentSignalTarget = -1;
  TsinceBCpolarisation = 0;
}
void cellTFR::TFR_bind_CC(const double& time, cellCC &cccell, space &l, AffinitySpace &shape) {
//    state = TFR_CCcontact;
    make_tfr_cc_link(cccell.index);

}
void cellTFR::make_tfr_cc_link(const long &index) {
    if (state != TFRnormal) { cerr<<"wrong state of tfr, already bound "<<state;exit(1); }
    state = TFR_CCcontact;
    // save the index of the CC on the lattice
    CC_nn_TFR[n_CC_nn_TFR] = index;
    // increase the number of used list members
    //cout << n_CC_nn_TFR << " ";
    ++n_CC_nn_TFR;
    tmp_blocked=0;
    //cout << n_CC_nn_TFR << endl;
    ///polarisation is not needed
//    set_changed_nn(1);
}

short cellTFR::get_n_boundCC() {
    /** @brief: returns n. of CC bound to Tfr
     **/
    for (short jj=0;jj<n_CC_nn_TFR;jj++) {
//        cerr<<jj<<"  ind  "<<CC_nn_TFR[jj]<<"\n";
    }
    return n_CC_nn_TFR;
}
void cellTFR::liberateCC(const long &index) {
  //  cerr << "in cellTFR::liberateCC() n_CC_nn_TFR=" << n_CC_nn_TFR<< endl;
   short where = -1;
   for (short w = 0; w < n_CC_nn_TFR; w++) {
     //     cerr << "CC_nn_TFR[" << w << "]=" << CC_nn_TFR[w] << ";  ";
     if (CC_nn_TFR[w] == index) {
       where = w;
     }
   }
//   cerr<<"tfr at ind"<< index << "is liberating cell with ind "<<CC_nn_TFR[0]<<endl;
   if (where==-1) {
       cerr<<" wrong liberate at index " <<index
                      <<" but index recorded:"<<CC_nn_TFR[0]<<endl;
       exit(1);
   } //else {cerr<<"tfr liberated wh ind:"<<CC_nn_TFR[where]<<endl;}
   if (where < n_CC_nn_TFR - 1) {
      CC_nn_TFR[where] = CC_nn_TFR[n_CC_nn_TFR - 1];
   }
   --n_CC_nn_TFR;
   // If this was the last contact the TFR had to CC, return to normal state.
   if (n_CC_nn_TFR == 0) { state = TFRnormal; } else {cerr<<"wrong update of nccnnTFR\n"; exit(1);}
   ///polarisation is not needed for the time being
   // If the signal target BCs is deliberated, a new target will have to be defined.
//   if (index == CurrentSignalTarget) { ini_polarisation_state(); }
//   set_changed_nn(1);
}
double cellTFR::set_cell_cycle_duration() {
   double arg;
   double cyclelength = 3.0 * meancycle;
   while (cyclelength <= 0. || cyclelength >= 2. * meancycle) {
      // This samples from a normal distribution with mean meancycle and width cyclewidth
      arg = (2. * drandom() - 1.);
      cyclelength = meancycle + sqrt(2.) * cyclewidth * meancycle * inverse_erf(arg);
   }
   return cyclelength;
}
void cellTFR::reset_cycle_times() {
   cell_cycle_duration = 0;
   cell_cycle_clock = 0;
}
void cellTFR::ask_enter_cell_cycle() {
   /* Takes a decision whether this cell shall enter the cell cycle.
    * This routine is called by TFR in state==TFRnormal exclusively.
    * If the cell cycle is enterred, the cycle_duration is set,
    * the clock is initialised and the state is set to state==TFRdivide.
    */
   // More mechanistic models of TFR division might be used here.

   // For the moment it is just homeostatic division with a constant rate:
   if (drandom() < proliferation) {
      cell_cycle_duration = set_cell_cycle_duration();
      cout<<cell_cycle_duration<<"::";
      cell_cycle_clock = 0;
      state = TFRdivide;
   }
}
bool cellTFR::progress_cell_cycle(double &dt) {
   /* Adds a time step to the cycle clock.
    * If the end of the cell cycle is reached,
    * true is returned.
    */
   bool endofcycle = false;
   // Add time step deltat to the internal clock:
   cell_cycle_clock += dt;
   if (cell_cycle_clock > cell_cycle_duration) { endofcycle = true; }
   return endofcycle;
}
long cellTFR::ask_mitosis(long * pp, space &l) {
   /* Because the same routine is used to determine a target position of a new cell,
    * <proliferation> is transported through, but ignored in the <find_mitosis_place(...)>
    * because the state of the TFR is always <TFRdivide> when this routine is called.
    * The place for the target cell is returned or -1 if none was found.
    */
   //
   //  if (volume==1 && target_volume==1) {
   return find_mitosis_place(proliferation,state == TFRdivide,max_distance2divide,pp,l);
   /* state==TFRdivide is used to force division irrespective of probabilities.
    * find_mitosis_place allows fpr two modes of division:
    * (i) Division based on a division rate. Here, a TFR in the state TFRnormal would call this
    *     routine and would then induce division with a probability as given by <proliferation>.
    * (ii) Division based on a cell cycle. Here, the decision for division is taken elsewhere.
    *      When the TFR is mature to divide (to do mitosis), the state is set to TFRdivide and
    *  only then this routine is called, which then imposes division.
    * Here, only (ii) was programmed and (i) is just included for symmetry with cellCB class.
    */
}


// ============================================================



// ============================================================
// ============================================================
// =================== FDC ======================================
// ============================================================
// ============================================================

cellFDC::cellFDC() {
   state = none;
   //  cerr<<"in cellFDC(): "<<ab_resolution<<"; ";
   //  cerr<<"in cellFDC(): "<<FDCmaxFrags<<"\n";
   if (FDCmaxFrags > FRAGMENT_STEP) {
      cerr << "FDCmaxFrags " << FDCmaxFrags << " larger than FRAGMENT_STEP "
           << FRAGMENT_STEP << " in cellFDC constructor. Abort.\n";
      exit(1);
   }

   // ag_fraction. = 0;

   // get memory for antigen_amount (double **)
   antigen_amount = new double*[FDCmaxFrags];
   for (int f = 0; f < FDCmaxFrags; f++) {
      antigen_amount[f] = new double[n_Antigen_max];
   }
   /* ##### Note that this way of initialisation requires that the number
    *       antigens provided in the parameter file is the highest number
    *   that will be achieved in the course of the reaction. Thus,
    *   all antigens have to be installed everywhere, just the
    *   antigen amount should be set zero in the beginning if the
    *   antigen is supposed to be added later.
    */

   // get memory for ic_amount (double ***)
   ic_amount = new double**[FDCmaxFrags];
   for (int f = 0; f < FDCmaxFrags; f++) {
      ic_amount[f] = new double*[ab_resolution + 1];
      for (int b = 0; b <= ab_resolution; b++) {
         ic_amount[f][b] = new double[n_Antigen_max];
      }
   }

   // initialise all this
   for (int f = 0; f < FDCmaxFrags; f++) {
      for (int a = 0; a < n_Antigen_max; a++) {
         antigen_amount[f][a] = 0.;
         for (int b = 0; b <= ab_resolution; b++) {
            ic_amount[f][b][a] = 0.;
         }
      }
      // ic_amount[f][5][0]=5.e-08; // in M
   }
}
// ============================================================

cellFDC::cellFDC(const cellFDC &x) : frag_cell(x) {
   ///Philippe danger
   antigen_amount = new double*[FDCmaxFrags];
   for (int f = 0; f < FDCmaxFrags; f++) {
      antigen_amount[f] = new double[n_Antigen_max];
   }
   ic_amount = new double**[FDCmaxFrags];
   for (int f = 0; f < FDCmaxFrags; f++) {
      ic_amount[f] = new double*[ab_resolution + 1];
      for (int b = 0; b <= ab_resolution; b++) {
         ic_amount[f][b] = new double[n_Antigen_max];
      }
   }
   operator =(x);
   for (int f = 0; f < FDCmaxFrags; f++) {
      for (int a = 0; a < n_Antigen_max; a++) {
         antigen_amount[f][a] = x.antigen_amount[f][a];
         for (int b = 0; b <= ab_resolution; b++) {
            ic_amount[f][b][a] = x.ic_amount[f][b][a];
         }
      }
   }
}
cellFDC&cellFDC::operator =(const cellFDC &x) {
   frag_cell::operator =(x);
   state = x.state;
   // cout<<"in cellFDC::operator=(const cellFDC): "<<ab_resolution<<"\n";
   for (int f = 0; f < FDCmaxFrags; f++) {
      for (int a = 0; a < n_Antigen_max; a++) {
         antigen_amount[f][a] = x.antigen_amount[f][a];
         for (int b = 0; b <= ab_resolution; b++) {
            ic_amount[f][b][a] = x.ic_amount[f][b][a];
         }
      }
   }
   return *this;
}
// ============================================================

cellFDC::~cellFDC() {
   // cout<<"in ~cellFDC() ...\n";
   delete[] antigen_amount;  // ###multiAg now this is 2D
   delete[] ic_amount;  // ### but this one is 2D !? ###multiAg now even 3D
}
// ============================================================

vector<double> cellFDC::ini_ag_fraction() {
   vector<double> v;
   v.reserve(100);
   v.clear();
   v.push_back(-1);
   return v;
}
double cellFDC::p_mksignal = 0.;
double cellFDC::p_mkCXCL13 = 0.;
double cellFDC::p_mkSEMA4D = 0.;
short cellFDC::vesicle = 1;
unsigned short cellFDC::use_antigen = 0;
double cellFDC::antigen_amount_per_FDC = 1.e+08;
double cellFDC::antigen_saturation = 1.;
double cellFDC::ag_threshold = 1.e-08; // Mol
double cellFDC::ic_k_on = 3600. * 1.e+06 * cellFDC::ag_threshold; // /(Mol hour)
double cellFDC::ic_k_off = 3600. * 1.e-03; // /hour
long cellFDC::ab_sign_errors = 0;
long cellFDC::ag_sign_errors = 0;
long cellFDC::ic_sign_errors = 0;
long cellFDC::ic_calculations = 0;
int cellFDC::DendriteLength = 0;
int cellFDC::FDCmaxFrags = 0;
int cellFDC::n_Antigen_max = 1;
int cellFDC::n_Antigen_dim_factor = 10;
vector<double> cellFDC::ag_fraction = cellFDC::ini_ag_fraction();
short cellFDC::ag_distribution_mode = 0;
short cellFDC::ag_detection_mode = 0;
bool cellFDC::BCreduceFDCag = true;

void cellFDC::set_statics(const Parameter &par, space &l, ofstream &ana) {
   if (par.Value.mksignal <= 0.) {
      p_mksignal = 0.;
   } else { p_mksignal = par.Value.mksignal * par.Value.deltat; }
   ana << "FDC differentiation signal production = " << p_mksignal << "\n";

   if (par.Value.mkCXCL13 <= 0.) { p_mkCXCL13 = 0.; } else {
      p_mkCXCL13 = par.Value.mkCXCL13                   // mol/(cell l hr)
                   * par.Value.deltat               // hr
                   * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15   // l
                   * par.N_A;                         // /mol
   }
   ana << "FDC CXCL13 production (molecules) = " << p_mkCXCL13 << "\n";

   if (par.Value.mk_SEMA4D <= 0.) { p_mkSEMA4D = 0.; } else {
      p_mkSEMA4D = par.Value.mk_SEMA4D                  // mol/(cell l hr)
                   * par.Value.deltat               // hr
                   * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15   // l
                   * par.N_A;                       // /mol
   }
   ana << "FDC SEMA4D production (molecules) = " << p_mkSEMA4D << "\n";

   vesicle = par.Value.FDCvesicle;
   ana << "FDC signal production mode = " << vesicle << "\n";

   // Number of fragments per FDC (note that FDC are not supposed to divide or grow):
   // dendrite length divided per lattice constant is the number of fragments per dendrite
   DendriteLength = int (par.Value.FDClength / par.Value.dx);
   ana << "FDC dendrites have a length of " << DendriteLength << " lattice constants.\n";
   // multiply this by the number of dendrites and add soma and one extra
   FDCmaxFrags = 2 * par.Value.DimSpace * DendriteLength + 2;
   cout << "FDCmaxFrags=" << FDCmaxFrags << "\n";
   // FDCmaxFrags can be used to define the dimension of FDC associated arrays.

   use_antigen = 0;

   // get the number of antigens to be used for start and the max dimension of the arrays
   int n_Antigen = 0;
   if (par.Value.use_sequence_space) {
      n_Antigen = par.Value.init_antigen_sequences;
   } else {
      n_Antigen = par.Value.APeakNumber;
   }
   n_Antigen_max = n_Antigen_dim_factor * n_Antigen;
   ana << "Generate FDC arrays for max " << n_Antigen_max << " antigens.\n"
       << "GC initialised with " << n_Antigen << " antigens at start.\n";

   // define ag_fraction from parameter data set:
   int a = 0;
   double sum = 0;
   ana << "Relative fraction of the different antigens from parameter file:\n";
   while (a < n_Antigen && par.Value.ag_fraction[a] >= 0) {
      ag_fraction[a] = par.Value.ag_fraction[a];
      sum += ag_fraction[a];
      ana << "ag_fraction[" << a << "]=" << ag_fraction[a] << "\n";
      ++a;
   }
   /* distribute the remaining ones:
    * 1-sum is to be distributed on n_Antigen-a Ag-types
    */
   ana
   << "Relative fraction of the different antigens equally distributed on the remaining Ags:\n";
   for (int i = a; i < n_Antigen; i++) {
      ag_fraction[i] = (1. - sum) / double (n_Antigen - a);
      ana << "ag_fraction[" << i << "]=" << ag_fraction[i] << "\n";
   }

   // save some modes of action
   ag_distribution_mode = par.Value.ag_distribution_mode;
   if (n_Antigen > 1) {
      if (ag_distribution_mode == 0) {
         ana << "Multiple antigens per FDC site mode:\n"
             << "  Attribute the " << n_Antigen
             << " antigens to each FDC site according to ag_fraction[]\n";
      } else {
         ana << "Single antigen per FDC site mode:\n"
             <<
         "  Attribute a single randomly chosen antigen with probability from ag_fraction[]\n";
      }
   }
   ag_detection_mode = par.Value.ag_detection_mode;
   if (n_Antigen > 1) {
      if (ag_detection_mode == 0) {
         ana
         << "Chose Ag for interaction with BCR at a FDC site according to highest affinity\n";
      } else {
         ana
         <<
         "Chose Ag for interaction with BCR at a FDC site according to highest availability\n";
      }
   }
   if (par.Value.ag_per_FDC >= 0) {
      use_antigen = 1;
      ag_threshold = par.Value.ag_threshold;
      ana << "Threshold Ag-concentration on FDC for CC activation = " << ag_threshold
          << " Mol.\n";
      antigen_amount_per_FDC = par.Value.ag_per_FDC;
      ana << "FDC-antigen is used and each FDC presents "
          << antigen_amount_per_FDC * ag_threshold << " Mol Ag.\n";
      antigen_saturation = par.Value.ag_saturation_FDC;
      ana << "Binding probability saturates at "
          << antigen_saturation * ag_threshold << " Mol Ag per fragment.\n";
      BCreduceFDCag = par.Value.BCreduceFDCag;
      if (BCreduceFDCag) {
	ana << "BC reduce Ag on FDC upon uptake.\n";
      } else {
	ana << "Ag is kept constant on FDCs.\n";
      }
   }

   if (par.Value.ic_k_on >= 0.) {
      // convert the value given per second to values per hour
      ic_k_on = par.Value.ic_k_on * 3600. * ag_threshold;
      ic_k_off = par.Value.ic_k_off * 3600.;
      ana << "Binding rates for Ab-Ag-complex: k_on=" << ic_k_on / ag_threshold
          << "/(M hr); and k_off=" << ic_k_off << "/hr.\n";
   } else {
      ic_k_on = 0.;
      ic_k_off = 0.;
   }
}
void cellFDC::set_statics(const double &time, const Parameter &par) {
   if (time >= par.Value.Start_Differentiation - 1.0e-09) {
      if (par.Value.mksignal < 1.e-08) { p_mksignal = 0.; } else {
         p_mksignal = par.Value.mksignal * par.Value.deltat;
         cout << "Signal production, ";
      }
   } else { p_mksignal = 0.; }
   // CXCL13 and SEMA4D are produced all the time! Thus only the first set_statics is used.
}
// ============================================================

void cellFDC::signal_production(const long &i, sigs &l) {
   // Produce signal for differentiation of CB to CC and for chemotaxis
   if (l.signal_use[sig_differ2CC] == 1) { signal_secretion(i,sig_differ2CC,vesicle,p_mksignal,l); }
   if (l.signal_use[CXCL13] == 1) { signal_secretion(i,CXCL13,vesicle,p_mkCXCL13,l); }
   if (l.signal_use[SEMA4D] == 1) { signal_secretion(i,SEMA4D,vesicle,p_mkSEMA4D,l); }
   // other signals may be included here
}
void cellFDC::add_antigen(int n_Antigen) {
   /* Attributes Ag to all <volume> fragments of the FDC
    * using ag_distribution_mode, n_Antigen and ag_fraction[]
    * Should be called at the beginning of a GCR or when all antigens
    * are thought to be exchanged.
    * ###
    * An extra routine is needed to add antigen to a pre-existing antigen distributions on FDCs.
    */
   double ag_per_frag = antigen_amount_per_FDC / double (volume);
   // ag_per_frag contains the total amount of ag to be distributed onto each fragment of the FDC
   if (ag_distribution_mode == 0) {
      /* take antigen_amount corresponding to ag_fraction on each Ag-presenting site:
       * In the parameter file ag_fraction may contain fractions of Ag that do not sum up to 1.
       * However, in cellFDC::set_statics(..) these were completed, such that all ag_fractions
       * are now having positive values.
       */
      for (int f = 0; f < volume; f++) {
         for (int a = 0; a < n_Antigen; a++) {
            antigen_amount[f][a] += ag_fraction[a] * ag_per_frag;
         }
      }
   } else {
      // if(ag_distribution_mode==1)
      // case of each fragment gets exactly one Ag
      // Go through all FDC fragments:
      for (int f = 0; f < volume; f++) {
         // make a random choice with weight from ag_fraction
         double rtmp = drandom();
         // Go through ag_fraction to determine the chosen Ag
         int a = 0;
         double sum = 0;
         while (rtmp > sum && a < n_Antigen) {
            sum += ag_fraction[a];
            if (rtmp > sum) { ++a; }
         }
         // a should contain the index of the Ag to be used
         antigen_amount[f][a] += ag_per_frag;
      }
   }
}
double cellFDC::get_total_antigen_at_site(int frag) {
  /* here n_Antigen_max is used instead of the real n_Antigen
   * which doesn't harm because it is only called for readout
   * and all variables are zero outside the used array part.
   */
  double t = 0.;
  for (int a = 0; a < n_Antigen_max; a++) {
    t += antigen_amount[frag][a];
  }
  return t;
}
double cellFDC::get_total_antigen() {
   double t = 0.;
   for (int f = 0; f < volume; f++) {
     t += get_total_antigen_at_site(f);
   }
   return t;
}
double cellFDC::get_total_antigen(int ag_index) {
   double t = 0.;
   for (int f = 0; f < volume; f++) {
      t += antigen_amount[f][ag_index];
   }
   return t;
}
void cellFDC::set_antigen_amount(double time, double dt) {
   cerr << "cellFDC::set_antigen_amount(double, double) was called and is not set.\n Abort.\n";
   exit(1);
   /*
    * for (int a=0; a<FDCmaxFrags; a++)
    * if (time>=144. && time<144.+dt) antigen_amount[a]=1.e-03/ag_threshold;
    */
}
int cellFDC::get_voidsites() {
  int nvoid = 0;
  for (int f = 0; f < volume; f++) {
    if (get_total_antigen_at_site(f) < 1.) { ++nvoid; }
  }
  return nvoid;
}
void cellFDC::add2histo_ag_per_site(int* free_ag_site, int* tot_ag_site, int resolution) {
  /* Adds the fraction of remaining antigen at each FDC-site to the histogram <ag_site>.
   * ag_site has resolution + 1 entries.
   */
  double initial_ag_per_site = antigen_amount_per_FDC / double (volume);
  /* Note that in principle it is possible that the ag_per_site differs for each
   * FDC, because the volume might be different and only the total amount is
   * determined in the parameter file.*/
  double free_ag, free_ag_frac, ics, tot_ag_frac;
  for (int f = 0; f < volume; f++) {
    // get the sum of all ag-types at each site
    free_ag = get_total_antigen_at_site(f);
    /* Note that when it comes to distinguishing ag-types, the normalization
     * to the initial amount of antigen of each type per site must be organized. 
     * ag_per_site is the sum of all ag-types that were used. */
    // total amount of ICs for all Ags together (antigens are not distinguished):
    ics = get_immune_complex(f) / cellFDC::ag_threshold;
    //cerr << "f=" << f << ": free_ag=" << free_ag << ", ICs=" << ics;
    // get the percentage of the original total ag per site
    free_ag_frac = free_ag / initial_ag_per_site;
    tot_ag_frac = (free_ag + ics) / initial_ag_per_site;
    //cerr << ", freefrac=" << free_ag_frac << ", totfrac=" << tot_ag_frac;
    // attribute this to the right bin
    int thisbin = int(free_ag_frac * resolution + 0.5);
    //cerr << ", freebin=" << thisbin;
    // the bins are organised centered around free_ag_frac*resolution
    ++free_ag_site[thisbin];
    thisbin = int(tot_ag_frac * resolution + 0.5);
    //cerr << ", totbin=" << thisbin << "\n";
    ++tot_ag_site[thisbin];
  }
}
double cellFDC::get_immune_complex(int fragment) {
  /* Returns the amount of immune complexes irrespective of the antigen,
   * i.e. summing all antigens at this site.
   */
  int length = AntibodyDyn::antibodies_resolution;
  double t = 0.;
  for (int b = 0; b <= length; b++) {
    for (int a = 0; a < n_Antigen_max; a++) {
      // see comment in get_total_antigen()
      t += ic_amount[fragment][b][a];
    }
  }
  return t;
}
double cellFDC::get_total_immune_complex() {
  /* This function is used exclusively for read-out (v2016-02-06).
   * This version is used with or without affinity bins (since v2019-02-08).
   */
   double t = 0.;
   for (int f = 0; f < volume; f++) {
     t += get_immune_complex(f);
   }
   return t;
}
double cellFDC::get_total_immune_complex(const int &ag_index) {
  /* this function is used exclusively for read-out (v2016-02-06, modified in v20190208)
   * this version is used with affinity bins 
   * this version is used to sum ICs for a particular Ag <ag_index>
   */
  int length = AntibodyDyn::antibodies_resolution;
  double t = 0.;
  for (int f = 0; f < volume; f++) {
    for (int b = 0; b <= length; b++) {
      t += ic_amount[f][b][ag_index];
    }
  }
  return t;
}
int cellFDC::get_fragment_index(const long &fragpos) {
   int fragfound = -1;
   int tmp = 0;
   while (fragfound == -1 && tmp < volume) {
      if (fragments[tmp] == fragpos) { fragfound = tmp; } else { ++tmp; }
   }
   if (fragfound == -1) {
      cerr << "Error in cellFDC::get_fragment_index(const long&). Abort.";
      exit(1);
   }
   // fragfound enthaelt den gesuchten Index auf der fragment-list!
   return fragfound;
}
int cellFDC::get_highest_amount_ag(const int &frag_index, AffinitySpace &AS) {
   // Returns the ag-index on fragment frag_index with the highest amount of ag.
   int a = -1;
   double amount = 0;
   int n_Antigen = AS.get_n_Antigen();
   if (n_Antigen > n_Antigen_max) {
      cerr << "ERROR in cellFDC::get_highest_amount_ag():\n"
           << "Number of antigens = " << n_Antigen
           << " is larger than the max number of antigens = " << n_Antigen_max << ".\n";
      exit(1);
   }
   for (int i = 0; i < n_Antigen; i++) {
      if ((antigen_amount[frag_index][i] > amount) && (antigen_amount[frag_index][i] >= 1.)) {
         amount = antigen_amount[frag_index][i];
         a = i;
      }
   }
   // returns -1 if none of the Ags is actually there.
   return a;
}
int cellFDC::get_highest_affinity_ag(const int &frag_index, const long &BCRposAS,
                                     AffinitySpace &AS) {
   /* Returns the ag-index of the Ag on fragment frag_index with the highest affinity to BCRposAS
    * As only the index for the Ag with highest affinity is required,
    * affinity_norm(.,.) can be used,
    * which ignores an eventual multiplication with an amplitude (as done in affinity(.,.)).
    */
   int a = -1;
   double best_aff = 0, aff = 0;
   int n_Antigen = AS.get_n_Antigen();
   if (n_Antigen > n_Antigen_max) {
      cerr << "ERROR in cellFDC::get_highest_affinity_ag():\n"
           << "Number of antigens = " << n_Antigen
           << " is larger than the max number of antigens = " << n_Antigen_max << ".\n";
      exit(1);
   }
   for (int i = 0; i < n_Antigen; i++) {
      // check for antigen being there first:
      if (antigen_amount[frag_index][i] >= 1.) {
         aff = AS.affinity(BCRposAS,AS.get_Antigen(i));
         if (aff > best_aff) {
            best_aff = aff;
            a = i;
         }
      }
   }
   // returns -1 if none of the Ags is actually there.
   return a;
}
short cellFDC::consume_ag(const long &frag_index, const int &ag_index) {
  /* Reduces the local (frag_index) amount of Ag (ag_index).
   * Decision is taken in dependence of antigen_saturation levels.
   * Success is returned as 1 or 0.
   */
  /*
   * cerr<<" ag_amount["<<ag_index<<"] at FDC_frag["<<frag_index<<"]="
   *    <<antigen_amount[frag_index][ag_index]
   *    <<"; ag-sat="<<antigen_saturation
   *    <<"; ratio="<<antigen_amount[frag_index][ag_index]/antigen_saturation<<"\n";
   */
  if ((antigen_amount[frag_index][ag_index] >= 1.)
      && (drandom() < antigen_amount[frag_index][ag_index] / antigen_saturation)) {
    /* Note that here the antigen_saturation limits is applied to each Ag-type.
     * As the total number of Ag-portions per type is less when the total Ag is distributed on
     * several Ag-types,
     * saturation is reached earlier than in the 1 Ag case. Just to know this. */
    // If binding response is positive, just use the antigen right here:
    if (BCreduceFDCag) {
      --antigen_amount[frag_index][ag_index];
    }
    // cerr<<"Bound it!\n";
    return 1;
  } else { 
    //cerr << "FDC-binding failed!\n";
    return 0; 
  }
}
int cellFDC::local_interaction_with_ag(const int &frag_index,
                                       const long &BCRpos_ss,
                                       AffinitySpace &shape) {
   /* Determines and returns the ag_index on the FDC fragment with position FDCposition on the
    * spatial lattice to the calling routine.
    * As in some settings the affinity between Ag and BCR enters the choice of the interacting Ag,
    * the AffinitySpace position of the BCR asking to mate is provided (BCRpos_ss) as well as the
    * AffinitySpace itself for affinity calculations.
    * ATTENTION: if (use_antigen==1) ag consumption has to be organised in the calling routine.
    */
   int ag_index = -1;
   // Distinguish different options:
   if (ag_distribution_mode == 1) {
      // only one ag was taken randomly with prob ~ ag_fraction[]
      /* As now a different Ag is positioned on every fragment, the position frag_index
       * has to be checked for the particular Ag in question. */
      // Determine Ag at this site -- as only one Ag is present here, this should return the right
      // one:
      ag_index = get_highest_amount_ag(frag_index,shape);
   } else {
      // if(ag_distribution_mode==0) -> ag distributed strictly as ag_fraction[]
      /* Multiple Ags are present at this site. */
      // Distinguish detection mode here:
      if (ag_detection_mode == 0) {
         // check for the Ag with highest affinity (irrespective of relative concentrations)
         ag_index = get_highest_affinity_ag(frag_index,BCRpos_ss,shape);
      } else {
         // if(ag_detection_mode==1) -> chose the Ag for interaction highest antigen_amount
         ag_index = get_highest_amount_ag(frag_index,shape);
      }
   }
   // returns -1 if none of the Ags is actually there.
   return ag_index;
}
short cellFDC::antigen_presence(const long &fragpos, int &ag_index) {
   /* multiAg: deprecated version for antigen detection and consumption
    * not called anymore from version v20160207
    */
   /* This function is called from cellCC::bind_antigen(...) and used for checking the presence of
    * the Ag.
    * fragpos is the spatial position of the FDC-fragment on the lattice.
    * ag_index is the index in the list of Ags in AffinitySpace.
    */
   // Distinguish the cases with or without consumption of Ag:
   if (use_antigen == 0) {
      if (ag_distribution_mode == 1) {
         cerr << "In cellFDC::antigen_presence(const long&):\n"
              << "called with use_antigen==0 and ag_distribution_mode==1.\n"
              << "Not yet programmed. Abort.\n";
         exit(1);
      }
      /* multiAg: As now a different Ag is positioned on every fragment, the position fragpos
       * has to be checked for the particular Ag in question.
       * Done. See local_interaction_with_ag(..)
       */
      // cout<<"no antigen-consumption: return 1 at FDC-fragment-position "<<fragpos<<".\n";
      return 1;
   } else {
      // cout<<"Ag-consumption is on:";
      // Es ist zunaechst das Fragment zu bestimmen, dass auf fragpos zeigt:
      int fragfound = get_fragment_index(fragpos);

      // if (antigen_amount[fragfound]>0.) return 1;
      // else return 0;
      /* Hier koennte man sich was komplizierteres ausdenken!!!
       * Problem, wie man eine absolute Skala einfuehrt. Man koennte
       * etwa antigen_amount=1. als 1 Quant einer Antigenmenge interpretieren,
       * dass bei der Bindung desselben verbraucht wird. Das ist eine Skala!
       * Man koennte dann sagen, dass mit 10 vorhandenen Quanten keine
       * Einschraenkung der Bindung vorliegt (Saturation: Habe dazu einen
       * Parameter antigen_saturation eingefuehrt!). Fuer kleinere
       * Mengen wird die Wahrscheinlichkeit der Bindung linear reduziert.
       * Das waere:
       */
      short consumption = consume_ag(fragfound,ag_index);
      return consumption;
   }
}
void cellFDC::mk_immune_complex(const double &d_t, AntibodyDyn &ABS, AffinitySpace &AS) {
   /* This version of chemical kinetics of Ab-Ag interactions is called when Ab affinity bins are
    * active.
    */
   int n_Antigen = AS.get_n_Antigen();
   if (n_Antigen > n_Antigen_max) {
      cerr << "ERROR in cellFDC::mk_immune_complex():\n"
           << "Number of antigens = " << n_Antigen
           << " is larger than the max number of antigens = " << n_Antigen_max << ".\n";
      exit(1);
   }
   // Chemical kinetics are done for each FDC fragment separately:
   for (int f = 0; f < volume; f++) {
      /* At each FDC fragment, there might be a multitude of different Ags.
       * Ab is treated in the thermic bath approximation. 
       * Therefore, no feedback onto Abs is considered.
       * This implies that the amount of Ab present to bind to any Ag 
       * is not changing in dependence on binding to other Ags. 
       * This implies that every Ag-type on this fragment can be calculated
       * independently with the same Ab affinity bin distribution.
       */
      for (int a = 0; a < n_Antigen; a++) {
         double d_ic[ABS.antibodies_resolution + 1];
         for (int i = 0; i <= ABS.antibodies_resolution; i++) {
            // Behebe Rundungsfehler um 1.e-08 (Diffusion von Ab kann lokal sehr klein sein):

            double ab = 1.e+15 * ABS.ab_bins[a].antibodies[i];     // antibodies[i] in
                                                                   // mol/micrometer^3, need
                                                                   // M=mol/l!
            // antigen_amount[f][a] kommt in Einheiten von ag_threshold=1.e-08 M;
            // ic_k_on definiert als =par.Value.ic_k_on*3600.*ag_threshold; // M/(M hr) = 1/hr
            // d.h. ic_k_on enthaelt einen Faktor ag_threshold
            // und das Produkt ic_k_on*antigen_amount[f][a] ist dann in M
            // ic_k_off definiert als par.Value.ic_k_off*3600., d.h. in 1/hr
            // ic_amount wird erst hier gefaellt, so dass das per Definition in M kommt
            // if (ab<1.e-08) ab=0.;
            // cout<<" ic_amount["<<f<<"]["<<i<<"]["<<a<<"]="<<ic_amount[f][i][a]<<"\n";
            /* ### This is using Euler scheme, which is not appropriate when aiming for exact
             * results.
             * Either go for a better scheme or use the steady state approximation.
             * The latter might also be faster.
             */
            d_ic[i] = d_t
                      * (ABS.ab_bins[a].k_on[i] * antigen_amount[f][a] * ab
                         - ABS.ab_bins[a].k_off[i]
                         * ic_amount[f][i][a]);
            /*      cout<<"fragment="<<f<<"; bin="<<i<<": dt="<<d_t<<", kon[bin]="<<kon[i]
             *  <<", koff[bin]="<<koff[i]
             *  <<", antigen_amount[frag]["<<a<<"]="<<antigen_amount[f][a]
             *  <<", ic_amount[frag][bin][ag]="<<ic_amount[f][i][a]
             *  <<", ab_bins[ag].antibodies[bin]="<<ab
             *  <<"--> d_ic[bin]="<<d_ic[i]<<"\n";
             */

            /* Fuer Untersuchungen:
             * if (l.sigsknot[fragments[f]].signal[antibody]>0. || ic_amount[f][a]!=0.) {
             * cout<<"frag "<<f<<": kon="<<ic_k_on<<"; koff="<<ic_k_off<<"; dt="<<d_t
             * <<"; ab="<<l.sigsknot[fragments[f]].signal[antibody]
             * <<"; ag="<<antigen_amount[f][a]
             * <<":
             * konterm="<<d_t*ABS.ab_bins[a].k_on[i]*
	     *            antigen_amount[f][a]*l.sigsknot[fragments[f]].signal[antibody]
             * <<"; koffterm="<<d_t*ABS.ab_bins[a].k_off[i]*ic_amount[f][i][a]
             * <<"; ic="<<ic_amount[f][i][a]<<"; d_ic="<<d_ic<<"\n";
             * }
             */
            // Eliminiere inkonsistente Berechnungen (treten auf wenn ic noch Null ist,
            // geringe Rundungen)
            /*
             * if (d_ic[i]>antibodies[i]) {
             * d_ic[i]=antibodies[i];
             * ++ab_sign_errors;
             * }
             */
            // cout.precision(20);
            if (d_ic[i] < -1. * ic_amount[f][i][a]) {
               // cout<<"d_ic["<<i<<"]="<<d_ic[i]<<";
               // ic_amount["<<f<<"]["<<i<<"]="<<ic_amount[f][i]<<"\n";
               d_ic[i] = -1. * ic_amount[f][i][a];
               ++ic_sign_errors;
               // cout<<"I was here";
            }
            ++ic_calculations;
         }    // end of for(i=0...) running over affinity bins

         // Prevent negative antigen values:
         double sum_dic = 0.;
         double new_dic = 0.;
         for (int i = 0; i <= ABS.antibodies_resolution; i++) {
            sum_dic += d_ic[i];
         }
         /*cout<<"frag="<<f<<": sum_dic="<<sum_dic<<", ag_threshold="<<cellFDC::ag_threshold
          * <<", antigen_amount["<<f<<"]["<<a<<"]="<<antigen_amount[f][a];
          */
         if (sum_dic / cellFDC::ag_threshold > antigen_amount[f][a]) {
            new_dic = antigen_amount[f][a] * cellFDC::ag_threshold;
            ++ag_sign_errors;
            for (int i = 0; i <= ABS.antibodies_resolution; i++) {
               d_ic[i] *= new_dic / sum_dic;
            }
            // cout<<", correct this!";
         }
         // cout<<"\n";

         // Aktualisiere ag, ab und ic:
         for (int i = 0; i <= ABS.antibodies_resolution; i++) {
            /* A reduction of antibodies by building of IC is not considered,
             * because Abs are treated in the thermic bath approximation.
             */
            // antibodies[i]-=d_ic[i];
            antigen_amount[f][a] -= d_ic[i] / cellFDC::ag_threshold; // in units of ag_threshold
            ic_amount[f][i][a] += d_ic[i];     // in M
            /*cout<<"Result f="<<f<<", i="<<i<<": antigen_amount[f][a]="<<antigen_amount[f][a]
             * <<", ic_amount[f][i][a]="<<ic_amount[f][i][a]<<"\n";
             */

            // Sicherheitsabfrage:
            if ((ic_amount[f][i][a] < 0.)
                || (ABS.ab_bins[a].antibodies[i] < 0.)
                || (antigen_amount[f][a] < 0.)) {
               if ((ic_amount[f][i][a] < -1.e-15)
                   //	    || antibodies[i]<-1.e-15
                   || (antigen_amount[f][a] < -1.e-08)) {
                  cout << "Negative values in mk_immune_complex!\n";
                  cout << "frag " << f << ": kon=" << ABS.ab_bins[a].k_on[i] << "; koff="
                       << ABS.ab_bins[a].k_off[i] << "; dt=" << d_t
                       << "; ab=" << ABS.ab_bins[a].antibodies[i]
                       << "; ag=" << antigen_amount[f][a]
                       << ": konterm=" << d_t * ABS.ab_bins[a].k_on[i]
                  * antigen_amount[f][a] * ABS.ab_bins[a].antibodies[i]
                       << "; koffterm=" << d_t * ABS.ab_bins[a].k_off[i] * ic_amount[f][i][a]
                       << "; ic=" << ic_amount[f][i][a] << "; d_ic=" << d_ic[i] << "\n";
                  // exit(1);
               }
               if (ic_amount[f][i][a] < 0.) { ic_amount[f][i][a] = 0.; }
               if (ABS.ab_bins[a].antibodies[i] < 0.) { ABS.ab_bins[a].antibodies[i] = 0.; }
               if (antigen_amount[f][a] < 0.) { antigen_amount[f][a] = 0.; }
            }
            /*cout<<"Result after check f="<<f
             * <<", i="<<i<<": antigen_amount[f][a]="<<antigen_amount[f][a]
             * <<", ic_amount[f][i][a]="<<ic_amount[f][i][a]<<"\n";
             */
         }
      }   // end of for(a=0...) running over Ags
   }  // end of for(f=0...) running over fragments
      /*cout<<"total ag="<<get_total_antigen()
       * <<", total_ic="<<get_total_immune_complex(length)
       * <<"\n";
       */
}
void cellFDC::mk_immune_complex(const double &d_t, sigs &l) {
   /* This version of chemical kinetics is called when no Ab affinity bins are being used.
    * Ab is taken from the signal lattice and only one antibody is considered.
    * All IC derived from this Ab are recollected in the first bin [b=0] of ic_amount.
    * ic_k_on is the same as ABS.ab_bins[Ag=a].k_on[bin=0] and could be used for different Ags.
    * However, ic_k_off is taken from the parameter file when Ab-affinity-bins are switched off.
    * In AntibodyDyn, k_offs are calculated in regimes. This cannot be applied here.
    * Therefore, this routine is only functional with a single Ag!
    */
   for (int f = 0; f < volume; f++) {
      // Behebe Rundungsfehler um 1.e-08 (Diffusion von Ab kann lokal sehr klein sein):
      double ab = l.sigsknot[fragments[f]].signal[antibody];
      if (ab < 1.e-08) { ab = 0.; }
      double d_ic = d_t * (ic_k_on * antigen_amount[f][0] * ab - ic_k_off * ic_amount[f][0][0]);
      /* Eliminiere inkonsistente Berechnungen 
       * (treten auf wenn ic noch Null ist, geringe Rundungen)
       */
      if (d_ic > l.sigsknot[fragments[f]].signal[antibody]) {
         d_ic = l.sigsknot[fragments[f]].signal[antibody];
         ++ab_sign_errors;
      }
      // Folgende Sicherheitsabfrage selten relevant:
      if (d_ic > antigen_amount[f][0]) {
         d_ic = antigen_amount[f][0];
         ++ag_sign_errors;
      }
      // Folgende Abfrage bisher noch nie relevant:
      if (d_ic < -1. * ic_amount[f][0][0]) {
         d_ic = -1. * ic_amount[f][0][0];
         ++ic_sign_errors;
      }
      ++ic_calculations;

      // Aktualisiere ag, ab und ic:
      /* In contrast to mk_immune_complex with Ab-affinity-bins, here there is a feedback onto Ab.
       * This is because here, Ab is not considered to be generated outside by many GCs but is
       * treated as produced even inside the GC. Thus, it could even have an inhomogeneous
       * distribution and needs to diffuse over the GC volume.
       */
      l.sigsknot[fragments[f]].signal[antibody] -= d_ic;
      antigen_amount[f][0] -= d_ic;
      ic_amount[f][0][0] += d_ic;

      // Sicherheitsabfrage:
      // if (antigen_amount[f][0]<0) { cout<<"\n ag:"<<antigen_amount[f][0]<<"\n"; exit(1); }
      if ((ic_amount[f][0][0] < 0.)
          || (l.sigsknot[fragments[f]].signal[antibody] < 0.)
          || (antigen_amount[f][0] < 0.)) {
         cout << "Negative values in mk_immune_complex!\n";
         cout << "frag " << f << ": kon=" << ic_k_on << "; koff=" << ic_k_off << "; dt=" << d_t
              << "; ab=" << l.sigsknot[fragments[f]].signal[antibody]
              << "; ag=" << antigen_amount[f][0]
              << ": konterm=" << d_t * ic_k_on * antigen_amount[f][0]
         * l.sigsknot[fragments[f]].signal[antibody]
              << "; koffterm=" << d_t * ic_k_off * ic_amount[f][0][0]
              << "; ic=" << ic_amount[f][0][0] << "; d_ic=" << d_ic << "\n";
         // exit(1);
         if (ic_amount[f][0][0] < 0.) { ic_amount[f][0][0] = 0.; }
         if (l.sigsknot[fragments[f]].signal[antibody] < 0.) {
            l.sigsknot[fragments[f]].signal[antibody] = 0.;
         }
         if (antigen_amount[f][0] < 0.) { antigen_amount[f][0] = 0.; }
      }
   }
}
// THIS FUNCTION IS NOT YET FUNCTIONAL (written May 2009)
// In February 2016 added some indices to make it compile and thinking of a single Ag only.
void cellFDC::mk_immune_complex(sigs &l) {
   cerr << "In mk_immune_complex(sigs&):\n"
        << "This function is not yet functional and should not be called.\n"
        << "Abort.\n";
   exit(1);
   for (int f = 0; f < volume; f++) {
      // Behebe Rundungsfehler um 1.e-08 (Diffusion von Ab kann lokal sehr klein sein):
      double ic = ic_amount[f][0][0];
      double ab = ic + l.sigsknot[fragments[f]].signal[antibody];
      if (ab < 1.e-08) { ab = 0.; }
      double ag = ic + antigen_amount[f][0];

      // Aktualisiere ag, ab und ic:
      ic_amount[f][0][0] = ic_k_on * ab * ag / ic_k_off;
      l.sigsknot[fragments[f]].signal[antibody] = ab - ic_amount[f][0][0];
      antigen_amount[f][0] = ag - ic_amount[f][0][0];

      // Sicherheitsabfrage:
      // if (antigen_amount[f][0]<0) { cout<<"\n ag:"<<antigen_amount[f][0]<<"\n"; exit(1); }
      if ((ic_amount[f][0][0] < 0.)
          || (l.sigsknot[fragments[f]].signal[antibody] < 0.)
          || (antigen_amount[f][0] < 0.)) {
         cout << "Negative values in mk_immune_complex!\n";
         cout << "frag " << f << ": kon=" << ic_k_on << "; koff=" << ic_k_off
              << "; ab=" << l.sigsknot[fragments[f]].signal[antibody]
              << "; ag=" << antigen_amount[f][0]
              << ": ss-term=" << ic_k_on * ab * ag / ic_k_off
              << "; ic=" << ic_amount[f][0][0] << "\n";
         //      exit(1);
         if (ic_amount[f][0][0] < 0.) { ic_amount[f][0][0] = 0.; }
         if (l.sigsknot[fragments[f]].signal[antibody] < 0.) {
            l.sigsknot[fragments[f]].signal[antibody] = 0.;
         }
         if (antigen_amount[f][0] < 0.) { antigen_amount[f][0] = 0.; }
      }
   }
}
// ============================================================

// ============================================================
// ============================================================
// ============================================================
// ===================== OUTPUT =================================
// ============================================================
// ============================================================

cellOUT::cellOUT() {
    //MS
    state = Outfree;
    out_to_die = 0;
    selfMutation = 0;
    nTFRcontacts = 0;
    tfr_clock = 0;
    tfr_index = -1;
   DEC205 = false;
   set_p_move();
}
cellOUT::~cellOUT() {
   // cout<<"in ~cellOUT() ,,,\n";
}
// ============================================================

double cellOUT::p_difu = 0.;
double cellOUT::p_difu_width = -1.;
double cellOUT::persistence = 0.;
short cellOUT::v_modi = 1;
short cellOUT::n_v_states = 1;
double cellOUT::v_slow_factor = 1.0;
double cellOUT::p_switch_v = 0.0;
double cellOUT::p_mk_ab = 0.;
short cellOUT::vesicle = 1;
double cellOUT::average_affinity = 0.;
double cellOUT::max_affinity = 0.;
short cellOUT::use_threshold = 1;
double cellOUT::initial_ab_affinity = -1;
bool cellOUT::exit2tz = false;
//MSchips
short cellOUT::TFR_CC_interaction_mode;

// ============================================================

cellOUT::cellOUT(const cellOUT &x) : frag_cell(x) {
   // cout<<"in cellOUT(const cellOUT&) ...\n";
   operator =(x);
   cerr << "ZZ.";
}
cellOUT&cellOUT::operator =(const cellCC &x) {
   cell::operator =(x);
   volume = 1;
   DEC205 = x.DEC205;
   fragments[0] = x.index;
   t_immobile[0] = 0;
   set_p_move();
   IgX.set_class(x.IgX);
   //MSchips
   state = Outfree;
   selfMutation = x.selfMutation;
   nTFRcontacts = x.nTFRcontacts;
   if ((TFR_CC_interaction_mode==19 || TFR_CC_interaction_mode==20 || TFR_CC_interaction_mode==18)
           && selfMutation && nTFRcontacts>0) {
       out_to_die=1;
   } else {out_to_die=0;}
   tfr_clock = 0;
   tfr_index = -1;

   // for (int i=0; i<3; i++) barycenter[i]=0.0;
   // !!! Achtung: die aufrufende Funktion muss die Koordinaten noch berechnen.
   // Siehe auch Erklaerungen in cellCB::operator=(cellCC&).
   return *this;
}
cellOUT&cellOUT::operator =(const cellCB &x) {
   cell::operator =(x);
   volume = 1;
   DEC205 = x.DEC205;
   fragments[0] = x.index;
   t_immobile[0] = 0;
   set_p_move();
   IgX.set_class(x.IgX);
   //MSchips
   state = Outfree;
   selfMutation = x.selfMutation;
   nTFRcontacts = x.nTFRcontacts;
   tfr_clock = 0;
   tfr_index = -1;
   return *this;
}
cellOUT&cellOUT::operator =(const cellOUT &x) {
   frag_cell::operator =(x);
   // Set those variables that are not part of frag_cell again!
   // This operator is used when writing to OUTlist and otherwise information is lost.
   DEC205 = x.DEC205;
   set_p_move();
   IgX.set_class(x.IgX);
   //MSchips
   selfMutation = x.selfMutation;
   nTFRcontacts = x.nTFRcontacts;
   state = x.state;
   tfr_clock = x.tfr_clock;
   tfr_index = x.tfr_index;
   out_to_die = x.out_to_die;
   return *this;
}
// ============================================================

void cellOUT::set_statics(const Parameter &par, space &l, ofstream &ana) {
   // bislang gleiche Diffusion fuer CC und OUT:
   if (par.Value.D_CC >= 0.) {
      p_difu = double (l.dim2) * par.Value.D_CC * par.Value.deltat
               / (par.Value.dx * par.Value.dx);
   } else {
      p_difu = 60. * par.Value.v_OUT * par.Value.deltat / par.Value.dx;
   }
   if (p_difu > 0.5) {
      cout << "OUTPUT-Motility = " << p_difu << " to large for dx and dt !!!\n";
      exit(1);
      // Siehe Kommentar bei Centroblasten!
   }
   p_difu_width = par.Value.v_OUT_width * p_difu;
   ana << "  OUT movement probability = " << p_difu << "\n";

   if (par.Value.OUT_persistence > 60. * par.Value.deltat) {
      persistence = 60. * par.Value.deltat / par.Value.OUT_persistence;
   } else {
      persistence = 1.;   // i.e. change polarity in every time step!
   }
   ana << "  OUTPUT polarity persistence = " << par.Value.OUT_persistence
       << " i.e. p_{change-polarity} = " << persistence << "\n";

   double tmpdt = par.Value.OUT_persistence;
   if (tmpdt < 60. * par.Value.deltat) {
      tmpdt = 60. * par.Value.deltat;
      ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
   }
   if (par.Value.D_CC < 0.) {
      ana << "  Expected effective diffusion constant for OUT(=CC) = "
          << 60. * par.Value.v_OUT * par.Value.v_OUT * tmpdt / double (l.dim2)
          << " microns^2/hr = "
          << par.Value.v_OUT * par.Value.v_OUT * tmpdt / double (l.dim2)
          << " microns^2/min\n";
   } else {
      ana << "  Expected mean cell velocity for OUTPUT = "
          << sqrt(double (l.dim2) * par.Value.D_CC / (60. * tmpdt))
          << " microns/min\n";
   }
   v_modi = par.Value.CC_v_modi;
   n_v_states = par.Value.CC_n_v_states;
   if (n_v_states > 1) {
      v_slow_factor = par.Value.v_CC_factor;
      p_switch_v = 60. * par.Value.deltat / par.Value.v_CC_switch_deltat;
      if (p_switch_v < 0.) { p_switch_v = 0.; }
   }

   // chemotaxis / exit direction
   exit2tz = par.Value.exit2tz;

   // production of soluble antibody:
   if (par.Value.mk_ab <= 0.) {
      p_mk_ab = 0.;
   } else { p_mk_ab = par.Value.mk_ab * par.Value.deltat / par.Value.ag_threshold; }

   vesicle = par.Value.FDCvesicle;
   ana << "OUT antibody production mode (=FDC-mode) = " << vesicle << "\n";

   initial_ab_affinity = par.Value.initial_ab_affinity;
   if (initial_ab_affinity < 0) {
      average_affinity = cellCB::average_seeder_affinity;
      max_affinity = cellCB::average_seeder_affinity;
   } else {
      average_affinity = initial_ab_affinity;
      max_affinity = initial_ab_affinity;
   }
   //if (par.Value.use_ab_dynamics > 0) { use_threshold = par.Value.use_ab_dynamics; }
   /* use_threshold is initialised as 1 (i.e. use average output quality)!
    * This implies that when use_ab_dynamics == 0, ab-feedback is still on!? Why?
    * In particular, when Abs are treated dynamically, this here should be switched off!
    * Furthermore, use_ab_dynamics was overwriting use_threshold in ss.cpp.
    * These inconsistencies were removed in v20180213.
    */
   use_threshold = par.Value.use_ab_dynamics;
   ana << "  Ag presenting Abs adapt to PC secreted Ab-affinity = " << use_threshold << "\n";

   //MSchips
   TFR_CC_interaction_mode = par.Value.TFR_mode;
}
void cellOUT::signal_production(const long &i, sigs &l) {
   // Produce antibodies
   signal_secretion(i,antibody,vesicle,p_mk_ab,l);
}
bool cellOUT::canbindTFR(const double &time) {
    bool try2find = true;
    //ASCs can bind TFR as long as they are in -- only one contact allowed // including those coming from CC
    if (nTFRcontacts>0 || (cellCC::TFR_bind_onlyself && !selfMutation)) {
        try2find = false;
    }
    return try2find;
}
void cellOUT::bind_TFR(const double& time, cellTFR &tcell, space &l, AffinitySpace &shape) {
//  if (state!=FDCselected) { // security check
//    cerr<<"in bind_tfr wrong state "<<state<<" num cont "<<nTFRcontacts<<" of cell_ss_pos="<<pos_ss<<"\n";
//    exit(1);
//  }
  state = OutTFRcontact;

  ++nTFRcontacts;

//  ccInd_TFRbound=index;

  tfr_clock = 0.;
  tfr_interaction_time = 0.05;//cellCC::get_tfr_interaction_time();
  tfr_index = tcell.index;
  cerr<<"out is binding "<<tfr_index<<"; is self? "<<selfMutation<<endl;

//  last_tfr_index = tcell.index;
//  last_tfr_id = tcell.id;
  tcell.make_tfr_cc_link(index);
//  tfr_tmp_index=-1;

}
short cellOUT::got_tfr_interaction(const double &time, const double &dt, cellTFR &tcell,
                             space &l, AffinitySpace &shape) {
    short staybound=1;

    tfr_clock+=dt;

    bool try2unbind = false;
        if (tfr_clock > tfr_interaction_time || tcell.to_delete) { try2unbind = true; }
    if (try2unbind) {
        if (selfMutation) {out_to_die=1;}
        state = Outfree;
        staybound = 0;
    }
    if (staybound==0) {
        cerr<<"out with ind "<<index<<" is getting liberated"<<endl;
        tcell.liberateCC(index);
    }
    return staybound;
}
short cellOUT::try_eject(long i, space &l, AffinitySpace &shape) {
    if (l.no_border(i) == 0) {
        if (out_to_die) {return 0;} //selfmutated with contacts with tfr should not produce abs
        else {
            // Fuege die Zelle zum output hinzu
            shape.add_cell(soutext,pos_ss);
            // sout counts all ever generated output cells => no change
            // total shall not be corrected, as it counts all living cells
            return 0;
        }
    }
    return 1;
}
void cellOUT::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
  if (p_difu_width > 0) {
    p_move = get_positive_sample_from_normal(p_difu, p_difu_width);
  }
}
short cellOUT::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   set_polarity_velocity(persistence,v_modi,p_switch_v,n_v_states,v_slow_factor,l,s);
   // if (polarity[2]>0) polarity[2]*=-1.0; ###
   // a repolarisation to south can be done always, because it only acts if polarity has changed:
   if (exit2tz) {
      set_south_polarity(l);
      // cout<<"polarity=("<<polarity[0]<<","<<polarity[1]<<","<<polarity[2]<<")\n";
   }
   if (writethis2track == polarisation) {
      double tmpi[l.dim];
      l.get_koord(index,tmpi);
      td.Write_movement(trackno,out,time,tmpi,polarity,writethis2track);
      // ## note that OUT is of type frag_cell sucht that barycenter (instead of tmpi) is declared
      // ## but it is not defined at any point !!! Reconsider when going to the
      // ## subcellular level.
      writethis2track = trackini;
   }
   return do_diffuse(out,li,l);
}
// ============================================================

// ============================================================
// ============================================================
// ============================================================
// ================== BETA-CELL ===============================
// ============================================================
// ============================================================

double cellbeta::p_pro = 0.;
double cellbeta::max_pro_distance = 0.;
double cellbeta::p_tension = 0.;
double cellbeta::p_grow = 0.;
double cellbeta::p_shrink = 0.;
double cellbeta::p_difu = 0.;
double cellbeta::diffusion_tolerance_min = 0.01;
double cellbeta::diffusion_tolerance_steepness = 0.01;
double cellbeta::elongation = 1.0;
double cellbeta::K_elongation = 2.0;
double cellbeta::smoothmove = 1.0;
double cellbeta::persistence = 0.0;
int cellbeta::target_volume = 1;
short cellbeta::v_modi = 1;
short cellbeta::n_v_states = 1;
double cellbeta::v_slow_factor = 1.0;
double cellbeta::p_switch_v = 0.0;
double cellbeta::max_adhesion = 0.0;

betaWerte cellbeta::p;
double cellbeta::glucose_rest;
double cellbeta::xi,cellbeta::xi_ER,cellbeta::xi_ERC,cellbeta::xi_S_ERC;
double cellbeta::J_K,cellbeta::J_Na,cellbeta::J_Ca,cellbeta::J_ions;
double cellbeta::I[N_currents];
double cellbeta::V_ER_0;
ode cellbeta::solver(cellbeta::N_equations);
double cellbeta::betadt;
long cellbeta::betandt;
long cellbeta::n_write;
suffix cellbeta::beta_index = "0000";

// ==============================================================================
// ==============================================================================

cellbeta::cellbeta() {
   // cerr<<"in cellbeta default constructor ...";
   // state=;
   get_adhesion(max_adhesion);

   method = RungeKutta_4th;
   prhs = &rhs;
   // cerr<<"prhs initialised.\n";

   y_n = new double[N_equations];
   y_n1 = new double[N_equations];
   y_n_old = new double[N_equations];
   rho = new double[betaWerte::N_beta_proteins];
   // cerr<<"y_n(1) memory reserved.\n";

   set_initial_values();
   set_p_move();
   step_count = 0;
   // cerr<<"initial values set.\n";
   // cerr<<"end of cellbeta default constructor.\n";
}
// ============================================================

cellbeta::cellbeta(const cellbeta &x) : frag_cell(x) {
   ///Philippe danger
   // cout<<"in cellbeta::cellbeta(cellbeta&) ...";
   y_n = new double[N_equations];
   y_n1 = new double[N_equations];
   y_n_old = new double[N_equations];
   rho = new double[betaWerte::N_beta_proteins];

   method = RungeKutta_4th;
   prhs = &rhs;  // Pointer on rhs(...)

   // set_initial_values(); // no set... for data are taken from the cell x
   // -> don't miss any variable!
   operator =(x);
   step_count = 0;
   set_p_move();
   // cout<<"end of cellbeta::cellbeta(cellbeta&)\n";
}
// ============================================================

cellbeta&cellbeta::operator =(const cellbeta &x) {
   // cout<<"in cellbeta::=(cellbeta&) ... ";
   frag_cell::operator =(x);
   // state=x.state;
   for (short i = 0; i < N_equations; i++) {
      y_n[i] = x.y_n[i];
      y_n1[i] = x.y_n1[i];
      y_n_old[i] = x.y_n_old[i];
   }
   for (short i = 0; i < betaWerte::N_beta_proteins; i++) {
      rho[i] = x.rho[i];
   }
   set_p_move();
   //  rho_gap=x.rho_gap;
   // rho_K_ATP=x.rho_K_ATP;
   // cout<<"end of cellbeta::=(cellbeta&)\n";
   return *this;
}
// ============================================================

void cellbeta::set_initial_values() {
   // Set initial values: (###3)
   // The following are set for the steady state situation without any stimulation:
   y_n[V] = p.V_0 / 1000.;  // in V while p.V_0 is given in mV
   // y_n[0]=-100;
   y_n[K] = p.K_0;
   y_n[Na] = p.Na_0;
   y_n[Ca] = p.Ca_0;
   y_n[g_K_ATP] = get_sigmoidal(p.s_h_K_ATP,p.glu_0,p.kappa_K_ATP);
   y_n[g_K_V] = get_sigmoidal(p.V_h_K_V,p.V_0,p.kappa_K_V);
   y_n[g_Na_V] = get_sigmoidal(p.V_h_Na_V,p.V_0,p.kappa_Na_V);
   y_n[g_fNa_V] = get_sigmoidal(p.V_h_fNa_V,p.V_0,p.kappa_fNa_V);
   y_n[g_Ca_L] = get_sigmoidal(p.V_h_Ca_L,p.V_0,p.kappa_Ca_L);
   y_n[g_Ca_T] = get_sigmoidal(p.V_h_Ca_T,p.V_0,p.kappa_Ca_T);
   y_n[h_K_V] = get_inactivation(p.W_h_K_V,p.V_0,p.lambda_K_V);
   y_n[h_Na_V] = get_inactivation(p.W_h_Na_V,p.V_0,p.lambda_Na_V);
   y_n[h_Ca_L] = get_inactivation(p.W_h_Ca_L,p.V_0,p.lambda_Ca_L);
   y_n[h_Ca_T] = get_inactivation(p.W_h_Ca_T,p.V_0,p.lambda_Ca_T);
   y_n[g_K_Ca] = get_sigmoidal(p.V_h_K_Ca,p.V_0,p.kappa_K_Ca);
   y_n[C_K_Ca] = get_dynamic_half(p.V_0,45.,30.);
   y_n[Ca_ER] = p.Ca_ER_0;
   y_n[V_ER] = V_ER_0;   // defined in set_parameter(..)
   y_n[IP3] = p.IP3_0;
   y_n[g_IP3] = get_Hill(p.Ca_0,p.C_IP3_act,p.n_IP3_act) * p.g_IP3_max;
   double C_IP3_inh_tmp = p.Cbar_IP3_inh * get_sigmoidal(p.P_IP3,p.IP3_0,p.kappa_IP3);
   y_n[h_IP3] = get_Hill(C_IP3_inh_tmp,p.Ca_0,p.n_IP3_inh);
   y_n[g_sK_Ca] = get_sigmoidal(p.C_sK_Ca,p.Ca_0,p.kappa_sK_Ca);
   y_n[glu] = p.glu_0;
   for (short j = 0; j <= 6; j++) {
      y_n[gap_K + j] = 0.;
      y_n[gap_Na + j] = 0.;
      y_n[gap_Ca + j] = 0.;
   }

   // for (int i=0; i<N_equations; i++) cout<<"y_n["<<i<<"]="<<y_n[i]<<"; ";

   // rho_gap=p.rho[betaWerte::gap];
   // rho_K_ATP=p.rho[betaWerte::K_ATP];
   for (short i = 0; i < betaWerte::N_beta_proteins; i++) {
      rho[i] = p.rho[i];
   }

   // in the sigmoidals all quantities are given in mV or mMol!

   /* Note that deviations from the steady state in the initial configuration
    * have to be given below explicitly! */

   // end of initial values.
}
void cellbeta::randomise_protein_expression() {
   if (p.randomisation_type == betaWerte::equal) {
      for (short i = 0; i < betaWerte::N_beta_proteins; i++) {
         // randomise only the proteins that are expressed anyway
         // if (p.rho[i]>0. && betaWerte::beta_proteins(i)==betaWerte::K_ATP)
         if (p.rho[i] > 0.) {
            rho[i] = p.rho[i]
                     * (drandom(2. * p.randomisation_range) + 1 - p.randomisation_range);
         } else {
            // keep the sequence of random numbers:
            double tmp = drandom(2. * p.randomisation_range);
            // avoid compiler warning "unused variable":
            tmp *= 2.;
         }
      }
   }
}
// ============================================================

void cellbeta::open_files() {
   char filename[30];
   strcpy(filename,"beta_a");
   strcat(filename,beta_index);
   strcat(filename,".out");
   // Initialise output file:
   result.open(filename);
   result << "! t[s]  glu  V[mV]  K[mMol] Na  Ca  "   // col 1-6
          << "g_K_ATP  g_K_V  g_Na_V  g_Ca_L  g_Ca_T  "  // col 7-11
          << "glu/glu0   K/K0  N/N0  Ca/Ca0  "    // col 12-15
          << "h_K_V  h_Na_V  h_Ca_L(1-H)  h_Ca_T  "  // col 16-19
          << "f_Ca  g_K_Ca  C_K_Ca  "             // col 20-22
          << "C_ER  C_ER/C_ER0  f_CaER  V_ER  "   // col 23-26
          << "\n";
   result << "!     IP3  IP3/IP30  g_IP3  h_IP3  "          // col 27-30
          << "1-H(NaK)  H(NCX)  H(PMCA)  H(KCa)  "    // col 31-34
          << "H(SERCA)  "                             // col 35
   //	<<"g_sK_Ca  "                                 // col 36
          << "\n";
   // result<<"! Used parameter values\n";

   strcpy(filename,"beta_b");
   strcat(filename,beta_index);
   strcat(filename,".out");
   result2.open(filename);
   result2 << "! t[s] g_sK_Ca gbar_gap g_fNa_V\n";

   strcpy(filename,"beta_i");
   strcat(filename,beta_index);
   strcat(filename,".out");
   currents.open(filename);
   currents << "! [pA] I_NaK, I_K_ATP, I_K_V, I_K_Ca, I_Na_V, "  // col 2-6
            << "I_NCX, I_PMCA, I_Ca_L, I_Ca_T, "         // col 7-10
            << "I_SERCA, I_IP3, "                        // col 11-12
            << "J_K, J_Na, J_Ca, J_ions, I_sK_Ca, "      // col 13-17
            << "sum(I_gap(K)), sum(I_gap(Na)), sum(I_gap(Ca)), "  // col 18-20
            << "I_fNa_V"                                 // col 21
            << "\n";

   strcpy(filename,"beta_r");
   strcat(filename,beta_index);
   strcat(filename,".out");
   currents_rho.open(filename);
   currents_rho << "! [pA/cell] rho* "
                << "(I_NaK, I_K_ATP, I_K_V, I_K_Ca, I_Na_V, "  // col 2-6
                << "I_NCX, I_PMCA, I_Ca_L, I_Ca_T, "    // col 7-10
                << "I_SERCA, I_IP3, "                   // col 11-12
                << "J_K, J_Na, J_Ca, J_ions, I_sK_Ca, "  // col 13-17
                << "sum(I_gap(K)/2D), sum(I_gap(Na)/2D), sum(I_gap(Ca)/2D), "  // col 18-20
                << "I_fNa_V"                                // col 21
                << ")\n";

   strcpy(filename,"beta_n");
   strcat(filename,beta_index);
   strcat(filename,".out");
   reversal.open(filename);
   reversal << "! t[s] V[mV]  Vbar_K  Vbar_Na  Vbar_Ca  Vbar_Ca_ER\n";

   strcpy(filename,"beta_g");
   strcat(filename,beta_index);
   strcat(filename,".out");
   gapfile.open(filename);
   gapfile << "! t[s] gap_K gap_Na gap_Ca gap_K_0..5 gap_Na_0..5 gap_Ca_0..5\n";

   addchar(beta_index);
}
// ============================================================

cellbeta::~cellbeta() {
   // cout<<"in ~cellbeta()...\n";
   // close cell specific files:
   if (LOCAL_FILES) {
      result.close();
      result2.close();
      currents.close();
      reversal.close();
      gapfile.close();
   }
}
// ============================================================
// ==============================================================================

void cellbeta::set_statics(const Parameter &par, space &l, ofstream &ana) {
   // save parameters:
   p = par.betaValue;
   // public glucose resting value:
   glucose_rest = p.glu_0;
   // make a surface density out of the number of gap-junctions per connection:
   p.rho[betaWerte::gap] *= (6. / (4. * pi * p.R_bc * p.R_bc));
   // cout<<"p.rho[betaWerte::gap]= "<<p.rho[betaWerte::gap]<<" \n";
   if (par.Value.show_mode == islet) {
      cout << "WARNING in cellbeta::set_statics(...): \n"
           <<
      "        There are 2 betacell radii: one in class betaWerte one in class cellbeta!!!\n";
   }

   // get target volume of betacells in number of lattice points
   if (l.dim
       == 2) {
      target_volume = int (3.1415 * pow(par.Value.BETA_radius,2.) / pow(l.dx,2.) + 0.5);
   } else { target_volume = int (4.1888 * pow(par.Value.BETA_radius,3.) / pow(l.dx,3.) + 0.5); }
   if (target_volume == 0) { target_volume = 1; }
   ana << "  betacell target volume = N = " << target_volume << "\n";
   // Maximaler Abstand fuer die betacell-Proliferation in Gittereinheiten
   max_pro_distance = par.Value.BETA_max_pro / par.Value.dx;

   ana << "Calculate action probabilities ... \n";
   // cell division
   p_pro = par.Value.BETA_proliferate * par.Value.deltat;
   ana << "  betacell proliferation = " << p_pro << "\n";
   p_grow = par.Value.BETA_grow * par.Value.deltat;
   ana << "  betacell growth = " << p_grow << "\n";
   p_shrink = par.Value.BETA_shrink * par.Value.deltat;
   ana << "  betacell shrink = " << p_shrink << "\n";
   // get the corrected proliferation rate including mitosis only (no growth)
   // p_proliferate=p_proliferate/(1.-growth_fraction);
   if (target_volume > 1) {
      // ignore this in the case of one cell = one lattice point
      p_pro = p_pro * p_grow / (p_grow - target_volume * p_pro);
   }
   ana << "   Probability of mitosis after completed growth = " << p_pro << "\n";
   // calculate growth-phase fraction
   // ana<<par.Value.BETA_proliferate<<"  "<<par.Value.BETA_grow<<"   "<<l.dx<<"
   //   "<<BETA_vol<<"\n";
   // double growth_fraction=0.5*BETA_vol*par.Value.BETA_proliferate/par.Value.BETA_grow;
   double growth_fraction = target_volume * par.Value.BETA_proliferate / par.Value.BETA_grow;
   if (target_volume == 1) { growth_fraction = 0.; }
   ana << "   Fraction of growth phase in the total cell cycle " << growth_fraction << "\n";
   if ((growth_fraction > 0.5) || (p_pro < 0.)) {
      cout << "beta-proliferation rate too large or beta-growth rate to small!\n";
      exit(1);
   }

   // Adhesion:
   // =========
   max_adhesion = par.Value.BETA_max_adhesion;

   // Motility:
   // =========
   // Toleranzparameter fuer die Abweichung von der Kugelform
   diffusion_tolerance_min = par.Value.BETA_distance_tolerance;
   ana << "  betacell diffusion_tolerance_min  = t_{\rm min} = " << diffusion_tolerance_min
       << "\n";
   /* get the exponential factor in the formula for the tolerance
    * for which the half-tolerance-deformation given in the par-file is respected: */
   diffusion_tolerance_steepness = -par.Value.BETA_half_tolerance_deformation
                                   / (log((1.0 - diffusion_tolerance_min)
                                          / (2.0 - diffusion_tolerance_min
                                             * (1.0 + diffusion_tolerance_min))));
   ana << "  betacell diffusion_tolerance_steepness = K_{1/2} = "
       << diffusion_tolerance_steepness << "\n";
   // elongation=par.Value.BETA_elongation*par.Value.BETA_radius/l.dx;
   // minimal elongation is 1 lattice constant!
   // if (elongation<1.) elongation=1.;
   elongation = par.Value.BETA_elongation;
   ana << "  betacell elongation = epsilon = " << elongation << "\n";
   K_elongation = par.Value.BETA_K_elongation;
   ana << "  betacell elongation for half reshaping force = K_epsilon = " << K_elongation << "\n";
   /*
    * if (par.Value.BETA_radius/l.dx<elongation) {
    * ana<<"!!!Warning!!! Elongation leads to barycenter shifts that are larger than r_BETA!\n";
    * cout<<"!!!Warning!!! Elongation leads to barycenter shifts that are larger than r_BETA!\n";
    * }
    */
   smoothmove = par.Value.BETA_smoothmove;
   ana << "  betacell smoothness of active walk = eta_max = " << smoothmove << "\n";
   if (par.Value.BETA_persistence > 60. * par.Value.deltat) {
      persistence = 60. * par.Value.deltat / par.Value.BETA_persistence;
   } else {
      persistence = 1.;   // i.e. change polarity in every time step!
   }
   ana << "  betacell polarity persistence = " << par.Value.BETA_persistence
       << " i.e. p_{change-polarity} = " << persistence << "\n";

   p_tension = 60. * par.Value.BETA_v_cytosol * par.Value.deltat / par.Value.dx;
   ana << "  betacell surface tension (tau): Fragment movement probability = " << p_tension
       << "\n";
   if (p_tension > 0.5) {
      cout << "betacell-Cytosol-Fragment-Movement-Probability (p="
           << p_tension
           << ") is too large for dx and dt !!!\n";
      exit(1);
      /* Hier wird die Wahrscheinlichkeit berechnet, mit der ein Centroblast
       * in irgendeine Richtung bewegt. Daher keine Normierung mit 1/(2*l.dim)! */
   }

   p_difu = 60. * par.Value.BETA_v * par.Value.deltat / par.Value.dx;
   double dfactor = 0.7;  // approximate value for use_D_correction==0
   if ((use_D_correction == 1)
       && (target_volume > 1)) { dfactor = (0.18 + (2.0 / double (target_volume))); }
   // Calculate smoothfactor for a spherical cell with target volume:
   double smoothfactor;
   if (l.dim == 2) {
      smoothfactor = smoothmove * par.Value.dx / (6.28 * par.Value.BETA_radius);
   } else {
      smoothfactor = smoothmove * par.Value.dx * par.Value.dx
                     / (12.56 * par.Value.BETA_radius * par.Value.BETA_radius);
   }
   if (p_difu * smoothfactor / (dfactor) > 1.0) {
      cout << "betacell-diffusion (p="
           << p_difu * smoothfactor / (dfactor)
           << ") is too large for dx and dt !!!\n";
      exit(1);
      /* Hier wird die Wahrscheinlichkeit berechnet, mit der ein Centroblast
       * in irgendeine Richtung diffundieren kann. Daher keine Normierung
       * mit 1/(2*l.dim)! */
   }
   ana << "  betacell movement probability = "
       << p_difu
       << ", with additional factors = "
       << p_difu * smoothfactor / dfactor
       << "\n";
   double tmpdt = par.Value.BETA_persistence;
   if (tmpdt < 60. * par.Value.deltat) {
      tmpdt = 60. * par.Value.deltat;
      ana << "  ! Expected value for diffusion/velocity scales with time resolution !\n";
   }
   ana << "  Expected effective diffusion constant for betacell = "
       << 60. * par.Value.BETA_v * par.Value.BETA_v * tmpdt / double (l.dim2)
       << " microns^2/hr = "
       << par.Value.BETA_v * par.Value.BETA_v * tmpdt / double (l.dim2)
       << " microns^2/min\n";
   v_modi = par.Value.BETA_v_modi;
   n_v_states = par.Value.BETA_n_v_states;
   if (n_v_states > 1) {
      v_slow_factor = par.Value.BETA_v_factor;
      p_switch_v = 60. * par.Value.deltat / par.Value.BETA_v_switch_deltat;
      if (p_switch_v < 0.) { p_switch_v = 0.; }
   }

   // end of motility.
   // ================

   // ===========================================
   // betacell electrophysiology:
   // get local time steps
   betadt = p.dt;
   if (betadt > 3600. * par.Value.deltat) {
      // beta-time-step cannot be larger than outer time-step!
      betadt = 3600. * par.Value.deltat;
      betandt = 1;
   } else {
      // the local timestep has to exactly fit into the outer timestep:
      betandt = int (3600. * par.Value.deltat / p.dt + 0.5);
      betadt = 3600. * par.Value.deltat / double (betandt);
   }
   ana << "Betacell timestep was " << p.dt << " sec. " << betadt << " sec is used.\n";
   cout << "Betacell timestep was " << p.dt << " sec.\n"
        << "hyphasma timestep was " << 3600. * par.Value.deltat << " sec.\n"
        << "betadt becomes " << betadt << " sec is used.\n"
        << "betandt becomes " << betandt << "\n";

   n_write = long (p.dt_output / betadt + 0.5);
   if (n_write == 0) { n_write = 1; }
   if (LOCAL_FILES) { ana << "Write beta every " << n_write << " betacell timestep\n"; }

   // xi, xi_ER, and xi_ERC are surface per volume (of cell and ER, respectively)
   xi = 3.0 / p.R_bc;
   xi_ERC = p.Sur_ER / ((4.0 * pi * pow(p.R_bc,3.0) / 3.0) - p.Vol_ER);
   // ## This is an approximation. One may think of improving.
   xi_S_ERC = p.Sur_ER / (4.0 * pi * pow(p.R_bc,2.0));
   if ((par.Value.show_mode == islet) && (p.Vol_ER == 0.0)) {
      xi_ER = 0.0;
      cout << "WARNING: xi_ER is set to zero because of zero ER-volume!\n";
      ana << "WARNING: xi_ER is set to zero because of zero ER-volume!\n";
   } else { xi_ER = p.Sur_ER / p.Vol_ER; }
   // cout<<xi<<","<<xi_ER<<","<<xi_ERC<<"\n";

   // Reversal potentials in V
   double Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER;
   double tmintmp = 3600. * par.Value.tmin;
   get_all_Nernst(tmintmp,p.K_0,p.Na_0,p.Ca_0,p.Ca_ER_0,
                  Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);
   // note that p.Value.tmin is given in hours, while get_all_Nernst expects seconds.
   // get_all_Nernst(p.t_0,p.K_0,p.Na_0,p.Ca_0,p.Ca_ER_0,
   //		  Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);

   // double Vbar_K=get_Nernst(p.K_0,p.K_ext,z_K);
   // double Vbar_Na=get_Nernst(p.Na_0,p.Na_ext,z_Na);
   // double Vbar_Ca=get_Nernst(p.Ca_0,p.Ca_ext,z_Ca)-p.Vbar_Ca_delta/1000.;
   // double Vbar_ER=get_Nernst(p.Ca_0,p.Ca_ER_0,z_Ca);

   ana << "Reversal potentials in resting state are:\n"
       << "Vbar_K=" << 1000. * Vbar_K << "mV; "
       << "Vbar_Na=" << 1000. * Vbar_Na << "mV; "
       << "Vbar_ER=" << 1000. * Vbar_ER << "mV; "
       << "Vbar_Ca=" << 1000. * Vbar_Ca << "mV; "
       << "corrected by " << p.Vbar_Ca_delta << "mV.\n";
   double V_0 = p.V_0 / 1000.;  // convert p.V_0 in mV to V_0 in V
   double V_0_mV = p.V_0;

   // calculate all currents for equilibrium values
   //     note that the V_0 are used in V in the brackets and in mV in sigmoidal(...)
   double I_NaK = p.Ihat_NaK
                  * (1.0 - get_Hill(p.K_0,p.H_NaK,p.n_NaK))
                  * get_Hill(p.Na_0,p.H2_NaK,p.n2_NaK);
   double I_K_ATP = (1.0 - get_sigmoidal(p.s_h_K_ATP,p.glu_0,p.kappa_K_ATP))
                    * p.gbar_K_ATP
                    * (V_0 - Vbar_K);
   double I_K_V = get_sigmoidal(p.V_h_K_V,p.V_0,p.kappa_K_V)
                  * get_inactivation(p.W_h_K_V,p.V_0,p.lambda_K_V)
                  * p.gbar_K_V * (V_0 - Vbar_K);
   double H_K_Ca = p.H_K_Ca;
   if (p.use_dynamic_H_K_Ca == 1) { H_K_Ca = get_dynamic_half(V_0_mV,45.,30.); }
   double I_K_Ca = p.gbar_K_Ca
                   * (V_0 - Vbar_K)
                   * get_Hill(p.Ca_0,H_K_Ca,p.n_K_Ca);
   if (p.use_voltage_gating_K_Ca == 1) {
      I_K_Ca *= get_sigmoidal(p.V_h_K_Ca,V_0_mV,p.kappa_K_Ca);
   }
   double I_sK_Ca = get_sigmoidal(p.C_sK_Ca,p.Ca_0,p.kappa_sK_Ca)
                    * p.gbar_sK_Ca * (V_0 - Vbar_K);
   double I_Na_V = get_sigmoidal(p.V_h_Na_V,p.V_0,p.kappa_Na_V)
                   * get_inactivation(p.W_h_Na_V,p.V_0,p.lambda_Na_V)
                   * p.gbar_Na_V
                   * (V_0 - Vbar_Na);
   double I_fNa_V = get_sigmoidal(p.V_h_fNa_V,p.V_0,p.kappa_fNa_V)
                    * p.gbar_fNa_V
                    * (V_0 - Vbar_Na);
   double I_NCX = p.Ihat_NCX
                  * get_Hill(p.Ca_0,p.H_NCX,p.n_NCX);
   double I_PMCA = p.Ihat_PMCA
                   * get_Hill(p.Ca_0,p.H_PMCA,p.n_PMCA);
   double I_Ca_L = get_sigmoidal(p.V_h_Ca_L,p.V_0,p.kappa_Ca_L)
                   * get_inactivation(p.W_h_Ca_L,p.V_0,p.lambda_Ca_L)
                   * p.gbar_Ca_L
                   * (V_0 - Vbar_Ca);
   I_Ca_L *= (1.0 - get_Hill(p.Ca_0,p.C_Ca_L,p.n_Ca_L));
   double I_Ca_T = get_sigmoidal(p.V_h_Ca_T,p.V_0,p.kappa_Ca_T)
                   * get_inactivation(p.W_h_Ca_T,p.V_0,p.lambda_Ca_T)
                   * p.gbar_Ca_T
                   * (V_0 - Vbar_Ca);

   double C_IP3_inh = p.Cbar_IP3_inh * get_sigmoidal(p.P_IP3,p.IP3_0,p.kappa_IP3);
   double I_SERCA = p.Ihat_SERCA * get_Hill(p.Ca_0,p.H_SERCA,p.n_SERCA);
   double I_IP3 = p.gbar_IP3 * p.g_IP3_max
                  * get_Hill(p.Ca_0,p.C_IP3_act,p.n_IP3_act)
                  * get_Hill(C_IP3_inh,p.Ca_0,p.n_IP3_inh);  // voltage missing here (see below)
   // double V_ER_0=p.V_ER_0/1000.; // V_ER_0 is now in V
   // Determine V_ER_0 from steady state condition:
   // This might be replaced by a fixed value if known later!!! ###
   if ((par.Value.show_mode == islet) && ((p.rho[betaWerte::IP3] == 0.0) || (I_IP3 == 0.0))) {
      V_ER_0 = 0.0;
      cout
      << "WARNING: V_ER_0 was set to zero because of vanishing IP3-current in steady state!\n";
      ana
      << "WARNING: V_ER_0 was set to zero because of vanishing IP3-current in steady state!\n";
   } else {
      V_ER_0 = V_0 - get_Nernst(p.Ca_0,p.Ca_ER_0,z_Ca)
               + (p.rho[betaWerte::SERCA] * I_SERCA) / (p.rho[betaWerte::IP3] * I_IP3);
      ana << "Calculated resting ER-potential is at " << 1000.0 * V_ER_0 << " mV.\n";
   }
   // +++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++++++
   ///*
   V_ER_0 = 0;
   if (par.Value.show_mode == islet) {
      cout << "WARNING: Set V_ER_0 by hand to the value " << 1000.0 * V_ER_0 << " mV.\n";
      cout << "         Change this in the code if wanted!\n";
      ana << "WARNING: Set V_ER_0 by hand to the value " << 1000.0 * V_ER_0 << " mV.\n";
      ana << "         Change this in the code if wanted!\n";
      // */
      I_IP3 *= (V_0 - V_ER_0 - Vbar_ER);
      // /*
      p.rho[betaWerte::SERCA] = (-1.0) * p.rho[betaWerte::IP3] * I_IP3 / I_SERCA;
      cout << "WARNING: As V_ER was set by hand the steady state has to be defined\n"
           << "         with the ER protein densities: set rho[betaWerte::SERCA] to "
           << p.rho[betaWerte::SERCA] << "/micron^2\n";
      ana << "WARNING: As V_ER was set by hand the steady state has to be defined\n"
          << "         with the ER protein densities: set rho[betaWerte::SERCA] to "
          << p.rho[betaWerte::SERCA] << "/micron^2\n";
   }
   // */
   // ++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++++++++

   /*
    * double n_Hill=1.0;
    * I_IP3    = p.Ihat_IP3
    * pow(get_Hill(p.Ca_0,p.H_IP3_act,n_Hill),p.n_IP3_act)
    * pow(get_Hill(p.H_IP3_inact,p.Ca_0,n_Hill),p.n_IP3_inact)
    * pow(get_Hill(p.IP3_0,p.H_IP3_IP3,n_Hill),p.n_IP3_IP3);
    */

   // calculate leakage currents J's
   J_K = -1.0 * (-2.0 * p.rho[betaWerte::NaK] * I_NaK
                 + p.rho[betaWerte::K_ATP] * I_K_ATP
                 + p.rho[betaWerte::K_V] * I_K_V
                 + p.rho[betaWerte::K_Ca] * I_K_Ca
                 + p.rho[betaWerte::sK_Ca] * I_sK_Ca);
   J_Na = -1.0 * (2.0 * p.rho[betaWerte::NaK] * I_NaK * p.alpha_NaK
                  + p.rho[betaWerte::NCX] * I_NCX * p.alpha_NCX
                  + p.rho[betaWerte::Na_V] * I_Na_V
                  + p.rho[betaWerte::fNa_V] * I_fNa_V);
   J_Ca = -1.0 * (p.rho[betaWerte::Ca_L] * I_Ca_L
                  + p.rho[betaWerte::Ca_T] * I_Ca_T
                  - z_Ca * p.rho[betaWerte::NCX] * I_NCX
                  + p.rho[betaWerte::PMCA] * I_PMCA
                  + p.rho[betaWerte::SERCA] * I_SERCA * xi_ERC / xi
                  + p.rho[betaWerte::IP3] * I_IP3 * xi_ERC / xi);

   if (p.set_leakage_zero == 1) {
      ana << "SS-leakage would be: J_K=" << J_K << ", J_Na=" << J_Na
          << ", J_Ca=" << J_Ca << ", J_ions=" << J_ions << " [pA/micron^2]\n";
      J_Ca = 0.0;
      J_Na = 0.0;
      J_K = 0.0;
   }
   /*
    * J_ions = -1.0*(p.rho[betaWerte::NaK]*I_NaK+
    * p.rho[betaWerte::K_ATP]*I_K_ATP+
    * p.rho[betaWerte::K_V]*I_K_V+
    * p.rho[betaWerte::K_Ca]*I_K_Ca+
    * p.rho[betaWerte::Na_V]*I_Na_V+
    * p.rho[betaWerte::fNa_V]*I_fNa_V+
    * p.rho[betaWerte::NCX]*I_NCX+
    * p.rho[betaWerte::PMCA]*I_PMCA+
    * p.rho[betaWerte::Ca_L]*I_Ca_L+
    * p.rho[betaWerte::Ca_T]*I_Ca_T+
    * p.rho[betaWerte::SERCA]*I_SERCA*xi_S_ERC+
    * p.rho[betaWerte::IP3]*I_IP3*xi_S_ERC+
    * J_K+J_Na+J_Ca);
    */
   J_ions = 0.0;  // per definition
   ana << "Leakage currents: J_K=" << J_K << ", J_Na=" << J_Na
       << ", J_Ca=" << J_Ca << ", J_ions=" << J_ions << " [pA/micron^2]\n";

   // Calculate the total calmodulin resting concentration needed to have f_Ca=0.01
   double cal_virtual = (p.Ca_0 + p.K_cal) * (99. - (p.buf_0 / (p.Ca_0 + p.K_buf)));
   ana << "Needed calmodulin for f_Ca(C_0)=0.01 is " << cal_virtual << " mMol\n";
}
// ============================================================
// ==============================================================================

void cellbeta::ini(const long &i, const long &li, const double &t, space &l) {
   /* Diese ini-Routine geht von einer bereits vordefinierten Zelle
    * aus, die in dem Objekt gespeichert ist. Es werden an dem Objekt
    * die Aenderungen vorgenommen, die es zu einem neuen Zellobjekt
    * machen. Ausserdem werden lattice und shapespace aktualisiert.
    * Achtung: state wird nicht gesetzt!
    */
   // Change lattice index in newBETA to i
   index = i;
   born_index = i;
   born_time = t;
   // initialize fragments (to empty) and add the cell center only:
   // (so fragment lists handed over in newBETA will be destroyed!!!)
   volume = 1;
   fragments[0] = i;
   l.get_koord(i,barycenter);
   for (short a = 0; a < l.dim; a++) {
      last_position[a] = barycenter[a];
   }
   // Actualize lattice point
   l.set_knot(i,BETA,li);
   // get initial polarity vector
   l.get_random_direction(polarity);
   // get initial radius
   get_radius(l.dim);
   // cout<<"polar-vector=("<<polarity[0]<<","<<polarity[1]<<")\n";
   if (LOCAL_FILES) { open_files(); }
}
void cellbeta::synchronise() {
   for (int i = 0; i < N_equations; i++) {
      y_n_old[i] = y_n[i];
   }
}
// ==============================================================================
// ==============================================================================

double cellbeta::get_2sigmoidal(double &t, double t_a, double t_b,
                                double rest, double factor,
                                double kappa_a, double kappa_b) {
   double f = 1. / ((1 + exp((t_a - t) / kappa_a)) * (1 + exp((t - t_b) / kappa_b)));  // weight
                                                                                       // function
   double g = rest * (1 - f) + factor * rest * f;  // mMol, external potassium dynamics
   return g;
}
double cellbeta::get_glucose(const long &i, sigs &s) {
   y_n[glu] = s.sigsknot[i].signal[glucose];
   // if a higher time resolution is needed, reuse sigs::get_2sigmoidal(...) in signals.x!
   return y_n[glu];
}
double cellbeta::get_IP3(double &t) {
   return p.IP3_0;
}
double cellbeta::get_tau_IP3(double &ip3) {
   double k_on = 4.6e+04;  // in 1/(mMol second)
   return 1.0 / (k_on * ip3);  // in seconds
}
double cellbeta::get_sigmoidal(double &half, double &x, double &kappa) {
   return 1.0 / (1.0 + exp((half - x) / kappa));
}
double cellbeta::get_inactivation(double &half, double &x, double &kappa) {
   if (p.use_inactivation == 0) { return 1.0; } else { return get_sigmoidal(x,half,kappa); }
}
double cellbeta::get_ca_buffer(double &b_0, double &c, double &dissociation) {
   return b_0 * dissociation / pow(c + dissociation,2.0);
}
double cellbeta::get_Ca_free_fraction(double &c) {
   // Calculate total calcium and free fraction:
   double c_total = c * (1.0 + (p.cal_0 / (c + p.K_cal)) + (p.buf_0 / (c + p.K_buf)));
   return c / c_total;
}
double cellbeta::get_Ca_ER_free_fraction(double &c) {
   // Calculate total calcium and free fraction:
   double c_total = c * (1.0 + (p.buf_ER_0 / (c + p.K_buf_ER)));
   return c / c_total;
}
double cellbeta::get_dynamic_half(double &x, double a, double b) {
   return exp((a - x) / b) / 1000.;  // value associated with mMol (no real units here)
}
double cellbeta::get_tau_V(double &x, double a, double b, double c, double Vx) {
   return c / (exp((x - Vx) / a) + exp((Vx - x) / b));
}
double cellbeta::get_tau_K_V_sherman88(double &x, double &VbarK) {
   // following sherman et al 1988 equation 2.3 for the dynamics of K,V channel
   // double c=0.0353; // sec;  c=60ms; lambda=1.7 in Sherman et al 1988
   double c = 2.0 * p.tau_K_V;  // use the parameter as half value
   // 0.060; // ms (fits well with range in Kelly 1991)
   double a = 0.065;  // V
   double b = 0.020;  // V
   return get_tau_V(x,a,b,c,VbarK);
}
double cellbeta::get_K_ext(double &t, double &k_0) {
   // This function is optionally called from rhs():
   // set the option there and fill values here!
   double t_a = 3.;  // s, time of external potassium increase
   double t_b = 53.;  // s, time of external potassium decrease
   double kappa_a = 0.2;  // s, width of external potassium increase
   double kappa_b = 0.5;  // s, width of external potassium decrease
   double f = 1. / ((1 + exp((t_a - t) / kappa_a)) * (1 + exp((t - t_b) / kappa_b)));  // weight
                                                                                       // function
   double k_1 = 13.0;  // mMol, high external potassium concentration
   double k = k_0 * (1 - f) + k_1 * f;  // mMol, external potassium dynamics
   // cout<< "t="<<t<<"s;  K_ext="<<k<<" mMol;  ";
   return k;
}
double cellbeta::get_Nernst(double &x, double &x_ext, double valence) {
   return Rydberg * p.T * log(x_ext / x) / (valence * Faraday);  // [V]
}
void cellbeta::get_all_Nernst(double &t, double &K, double &Na, double &Ca, double &Ca_ER,
                              double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER) {
   double K_ext = p.K_ext;
   // +++++++++++++++++++++++++ OPTION ++++++++++++++++++++
   // Set here alternative values for external potassium:
   // K_ext=get_K_ext(t,p.K_ext);
   // +++++++++++++++++++++ end OPTION ++++++++++++++++++++

   if (p.use_Nernst == 1) {
      Vbar_K = get_Nernst(K,K_ext,z_K);   // in V
      Vbar_Na = get_Nernst(Na,p.Na_ext,z_Na);
      Vbar_Ca = get_Nernst(Ca,p.Ca_ext,z_Ca) - p.Vbar_Ca_delta / 1000.;
      Vbar_ER = get_Nernst(Ca,Ca_ER,z_Ca);
   } else {
      Vbar_K = get_Nernst(p.K_0,K_ext,z_K);   // in V
      Vbar_Na = get_Nernst(p.Na_0,p.Na_ext,z_Na);
      Vbar_Ca = get_Nernst(p.Ca_0,p.Ca_ext,z_Ca) - p.Vbar_Ca_delta / 1000.;
      Vbar_ER = get_Nernst(p.Ca_0,p.Ca_ER_0,z_Ca);
   }
}
void cellbeta::insert(double &ion, double &voltage, double I_load, double valence_sign) {
   short showit = false;
   // I_load is <0 if positive ions enter the cell
   if (showit == true) { cout << "ion_before=" << ion << "  "; }
   ion -= I_load / (abs(int (valence_sign)) * MFaraday * (4.0 * pi * pow(p.R_bc,3.0) / 3.0));
   if (showit == true) { cout << "ion_after=" << ion << "\n"; }
   // ion-=I_load*xi/(abs(int(valence_sign))*MFaraday);

   if (showit == true) { cout << "v_before=" << 1000. * voltage << "  "; }
   double sign = -1.0;
   if (valence_sign > 0) { sign = 1.0; }
   voltage -= sign * I_load / (p.C_m * (4.0 * pi * pow(p.R_bc,2.0)));
   // I_load in pA, then voltage in V
   // voltage-=sign*I_load/p.C_m;
   // I_load in pA/micron^2, then voltage in V
   if (showit == true) {
      cout << "v_after=" << 1000. * voltage << "\n";
      char cc;
      cin >> cc;
   }
}
// ==============================================================================
// ==============================================================================

void cellbeta::get_currents(double * y,
                            double &Vbar_K, double &Vbar_Na, double &Vbar_Ca, double &Vbar_ER,
                            double &C_IP3_inh) {
   I[NaK] = p.Ihat_NaK * (1.0 - get_Hill(y[K],p.H_NaK,p.n_NaK))
            * get_Hill(y[Na],p.H2_NaK,p.n2_NaK);
   I[K_ATP] = (1.0 - y[g_K_ATP]) * p.gbar_K_ATP * (y[V] - Vbar_K);
   I[K_V] = y[g_K_V] * y[h_K_V] * p.gbar_K_V * (y[V] - Vbar_K);
   double H_K_Ca = p.H_K_Ca;
   // cerr<<t<<": H_K_Ca="<<H_K_Ca;
   if (p.use_dynamic_H_K_Ca == 1) { H_K_Ca = y[C_K_Ca]; }
   // cerr<<"; H_K_Ca(dyn)="<<H_K_Ca<<"\n";
   I[K_Ca] = p.gbar_K_Ca * (y[V] - Vbar_K) * get_Hill(y[Ca],H_K_Ca,p.n_K_Ca);
   // cerr<<"; I[K_Ca]="<<I[K_Ca];
   if (p.use_voltage_gating_K_Ca == 1) { I[K_Ca] *= y[g_K_Ca]; }
   // cerr<<"; I[K_Ca](with V)="<<I[K_Ca]<<"\n";
   I[sK_Ca] = y[g_sK_Ca] * p.gbar_sK_Ca * (y[V] - Vbar_K);
   I[Na_V] = y[g_Na_V] * y[h_Na_V] * p.gbar_Na_V * (y[V] - Vbar_Na);
   I[fNa_V] = y[g_fNa_V] * p.gbar_fNa_V * (y[V] - Vbar_Na);
   I[NCX] = p.Ihat_NCX * get_Hill(y[Ca],p.H_NCX,p.n_NCX);
   I[PMCA] = p.Ihat_PMCA * get_Hill(y[Ca],p.H_PMCA,p.n_PMCA);
   I[Ca_L] = y[g_Ca_L] * y[h_Ca_L] * p.gbar_Ca_L * (y[V] - Vbar_Ca);
   I[Ca_L] *= (1.0 - get_Hill(y[Ca],p.C_Ca_L,p.n_Ca_L));
   // cout<<t<<" Ca,L inact="<<get_Hill(Ca,p.C_Ca_L,p.n_Ca_L)<<"\n";
   I[Ca_T] = y[g_Ca_T] * y[h_Ca_T] * p.gbar_Ca_T * (y[V] - Vbar_Ca);
   // cerr<<I[Ca_T]<<"\n";

   C_IP3_inh = p.Cbar_IP3_inh * get_sigmoidal(p.P_IP3,y[IP3],p.kappa_IP3);
   I[SERCA] = p.Ihat_SERCA * get_Hill(y[Ca],p.H_SERCA,p.n_SERCA);
   I[cIP3] = p.gbar_IP3 * y[g_IP3] * y[h_IP3] * (y[V] - y[V_ER] - Vbar_ER);
   /*
    * double n_Hill=1.0;
    * I[cIP3]   = p.Ihat_IP3  * pow(get_Hill(y[Ca],p.H_IP3_act,n_Hill),p.n_IP3_act)
    * pow(get_Hill(p.H_IP3_inact,y[Ca],n_Hill),p.n_IP3_inact)
    * pow(get_Hill(y[IP3],p.H_IP3_IP3,n_Hill),p.n_IP3_IP3);
    */
}
void cellbeta::get_current_factors(double t) {
   // ++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++
   // Here it is possible to block parts of the currents
   // Specify behind t: time of block, time of release, resting factor (normally 1),
   // target activity (in remaining % of rest factor), steepness of block and release.
   // I[K_ATP]*=get_2sigmoidal(t,20.,45.,1.,0.5,0.2,0.2);
   // I[K_ATP]*=get_2sigmoidal(t,50.,75.,1.,0.4,0.2,0.2);
   // I[K_ATP]*=get_2sigmoidal(t,80.,105.,1.,0.3,0.2,0.2);
   // I[Ca_T]*=get_2sigmoidal(t,3.,80.,1.,0.0,0.2,0.2);
   // I[Ca_L]*=get_2sigmoidal(t,3.,80.,1.,0.2,0.2,0.2);
   // I[PMCA]*=get_2sigmoidal(t,3.,80.,1.,0.7,0.2,0.2);
   // I[NCX]*=get_2sigmoidal(t,3.,80.,1.,0.7,0.2,0.2);
   // I[NaK]*=get_2sigmoidal(t,3.,80.,1.,0.9,0.2,0.2);
   // I[Na_V]*=get_2sigmoidal(t,3.,80.,1.,1.4,0.2,0.2);
   // I[fNa_V]*=get_2sigmoidal(t,20.,120.,1.,2.0,0.2,0.2);
   // I[K_V]*=get_2sigmoidal(t,3.,80.,1.,0.0,0.2,0.2);
   // I[K_Ca]*=get_2sigmoidal(t,3.,80.,1.,0.0,0.2,0.2);
   // I[sK_Ca]*=get_2sigmoidal(t,3.,80.,1.,0.5,0.2,0.2);
   // +++++++++++++++ end OPTION ++++++++++++++++++++++++++++
}
void cellbeta::rhs(double t, double * y, double * derivative) {
   short showit = false;

   // ==============================================
   // Quantities
   // double glucose=get_glucose(t);
   /* ### Note:
    * With a diffusion constant of 6300 micron^2/min the time scale of glucose
    * changes is one second. This corresponds to the external time resolution
    * which is used by hyphasma. On the cellular level the electrophysiology
    * (which calls this routine here) is more on the msecond scale. It is, thus,
    * justified to use the glucose level in the approximation that glucose
    * remains constant during the calculation of the electrophysiology.
    * However, once the diffusion constant is chosen larger or the
    * time resolution of the electrophysiology dynamics is reduced or
    * the time resolution of the external dynamics is reduced, this
    * has to be reconsidered!
    *
    * ATTENTION WHEN CHANGING!
    * Don't forget to also change cellman::mk_single_beta_files when adding output!
    */
   // double IP3=get_IP3(t);

   double x_cal = get_ca_buffer(p.cal_0,y[Ca],p.K_cal);
   double x_buf = get_ca_buffer(p.buf_0,y[Ca],p.K_buf);
   double x_ER = get_ca_buffer(p.buf_ER_0,y[Ca_ER],p.K_buf_ER);
   // cout<<x_cal<<"\n";

   double V_mV = 1000. * y[V];  // sigmoidals use voltage in mV while V is given in Volts

   // ==============================================
   // Reaction equations:
   /* Returns all N rhs in dydt */

   double Vbar_K, Vbar_Na, Vbar_Ca, Vbar_ER;
   get_all_Nernst(t,y[K],y[Na],y[Ca],y[Ca_ER],Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);
   // cout << "Vbar_K=" << 1000.*Vbar_K<<" mV; \n";

   double C_IP3_inh = 0;
   get_currents(y,Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER,C_IP3_inh);
   get_current_factors(t);

   if (showit == true) {
      // all currents in pA/micron^2
      cout << "time=" << t << ": [pA, pA/micron^2]\n";
      cout << "I_NaK=" << I[NaK] << ", " << p.rho[betaWerte::NaK] * I[NaK] << "\n";
      cout << "I_K_ATP=" << I[K_ATP] << ", " << p.rho[betaWerte::K_ATP] * I[K_ATP] << "\n";
      cout << "I_K_V=" << I[K_V] << ", " << p.rho[betaWerte::K_V] * I[K_V] << "\n";
      cout << "I_K_Ca=" << I[K_Ca] << ", " << p.rho[betaWerte::K_Ca] * I[K_Ca] << "\n";
      cout << "I_sK_Ca=" << I[sK_Ca] << ", " << p.rho[betaWerte::sK_Ca] * I[sK_Ca] << "\n";
      cout << "I_Na_V=" << I[Na_V] << ", " << p.rho[betaWerte::Na_V] * I[Na_V] << "\n";
      cout << "I_fNa_V=" << I[fNa_V] << ", " << p.rho[betaWerte::fNa_V] * I[fNa_V] << "\n";
      cout << "I_NCX=" << I[NCX] << ", " << p.rho[betaWerte::NCX] * I[NCX] << "\n";
      cout << "I_PMCA=" << I[PMCA] << ", " << p.rho[betaWerte::PMCA] * I[PMCA] << "\n";
      cout << "I_Ca_L=" << I[Ca_L] << ", " << p.rho[betaWerte::Ca_L] * I[Ca_L] << "\n";
      cout << "I_Ca_T=" << I[Ca_T] << ", " << p.rho[betaWerte::Ca_T] * I[Ca_T] << "\n";
   }

   // cerr<<"y[gapK]="<<y[gap_K]<<"\n";
   // cout<<"p.rho[betaWerte::K_ATP]= "<<p.rho[betaWerte::K_ATP]<<"\n";

   // membrane potential V
   // ++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++
   // Voltage-clamp
   // if (t>3. && t<=10.) { V=-0.04; derivative[0]=0.; } else
   // ++++++++++++++++++++ end OPTION ++++++++++++++++++++++++++
   derivative[V] = -1.0 * (p.rho[betaWerte::NaK] * I[NaK]
                           + p.rho[betaWerte::K_ATP] * I[K_ATP]
                           + p.rho[betaWerte::K_V] * I[K_V]
                           + p.rho[betaWerte::K_Ca] * I[K_Ca]
                           + p.rho[betaWerte::sK_Ca] * I[sK_Ca]
                           + p.rho[betaWerte::Na_V] * I[Na_V]
                           + p.rho[betaWerte::fNa_V] * I[fNa_V]
                           + p.rho[betaWerte::NCX] * I[NCX]
                           + p.rho[betaWerte::PMCA] * I[PMCA]
                           + p.rho[betaWerte::Ca_L] * I[Ca_L]
                           + p.rho[betaWerte::Ca_T] * I[Ca_T]
                           + p.rho[betaWerte::SERCA] * I[SERCA] * xi_S_ERC
                           + p.rho[betaWerte::IP3] * I[cIP3] * xi_S_ERC
                           + y[gap_K] + y[gap_Na] + y[gap_Ca]
                           + J_K + J_Na + J_Ca + J_ions) / p.C_m;

   // K potassium concentration
   derivative[K] = -1.0 * xi * (p.rho[betaWerte::K_ATP] * I[K_ATP]
                                - 2.0 * p.rho[betaWerte::NaK] * I[NaK]
                                + p.rho[betaWerte::K_V] * I[K_V]
                                + p.rho[betaWerte::K_Ca] * I[K_Ca]
                                + p.rho[betaWerte::sK_Ca] * I[sK_Ca]
                                + y[gap_K]
                                + J_K) / MFaraday;

   // Na sodium concentration
   derivative[Na] = -1.0 * xi * (p.rho[betaWerte::Na_V] * I[Na_V]
                                 + p.rho[betaWerte::fNa_V] * I[fNa_V]
                                 + p.rho[betaWerte::NaK] * I[NaK] * 2.0 * p.alpha_NaK
                                 + p.rho[betaWerte::NCX] * I[NCX] * p.alpha_NCX
                                 + y[gap_Na]
                                 + J_Na) / MFaraday;

   // Ca calcium concentration (intracellular)
   derivative[Ca] = -1.0 * xi * (p.rho[betaWerte::Ca_L] * I[Ca_L]
                                 + p.rho[betaWerte::Ca_T] * I[Ca_T]
                                 - p.rho[betaWerte::NCX] * I[NCX] * z_Ca
                                 + p.rho[betaWerte::PMCA] * I[PMCA]
                                 + p.rho[betaWerte::SERCA] * I[SERCA] * xi_ERC / xi
                                 + p.rho[betaWerte::IP3] * I[cIP3] * xi_ERC / xi
                                 + y[gap_Ca]
                                 + J_Ca) / (z_Ca * MFaraday * (1.0 + x_cal + x_buf));

   // open probability of ATP-sensitive K+ channel
   derivative[g_K_ATP]
      = (get_sigmoidal(p.s_h_K_ATP,y[glu],p.kappa_K_ATP) - y[g_K_ATP]) / p.tau_K_ATP;
   // open probability of delayed rectifier K+ channel, V-gated
   double tau_K_V = p.tau_K_V;
   if (p.use_dynamic_tau_K_V == 1) { tau_K_V = get_tau_K_V_sherman88(y[V],Vbar_K); }
   // cout<<"V="<<y[V]<<" -> tau_K_V="<<tau_K_V<<"\n";
   derivative[g_K_V] = (get_sigmoidal(p.V_h_K_V,V_mV,p.kappa_K_V) - y[g_K_V]) / tau_K_V;

   // open probability of V-gated Na+ channel
   double tau_Na_V = p.tau_Na_V;
   if (p.use_dynamic_tau_Na_V == 1) { tau_Na_V = get_tau_V(V_mV,40.,50.,0.0115,-70.); }
   // cout<<"V="<<y[V]<<" -> tau_Na_V="<<tau_Na_V<<"\n";
   derivative[g_Na_V] = (get_sigmoidal(p.V_h_Na_V,V_mV,p.kappa_Na_V) - y[g_Na_V]) / tau_Na_V;
   // open probability of non-inactivating V-gated Na+ channel
   double tau_fNa_V = p.tau_fNa_V;
   if (p.use_dynamic_tau_fNa_V == 1) { tau_fNa_V = get_tau_V(V_mV,40.,50.,0.0115,-70.); }
   // cout<<"V="<<V<<" -> tau_fNa_V="<<tau_fNa_V<<"\n";
   derivative[g_fNa_V] = (get_sigmoidal(p.V_h_fNa_V,V_mV,p.kappa_fNa_V) - y[g_fNa_V]) / tau_fNa_V;
   // open probability of V-gated Ca2+ L-type channel
   derivative[g_Ca_L] = (get_sigmoidal(p.V_h_Ca_L,V_mV,p.kappa_Ca_L) - y[g_Ca_L]) / p.tau_Ca_L;
   // open probability of V-gated Ca2+ T-type channel
   derivative[g_Ca_T] = (get_sigmoidal(p.V_h_Ca_T,V_mV,p.kappa_Ca_T) - y[g_Ca_T]) / p.tau_Ca_T;
   // open probability inactivation of delayed rectifier K+ channel, V-gated
   derivative[h_K_V] = (get_inactivation(p.W_h_K_V,V_mV,p.lambda_K_V) - y[h_K_V]) / p.theta_K_V;
   // open probability of V-gated Na+ channel
   derivative[h_Na_V]
      = (get_inactivation(p.W_h_Na_V,V_mV,p.lambda_Na_V) - y[h_Na_V]) / p.theta_Na_V;
   // open probability of V-gated Ca2+ L-type channel
   derivative[h_Ca_L]
      = (get_inactivation(p.W_h_Ca_L,V_mV,p.lambda_Ca_L) - y[h_Ca_L]) / p.theta_Ca_L;
   // open probability of V-gated Ca2+ T-type channel
   derivative[h_Ca_T]
      = (get_inactivation(p.W_h_Ca_T,V_mV,p.lambda_Ca_T) - y[h_Ca_T]) / p.theta_Ca_T;
   // open probability of Ca2+ and V-gated potassium channel K,Ca
   derivative[g_K_Ca] = (get_sigmoidal(p.V_h_K_Ca,V_mV,p.kappa_K_Ca) - y[g_K_Ca]) / p.tau_K_Ca;
   // feedback of voltage on calcium half opening concentration for K,Ca
   derivative[C_K_Ca] = (get_dynamic_half(V_mV,45.,30.) - y[C_K_Ca]) / p.tau_K_Ca;
   // calcium concentration in ER
   derivative[Ca_ER] = xi_ER * (p.rho[betaWerte::SERCA] * I[SERCA]
                                + p.rho[betaWerte::IP3] * I[cIP3])
                       / (z_Ca * MFaraday * (1.0 + x_ER));
   // potential in ER
   derivative[V_ER] = (p.rho[betaWerte::SERCA] * I[SERCA] + p.rho[betaWerte::IP3] * I[cIP3])
                      / (p.C_m);
   // IP3
   if (p.use_dynamic_IP3 == 1) {
      derivative[IP3] = p.k_IP3_plus * y[g_IP3] * y[h_IP3] * get_Hill(y[Ca],p.C_P,p.n_P)
                        - p.k_IP3_minus * y[IP3];
   } else {
      derivative[IP3] = 0.0;
   }
   // g_IP3
   double tau_IP3 = p.tau_IP3;
   if (p.use_dynamic_tau_IP3 == 1) { tau_IP3 = get_tau_IP3(y[IP3]); }
   derivative[g_IP3]
      = (p.g_IP3_max * get_Hill(y[Ca],p.C_IP3_act,p.n_IP3_act) - y[g_IP3]) / tau_IP3;
   // h_IP3
   derivative[h_IP3] = (get_Hill(C_IP3_inh,y[Ca],p.n_IP3_inh) - y[h_IP3]) / p.theta_IP3;
   // g_sK_Ca
   derivative[g_sK_Ca] = (get_sigmoidal(p.C_sK_Ca,y[Ca],p.kappa_sK_Ca) - y[g_sK_Ca]) / p.tau_sK_Ca;
   // glucose
   derivative[glu] = 0.;
   // gap-junction are calculated outside the cell-object! No changes here:
   derivative[gap_K] = 0.;
   derivative[gap_Na] = 0.;
   derivative[gap_Ca] = 0.;

   if (showit == true) {
      cout << "d/dt of X:\n";
      cout << "V=" << derivative[V] << "\n"
           << "K=" << derivative[K] << "\n"
           << "Na=" << derivative[Na] << "\n"
           << "Ca=" << derivative[Ca] << "\n"
           << "glu=" << derivative[glu] << "\n"
           << "g_K_ATP=" << derivative[g_K_ATP] << "\n"
           << "g_K_V=" << derivative[g_K_V] << "\n"
           << "g_K_Ca=" << derivative[g_K_Ca] << "\n"
           << "C_K_Ca=" << derivative[C_K_Ca] << "\n"
           << "g_sK_Ca=" << derivative[g_sK_Ca] << "\n"
           << "g_Na_V=" << derivative[g_Na_V] << "\n"
           << "g_fNa_V=" << derivative[g_fNa_V] << "\n"
           << "g_Ca_L=" << derivative[g_Ca_L] << "\n"
           << "g_Ca_T=" << derivative[g_Ca_T] << "\n"
           << "h_K_V=" << derivative[h_K_V] << "\n"
           << "h_Na_V=" << derivative[h_Na_V] << "\n"
           << "h_Ca_L=" << derivative[h_Ca_L] << "\n"
           << "h_Ca_T=" << derivative[h_Ca_T] << "\n"
           << "Ca_ER=" << derivative[Ca_ER] << "\n"
           << "V_ER=" << derivative[V_ER] << "\n"
           << "IP3=" << derivative[IP3] << "\n"
           << "g_IP3=" << derivative[g_IP3] << "\n"
           << "h_IP3=" << derivative[h_IP3] << "\n"
           << "\n";
      char cc;
      cin >> cc;
   }
}
void cellbeta::get_gap_junction(dynarray<cellbeta> &bl, space &l) {
   // cerr<<"in get_gap ... ";
   if ((volume == 1) && (target_volume == 1)) {
      // initialise the variables for all ions with the own potential
      y_n[gap_K] = 0.;   // K
      y_n[gap_Na] = 0.;   // Na
      y_n[gap_Ca] = 0.;   // Ca
      // go through all neighbours and sum the contributions:
      for (short i = 0; i < l.dim2; i++) {
         long j = l.knot[index].near_n[i];
         // cerr<<"index="<<index<<"  j="<<j<<"  ";
         if ((j != -1) && (l.cellknot[j].cell == BETA)) {
            long li = l.cellknot[j].listi;
            // cerr<<"li="<<li<<"  ";
            y_n[gap_K] += y_n[V] - bl[li].y_n_old[V] - get_Nernst(y_n[K],bl[li].y_n_old[K],z_K);
            y_n[gap_Na] += y_n[V] - bl[li].y_n_old[V] - get_Nernst(y_n[Na],
                                                                   bl[li].y_n_old[Na],
                                                                   z_Na);
            y_n[gap_Ca] += y_n[V] - bl[li].y_n_old[V] - get_Nernst(y_n[Ca],
                                                                   bl[li].y_n_old[Ca],
                                                                   z_Ca);
         }
      }
      // multiply by the factor which is identical for all ions:
      // double factor=p.gbar_gap*p.rho[betaWerte::gap]/l.dim2;
      double factor = p.gbar_gap * p.rho[betaWerte::gap] / 6.;
      y_n[gap_K] *= factor;
      y_n[gap_Na] *= factor;
      y_n[gap_Ca] *= factor;
      // cerr<<y_n[gap_K]<<"  "<<y_n[gap_Na]<<"  "<<y_n[gap_Ca]<<"\n";
   } else {
      cout << "WARNING! A multi-element object is not defined in cellbeta::get_gap_junction().\n"
           << "         Set all gap-junctions currents to zero.\n";
      y_n[gap_K] = 0.;
      y_n[gap_Na] = 0.;
      y_n[gap_Ca] = 0.;
   }
}
void cellbeta::electrophysiology(double thr, double dthr, sigs &s,
                                 dynarray<cellbeta> &bl, space &l) {
   double t = 3600. * (thr - dthr);
   double tsend;
   // double dt=3600.*dthr;
   // cerr<<"betandt="<<betandt<<"; betadt="<<betadt<<" sec; t="<<t<<" sec.\n";
   for (int n = 0; n < betandt; n++) {
      // cerr<<"thr="<<thr<<" hr; dthr="<<dthr<<" t="<<t<<" s; dt="<<betadt<<" s\n";
      ++step_count;
      // Calculate the gap-junction explicitly:
      //     It is not possible to get the neighbours in the static rhs-routine.
      //     Making rhs non-static does not help
      //     for the neighbour intermediate time points not being stored.
      // ==> Take care that the solver does not make extra time steps!!!
      get_gap_junction(bl,l);   // non-static function defining additional y_n
      // if (index==12) cerr<<"t="<<t<<": gap_K="<<y_n[gap_K]<<" Na="<<y_n[gap_Na]<<"
      // Ca="<<y_n[gap_Ca]<<"\n";
      tsend = t;
      solver.ode_step(y_n, y_n1, N_equations, tsend, betadt, p.dy, prhs, method);
      t += betadt;
      if ((step_count == n_write) && LOCAL_FILES) { show_all(t); step_count = 0; }
   }
}
// ==============================================================================
// ==============================================================================

void cellbeta::get_new_state(const long &i, sigs &s) {
   // Get local glucose concentration
   if (s.signal_use[glucose] == 1) { get_glucose(i,s); }
}
// ============================================================

long cellbeta::ask_mitosis(long * pp, space &l) {
   if ((volume == 1) && (target_volume == 1)) {
      return find_mitosis_place(p_pro,false,max_pro_distance,pp,l);
   } else {
      // case of more fragment object
      // probabilistic decision if proliferation is done
      if ((volume > 0.9 * target_volume) && (drandom() < p_pro)) { return 0; }
      // Proliferation is allowed if the total volume is near the target-volume of a cell!
      return 1;
   }
}
// ============================================================
void cellbeta::set_p_move() {
  /* p_move is a property of the mother class cell.
   * p_move controls the probability of movement of the cell in question.
   * Here, p_move is set for each cell separately.
   * Either the reference value is used or a random variation around the
   * reference value.
   * The reference value is a static property of cellCB and is set with
   * parameters from the parameter file in the cellCB constructor.
   */
  p_move = p_difu;
}
double cellbeta::move(const long &li, space &l, sigs &s, TRACK &td, double &time) {
   // call polarisation
   set_polarity_velocity(persistence,v_modi,p_switch_v,n_v_states,v_slow_factor,l,s);
   if (writethis2track == polarisation) {
      td.Write_movement(trackno,BETA,time,barycenter,polarity,writethis2track);
      writethis2track = trackini;
   }
   return fragmove(BETA,li,diffusion_tolerance_min,diffusion_tolerance_steepness,
                   elongation,smoothmove,p_tension,K_elongation,l);
}
// ==============================================================================
// ==============================================================================

void cellbeta::show(double t, double * y_n, ofstream &file, ofstream &file2) {
   // cout<<t<<"   "<<get_Nernst(y_n[3],p.Ca_ext,z_Ca)-p.Vbar_Ca_delta/1000.<<"\n";
   file << t << "  "
        << y_n[glu] << "  "  // glucose
        << y_n[V] * 1000.0 << "  ";  // voltage is converted from V to mV
   for (int i = K; i <= Ca; ++i) {
      file << y_n[i] << "  ";
   }
   file << 1 - y_n[g_K_ATP] << "  ";
   for (int i = g_K_V; i <= g_Ca_T; ++i) {
      file << y_n[i] << "  ";
   }
   // add some analysis if necessary...
   file << y_n[glu] / p.glu_0 << "  "  // glucose
        << y_n[K] / p.K_0 << "  "
        << y_n[Na] / p.Na_0 << "  "
        << y_n[Ca] / p.Ca_0 << "  ";
   for (int i = h_K_V; i <= h_Na_V; ++i) {
      file << y_n[i] << "  ";
   }
   file << y_n[h_Ca_L] * (1.0 - get_Hill(y_n[Ca],p.C_Ca_L,p.n_Ca_L)) << "  "
        << y_n[h_Ca_T] << "  "
        << get_Ca_free_fraction(y_n[Ca]) << "  ";
   if (p.use_voltage_gating_K_Ca == 1) {
      file << y_n[g_K_Ca] << "  ";
   } else {
      file << "1  ";
   }
   if (p.use_dynamic_H_K_Ca == 1) {
      file << y_n[C_K_Ca] << "  ";
   } else {
      file << p.H_K_Ca << "  ";
   }
   file << y_n[Ca_ER] << "  "
        << y_n[Ca_ER] / p.Ca_ER_0 << "  "
        << get_Ca_ER_free_fraction(y_n[Ca_ER]) << "  "
        << y_n[V_ER] * 1000.0 << "  "
        << y_n[IP3] << "  ";
   if (p.IP3_0 > 0) { file << y_n[IP3] / p.IP3_0 << "  "; } else { file << "0  "; }
   for (int i = g_IP3; i <= h_IP3; ++i) {
      file << y_n[i] << "  ";
   }
   file << I[NaK] / p.Ihat_NaK << "  "  // = (1-H(...))
        << I[NCX] / p.Ihat_NCX << "  "
        << I[PMCA] / p.Ihat_PMCA << "  ";
   if (p.use_dynamic_H_K_Ca == 1) {
      file << get_Hill(y_n[Ca],y_n[C_K_Ca],p.n_K_Ca) << "  ";
   } else {
      file << get_Hill(y_n[Ca],p.H_K_Ca,p.n_K_Ca) << "  ";
   }
   if (p.Ihat_SERCA > 0) { file << I[SERCA] / p.Ihat_SERCA << "  "; } else { file << "0  "; }
   file << "\n";

   file2 << t << "  "
         << y_n[g_sK_Ca] << "  "
         << p.gbar_gap << "  "
         << y_n[g_fNa_V] << "  "
         << "\n";
}
void cellbeta::show_rev(double t, double * v, ofstream &file) {
   file << t << "  "
        << v[V] * 1000.0 << "  ";  // voltage is converted from V to mV

   // Reversal potentials in V
   double Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER;
   get_all_Nernst(t,v[K],v[Na],v[Ca],v[Ca_ER],Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);

   file << 1000. * Vbar_K << "  "
        << 1000. * Vbar_Na << "  "
        << 1000. * Vbar_Ca << "  "
        << 1000. * Vbar_ER
        << "\n";
}
void cellbeta::show_c(double t, ofstream &file) {
   file << t << "   "
        << I[NaK] << "   " << I[K_ATP] << "   " << I[K_V] << "   " << I[K_Ca] << "   "
        << I[Na_V] << "   " << I[NCX] << "   " << I[PMCA] << "   " << I[Ca_L] << "   "
        << I[Ca_T] << "   "
        << I[SERCA] << "   " << I[cIP3] << "   "
        << J_K << "   " << J_Na << "   " << J_Ca << "   " << J_ions << "   "
        << I[sK_Ca] << "   ";
   if (p.rho[betaWerte::gap] > 0) {
      file << y_n[gap_K] * double (lattice_dim) / p.rho[betaWerte::gap] << "  "   // sum(I_gap(K))
           << y_n[gap_Na] * double (lattice_dim) / p.rho[betaWerte::gap] << "  "   // sum(I_gap(Na))
           << y_n[gap_Ca] * double (lattice_dim) / p.rho[betaWerte::gap] << "  ";   // sum(I_gap(Ca))
   } else { file << "0  0  0  "; }
   file << I[fNa_V] << "  ";
   file << "\n";
}
void cellbeta::show_cr(double t, ofstream &file) {
   double surface = 4.0 * pi * pow(p.R_bc,2.0);
   file << t << "   "
        << p.rho[betaWerte::NaK] * I[NaK] * surface << "   "
        << p.rho[betaWerte::K_ATP] * I[K_ATP] * surface << "   "
        << p.rho[betaWerte::K_V] * I[K_V] * surface << "   "
        << p.rho[betaWerte::K_Ca] * I[K_Ca] * surface << "   "
        << p.rho[betaWerte::Na_V] * I[Na_V] * surface << "   "
        << p.rho[betaWerte::NCX] * I[NCX] * surface << "   "
        << p.rho[betaWerte::PMCA] * I[PMCA] * surface << "   "
        << p.rho[betaWerte::Ca_L] * I[Ca_L] * surface << "   "
        << p.rho[betaWerte::Ca_T] * I[Ca_T] * surface << "   "
        << p.rho[betaWerte::SERCA] * I[SERCA] * p.Sur_ER << "   "
        << p.rho[betaWerte::IP3] * I[cIP3] * p.Sur_ER << "   "
        << J_K * surface << "   "
        << J_Na * surface << "   "
        << J_Ca * surface << "   "
        << J_ions * surface << "   "
        << p.rho[betaWerte::sK_Ca] * I[sK_Ca] * surface << "   "
        << y_n[gap_K] * surface << "   "
        << y_n[gap_Na] * surface << "   "
        << y_n[gap_Ca] * surface << "   "
        << p.rho[betaWerte::fNa_V] * I[fNa_V] * surface
        << "\n";
}
void cellbeta::show_gap(double t, double * y, ofstream &file) {
   file << t << "   "
        << y[gap_K] << "   "
        << y[gap_Na] << "   "
        << y[gap_Ca] << "   ";
   for (int i = gap_K_0; i <= gap_K_5; i++) {
      file << y[i] << "   ";
   }
   for (int i = gap_Na_0; i <= gap_Na_5; i++) {
      file << y[i] << "   ";
   }
   for (int i = gap_Ca_0; i <= gap_Ca_5; i++) {
      file << y[i] << "   ";
   }
   file << "\n";
}
void cellbeta::show_all(double t) {
   double Vbar_K, Vbar_Na, Vbar_Ca, Vbar_ER;
   get_all_Nernst(t,y_n[K],y_n[Na],y_n[Ca],y_n[Ca_ER],Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);
   double C_IP3_inh = 0;
   get_currents(y_n,Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER,C_IP3_inh);  // defines I[]
   get_current_factors(t);
   show(t,y_n,result,result2);
   show_c(t,currents);
   show_cr(t,currents_rho);
   show_rev(t,y_n,reversal);
   show_gap(t,y_n,gapfile);
}
void cellbeta::show_all(double t,
                        ofstream &f_a,ofstream &f_b,ofstream &f_i,
                        ofstream &f_r,ofstream &f_n,ofstream &f_g) {
   double Vbar_K, Vbar_Na, Vbar_Ca, Vbar_ER;
   get_all_Nernst(t,y_n[K],y_n[Na],y_n[Ca],y_n[Ca_ER],Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER);
   double C_IP3_inh = 0;
   get_currents(y_n,Vbar_K,Vbar_Na,Vbar_Ca,Vbar_ER,C_IP3_inh);  // defines I[]
   get_current_factors(t);
   f_a << index << "   ";
   f_b << index << "   ";
   f_i << index << "   ";
   f_r << index << "   ";
   f_n << index << "   ";
   f_g << index << "   ";
   show(t,y_n,f_a,f_b);
   show_c(t,f_i);
   show_cr(t,f_r);
   show_rev(t,y_n,f_n);
   show_gap(t,y_n,f_g);
}
// ============================================================
// ============================================================
// ============================================================
// ===========================================================
// ======= Vergleichsoperatoren ==============================
// ============================================================

char operator ==(const cellCB &a, const cellCB &b) {
   char back = 1;
   if ((a.index != b.index) || (a.state != b.state) || (a.pos_ss != b.pos_ss)) {
      back = 0;
   }
   return back;
}
// ============================================================

char operator !=(const cellCB &a, const cellCB &b) {
   return !((a == b) == 1);
}
// ============================================================
