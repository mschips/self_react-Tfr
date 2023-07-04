#include "cellman.h"
#include <math.h>
#include <string.h>

// Zahl der Punkte auf der v-Achse fuer die v-Verteilung der diffusiven Bewegung
int cellman::v_resolution = 100;
// intervall of velocities in microns/min (used for time intervall averaged velocities)
double cellman::delta_v = 2.;
bool cellman::use_antibody_bins = false;
double cellman::mutation_start_time = 0.;
int cellman::mutation_frequency[max_mutation_bin] = { 0 };
double cellman::dec205_max_antigen = 10000.;
double cellman::smooth_stopBCinflux = -1;
bool cellman::do_smooth_stopBCinflux = false;

// Number of cells of different celltypes to be tracked
// If tALL<=0 the individual values are used
// If tALL>0 only the total number of tracked cells is fixed irrespective of their type
int cellman::tALL = 0;
int cellman::tCB = 0;
int cellman::tCC = 0;
int cellman::tOUT = 0;
int cellman::tTC = 0;
int cellman::tBETA = 0;

void cellman::set_tracking(Parameter &p) {
   v_resolution = p.Value.v_resolution;
   delta_v = p.Value.delta_v;
   tALL = p.Value.tALL;
   tCB = p.Value.tCB;
   tCC = p.Value.tCC;
   tOUT = p.Value.tOUT;
   tTC = p.Value.tTC;
   tBETA = p.Value.tBETA;
}
// *************************************
// *************************************

cellman::cellman()
   : CB_list(1000, 1, 0),
   kinetics(),
   CC_list(1000, 1, 0),
   TC_list(200, 1, 0),
   FDC_list(100, 1, 0),
   OUT_list(100, 1, 0),
   STROMA_list(100, 1, 0),
   BETA_list(100, 1, 0),
   TFR_list(200, 1, 0), //msc
   solver(10) {
   // Array for CB velocities and initialization
   velocity = new long[v_resolution + 1];
   delta_velocity = new long[v_resolution + 1];
   for (int t = 0; t <= v_resolution; t++) {
      velocity[t] = 0;
      delta_velocity[t] = 0;
   }
   // trackdatainit() is not called as it is supposed
   // that the other constructor cellman(par,l,s,shape,ana) is used.
   photoactivation_rmin = new long[3];
   photoactivation_rmax = new long[3];
   fateTRACK fT;
   fT.mk_fate_track_file();
   fT.mk_fate_track_legend();
}
cellman::cellman(Parameter &par, space &l, sigs &s, AffinitySpace &shape, ofstream &ana)
   : CB_list(par.Value.CB_Narray, 1, 0),
   kinetics(par.Value.deltat),
   CC_list(par.Value.CC_Narray, 1, 0),
   TC_list(par.Value.TC_Narray, 1, 0),
   FDC_list(par.Value.FDC_Narray, 1, 0),
   OUT_list(par.Value.OUT_Narray, 1, 0),
   STROMA_list(par.Value.STROMA_Narray, 1, 0),
   BETA_list(par.Value.BETA_Narray, 1, 0),
   TFR_list(par.Value.TFR_Narray, 1, 0), //msc
   solver(par.Value.BETA_Narray * cellbeta::N_equations) {
   cout << "Start cellman constructor ...\n";
   // Setze Anzeige-Variablen und Schalter
   checkit = par.Value.safety_checks;
   show_Ki67 = par.Value.show_Ki67;
   show_mode = par.Value.show_mode;
   // Ende

   tcounter = 0;

   ana << "Get time parameters ...\n";
   time = par.Value.tmin;
   last2days = par.Value.tmax - 48.;
   t_dark_end = par.Value.tmax + 10.;
   t_1st_under100 = -1.;
   //MS --> this should be coming from a parameter value
   first_record = 3.0*24;
   time_window = 6.0;
   CC2CBratio = -1.;
   mutation_start_time = par.Value.Start_Mutation;
   dt = par.Value.deltat;
   p_macrophage = par.Value.p_macrophage * par.Value.deltat;
   ana << "  Probability of macrophagocytosis of necrotic cells = " << p_macrophage << "\n";

   // Array for CB velocities and initialization
   velocity = new long[v_resolution + 1];
   delta_velocity = new long[v_resolution + 1];
   delta_velocity_euklid = new long[v_resolution + 1];
   for (int t = 0; t <= v_resolution; t++) {
      velocity[t] = 0;
      delta_velocity[t] = 0;
      delta_velocity_euklid[t] = 0;
   }
   allow_exchange = par.Value.allow_exchange;

   set_tracking(par);
   track_mutations = true;   // ### make parameter in setparam

   tmx.set_parameter(par.Value);

   photoactivation_rmin = new long[l.dim];
   photoactivation_rmax = new long[l.dim];
   if (par.Value.photoactivation) {
      photoactivation = true;
      photoactivation_t0 = par.Value.photoactivation_t0;
      /*     cout<<par.Value.dx<<" r0=("<<par.Value.photoactivation_x0<<","
       *  <<par.Value.photoactivation_y0<<"), dr=("
       *  <<par.Value.photoactivation_delta_x<<","
       *  <<par.Value.photoactivation_delta_y<<")\n";*/
      photoactivation_rmin[0] = long (par.Value.photoactivation_x0 / par.Value.dx);
      photoactivation_rmax[0]
         = photoactivation_rmin[0]
           + long (par.Value.photoactivation_delta_x / par.Value.dx + 0.5) - 1;
      photoactivation_rmin[1] = long (par.Value.photoactivation_y0 / par.Value.dx);
      photoactivation_rmax[1]
         = photoactivation_rmin[1]
           + long (par.Value.photoactivation_delta_y / par.Value.dx + 0.5) - 1;
      if (l.dim == 3) {
         photoactivation_rmin[2] = long (par.Value.photoactivation_z0 / par.Value.dx);
         photoactivation_rmax[2]
            = photoactivation_rmin[2]
              + long (par.Value.photoactivation_delta_z / par.Value.dx + 0.5) - 1;
      }
      ana << "Set photoactivation area to rmin=(" << photoactivation_rmin[0] << ","
          << photoactivation_rmin[1];
      if (l.dim == 3) {
         ana << "," << photoactivation_rmin[2];
      }
      ana << "); rmax=(" << photoactivation_rmax[0] << "," << photoactivation_rmax[1];
      if (l.dim == 3) {
         ana << "," << photoactivation_rmax[2];
      }
      ana << ");\n";
   } else {
      for (short j = 0; j < l.dim; ++j) {
         photoactivation = false;
         photoactivation_rmin[j] = 0;
         photoactivation_rmax[j] = 0;
      }
   }

   //MS
   TFR_mode = par.Value.TFR_mode;
   tfr_depletion_time = par.Value.tfr_depletion_time;
   exp_stop_apo = par.Value.exp_stop_apo;
   exp_CCL3KO = par.Value.exp_CCL3KO;
   CCL3KO = par.Value.CCL3KO;
   //MSchips -- tfr/selfmutation parameter
   /* -1: No SelfReactivity
      0: With probability p_selfMut at mutation
      1: With probability p_selfMut if d(ag,cell)>2
      2: With probability p_selfMut if aff(daughter)<aff(mother) */
   SelfReact_mode = par.Value.SelfReact_mode;
   if (SelfReact_mode >= 0) {
       p_selfMut = par.Value.selfMut_prob;
       cerr <<"P SELF "<<p_selfMut<<"\n";
       if (par.Value.redemption_prob>=1) {
           p_redeem = (p_selfMut/par.Value.redemption_prob);
           cerr <<"P redeem "<<p_redeem<<"\n";
       } else if (/*TFR_mode==10 || TFR_mode==11 ||*/ TFR_mode==12) { //7 trogocytosis
           p_redeem = (1-p_selfMut);
       } else {p_redeem = par.Value.redemption_prob;} //for now use value -3 to avoid inheritance
   } else {
       p_selfMut = -1;
       p_redeem = -1;
   }
   cerr<<"THIS IS THE SelfReact MODEL "<<SelfReact_mode
      <<"\nwith this p_selfMut "<<p_selfMut
     <<"\nand with this p_redeem "<<p_redeem
    <<"\nand THIS IS THE TFR MODEL "<<TFR_mode<<endl;
   asc_got_contact=0;
   newly_generated_redeemed=0;
   newly_generated_sCBs=0;

   def_DEC205 = par.Value.def_DEC205;
   def_DEC205_t0 = par.Value.def_DEC205_t0;
   p_DEC205 = par.Value.p_DEC205;
   inject_antiDEC205OVA = par.Value.inject_antiDEC205OVA;
   inject_antiDEC205OVA_t0 = par.Value.inject_antiDEC205OVA_t0;
   DEC205_ova_activity = par.Value.inject_antiDEC205OVA_t0 + par.Value.antiDEC205OVA_tend;
   TC_factor_dec205ova = par.Value.TC_factor_dec205ova;
   DEC205_induce_differentiation = par.Value.DEC205_induce_CBdifferentiation;
   DEC205_forces_output = par.Value.DEC205_forces_output;
   retain_DEC205_ag = par.Value.retain_DEC205_ag;
   p_CB2OUT = par.Value.CB2OUT_prob;

   // parameters for antigen retention
   retain_ag = par.Value.retain_ag;
   divide_ag_asymmetric = par.Value.divide_ag_asymmetric;

   if (par.Value.antibodies_resolution > 0) {
      use_antibody_bins = true;
   }

   // BrdU staining:
   // tfirst saves the first time of BrdU injection (used to initiate read-out)
   tfirst_inject_BrdU = par.Value.t_inject_BrdU;
   do_inject_BrdU = true;
   if (tfirst_inject_BrdU < 0) { do_inject_BrdU = false; }
   // tnow is the last performed BrdU injection (used for exponential degradation)
   tnow_inject_BrdU = -1;
   // tnext is the next BrdU injection to be done (criterion in time_step)
   tnext_inject_BrdU = tfirst_inject_BrdU;
   // eventually use a time interval for multiple injections
   deltat_inject_BrdU = par.Value.deltat_inject_BrdU;
   n_inject_BrdU = par.Value.n_inject_BrdU;
   BrdU_detection_threshold = par.Value.BrdU_detection_threshold;
   if (do_inject_BrdU) {
     ofstream brdu_cell_t;
     brdu_cell_t.open("brdu_t.out");
     brdu_cell_t << "! Time since injection : Time : Generation : BrdU bin (centre) "
		 << ": total # of cells\n";
     brdu_cell_t.close();
     ofstream brduphases;
     brduphases.open("cellcycle_phases_brdu.out", ofstream::app);
     brduphases << "! Time since injection : Time: #BC with BrdU>threshold "
		<< ": [G1,G0,S,G2,M,CCnotyetselected] : fractions\n";
     brduphases.close();
     ana << "Do inject BrdU at t = " << tnext_inject_BrdU << " hours\n"; 
   }
   else { 
     ana << "No BrdU injections.\n"; 
   }
   ana << "... done\n";
   cout << "... time parameters read\n";

   // ==============================================================

   ana << "Calculate action probabilities ... \n";

   for (int n = 0; n < 4; n++) {
      pp[n] = 0;
   }
   // Zaehler setzen
   CC_total = 0;
   CB_average = 0;
   CB_variance = 0;
   CB_end = 0;
   n_outs = 0;
   n_outs_with_recycling = 0;
   n_muts = 0;
   n_recmuts = 0;
   n_recycling_events = 0;
   n_recycling_events_last = 0;
   for (int n = 0; n < max_mutation_bin; n++) {
      mutation_frequency[n] = 0;
   }

   // define static variables in cell classes:
   /*multiAg: shifted this to hyphasma-main in order to set some variables before the call
    * of the cell-specific constructors when the cellman-class is constructed. */
   /*
    * cell::set_statics(par, ana);
    * immunoglobulin_class::load_matrix(par, ana); // ### didn't check whether shifting this one
    * makes problems
    * cellCB::set_statics(par, l, ana);
    * cellCC::set_statics(par, l, ana);
    * cellTC::set_statics(par, l, ana);
    * cellOUT::set_statics(par, l, ana);
    * cellFDC::set_statics(par, l, ana);
    * cellbeta::set_statics(par, l, ana);
    */
   set_pars(par, l);   // sets parameters dependent on time (can be repeated)
   // CXCL12-production is done without explicitly introducing the producing stroma cells:
   if (par.Value.mkCXCL12 <= 0.) {
      p_mkCXCL12 = 0.;
   } else {
       p_mkCXCL12 = par.Value.mkCXCL12                                       // mol/(cell l hr)
	 * par.Value.deltat                                       // hr
	 * par.Value.dx * par.Value.dx * par.Value.dx * 1.e-15    // l
	 * par.N_A;                                               // /mol
       // one may think of a slice of thickness <lattice-constant> in 2D, thus valid for 2 and 3D.
   }
   ana << "STROMA CXCL12 production (molecules) = " << p_mkCXCL12 << "\n";

   // calculated the maximum of possible antigen collection ...
   if (cellCC::test_delay > 0) {
      dec205_max_antigen = cellCC::collectFDCperiod / (double (cellCC::test_delay) * dt);
   } else if ((cellCC::test_delay == -1) && (cellCC::ICAM_delay > 0)) {
      dec205_max_antigen = cellCC::collectFDCperiod / (double (cellCC::ICAM_delay) * dt);
   } else {
      dec205_max_antigen = cellCC::collectFDCperiod / dt;
   }
   // ... and multiply this with 3 according to Victora et al 2010
   // dec205_max_antigen*=3.0;
   // ... I removed this multiplication because the reality of antigen collection
   //     already induces a factor of three difference as seen in Victora et al 2010.
   /*
    * cout<<"dt="<<dt
    *  <<"; collect-dt="<<cellCC::collectFDCperiod
    *  <<"; test-dt="<<cellCC::test_delay
    *  <<"; ICAM-dt="<<cellCC::ICAM_delay
    *  <<"; max-ag="<<dec205_max_antigen<<"\n";
    * char cc; cin>>cc;
    */
   // for the usage in the case of DEC205 targeting.

   ana << "... done.\n";
   cout << "... probabilities defined\n";

   // =================================================================

   // Initialize Cell-lists
   ana << "Initialize the cell-lists ...";
   CB_list.setstep(1000);
   CC_list.setstep(1000);
   TC_list.setstep(200);
   TFR_list.setstep(200); //msc
   FDC_list.setstep(100);
   OUT_list.setstep(100);
   BETA_list.setstep(100);
   BC_integral = 0;
   ana << "done\n";
   cout << "... cell lists initialized\n";

   // create stroma cells:
   if (p_mkCXCL12 > 0.) {
      // make 30 (2D) or 300 (3D)
      int stroma_max = 30;
      if (l.dim == 3) {
         stroma_max = 300;
      }
      ana << "Put " << stroma_max << " stromal cells ...";
      cout << "Put " << stroma_max << " stromal cells ...";
      // ++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++
      // if true stroma are distributed in DZ
      // if false stroma are set at the border of the DZ to T zone
      // both only happens if stroma cells produce CXCL12 (see if)
      bool put_in_DZ = true;
      if (put_in_DZ) {
         ana << " ... distribute STROMA in DZ ... ";
         cout << "STROMA are distributed in the whole DZ.\n";
      } else {
         ana << " ... distribute STROMA on the border of the DZ to the T zone ... ";
         cout << "STROMA are distributed on the DZ-TZ boundary.\n";
      }
      // ++++++++++++++++ end OPTION +++++++++++++++++++++++++++++
      long k[l.dim];
      while (STROMA_list.benutzt() < stroma_max) {
         // randomly chose an x-component
         for (int i = 0; i < l.dim - 1; i++) {
            k[i] = irandom(l.prodimvec[i]);
         }
         if (put_in_DZ) {
            // Put anywhere in DZ
            k[l.dim - 1] = irandom(int (par.Value.FDCnetwork * (l.prodimvec[l.dim - 1])));
            // Now this is a value in the LZ --> shift to DZ
            k[l.dim - 1] += int ((l.prodimvec[l.dim - 1] / 2.));  
	    // ... resulting value is at least radius
            // Shift a little into the DZ:
            if (k[l.dim - 1]
                > int ((double (l.prodimvec[l.dim - 1]) / 2.)
                       + (0.1 * double (l.prodimvec[l.dim - 1])))) {
               long getpos = l.Index(k);
               if (l.knot[getpos].status != external) {
                  STROMA_list.add(getpos);
               }
            }
         } else {
            // put at south border of DZ
            k[l.dim - 1] = int (par.Value.FDCnetwork * l.prodimvec[l.dim - 1] + 0.5);
            long getpos = l.Index(k);
            short found = 0;      // 0:not found, 1:found target, 2:external
            while (found == 0) {
               // go down until a neighbour is external or
               if (l.knot[getpos].status == external) {
                  found = 2;
               } else if ((l.knot[getpos].near_n[2 * l.dim - 1] == -1)
                          || (l.knot[l.knot[getpos].near_n[2 * l.dim - 1]].status
                              == external)) {
                  found = 1;
               } else {
                  getpos = l.knot[getpos].near_n[2 * l.dim - 1];
               }
            }
            if (found == 1) {
               STROMA_list.add(getpos);
            }
         }
      }
      ana << " done.\n";
   }

   // create the FDCs (in the lower sphere), i.e. k[3]<prodim/2
   ana << "Put " << par.Value.FDCnumber << " FDCs ...";
   cout << "Put " << par.Value.FDCnumber << " FDCs ...";
   FDCtransparent = par.Value.FDCtransparent;
   int nFDC = 0;
   long getpos;
   short getrandompos = 0;
   long k[l.dim];
   while (nFDC < par.Value.FDCnumber) {
      if (getrandompos == 1) {
         getpos = -1;
      } else {
         getpos = par.Value.posFDC[nFDC];
      }
      if (getpos == -1) {
         getrandompos = 1;     // Now, only random FDCs are used
         for (int i = 0; i < l.dim - 1; i++) {
            k[i] = irandom(l.prodimvec[i]);
         }
         k[l.dim - 1] = irandom(int (par.Value.FDCnetwork * (l.prodimvec[l.dim - 1] + 1)));
         getpos = l.Index(k);
      }
      if (put_FDC(getpos, l, shape, par.Value.FDCtransparent) == 0) {
         if (getrandompos == 0) {
            ana << " Fix:";
         } else {
            ana << " Random:";
         }
         ana << " " << getpos << " ";
         ++nFDC;
      } else {
         getrandompos = 1;
      }
   }

   // Go through the lattice and define GCpointnumber, FDCpointnumber
   GCpointnumber = 0;
   FDCpointnumber = 0;
   for (int n = 0; n < l.pointnumber; n++) {
      if (l.knot[n].status != external) {
         ++GCpointnumber;
      }
      if (l.cellknot[n].cell == FDC) {
         ++FDCpointnumber;
      }
   }
   TestedFDCsiteVoidOfAg = 0;
   ana << " done.\n";
   ana << "Total number of lattice points for GC = " << GCpointnumber << "\n";
   ana << "Total number of lattice points occupied by FDC = " << FDCpointnumber << "\n";
   cout << "... FDCs created\n";

   if (par.Value.total_blast2 > 0) {
      ana << "Put " << par.Value.total_blast2 << " naive BC (blast2) ...";
      int nblast2 = 0;
      long getpos;
      long k[l.dim];
      while (nblast2 < par.Value.total_blast2) {
         for (int i = 0; i < l.dim - 1; i++) {
            k[i] = irandom(l.prodimvec[i]);
         }
         k[l.dim - 1] = irandom(int (par.Value.FDCnetwork * (l.prodimvec[l.dim - 1] + 1)));
         getpos = l.Index(k);
         // cout<<"("<<k[0]<<","<<k[1]<<","<<k[2]<<"):"<<l.cellknot[getpos].cell<<" -- ";
         if (l.cellknot[getpos].cell == nocell) {
            l.set_knot(getpos, blast2, nblast2);
            ++nblast2;
            // if (int(n/100)*100==n) cout<<"... "<<n;
         }
      }
      ana << " done.\n";
      cout << "... naive BC distributed on lattice\n";
   }

   // put seeder cells on the lattice
   ana << "Put " << par.Value.totalB << " seeder cells ...";
   // cout << "Put " << par.Value.totalB << " seeder cells ...";
   n_founder = 0;   // this counter is further used in cellman::BC influx()
   getrandompos = 0;
   // Variable for seeder cells
   // Determine the number of divisions (if appropriate):
   determine_founder_times_of_division(par.Value);
   // determine the new seeder cells
   // do this even when influx of new BC is used later on
   while (n_founder < par.Value.totalB) {
     cellCB seederCB;
      // cout<<"Seeder cell "<<n<<": ";
      if (getrandompos == 2) {
         long getposnew = l.next_border(getpos);
         getpos = -1;
         for (int alpha = 0; alpha < l.dim2; alpha++) {
            if (l.cellknot[l.knot[getposnew].near_n[alpha]].cell == nocell) {
               getpos = l.knot[getposnew].near_n[alpha];
            }
         }
      } else if ((getrandompos == 1) || (n_founder >= MAXDIM)) {
         getpos = -1;
      } else {
         getpos = par.Value.posCB[n_founder];
      }
      if (getpos == -1) {
         for (int i = 0; i < l.dim; i++) {
            k[i] = irandom(l.prodimvec[i]);
         }
         // k[1]=irandom((l.prodimvec[1]+1)/4);
         // k[l.dim-1]=irandom(int((l.prodimvec[l.dim-1]+1)/2));
         getpos = l.Index(k);
      }
      // use the list of seeder cells
      // if more cells are chosen than the seeder cell list is long
      // use random seeder cells within the same list
      seederCB.pos_ss = shape.get_Seeder(int (n_founder));
      // cout<<seederCB.pos_ss<<"\n";
      seederCB.make_CB_new();
      if (cellCB::fixed_number_of_divisions()) {
	if (par.Value.fixed_time_of_divisions_mode == 2) {
	// +++++++++++++ OPTION +++++++++++++++++++++++++++++++++++
	// define whatever you want
	/* The following was used in the Figure S1 of Cell Reports 2012.
	 * Refer to hyphasma12.03.2 for the original simulation.
	 * Parameter file is brdu-ag006.par.
	 */
	double xdivp=drandom(2.0*par.Value.CB_fixed_times_of_divisions);
	double xndiv = int (xdivp);
	xdivp -= (double (xndiv));    // prob for rest
	if (drandom()<xdivp) ++xndiv;
	seederCB.n_divisions2do = xndiv;
	// +++++++++ end OPTION +++++++++++++++++++++++++++++++++++
	} else {
	  seederCB.n_divisions2do = get_founder_times_of_division();
	}
	if (seederCB.n_divisions2do == 0) {
	  seederCB.state = cb_stop_dividing;
	}
      }
      if (s.signal_use[CXCL12]) {
         seederCB.responsive2signal[CXCL12] = true;
      }
      seederCB.preload_with_ag();
      /* If DEC205 is used, attribute it with the respective probability to the new BC;
       * only do so if the time of DEC205 attribution is set negative;
       * otherwise it is done in a shot at this time (cellman::attribute_DEC205();). 
       */
      if (def_DEC205 && def_DEC205_t0 < 0) {
	seederCB.attribute_DEC205(p_DEC205);
      }
      if (put_CB(getpos, seederCB, l, shape) != -1) {
         if ((n_founder < MAXDIM) && (par.Value.posCB[n_founder] != -1) && (getrandompos == 0)) {
            ana << " Fix:";
         }
         // adapt the average seeder cell affinity
         cellCB::average_seeder_affinity
            = (n_founder * cellCB::average_seeder_affinity
               + shape.best_affinity_norm(seederCB.pos_ss))
              / (n_founder + 1);
         // cout<<"cellCB::average_seeder_affinity="<<cellCB::average_seeder_affinity<<"\n";
         ana << " " << getpos << " ";
         getrandompos = 0;
         ++n_founder;
      } else {
         // OPTION +++
         // chose 1 if random position to be used
         // chose 2 if the next nearest position is to be used
         // getrandompos=1;
         getrandompos = 2;
         // end OPTION +++
      }
   }
   if (track_mutations) {
      // write the new cell as founder cell to the brainbow class
      // parameters: (double time, bool founder, bool birth, long mother_index, long ss_position)
      for (long i = 0; i < CB_list.benutzt(); i++) {
         CB_list[i].brainbow_index
            = trackmutations.write(time, true, true, -1, CB_list[i].pos_ss);
      }
   }
   // fate tracker:
   if (par.Value.write_trackfate > 0) {
     fateTRACK fT;
     fT.mk_fate_track_file();
     fT.mk_fate_track_legend();
   }
   // Fix the probability of new founder cell influx
   p_BCinflux = par.Value.newBCinflux_rate * par.Value.deltat;
   t_stopBCinflux = par.Value.newBCinflux_stop;
   smooth_stopBCinflux = par.Value.smooth_stopBCinflux;
   do_smooth_stopBCinflux = false;
   if (smooth_stopBCinflux > 0) {
      do_smooth_stopBCinflux = true;
   }

   // Einschraenkung der seeder cells auf FDC-network?
   /*
    * while (n<par.Value.totalB) {
    * for (int i=0; i<l.dim; i++) k[i]=irandom(l.prodimvec[i]);
    * k[l.dim-1]=irandom(int((l.prodimvec[l.dim-1]+1)/2));
    * cout<<"k=("<<k[0]<<" "<<k[1]<<") Index="<<l.Index(k)<<"\n";
    * if ( put_CB(Index(k),shape.get_Seeder(int(n)),0,0,shape)==0 ) n++;
    * }
    */
   // Set cellOUT-parameters if needed:
   if ((par.Value.use_ab_dynamics > 0) && (par.Value.initial_ab_affinity < 0)) {
      cellOUT::average_affinity = cellCB::average_seeder_affinity;
      // #### this is wrong because max_affinity shoud be set to the best affinity 
      // within the seeder cells
      cellOUT::max_affinity = cellCB::average_seeder_affinity;
      // #### better is the following (activate later, don't want to change too many things at a
      // time):
      // double aff=shape.best_affinity_norm(newCB.pos_ss);
      // if (aff>cellOUT::max_affinity) cellOUT::max_affinity=aff;
      // #### note that this has to be repeated for each new cell, so shift it above!
   }

   // Add extrafollicular output cells:
   if (par.Value.add_extrafollicular_PC > 0) {
     for (int n = 0; n < par.Value.add_extrafollicular_PC; n++) 
       { add_extrafollicular_PC(shape, par.Value.pos_extrafollicular_PC); }
   }
   ana << " done.\n";
   cout << "... seeder cells incorporated\n";

   /* If BrdU is set at the beginning of the simulation,
      this run is meant to do an in vitro experiment with all BCs dividing
      and stained by BrdU (used in Cell Reports 2012).
   */
   if (tfirst_inject_BrdU < 0) { label_BrdU(true); }
   

   // +++ OPTION restrict TC position in space to the upper part of the volume:
   addTC(par.Value.totalTC, l, 2.0, shape, ana);   // restrict to upper 1/"2" space part
   // +++ end OPTION
   /*
    * ana << "Put " << par.Value.totalTC << " T cells ...";
    * n=0;
    * cellTC seederTC;   // Variable for T cells
    * while (n<par.Value.totalTC) {
    * for (int i=0; i<l.dim; i++) k[i]=irandom(l.prodimvec[i]);
    * // +++ OPTION restrict TC position in space to the upper part of the volume:
    * k[l.dim-1]=irandom(int((l.prodimvec[l.dim-1]+1)/3)); // ##### remove +1 #####
    * // +++ end OPTION
    * getpos=l.Index(k);
    * // Assume TC to be specific for the Antigen!
    * seederTC.pos_ss=shape.get_Antigen();
    * if (put_TC(getpos,seederTC,l,shape)!=-1) {
    *  ana<<" "<<getpos<<" ";
    * ++n;
    * }
    * }
    * ana << " done.\n";
    */
   cout << "... T cells incorporated\n";
   //MSchips
   addTFR(par.Value.TFR_mode,par.Value.num_TFR, l , 2.0 , shape , ana);
   cout << par.Value.num_TFR << " ... TFR cells incorporated\n";
   writeZPositions(time,l);

   // initialise the counters of class-specific output cells
   for (int i = 0; i <= nIg_classes; i++) {
     integrate_out[i] = 0;
     integrate_out_old[i] = 0;
   }

   // ===========================================================
   // CREATE BETA-CELLS

   // cerr<<"Start BETA ...\n";
   // put seeder cells on the lattice
   ana << "Put " << par.Value.BETA_Nini << " beta-cells ...";
   // cout << "Put " << par.Value.BETA_Nini << " beta-cells ...";
   int nbeta = 0;

   getrandompos = 0;
   // Variable for seeder cells
   // cerr<<"before seederBETA...\n";
   // cerr<<"first BETA call done.\n";
   if (par.Value.BETA_Nini >= l.pointnumber) {
      cout << "... fill the lattice with beta-cells ";
      ana << " -> Fill the lattice with beta-cells ... ";
      for (long i = 0; i < l.pointnumber; i++) {
	 cellbeta seederBETA;
         if (cellbeta::p.randomise_beta_proteins) {
            // randomise the protein expression level in the cell
            // to restrict this to specific proteins do changes
            // in cellbeta::randomise_protein_expression() in cellthis.C
            seederBETA.randomise_protein_expression();
            for (short j = 0; j < betaWerte::N_beta_proteins; j++) {
               ana << "rho[" << j << "]= " << seederBETA.rho[j] << "; ";
               cout << "rho[" << j << "]= " << seederBETA.rho[j] << "; ";
            }
            ana << "\n";
            cout << "\n";
         }

         put_BETA(i, seederBETA, l);
         // for (int ii=0; ii<cellbeta::N_equations; ii++)
         // cout<<"y_n["<<ii<<"]="<<seederBETA.y_n[ii]<<"; ";
         // char cc; cin>>cc;
      }
      cout << "\n";
   } else {
      while (nbeta < par.Value.BETA_Nini) {
         cellbeta seederBETA;
         // cout<<"Seeder cell "<<nbeta<<": ";
         if (getrandompos == 2) {
            long getposnew = l.next_border(getpos);
            getpos = -1;
            for (int alpha = 0; alpha < l.dim2; alpha++) {
               if (l.cellknot[l.knot[getposnew].near_n[alpha]].cell == nocell) {
                  getpos = l.knot[getposnew].near_n[alpha];
               }
            }
         } else if ((getrandompos == 1) || (nbeta >= MAXDIM)) {
            getpos = -1;
         } else {
            getpos = par.Value.BETA_pos[nbeta];
         }
         if (getpos == -1) {
            for (int i = 0; i < l.dim; i++) {
               k[i] = irandom(l.prodimvec[i]);
            }
            getpos = l.Index(k);
         }
         if (cellbeta::p.randomise_beta_proteins) {
            seederBETA.randomise_protein_expression();
            for (short j = 0; j < betaWerte::N_beta_proteins; j++) {
               ana << "rho[" << j << "]= " << seederBETA.rho[j] << "; ";
               cout << "rho[" << j << "]= " << seederBETA.rho[j] << "; ";
            }
            ana << "\n";
            cout << "\n";
         }
         if (put_BETA(getpos, seederBETA, l) != -1) {
            if ((nbeta < MAXDIM) && (par.Value.BETA_pos[nbeta] != -1) && (getrandompos == 0)) {
               ana << " Fix:";
            }
            ana << " " << getpos << " ";
            getrandompos = 0;
            ++nbeta;
         } else {
            // OPTION +++
            // chose 1 if random position to be used
            // chose 2 if the next nearest position is to be used
            // getrandompos=1;
            getrandompos = 2;
            // end OPTION +++
         }
      }
   }
   ana << " done.\n";
   if (par.Value.show_mode == islet) {
      cout << "... initial beta-cells incorporated\n";
   }
   // for global Runge-Kutta only:
   method = RungeKutta_4th;
   prhs = &rhs;
   if (par.Value.show_mode == islet) {
      if (cellbeta::FULL_RUNGE) {
         cout << "Use consistent implicit solution of the betacell ode-system.\n";
         ana << "Use consistent implicit solution of the betacell ode-system.\n";
      } else {
         cout << "WARNING: Use inconsistent implicit solution of the betacell ode-system.\n"
              << "         with gap-junctions treated explicitly.\n";
         ana << "WARNING: Use inconsistent implicit solution of the betacell ode-system.\n"
             << "         with gap-junctions treated explicitly.\n";
      }
   }

   // ===========================================================

   // save the new variables as standard
   // ana<<"Actualize changes ...";
   // for (int n=0; n<l.pointnumber; n++) l.knot[n].actualize();
   // ana << " done.\n";

   outputfiles = par.Value.outputfiles;
   ignore_apoptotic_CC = par.Value.ignore_apoptotic_CC;

   ana << "Write initial distribution to files ...\n";
   char ini[5] = "_ini";
   xfiles(ini, l);
   shape.write_gcbc_hamming(0);
   shape.write_gcbc_affinity(0);
   if (outputfiles < 2) {
      s.write_files(ini, true);
      zone_files(ini, l);
   }
   // cerr<<"Wrote files!\n";

   movie.open("movie");
   movie << "nice -19 gifsicle --delay 10 --scale " << 750. / double (l.prodimvec[1])
         << " $1\"_ini.gif\" ";
   if (show_mode != islet) {
      // Note that the scale is defined with respect to the size of the 1-direction
      // because the 1-direction is shown in 2D (vertical) and 3D (horizontal axis)!
      // movie<<"nice -19 gifsicle --delay 10 --scale 1 $1\"_ini.gif\" ";
      // mkgifs.open("mkgifs");
      // mkgifs<<"ppmtogif xy_ini.ppm > xy_ini.gif\n";

      dark_zone.open("zone_da.out");
      dark_zone << "! Ratios and absolute numbers of cell types in the dark zone:\n";
      dark_zone << "! t #CB/#all #notrecycledCB/#all #CC/#all #CB/#CC #nrCB/#CC "
                << "#all #CB #CC #nrCB\n";
      light_zone.open("zone_li.out");
      light_zone << "! Ratios and absolute numbers of cell types in the light zone:\n";
      light_zone << "! t #CB/#all #notrecycledCB/#all #CC/#all #CB/#CC #nrCB/#CC "
                 << "#all #CB #CC #nrCB\n";

      xsums[CB].open("xsumcb.out");
      xsums[CB] << "! Summe aller Centroblasten\n";
      xsums[CB] << "! Time : all : not recycled : total#recyclings : new#recyclings :"
                << " fraction of recycled : current fraction of recycled :"
                << " (recycled once) : (recycled often)\n";
      xsums[CC].open("xsumcc.out");
      xsums[CC] << "! Summe aller Centrocyten\n";
      xsums[CC] << "! Time : all : not selected : contact2FDC : FDC-selected : "
                << "contact2TFH : selected : apoptotic\n";
      xsums[out].open("xsumout.out");
      xsums[out] << "! Sum of all output cells in the GC area\n"
                 << "! time : gc-output : in DZ : in LZ : fraction of GC cells\n";
      xvolume.open("xvolume.out");
      xvolume
         << "! cellsum time : notempty : notempty-FDC : CB+CC : frac notempty :"
	 << " frac notempty-FDC : frac notempty-FDC-out : CB+CC+out : BC integral\n";
      xsumBCig.open("xsumbc_ig.out");
      xsumBCig << "! time : CB+CC : DZ : LZ : DZ/LZ: IgM : %IgM : IgM DZ : IgM LZ : IgM DZ/LZ : "
               << "IgG... : IgE... : IgA...\n";
      xsumBCOig.open("xsumbco_ig.out");
      xsumBCOig
         <<
         "! time : CB+CC+OUT : DZ : LZ : DZ/LZ : IgM : %IgM : IgM DZ : IgM LZ : IgM DZ/LZ : "
         << "IgG ... : IgE... : IgA...\n";
      cum_vol.open("cum_vol.txt", std::ios_base::app);
      xapoig.open("xapo_ig.out");
      xapoig << "! time : all apoptotic on lattice : % of all BC (CB+CC+OUT) "
             << ": apoptotic IgM : % of all IgM BC : ... classes ...\n";
      xapoigdz.open("xapo_ig_dz.out");
      xapoigdz << "! time : all apoptotic in DZ : % of all DZ-BC (CB+CC+OUT) "
               << ": apoptotic IgM : % of all IgM BC : ... classes ...\n";
      xapoiglz.open("xapo_ig_lz.out");
      xapoiglz << "! time : all apoptotic in LZ : % of all LZ-BC (CB+CC+OUT) "
               << ": apoptotic IgM : % of all IgM BC : ... classes ...\n";
      xaffig.open("xaff_ig.out");
      xaffig << "! time : all CB : all CC : all OUT : all CB+CC : all CB+CC+OUT "
             << ": IgM CB : IgM CC : IgM OUT : IgM CB+CC : IgM CB+CC+OUT : ... classes ...\n";

      log_out_aff.open("xlogouts.out");
      log_out_aff << "! Protokol of all produced output-cells and affinity to antigen\n";

      log_bc_aff.open("xlogbcaff.out");
      log_bc_aff << "! time : list of affinities to all Ags at this time "
                 << "(best : Ag0 : ... : AgN)\n";

      ag_presentation.open("antigen.out");
      ag_presentation << "! Time course of antigen presentation on FDC\n";
      ag_presentation
         <<
         "! time & free on FDC & ic on FDC & ab on lattice & free (Mol) & ic (Mol) & ab (Mol)\n";

      ofstream cycle_phases("cellcycle_phases.out");
      cycle_phases << "! time[hr] number of CBs (+selected CCs in delay mode)"
		   << "[normal,differentiate,G1,G0,S,G2,M,divide,stop_dividing]\n";
      cycle_phases.close();
      ofstream cycle_phases_zones("cellcycle_phases_zones.out");
      cycle_phases_zones << "! time[hr] total number of dividing cells,"
			 << "[normal,differentiate,G1,G0,S,G2,M,divide,stop_dividing] "
			 <<": in DZ : in LZ : fraction of DZ of all in this phase\n";
      cycle_phases_zones.close();

      ofstream K("TFHsignalintensityK.out");
      K << "! time[hr] : mean : sd : n\n";
      K.close();

      ofstream antigens_out("antigens.out");
      antigens_out
         << "! Free antigen presented on FDCs for all antigen types separately (columns)"
	 << " [# of ag_portions] \n";
      antigens_out.close();

      ofstream ics_out("ics.out");
      ics_out << "! Immune complexes on FDCs for all antigen types separately (columns) [Mol] \n";
      ics_out.close();

      ofstream icsfrac_out("ics_frac.out");
      icsfrac_out << "! Fraction of immune complex versus total antigen"
		  <<" for each antigen type separately (columns)\n";
      icsfrac_out.close();

      cum_ag.open("cum_ag.txt", std::ios_base::app);

      xtc.open("tc.out");
      xtc << "! t  .  #tc  .  #tc0bc  .  #tc1bc  .  #tc2bc  ....  #tc6bc\n";

      xdec205.open("dec205.out");
      xdec205 << "! t . #dec+ . #dec- . #dec+/#dec+- . "
              << "#dec+dz . #dec-dz . #dec+dz/#dec+-dz . "
              << "#dec+lz . #dec-lz . #dec+lz/#dec+-lz . "
              << "#dec+CB . #dec-CB . #dec+CB/#dec+-CB . "
              << "#dec+CC . #dec-CC . #dec+CC/#dec+-CC . "
              << "#dec+CCsel . #dec-CCsel . #dec+CCsel/#dec+-CCsel . "
              << "#dec+CCapo . #dec-CCapo . #dec+CCapo/#dec+-CCapo . "
              << "#dec+CCapo-dz . #dec-CCapo-dz\n";
      ofstream MHCdeficient_out("MHCdeficient.out");
      MHCdeficient_out << "! time : #BC : #CB : #CC : %BC : %CB : %CC "
		       << "of MHC-deficient in DEC205+\n";
      MHCdeficient_out.close();

      aff_cb_lz.open("aff_cb_lz.out");
      aff_cb_lz << "! t . LZ-CB-bin0 . LZ-CB-DEC+-bin0 . ... \n";
      aff_cc_lz.open("aff_cc_lz.out");
      aff_cc_lz << "! t . LZ-CC-bin0 . LZ-CC-DEC+-bin0 . ... \n";
      aff_out_lz.open("aff_out_lz.out");
      aff_out_lz << "! t . LZ-OUT-bin0 ... \n";
      aff_cb_dz.open("aff_cb_dz.out");
      aff_cb_dz << "! t . DZ-CB-bin0 . DZ-CB-DEC+-bin0 . ... \n";
      aff_cc_dz.open("aff_cc_dz.out");
      aff_cc_dz << "! t . DZ-CC-bin0 . DZ-CC-DEC+-bin0 . ... \n";
      aff_out_dz.open("aff_out_dz.out");
      aff_out_dz << "! t . DZ-OUT-bin0 ... \n";
   }

   if ((show_mode != islet) && (outputfiles == 0)) {
      // Write for movement analysis only
      movements.open("movement.out");
      movements << "! Bewegungsschritte von Zellen pro Zeitschritt dt=" << dt << " min\n"
                << "! in Einheiten von dx=" << l.dx << " microns\n"
                << "! v_center=moveOF*dx/dt in Einheiten von microns/min\n"
                << "! pos in Einheiten von microns\n"
                << "! D in Einheiten von microns^2/min\n"
                << "! t : sqrt(t) : cellvol : moveOFbarycenter : v_center : pos(born_pos) : "
                << "Diffconstant : rel.coords\n";
      move_vdt.open("move_vdt.out");
      // move_vdt<<"! Bewegungsschritte von Zellen pro Zeitintervall="<<CB_list[0].deltat_v<<"
      // min\n"
      move_vdt << "! Bewegungsschritte von Zellen pro Zeitintervall=" << cell::deltat_v
               << " min\n"
               << "! in Einheiten von dx=" << l.dx << " microns\n"
               << "! v_center=moveOF*dx/dt in Einheiten von microns/min\n"
               << "! pos in Einheiten von microns\n"
               << "! D in Einheiten von microns^2/min\n"
               << "! t : sqrt(t) : cellvol : moveOFbarycenter : v_center_real : "
	       << "v_center : pos(born_pos) : Diffconstant : rel.coords\n";
      axis.open("axis.out");
      axis << "! The ratio of long to short axis of the cell\n"
         // <<"! measured in regular time intervalls of "<<CB_list[0].deltat_v<<" min\n"
           << "! measured in regular time intervalls of " << cell::deltat_v << " min\n"
           << "!\n"
           << "! time [min] :: long axis [microns] :: short axis [microns] :: long2short_axis\n";
      polarities.open("polarities.out");
      polarities << "! times of polarity changes\n"
		 << "! cell#   time[min]   0\n";
      corrs.open("corrs.out");
   }
   //
   if (outputfiles == 0) {
      prolog.open("prolog.out");
      prolog << "! time  #prol_nn  #prol_diagn  #prol_othern  failures  distance\n";
      /*
       * cbquality.open("cbqual.out");
       * cbquality<<"! Quality of all centroblasts on the lattice\n";
       * cbquality<<"time   index   ss-index   best-affinity\n";
       */
   }
   //
   if (show_mode != islet) {
     cbhighaff.open("cbha.out");
     cbhighaff << "! Quality of all, non-recycled, recycled CBs in %\n";
     cbhighaff << "time   all   non-recycled   recycled\n";
     
     mutation_out.open("mutation.out");
     mutation_out << "! Number of mutations in produced output cells\n";
     mutation_out << "time  #outs  total#muts  #mutsonrecs  c3/#outs  "
		  << "c4/#outs  #outswithrecycle  c4/#outswithrecycle\n";
     
     ag_collected.open("ag_collected.out");
     ag_collected
       << "! Average number of collected antigen portions by FDCselected and TCcontact CC (1).\n"
       << "!                or in addition of selected and apoptosis CC (2).\n"
       << "! time    N(1)  (1)  sd(1)     N(2)  (2)  sd(2)\n";
     
     ndivtime.open("ndivtime.out");
     ndivtime << "! Number of divisions attributed to selected BCs: time.average.sd\n";
     //MSchips
     selfBC_ndivtime.open("selfBC_ndivtime.out");
     selfBC_ndivtime << "! Number of divisions attributed to selected BCs: time.average.sd\n";

     Selfmutation_time.open("Selfmutation_time.out");
     Selfmutation_time << "! Mutation probability attributed to SELF BCs: time . averag . sd\n"
                   << "! c2+3: version based on selected BCs since last writing.\n"
                   << "! c4+5: version based on all current CBs.\n";

     freq_self.open("freq_tfrCont_self.out");
     freq_self << "! Tfr:CC cont / #CC (self) - time : avg : sd\n";
     freq_nonself.open("freq_tfrCont_nonself.out");
     freq_nonself << "! Tfr:CC cont / #CC (nonself) - time : avg : sd\n";

     mutation_time.open("mutation_time.out");
     mutation_time << "! Mutation probability attributed to BCs: time . averag . sd\n"
		   << "! c2+3: version based on selected BCs since last writing.\n"
		   << "! c4+5: version based on all current CBs.\n";
     // add a file for a histogram of GC-BC mutation numbers
     ofstream mutation_histo("mutation_histo.out");
     mutation_histo 
       << "! time[days] : time[hr] : mutation # : # of GC-BC with this mutation #\n";
     mutation_histo.close();
     write_mutations(0);
     ofstream igoutput("xsumout_ig.out");
     igoutput << "! time[hr] : time[days] : total out:since last : IgM[total:sincelast] :"
	      << "IgG[] : IgE[] : IgA[]\n";
     igoutput.close();
     ofstream foxofile("foxoCC.out");
     foxofile << "! time[hr] : time[days] "
	      << ": [CC : FoxO 0<=0.25 : 0.25<0.5 : 0.5<0.75 : 0.75<1 : 1<=] "
	      << ": for [CC : CC-apo : unselected+contact : FDCselect+TCcontact]\n";
     foxofile.close();
     ofstream FDCfailure("FDCagfailure.out");
     FDCfailure << "! time[hr] : time[days] : N_EventsOfVoidAgLastHr : this/N_CC : "
		<< "NofAgVoidFDCsites : TotalFDCsites : FractionOfAgVoidFDCsites \n";
     FDCfailure.close();
   }

   if ((show_mode == islet) && (cellbeta::LOCAL_FILES == false)) {
      // Initialise output file:
      beta_a.open("beta_a.out");
      beta_b.open("beta_b.out");
      beta_i.open("beta_i.out");
      beta_r.open("beta_r.out");
      beta_n.open("beta_n.out");
      beta_g.open("beta_g.out");
   }

   int objno = tALL;
   if (objno <= 0) {
      objno = tCB + tCC + tOUT + tTC + tBETA;
   }
   if (photoactivation == 0) {
      trackdata.init(l.dx,
                     dt,
                     l.dim,
                     l.prodimvec,
                     objno,
                     par.Value.trackfrom,
                     par.Value.trackuntil,
                     par.Value.track_delta_t,
                     par.Value.v_resolution,
                     par.Value.delta_v,
                     par.Value.alpha_resolution,
                     par.Value.delta_alpha,
                     par.Value.s_resolution,
                     par.Value.delta_s);
   } else {
      // if photoactivation==1
      if (objno > 0) {
         cout << "WARNING: Photoactivation is active and does not allow "
              << "for independent cell tracking.\n"
              << "Cell tracking object numbers are ignored.\n\n";
      }
      objno = (int (par.Value.photoactivation_delta_x / par.Value.dx) + 1);
      objno *= (int (par.Value.photoactivation_delta_y / par.Value.dx) + 1);
      if (l.dim == 3) {
         objno *= (int (par.Value.photoactivation_delta_z / par.Value.dx) + 1);
      }
      ana << "Initialise TRACK with " << objno << " objects.\n";
      trackdata.init(l.dx,
                     dt,
                     l.dim,
                     l.prodimvec,
                     objno,
                     par.Value.trackfrom,
                     par.Value.trackuntil,
                     par.Value.track_delta_t,
                     par.Value.v_resolution,
                     par.Value.delta_v,
                     par.Value.alpha_resolution,
                     par.Value.delta_alpha,
                     par.Value.s_resolution,
                     par.Value.delta_s);
      if (par.Value.trackfrom < par.Value.photoactivation_t0) {
         cout << "WARNING: Tracking starts before time of photoactivation.\n";
      }
   }

   ana << "done.\n";
   // cerr<<"... end of cellman::constructor.\n";
}
cellman::~cellman() {
   /*
    * cout<<"In ~cellman() deliberate memory [ ";
    * delete[] velocity;
    * cout<<"velocity ";
    *
    * CB_list.setcontrol(0);
    * for (long j=0; j<CB_list.lang(); ++j) CB_list[j].deliberate_memory();
    * cout<<"fragments of CB_ ";
    * FDC_list.setcontrol(0);
    * for (long j=0; j<FDC_list.lang(); ++j) FDC_list[j].deliberate_memory();
    * cout<<"FDC_ ";
    * OUT_list.setcontrol(0);
    * for (long j=0; j<OUT_list.lang(); ++j) OUT_list[j].deliberate_memory();
    * cout<<"OUT_list ] ";
    * cout<<"done.\n";
    */
   /*
    * CB_list.~dynarray();
    * CC_list.~dynarray();
    * FDC_list.~dynarray();
    * OUT_list.~dynarray();
    */
}
void cellman::close_files() {
   cout << "Close files in cellman ... ";
   if (show_mode != islet) {
      dark_zone.close();
      light_zone.close();
      // for (short j=0; j<xlogs; j++) xsums[j].close();
      xsums[CB].close();
      xsums[CC].close();
      xsums[out].close();
      xvolume.close();
      xsumBCig.close();
      xsumBCOig.close();
      xapoig.close();
      xapoigdz.close();
      xapoiglz.close();
      xaffig.close();
      cum_vol.close();
      log_out_aff.close();
      log_bc_aff.close();
      ag_presentation.close();
      cum_ag.close();
      xtc.close();
      xdec205.close();
      aff_cb_lz.close();
      aff_cc_lz.close();
      aff_out_lz.close();
      aff_cb_dz.close();
      aff_cc_dz.close();
      aff_out_dz.close();

      if (outputfiles == 0) {
         movements.close();
         move_vdt.close();
         axis.close();
         polarities.close();
         corrs.close();
      }
      mutation_out.close();
      cbhighaff.close();
      ag_collected.close();
      ndivtime.close();
      mutation_time.close();
      //MSchips
      selfBC_ndivtime.close();
      Selfmutation_time.close();
      freq_nonself.close();
      freq_self.close();
   }
   if (outputfiles == 0) {
      prolog.close();
      // cbquality.close();
   }
   if ((show_mode == islet) && (cellbeta::LOCAL_FILES == false)) {
      beta_a.close();
      beta_b.close();
      beta_i.close();
      beta_r.close();
      beta_n.close();
      beta_g.close();
   }
   cout << "done.\n";
}
void cellman::determine_founder_times_of_division(Werte &p) {
   // determines the values of the static variables founder_ndiv and p_founder_ndiv_plus
   // By default the number of division of founder cells is not fixed but determined
   // by the duration of the expansion phase:
   founder_ndiv = 0;
   p_founder_ndiv_plus = 0;
   // If a fixed number of divisions of founder cells is wanted there are different options:
   if (cellCB::fixed_number_of_divisions()) {
      if (p.fixed_time_of_divisions_mode == 0) {
         // calculate the number of divisions from the duration of the monoclonal expansion phase
         p_founder_ndiv_plus = (p.Start_Differentiation - p.tmin)
                               / cellCB::total_cell_cycle_duration();
      } else if (p.fixed_time_of_divisions_mode == 1) {
         // use the same number as after monoclonal expansion
         p_founder_ndiv_plus = p.CB_fixed_times_of_divisions;
      } else if (p.fixed_time_of_divisions_mode == 2) {
         // +++++++++++++ OPTION +++++++++++++++++++++++++++++++++++
         // define whatever you want
	/* The following was used in the Figure S1 of Cell Reports 2012.
	 * Refer to hyphasma12.03.2 for the original simulation.
	 * Parameter file is brdu-ag006.par.
	 */
         p_founder_ndiv_plus=drandom(2.0*p.CB_fixed_times_of_divisions);
	 //
         // p_founder_ndiv_plus=drandom(2.0*(p.CB_fixed_times_of_divisions-1.))+1.0;
         // p_founder_ndiv_plus = 8.0;
         // +++++++++ end OPTION +++++++++++++++++++++++++++++++++++
      } else if (p.fixed_time_of_divisions_mode == 3) {
         // use the value provided in the parameter file
         p_founder_ndiv_plus = p.CB_fixed_times_of_divisions_in_expansion;
      }
      founder_ndiv = int (p_founder_ndiv_plus);
      p_founder_ndiv_plus -= (double (founder_ndiv));    // prob for rest
   }
}
void cellman::add_extrafollicular_PC(AffinitySpace& shape, long pos_ss) {
  /* This routine adds an output cell at AS-position pos_ss.
   * The cell directly enters the soutext compartment, and will
   * subsequently differentiate to soutextproduce cells.
   * No adaptation of effective threshold variables like cellOUT::average_affinity is done.
   */
  shape.add_cell(soutext, pos_ss);
  /* "soutext" counted all output cells that left the GC so far. Now cells in "soutext"
   * have to be interpreted as the sum of GC output and extrafollicular PC.
   * As "sout" counts all GC-generated output cells, this is not added to "sout".
   * "total" counts all living cells in or from the GC and is not incremented.
   * This might have implications, when differences of "total" and "soutext" are used.
   * This has to be double checked. Could not identify any mistake (MMH 2018-09-20).
   */
}  
int cellman::get_founder_times_of_division() {
   // This routine uses predetermined variables and returns a number of divisions to do
   // ### This might become a property of the cells (cellthis.*) because there
   //     an analogous thing is done for cells differentiating from CC to CB in the
   //     equality operator.
   if (drandom() < p_founder_ndiv_plus) {
      return founder_ndiv + 1;
   }
   return founder_ndiv;
}
void cellman::get_average_and_sd(int * values, int &n, double &average, double &sd) {
   average = 0.;
   sd = 0.;
   if (n > 1) {
      for (int i = 0; i < n; i++) {
         average += double (values[i]);
      }
      average /= double (n);
      for (int i = 0; i < n; i++) {
         sd += (values[i] - average) * (values[i] - average);
      }
      sd /= double (n - 1);
      sd = sqrt(sd);
   } else if (n == 1) {
      average = values[0];
   }
}
void cellman::make_photoactivation(space &l) {
   cout << "Photoactivation of B cells ...\n";
   long noCB = CB_list.benutzt();
   long noCC = CC_list.benutzt();
   long noOUT = OUT_list.benutzt();
   long a = 0;
   long nCB = 0, nCC = 0, nOUT = 0;
   double r[l.dim];
   if (def_DEC205) {
      cout << " ... activate DEC205 positive B cells ...\n";
      /*
       * for (long i=0; i<noCB; ++i) {
       * if (CB_list[i].DEC205) {
       *  l.get_koord(CB_list[i].index,r);
       *  CB_list[i].trackit=true;
       *  CB_list[i].trackno=a;
       *  double l2saxis=CB_list[i].get_long2short_axis(i,CB,l);
       *  double elongat=CB_list[i].get_elongation(i,CB,l);
       *  trackdata.Write_movement(a,CB,time,r,CB_list[i].polarity,elongat,l2saxis,trackini);
       * ++nCB;
       * ++a;
       * }
       * }
       * cout<<"Photoactivation of "<<nCB<<" centroblasts.\n";
       */
      //
      for (long i = 0; i < noCC; ++i) {
         if (CC_list[i].DEC205) {
            if (CC_list[i].state != apoptosis) {
               l.get_koord(CC_list[i].index, r);
               CC_list[i].trackit = true;
               CC_list[i].trackno = a;
               trackdata.Write_movement(a, CC, time, r, CC_list[i].polarity, trackini);
               ++nCC;
               ++a;
            } else {
               cout << "*";
            }
         }
      }
      cout << "Photoactivation of " << nCC << " centrocytes.\n";
   } else {
      // photoactivate in a space regime
      cout << " ... activate B cells in a specified space regime ...\n";
      // go through all centroblasts
      for (long i = 0; i < noCB; ++i) {
         // get the coordinates of the cell i
         l.get_koord(CB_list[i].index, r);
         /*cout<<"r=("<<r[0]<<","<<r[1]<<"); rmin=("
          * <<photoactivation_rmin[0]<<","<<photoactivation_rmin[1]<<"); rmax=("
          * <<photoactivation_rmax[0]<<","<<photoactivation_rmax[1]<<");\n";*/
         // check whether the cell is in the photoactivation area
         bool itsin = true;
         for (short j = 0; j < l.dim; ++j) {
            if ((r[j] < photoactivation_rmin[j]) || (r[j] > photoactivation_rmax[j])) {
               itsin = false;
            }
         }
         // if yes mark this cell for tracking
         if (itsin) {
            CB_list[i].trackit = true;
            CB_list[i].trackno = a;
            double l2saxis = CB_list[i].get_long2short_axis(i, CB, l);
            double elongat = CB_list[i].get_elongation(i, CB, l);
            trackdata.Write_movement(a,
                                     CB,
                                     time,
                                     r,
                                     CB_list[i].polarity,
                                     elongat,
                                     l2saxis,
                                     trackini);
            ++nCB;
            ++a;
         }
      }
      cout << "Photoactivation of " << nCB << " centroblasts.\n";
      //
      for (long i = 0; i < noCC; ++i) {
         // get the coordinates of the cell i
         l.get_koord(CC_list[i].index, r);
         // check whether the cell is in the photoactivation area
         bool itsin = true;
         for (short j = 0; j < l.dim; ++j) {
            if ((r[j] < photoactivation_rmin[j]) || (r[j] > photoactivation_rmax[j])) {
               itsin = false;
            }
         }
         // if yes mark this cell for tracking
         // +++++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++
         // here one can chose to include apoptotic CC in photoactivation or not
         if (itsin) {
            // && CC_list[i].state!=apoptosis) { // || ignore_apoptotic_CC==false)) {
            // if the back part of the line above is active this option
            // is set by the parameter "Ignore apoptotic ..."
            // which then impacts on cell counts AND photoactivation !
            // +++++++++++++++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++
            CC_list[i].trackit = true;
            //	cout<<CC_list[i].state<<", ";
            CC_list[i].trackno = a;
            trackdata.Write_movement(a, CC, time, r, CC_list[i].polarity, trackini);
            ++nCC;
            ++a;
         } else if (itsin) {
            cout << "*";
         }
      }
      cout << "Photoactivation of " << nCC << " centrocytes.\n";
      //
      for (long i = 0; i < noOUT; ++i) {
         // get the coordinates of the cell i
         l.get_koord(OUT_list[i].index, r);
         // check whether the cell is in the photoactivation area
         bool itsin = true;
         for (short j = 0; j < l.dim; ++j) {
            if ((r[j] < photoactivation_rmin[j]) || (r[j] > photoactivation_rmax[j])) {
               itsin = false;
            }
         }
         // if yes mark this cell for tracking
         if (itsin) {
            OUT_list[i].trackit = true;
            OUT_list[i].trackno = a;
            double l2saxis = OUT_list[i].get_long2short_axis(i, out, l);
            double elongat = OUT_list[i].get_elongation(i, out, l);
            // cout<<"OUT-write-movement: elongat="<<elongat<<"\n";
            trackdata.Write_movement(a,
                                     out,
                                     time,
                                     r,
                                     OUT_list[i].polarity,
                                     elongat,
                                     l2saxis,
                                     trackini);
            ++nOUT;
            ++a;
         }
      }
      cout << "Photoactivation of " << nOUT << " output cells.\n";
   }   //
       // reset the number of tracked objects correspondingly:
   trackdata.set_N_OBJECTS(a);
   cout << "Photoactivation of " << a << " B cells in total.\n";
}
void cellman::make_fluorescent(space &l) {
   // Attributes trackit and trackno to the cells that are to be tracked
   // OPTION: The way to chose the cells can be defined differently

   // Start this only, if the timesteps are chosen accordingly.
   // Otherwise send out an error message!
   if (dt > TRACK::DELTA_T) {
      cout << "\nWARNING!\n"
           << "  dt_simulation=" << dt << " is larger than dt_tracking=" << TRACK::DELTA_T
           << " hours \n"
           << "  No meaningful tracking possible and prone to errors in class TRACK.\n"
           << "  Abort tracking.\n\n";
   } else {
      long noCB = CB_list.benutzt();
      long noCC = CC_list.benutzt();
      long noOUT = OUT_list.benutzt();
      long noTC = TC_list.benutzt();
      long noBETA = BETA_list.benutzt();
      // FDC are ignored
      int ctCB = 0, ctCC = 0, ctOUT = 0, ctTC = 0, ctBETA = 0;
      int a;

      if (tALL > 0) {
         // Version 1: Take the sum of all BC and TC and BETAcells,
         //            chose a random number up to the total cell number,
         //            attribute to the cell (cell-type is defined by the number).
         long noCELL = noCB + noCC + noOUT + noTC + noBETA;
         int b;
         a = 0;
         while (a < tALL) {
            b = irandom(noCELL);
            if (b < noCB) {
               // Find the CB
               if (CB_list[b].trackit == false) {
                  CB_list[b].trackit = true;
                  CB_list[b].trackno = a;
                  double r[l.dim];
                  l.get_koord(CB_list[b].index, r);
                  double l2saxis = CB_list[b].get_long2short_axis(b, CB, l);
                  double elongat = CB_list[b].get_elongation(b, CB, l);
                  trackdata.Write_movement(a,
                                           CB,
                                           time,
                                           r,
                                           CB_list[b].polarity,
                                           elongat,
                                           l2saxis,
                                           trackini);
                  ++ctCB;
                  ++a;
               }
            } else {
               b -= noCB;
               if (b < noCC) {
                  // Find the CC
                  if (CC_list[b].trackit == false) {
                     CC_list[b].trackit = true;
                     CC_list[b].fdc_clock = 0.;
                     CC_list[b].trackno = a;
                     double r[l.dim];
                     l.get_koord(CC_list[b].index, r);
                     trackdata.Write_movement(a, CC, time, r, CC_list[b].polarity, trackini);
                     ++ctCC;
                     ++a;
                  }
               } else {
                  b -= noCC;
                  if (b < noOUT) {
                     // Find the output cell
                     if (OUT_list[b].trackit == false) {
                        OUT_list[b].trackit = true;
                        OUT_list[b].trackno = a;
                        double r[l.dim];
                        l.get_koord(OUT_list[b].index, r);
                        trackdata.Write_movement(a,
                                                 out,
                                                 time,
                                                 r,
                                                 OUT_list[b].polarity,
                                                 trackini);
                        ++ctOUT;
                        ++a;
                     }
                  } else {
                     b -= noOUT;
                     if (b < noBETA) {
                        // Find the beta-cells
                        if (BETA_list[b].trackit == false) {
                           BETA_list[b].trackit = true;
                           BETA_list[b].trackno = a;
                           double r[l.dim];
                           l.get_koord(BETA_list[b].index, r);
                           trackdata.Write_movement(a,
                                                    BETA,
                                                    time,
                                                    r,
                                                    BETA_list[b].polarity,
                                                    trackini);
                           ++ctBETA;
                           ++a;
                        }
                     } else {
                        b -= noBETA;
                        // Find the TC
                        if (TC_list[b].trackit == false) {
                           TC_list[b].trackit = true;
                           TC_list[b].trackno = a;
                           double r[l.dim];
                           l.get_koord(TC_list[b].index, r);
                           trackdata.Write_movement(a,
                                                    TC,
                                                    time,
                                                    r,
                                                    TC_list[b].polarity,
                                                    trackini);
                           ++ctTC;
                           ++a;
                        }
                     }
                  }
               }
            }
         }
      } else {
         // if tALL<=0
         // Version 2: Chose the cells within a specified cell-type only
         ctCB = tCB;
         ctCC = tCC;
         ctOUT = tOUT;
         ctTC = tTC;
         ctBETA = tBETA;
         int b;
         a = 0;
         if (ctCB > noCB) {
            cout << "\nWARNING!!!\n"
                 << "  Not sufficient CB for tracking. Track all CB around.\n"
                 << "end of WARNING.\n";
            ctCB = noCB;
         }
         while (a < ctCB) {
            b = irandom(noCB);
            if (CB_list[b].trackit == false) {
               CB_list[b].trackit = true;
               CB_list[b].trackno = a;
               double r[l.dim];
               l.get_koord(CB_list[b].index, r);
               double l2saxis = CB_list[b].get_long2short_axis(b, CB, l);
               double elongat = CB_list[b].get_elongation(b, CB, l);
               trackdata.Write_movement(a,
                                        CB,
                                        time,
                                        r,
                                        CB_list[b].polarity,
                                        elongat,
                                        l2saxis,
                                        trackini);
               ++a;
            }
         }
         if (ctCC > noCC) {
            cout << "\nWARNING!!!\n"
                 << "  Not sufficient CC for tracking. Track all CC around.\n"
                 << "end of WARNING.\n";
            ctCC = noCC;
         }
         while (a < ctCB + ctCC) {
            b = irandom(noCC);
            // +++++++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++
            // Here, a subset of CC might be chosen for tracking
            if (CC_list[b].trackit == false) {
               // if (CC_list[b].trackit==false && (CC_list[b].state==unselected ||
               // CC_list[b].state==contact)) {
               // +++++++++++++++++++++++++ end OPTION ++++++++++++++++++++++++++++++++
               CC_list[b].trackit = true;
               CC_list[b].fdc_clock = 0.;
               CC_list[b].trackno = a;
               double r[l.dim];
               l.get_koord(CC_list[b].index, r);
               trackdata.Write_movement(a, CC, time, r, CC_list[b].polarity, 1, 1, trackini);
               ++a;
            }
         }
         if (ctOUT > noOUT) {
            cout << "\nWARNING!!!\n"
                 << "  Not sufficient OUT for tracking. Track all OUT around.\n"
                 << "end of WARNING.\n";
            ctOUT = noOUT;
         }
         while (a < ctCB + ctCC + ctOUT) {
            b = irandom(noOUT);
            if (OUT_list[b].trackit == false) {
               OUT_list[b].trackit = true;
               OUT_list[b].trackno = a;
               double r[l.dim];
               l.get_koord(OUT_list[b].index, r);
               trackdata.Write_movement(a, out, time, r, OUT_list[b].polarity, 1, 1, trackini);
               ++a;
            }
         }
         if (ctTC > noTC) {
            cout << "\nWARNING!!!\n"
                 << "  Not sufficient TC for tracking. Track all TC around.\n"
                 << "end of WARNING.\n";
            ctTC = noTC;
         }
         while (a < ctCB + ctCC + ctOUT + ctTC) {
            b = irandom(noTC);
            if (TC_list[b].trackit == false) {
               TC_list[b].trackit = true;
               TC_list[b].trackno = a;
               double r[l.dim];
               l.get_koord(TC_list[b].index, r);
               trackdata.Write_movement(a, TC, time, r, TC_list[b].polarity, 1, 1, trackini);
               ++a;
            }
         }
         if (ctBETA > noBETA) {
            cout << "\nWARNING!!!\n"
                 << "  Not sufficient betacells for tracking. Track all betacells around.\n"
                 << "end of WARNING.\n";
            ctBETA = noBETA;
         }
         while (a < ctBETA) {
            b = irandom(noBETA);
            if (BETA_list[b].trackit == false) {
               BETA_list[b].trackit = true;
               BETA_list[b].trackno = a;
               double r[l.dim];
               l.get_koord(BETA_list[b].index, r);
               double l2saxis = BETA_list[b].get_long2short_axis(b, CB, l);
               double elongat = BETA_list[b].get_elongation(b, CB, l);
               trackdata.Write_movement(a,
                                        BETA,
                                        time,
                                        r,
                                        BETA_list[b].polarity,
                                        elongat,
                                        l2saxis,
                                        trackini);
               ++a;
            }
         }
      }

      double rtmp[l.dim];
      for (int i = 0; i < l.dim; i++) {
         rtmp[i] = 0.;
      }
      int u = tALL;
      if (u <= 0) {
         u = tCB + tCC + tOUT + tTC + tBETA;
      }
      for (int i = a; i < u; i++) {
         trackdata.Write_movement(i, nocell, TRACK::TRACKFROM, rtmp, rtmp, trackini);
      }

      // Declare what is tracked on the screen and in the analysis-file
      cout << "\nAt time " << time - dt << " hours: Start of cell tracking ...\n";
      cout << "   Tracking " << ctCB << " CB, " << ctCC << " CC, " << ctOUT << " OUT, " << ctTC
           << " TC, " << ctBETA
           << " betacells.\n";
      if (tALL > 0) {
         cout << "   These numbers were randomly chosen.\n";
      } else {
         cout << "   These numbers were fixed in the parameter file.\n";
      }
   }   // end else (no timestep error)
}
void cellman::stop_fluorescent() {
   long n;
   for (n = 0; n < CB_list.benutzt(); n++) {
      CB_list[n].trackit = false;
      CB_list[n].trackno = -1;
   }
   for (n = 0; n < CC_list.benutzt(); n++) {
      CC_list[n].trackit = false;
      CC_list[n].fdc_clock = 0.;
      CC_list[n].trackno = -1;
   }
   for (n = 0; n < OUT_list.benutzt(); n++) {
      OUT_list[n].trackit = false;
      OUT_list[n].trackno = -1;
   }
   for (n = 0; n < TC_list.benutzt(); n++) {
      TC_list[n].trackit = false;
      TC_list[n].trackno = -1;
   }
   for (n = 0; n < BETA_list.benutzt(); n++) {
      BETA_list[n].trackit = false;
      BETA_list[n].trackno = -1;
   }

   cout << "   ... stop tracking at time " << time << " hours.\n";
}
void cellman::recombine_cells(double t) {
  long noCB = CB_list.benutzt();
  long noCC = CC_list.benutzt();
  //double p_ideal = tmx.get_recombination_probability(t); // deprecated
  //==================================================================================
  //========================= OPTION =================================================
  // Work on CBs:
  long countdefsCB = 0, countdecsCB = 0;
  for (int n = 0; n < noCB; n++) { 
    if (CB_list[n].MHCdeficient) { ++countdefsCB; } 
    if (CB_list[n].DEC205) { ++countdecsCB; }
  }
  double prob = tmx.get_recombination_probability(t, countdefsCB, countdecsCB);
  long countCB = 0;
  for (int n = 0; n < noCB; n++) {
    if (CB_list[n].DEC205 && not(CB_list[n].MHCdeficient) && drandom() < prob) { 
      CB_list[n].MHCdeficient = true; 
      ++countCB;
    }
  }
  /*
  cout << "Calculated recombination probability is " << 100*p_ideal << "%\n";
  cout << "Did recombine " << countCB << " of " << countdecsCB-countdefsCB 
       << " not recombined CBs in a total of " << countdecsCB << " Cre+ CBs,\n";
  cout << "This corresponds to ";
  if (countdecs > 0) { cout << 100*double(countCB)/double(countdecsCB); } else { cout << "0"; }
  cout << "% of Cre+ CBs and\n";
  if (countdecsCB-countdefsCB > 0) 
    { cout << 100*double(countCB)/double(countdecsCB-countdefsCB); } else { cout << "0"; }
  cout << "% of not-recombined CBs,\n";
  cout << "which compares to a rescaled probability of " << 100*prob << "%\n";
  */
  // Do the same for CC:
  long countCC = 0;
  long countdefsCC = 0, countdecsCC = 0;
  for (int n = 0; n < noCC; n++) { 
    if (CC_list[n].MHCdeficient) { ++countdefsCC; } 
    if (CC_list[n].DEC205) { ++countdecsCC; }
  }
  prob = tmx.get_recombination_probability(t, countdefsCC, countdecsCC);
  for (int n = 0; n < noCC; n++) {
    if (CC_list[n].DEC205 && not(CC_list[n].MHCdeficient) && drandom() < prob) { 
      CC_list[n].MHCdeficient = true; 
      ++countCC;
    }
  }
  /*
  cout << "Did recombine " << countCC << " of " << countdecsCC-countdefsCC 
       << " not recombined CCs in a total of " << countdecsCC << " Cre+ CCs,\n";
  cout << "which compares to a rescaled probability of " << prob << "\n";
  cout << "This corresponds to ";
  if (countdecsCC > 0) { cout << 100*double(countCC)/double(countdecsCC); } else { cout << "0"; }
  cout << "% of Cre+ CCs and\n";
  if (countdecsCC-countdefsCC > 0) 
    { cout << 100*double(countCC)/double(countdecsCC-countdefsCC); } else { cout << "0"; }
  cout << "% of not-recombined CCs,\n";
  cout << "which compares to a rescaled probability of " << 100*prob << "%\n";
  */
  long countCre = countdecsCB+countdecsCC;
  cout << "Recombined " << countCB+countCC << " of " 
       << countCre << " Cre+ B cells; ";
  double recombpercent = 0.;
  if (countCre > 0) 
    { recombpercent = double(countdefsCB+countdefsCC+countCB+countCC) / double(countCre); }
  cout << 100.*recombpercent  << "% Cre+ B cells recombined.\n";
  /* #### DEC205 was used as synonym for Cre+. This might be replaced later
   * by an extra cell variable. Not needed for now (8/2018). */
  //===================== end OPTION =================================================
  //==================================================================================
}
void cellman::attribute_DEC205() {
   cout << "Attribute DEC205 to B cells ...\n";
   long noCB = CB_list.benutzt();
   long noCC = CC_list.benutzt();
   long nCB = 0, nCC = 0;
   // go through all centroblasts
   for (long i = 0; i < noCB; ++i) {
      // attribute DEC205 positive with probability p_DEC205
      if (drandom(1.) < p_DEC205) {
         CB_list[i].DEC205 = true;
         ++nCB;
      }
   }
   cout << "Attributed DEC205 to " << nCB << " of " << noCB << " centroblasts.\n";
   //
   // go through all centrocytes
   for (long i = 0; i < noCC; ++i) {
      // attribute DEC205 positive with probability p_DEC205
      if (drandom(1.) < p_DEC205) {
         CC_list[i].DEC205 = true;
         ++nCC;
      }
   }
   cout << "Attributed DEC205 to " << nCC << " of " << noCC << " centrocytes.\n"
        << "Attributed DEC205 to " << nCB + nCC << " of " << noCB + noCC 
	<< " B cells in total.\n";
}
void cellman::label_BrdU(bool stainall) {
  /* Attributes BrdU to all CBs in S-phase now.
     Attributes BrdU to all CCs after Tfh selection in virtual S-phase now.
  */
  tnow_inject_BrdU = tnext_inject_BrdU;
  double BrdUmean = 100.;
  double BrdUwidth = 0.1 * BrdUmean;
  for (long i = 0; i < CB_list.benutzt(); ++i) {
    if (stainall || CB_list[i].state == cb_S) {
      CB_list[i].BrdU = get_positive_sample_from_normal(BrdUmean, BrdUwidth);
    }
  }
  /* stainall == true is set for in vitro experiments with CBs only, thus, the following
   * is only relevant for tfirst_injection_BrdU > 0 (which implies stainall == false): 
   */
  if (not(stainall)) {
    for (long i = 0; i < CC_list.benutzt(); ++i) {
      if (CC_list[i].state == selected) {
	// get virtual cell cycle phase of CCs:
	centroblasts phase = cellCB::get_virtual_cell_cycle_phase(CC_list[i].selected_clock);
	if (phase == cb_S) {
	  CC_list[i].BrdU = get_positive_sample_from_normal(BrdUmean, BrdUwidth);
	}
      }
    }
    /* If BrdU is repeatedly injected shift t_inject_BrdU to the next injection time point.
     * This is only done if the number of injections n_inject_BrdU was not reached.
     */
    if (deltat_inject_BrdU > 0 && 
	int(0.5 + (tnext_inject_BrdU - tfirst_inject_BrdU) / deltat_inject_BrdU) 
	< n_inject_BrdU - 1 ) 
      { tnext_inject_BrdU += deltat_inject_BrdU; }
  }
}
void cellman::do_inject_antiDEC205OVA(space &l, AffinitySpace &shape, ofstream &ana) {
  cout << "Inject anti-DEC205-OVA and activate DEC+/+ B cells ...\n";
  long noCB = CB_list.benutzt();
  long noCC = CC_list.benutzt();
  long noTC = TC_list.benutzt();
  long nCB = 0, nCC = 0;
  
  // go through all centroblasts
  // +++++DEC+++++
  if (inject_antiDEC205OVA_t0 == DEC205_ova_activity) {
    // then there is no time window for activation of CB so that they are activated here:
    for (long i = 0; i < noCB; ++i) {
      // activate all DEC205 positive cells
      if (CB_list[i].DEC205) {
	CB_list[i].DEC205_ova = true;
	++nCB;
      }
    }
  }
  // +++++DEC+++++
  
  // go through all centrocytes
  for (long i = 0; i < noCC; ++i) {
    // activate all DEC205 positive cells
    // possible states are
    // enum centrocytes {unselected,contact,FDCselected,TCcontact,selected,apoptosis};
    // only set DEC205_ova if neither selected nor apoptosis
    if ( CC_list[i].DEC205 && 
	 ((CC_list[i].state == unselected) || 
	  (CC_list[i].state == contact) || 
	  (CC_list[i].state == FDCselected) || 
	  (CC_list[i].state == TCcontact)) ) {
      CC_list[i].DEC205_ova = true;
      // reset state as now TC selection is FDC-independent:
      if (CC_list[i].state == contact) {
	shape.rem_cell(sCCcontact, CC_list[i].pos_ss);
	shape.add_cell(sCCunselected, CC_list[i].pos_ss);
	CC_list[i].state = unselected;
      }
      if (CC_list[i].state == FDCselected || CC_list[i].state == TCcontact) {
	// Adapt the attributed tc_search_duration
	CC_list[i].attribute_tc_search_duration();
	/* I cannot remember a reason to reset the aDEC205OVA targeted CCs to
	 * the state <unselected>. This has secondary effects since v20180120,
	 * where (tc_signal > FDCselected_clock) can happen, leading to
	 * a (TfhSignalSpeed > 1).
	 * Can this be savely removed? MMH 21.01.2018
	 * I prefer to do so now, because I also see other problems associated
	 * with simply resetting the state. There are other variables which
	 * have to be adapted, like counters etc. Might well induce artefacts.
	 */
	//shape.rem_cell(sCCFDCselected, CC_list[i].pos_ss);
	//shape.add_cell(sCCunselected, CC_list[i].pos_ss);
	//CC_list[i].state = unselected;
	/* Note that if the following lines are activated also CC in contact
	 * to TC would be reset to the unselected state. This might be wished,
	 * however, if activated, the CC have to be actively released from
	 * the TC using cellTC::liberateCC(...). Otherwise this leads to
	 * miscounting of the number of connected CC in the TCs. */
	// shape.rem_cell(sCCTCcontact,CC_list[i].pos_ss);
	// shape.add_cell(sCCunselected,CC_list[i].pos_ss);
	// CC_list[i].state=unselected;
      }
      ++nCC;
    }
  }
  // set changed_nn to 1 in order to force recalculation of all TC-polarities
  for (long i = 0; i < noTC; ++i) {
    TC_list[i].set_changed_nn(1);
  }
  cout << "Activated " << nCB + nCC << " DEC205+/+ B cells.\n";
  
  // Now change number of TC if required:
  if (TC_factor_dec205ova != 1.0) {
    long ntc = TC_list.benutzt();    // number of TC now
    long target_ntc = long (double (ntc) * TC_factor_dec205ova + 0.5);
    if (TC_factor_dec205ova < 1.0) {
      // reduce number of TC
      // remove (ntc-target_ntc) TC:
      cout << "\n\nERROR! A reduction of TC upon DEC205OVA stimulation "
	"is not yet programmed.\n\n\n";
      ana  << "\n\nERROR! A reduction of TC upon DEC205OVA stimulation "
	"is not yet programmed.\n\n\n";
      // ###
      // The problem with TC deletion is twofold:
      // 1) Deleted TC might be in contact to CC, thus, the state of CC has to be reset to
      // FDCselected.
      // 2) A changed TC listindex might induce trouble to the bound CC. But actually the space
      // index
      //    is attributed to bound CC so that this might be find. Check before programming it.
      /*
       * long todel;
       * short deleted;
       * while (ntc>target_ntc) {
       * todel=irandom(ntc);
       * deleted=del_TC(TC_list[todel].index,todel,l);
       * // in del_TC the CC_list has to be available and bound CC to be treated appropriately!
       *###
       * if (deleted==0) --ntc;
       * }
       */
    } else {
      // increase number of TC
      long add_nTC = target_ntc - ntc;
      // add add_nTC TC:
      addTC(add_nTC, l, 2.0, shape, ana);     // restrict to upper 1/"2" space part
      cout << "Added " << add_nTC << " TC in response to DEC205OVA stimulation.\n\n";
    }
  }
}
// =======================================================
// Fragment checks =======================================
// =======================================================

// ### diese Routinen sind nicht allgemein (Typ-unabhaengig)

void cellman::fragment_consistency(space &l) {
   /*
    * Checks if all fragments on the fragment lists of the cells
    * point to lattice point that point back (listi) to the same cell.
    */
   cout << "Run fragment consistency check of cell-list and lattice ... ";
   for (int n = 0; n < CB_list.benutzt(); n++) {
      for (int nn = 0; nn < CB_list[n].volume; nn++) {
         if (l.cellknot[CB_list[n].fragments[nn]].listi != n) {
            cout << "\nCB_list index=" << n << " while l.cellknot[CB_list[" << n
                 << "].fragments[" << nn
                 << "]].listi=" << l.cellknot[CB_list[n].fragments[nn]].listi << "\n";
            exit(1);
         }
      }
   }
   for (int n = 0; n < BETA_list.benutzt(); n++) {
      for (int nn = 0; nn < BETA_list[n].volume; nn++) {
         if (l.cellknot[BETA_list[n].fragments[nn]].listi != n) {
            cout << "\nBETA_list index=" << n << " while l.cellknot[BETA_list[" << n
                 << "].fragments[" << nn
                 << "]].listi=" << l.cellknot[BETA_list[n].fragments[nn]].listi << "\n";
            exit(1);
         }
      }
   }
   cout << "done.\n";
}
void cellman::check_listi(space &l) {
   // Checks if lattice points containing CBs BETAs
   // are contained in the fragment-lists of CB/BETA listi
   cout << "Check of l.knot[].listi-values for subelement cells ...";
   for (long f = 0; f < l.pointnumber; f++) {
      if ((f != -1) && (l.cellknot[f].listi != -1) && (l.cellknot[f].cell == CB)) {
         // check if f is contained in the fragments off listi
         // if (l.cellknot[f].listi==12) cout<<f<<" ";
         if (CB_list[l.cellknot[f].listi].where_fragment(f) == -1) {
            cout << "in proliferate: point " << f << " has listi = " << l.cellknot[f].listi
                 << " but is not saved as fragment.\n";
            exit(1);
         }
      }
      if ((f != -1) && (l.cellknot[f].listi != -1) && (l.cellknot[f].cell == BETA)) {
         // check if f is contained in the fragments off listi
         // if (l.cellknot[f].listi==12) cout<<f<<" ";
         if (BETA_list[l.cellknot[f].listi].where_fragment(f) == -1) {
            cout << "in proliferate: point " << f << " has listi = " << l.cellknot[f].listi
                 << " but is not saved as fragment.\n";
            exit(1);
         }
      }
   }
   cout << " done.\n";
}
/* Controls that the number of bound CC to any TC is staying
 *   below the maximum possible number of 2*space-dimension.
 *   In a more sophisticated version one might consider to
 *   check that the number of bound CC really fits the
 *   currently bound number of CC. */
/*
 * void cellman::check_TC_neighbours(space& l) {
 * for (int n=0; n<TC_list.benutzt(); n++)
 *  if (TC_list[n].n_CC_nn>2*l.dim) {
 *    cerr<<"ERROR!\n"
 *        <<"Number of CC bound to TC "<<TC_list[n].n_CC_nn
 *        <<" larger than "<<2*l.dim
 *        <<"for TC with list-index "<<n
 *        <<" and space index "<<TC_list[n].index<<".\n\n";
 *    exit(1);
 *  }
 * }
 */

// =======================================================
// All cells =============================================
// =======================================================
// === Those actions that concern lists of any type ======
// =======================================================

//MSc
long cellman::try2exchange_cells(long &i,
                                 long &li,
                                 states my_type,
                                 double * pol,
                                 dynarray<long> &redolist,
                                 space &l) {
   /* Exchanges cells if a movement was initiated but suppressed for lack of space
    * and if the cells are represented by single nodes. The exchange is only performed
    * if the polarities of both cells point to each other.
    * The routine can be called by any kind of cell because all necessary
    * properties are used as argument.
    */
   // cerr<<"ENTER try2exchange ... time="<<time<<"; i="<<i<<"; li="<<li<<"; my_type="<<my_type
   //  <<"; pol=("<<pol[0]<<","<<pol[1]<<")\n";
   // Check the requirements for an exchange (volume==1, type is in CB,CC,TC,out)
   // ... find the neighbour position in direction of polarity:
   long nn_index = l.get_nn_directed2(i, pol);
   // ... check cell-type of this neighbour
   states nn_type = l.cellknot[nn_index].cell;
   // cerr<<"Neighbour in direction pol has index = "<<nn_index
   //  <<" and type "<<nn_type
   //  <<" and position ("<<l.knot[i].x[0]<<","<<l.knot[i].x[1]<<").\n";
   // ##### What if this neighbour is empty, which might happen
   // because the blocking cell went away???
   // Then the cell could simply walk! #####end
   if ((nn_index != -1)
       && ((nn_type == CB) || (nn_type == CC) || (nn_type == TC) || (nn_type == TFR)
           || (nn_type == out) || (nn_type == BETA))) {
      // ### This routine would greatly benefit from a using a single list for all cell-types
      // #### Define a reference to a list !? and use it instead of having a list of ifs.
      // ... check volume of this neighbour and polarity
      bool volume_fine = false;
      bool state_fine = true;
      long nn_li = l.cellknot[nn_index].listi;
      double nn_pol[l.dim];
      if (nn_type == CB) {
         if (CB_list[nn_li].volume == 1) {
            volume_fine = true;
         }
         for (short j = 0; j < l.dim; j++) {
            nn_pol[j] = CB_list[nn_li].polarity[j];
         }
         // cerr<<" and volume="<<CB_list[nn_li].volume;
      } else if (nn_type == CC) {
         if (CC_list[nn_li].volume == 1) {
            volume_fine = true;
         }
         if ((CC_list[nn_li].state == contact) || (CC_list[nn_li].state == TCcontact)
                 || (CC_list[nn_li].state == TFRcontact) || (CC_list[nn_li].tfr_tmp_index!=-1) ) {
            state_fine = false;
         }
         for (short j = 0; j < l.dim; j++) {
            nn_pol[j] = CC_list[nn_li].polarity[j];
         }
         // cout<<" and volume="<<CC_list[nn_li].volume;
      } else if (nn_type == TC) {
         if (TC_list[nn_li].volume == 1) {
            volume_fine = true;
         }
         if (TC_list[nn_li].state == TC_CCcontact) {
            state_fine = false;
         }
         for (short j = 0; j < l.dim; j++) {
            nn_pol[j] = TC_list[nn_li].polarity[j];
         }
         // cout<<" and volume="<<TC_list[nn_li].volume;
      } else if (nn_type == TFR) {
          if (TFR_list[nn_li].volume == 1) {
             volume_fine = true;
          }
          if (TFR_list[nn_li].state == TFR_CCcontact || TFR_list[nn_li].tmp_blocked) {
             state_fine = false;
          }
          for (short j = 0; j < l.dim; j++) {
             nn_pol[j] = TFR_list[nn_li].polarity[j];
          }
      } else if (nn_type == out) {
         if (OUT_list[nn_li].volume == 1) {
            volume_fine = true;
         }
         if (OUT_list[nn_li].state == OutTFRcontact) {
             state_fine = false;
         }
         for (short j = 0; j < l.dim; j++) {
            nn_pol[j] = OUT_list[nn_li].polarity[j];
         }
         // cout<<" and volume="<<OUT_list[nn_li].volume;
      } else if (nn_type == BETA) {
         if (BETA_list[nn_li].volume == 1) {
            volume_fine = true;
         }
         for (short j = 0; j < l.dim; j++) {
            nn_pol[j] = BETA_list[nn_li].polarity[j];
         }
         // cout<<" and volume="<<CB_list[nn_li].volume;
      }

      // cout<<"  and list-index="<<nn_li<<"; pol=("<<nn_pol[0]<<","<<nn_pol[1]<<")\n";

      if (volume_fine && state_fine) {
         // check for the polarity now:
         // cout<<"Volume is fine.\n";

         // long nnn_index=l.get_nn_directed2(nn_index,nn_pol);
         // if (nnn_index==i) {
         // Alternatively one may ask whether the scalarproduct of both polarities is negative
         if (l.get_scalarproduct(pol, nn_pol) <= 0) {
            // cout<<"Polarities point to each other.\n";
            /* Do it with 0.5 probability only because the move is effectively two l.dx!
             * However, if both cells were contact inhibited, the exchange might be
             * performed anyway.
             */
            if ((drandom() < 0.5) || (redolist.find(nn_index) >= 0)) {
               /*cout<<"Do exchange ... \n";
                * cout<<"  i="<<i<<"; nn_index="<<nn_index<<"\n";
                * cout<<"  li="<<li<<"; nn_li="<<nn_li<<"\n"; */
               if ((i < 0) || (nn_index < 0) || (li < 0) || (nn_li < 0)) {
                  exit(1);
               }
               // perform exchange
               // move calling cell
               l.clear_knot(i);
               l.clear_knot(nn_index);

               l.set_knot(nn_index, my_type, li);
               if (my_type == CB) {
                  CB_list[li].index = nn_index;
                  CB_list[li].fragments[0] = nn_index;
                  // As the barycenter has changed it has to be actualised:
                  CB_list[li].get_barycenter(l);
                  // movement by one lattice constant is added
                  CB_list[li].add_moves += 1.;
                  // When fragmove was called n_immobile was set up by 1 because of no movement,
                  // now the movement was still performed in the same time step, thus, this has
                  // to be reversed.
                  --CB_list[li].n_immobile;
                  // Migration analysis is performed in the calling routine
                  // on the basis of the returned value.
               } else if (my_type == CC) {
                  CC_list[li].index = nn_index;
               } else if (my_type == TC) {
                  TC_list[li].index = nn_index;
               } else if (my_type == TFR) {
                   TFR_list[li].index = nn_index;
               } else if (my_type == out) {
                  OUT_list[li].index = nn_index;
               } else if (my_type == BETA) {
                  BETA_list[li].index = nn_index;
                  BETA_list[li].fragments[0] = nn_index;
                  // As the barycenter has changed it has to be actualised:
                  BETA_list[li].get_barycenter(l);
                  // movement by one lattice constant is added
                  BETA_list[li].add_moves += 1.;
                  // When fragmove was called n_immobile was set up by 1 because of no movement,
                  // now the movement was still performed in the same time step, thus, this has
                  // to be reversed.
                  --BETA_list[li].n_immobile;
                  // Migration analysis is performed in the calling routine
                  // on the basis of the returned value.
               }

               l.set_knot(i, nn_type, nn_li);
               if (nn_type == CB) {
                  CB_list[nn_li].index = i;
                  CB_list[nn_li].fragments[0] = i;
                  // As the barycenter has changed it has to be actualised:
                  CB_list[nn_li].get_barycenter(l);
                  // movement by one lattice constant is added
                  CB_list[nn_li].add_moves += 1.;
                  // n_immobile cannot be reversed because it is not known
                  // whether this cell was already moved in the same time step or not.
                  // #### ignored error for the moment (n_immobile is not important)
                  if (CB_list[nn_li].trackit) {
                     // Also write to the migration analysis file
                     // (elongation parameters need not to be calculated for volume==1)
                     trackdata.Write_movement(CB_list[nn_li].trackno,
                                              CB,
                                              time,
                                              CB_list[nn_li].barycenter,
                                              CB_list[nn_li].polarity,
                                              movement);
                  }
               } else if (nn_type == CC) {
                  CC_list[nn_li].index = i;
                  if (CC_list[nn_li].trackit) {
                     double tmpi[l.dim];
                     l.get_koord(CC_list[nn_li].index, tmpi);
                     trackdata.Write_movement(
                        CC_list[nn_li].trackno,
                        CC,
                        time,
                        tmpi,
                        CC_list[nn_li].polarity,
                        movement);
                  }
               } else if (nn_type == TC) {
                  TC_list[nn_li].index = i;
                  if (TC_list[nn_li].trackit) {
                     double tmpi[l.dim];
                     l.get_koord(TC_list[nn_li].index, tmpi);
                     trackdata.Write_movement(
                        TC_list[nn_li].trackno,
                        TC,
                        time,
                        tmpi,
                        TC_list[nn_li].polarity,
                        movement);
                  }
               } else if (nn_type == TFR) {
                   TFR_list[nn_li].index = i;
                   if (TFR_list[nn_li].trackit) {
                      double tmpi[l.dim];
                      l.get_koord(TFR_list[nn_li].index, tmpi);
                      trackdata.Write_movement(
                         TFR_list[nn_li].trackno,
                         TFR,
                         time,
                         tmpi,
                         TFR_list[nn_li].polarity,
                         movement);
                   }
               } else if (nn_type == out) {
                  OUT_list[nn_li].index = i;
                  if (OUT_list[nn_li].trackit) {
                     double tmpi[l.dim];
                     l.get_koord(OUT_list[nn_li].index, tmpi);
                     trackdata.Write_movement(
                        OUT_list[nn_li].trackno,
                        out,
                        time,
                        tmpi,
                        OUT_list[nn_li].polarity,
                        movement);
                  }
               } else if (nn_type == BETA) {
                  BETA_list[nn_li].index = i;
                  BETA_list[nn_li].fragments[0] = i;
                  // As the barycenter has changed it has to be actualised:
                  BETA_list[nn_li].get_barycenter(l);
                  // movement by one lattice constant is added
                  BETA_list[nn_li].add_moves += 1.;
                  // n_immobile cannot be reversed because it is not known
                  // whether this cell was already moved in the same time step or not.
                  // #### ignored error for the moment (n_immobile is not important)
                  if (BETA_list[nn_li].trackit) {
                     // Also write to the migration analysis file
                     // (elongation parameters need not to be calculated for volume==1)
                     trackdata.Write_movement(BETA_list[nn_li].trackno,
                                              BETA,
                                              time,
                                              BETA_list[nn_li].barycenter,
                                              BETA_list[nn_li].polarity,
                                              movement);
                  }
               }
               // cout<<"\n";
               return nn_index;
            }
         }
      }
   }
   // cout<<"End of try2exchange: did not!\n\n";
   return -1;
}
//MSc
long cellman::retry_movement(dynarray<long> &redolist, space &l) {
   // Any cell is only on the list if its volume is 1, thus, no need to check.
   // Get the properties of the cell at this position:
   long i = redolist[0];
   states my_type = l.cellknot[i].cell;
   long my_li = l.cellknot[i].listi;
   long partner = -1;   // position of the partner before movement
   // and call the exchange routine
   if (my_type == CB) {
      partner = try2exchange_cells(i, my_li, my_type, CB_list[my_li].polarity, redolist, l);
   } else if ((my_type == CC) && (CC_list[my_li].state != contact)
              && (CC_list[my_li].state != TCcontact) && (CC_list[my_li].state != TFRcontact)
              && (CC_list[my_li].tfr_tmp_index==-1)) {
      partner = try2exchange_cells(i, my_li, my_type, CC_list[my_li].polarity, redolist, l);
   } else if ((my_type == TC) && (TC_list[my_li].state != TC_CCcontact)) {
      partner = try2exchange_cells(i, my_li, my_type, TC_list[my_li].polarity, redolist, l);
   } else if ((my_type == TFR) && (TFR_list[my_li].state != TFR_CCcontact) && !TFR_list[my_li].tmp_blocked) {
       partner = try2exchange_cells(i, my_li, my_type, TFR_list[my_li].polarity, redolist, l);
   } else if (my_type == out && OUT_list[my_li].state!=OutTFRcontact) {
      partner = try2exchange_cells(i, my_li, my_type, OUT_list[my_li].polarity, redolist, l);
   } else if (my_type == BETA) {
      partner = try2exchange_cells(i, my_li, my_type, BETA_list[my_li].polarity, redolist, l);
   }

   // if no exchange was performed return -1:
   return partner;
}
// =======================================================
// Centroblasts ==========================================
// =======================================================
// === Only those actions that concern the CB_list =======
// === at a whole stay in cellman.x, more in cells.x =====
// =======================================================

long cellman::put_CB(long i, cellCB &newCB, space &l, AffinitySpace &shape) {
   /* This routine is to be understood as construction of a
    * first CB (one may think of naive CBs) which is smaller
    * compared to its asymptotic final state. put_CB cannot
    * be used to restore large CBs after their deletion.
    * This is done by the shift of single fragments following
    * Pott's ideas.
    * ==> no generalization to CB with large radii is necessary!
    * Returns the CB_list-index of the new cell.
    */
   if (i < 0) {
      cout << "Negative Index in put_CB!!!\n";
      exit(1);
   }
   // cout<<"newCB.switch="<<newCB.time_of_cycle_state_switch<<", ";
   if (l.cellknot[i].cell != nocell) {
      return -1;
   }
   // Add new cell on the list
   newCB.get_radius(l.dim);   // ### is repeated in ini() below -> remove at one place
   // cout<<"Add CB at position i="<<i<<"; ";
   long li = CB_list.add(newCB);
   // cout<<"list-index is li="<<li<<"; ";
   // cout<<"after add:newCB.switch="<<CB_list[li].time_of_cycle_state_switch<<", ";
   // Initialize it as new cell:
   CB_list[li].ini(i, li, time, l, shape);
   // cout<<"after ini:newCB.switch="<<CB_list[li].time_of_cycle_state_switch<<"; ... ; ";
   // cout<<"initialized the new cell!\n";
   // That's it.
   return li;
}
/*
 * long int cellman::find_CB(long i) {
 * long int found=-1;
 * long int n=0;
 * long int bis=CB_list.benutzt();
 * while (n<bis && found==-1) {
 *    if (CB_list[n].index==i) { found=n; }
 *    else { n++; }
 * }
 * if (found==-1) { cout<<"Error: Index "<<i<<" not found in find_CB!"; exit(1); }
 * return found;
 * }
 */

short int cellman::del_CB(long int i, long li, space &l, AffinitySpace &shape) {
   if (l.cellknot[i].cell != CB) {
      cout << "Error in del_CB: Zelltyp ist = " << l.cellknot[i].cell << " !!!\n";
      exit(1);
   }
   // Delete all fragments from the lattice
   for (int a = 0; a < CB_list[li].volume; a++) {
      l.clear_knot(CB_list[li].fragments[a]);
      // l.cellknot[CB_list[li].fragments[a]].cell=nocell;
      // l.cellknot[CB_list[li].fragments[a]].listi=-1;
   }
   shape.rem_cell(sCB, CB_list[li].pos_ss);
   shape.rem_cell(total, CB_list[li].pos_ss);
   // Delete cell on position li on the list
   CB_list.erase_jump(li);
   // Correct on the lattice for the shift of last cell on the list to position li
   // Don't if the deleted cell is the last cell on the list
   if (li < CB_list.benutzt()) {
      for (int a = 0; a < CB_list[li].volume; a++) {
         l.cellknot[CB_list[li].fragments[a]].listi = li;
      }
   }
   return 0;
}
short cellman::macrophagocyte_CB(long i, long li, space &l, AffinitySpace &shape) {
   if (drandom() < p_macrophage) {
      del_CB(i, li, l, shape);
      return 1;
   }
   return 0;
}
short int cellman::differ2CC_CB(long i, long li, space &l, AffinitySpace &shape) {
   // Uebernehme die aktuellen Daten from cell at position li
   cellCC ctmp;
   ctmp = CB_list[li];
   ctmp.responsive2signal[CXCL12] = false;
   ctmp.responsive2signal[CXCL13] = true;
   if (ctmp.IgX.Ig_class == IgE) {
      if (ctmp.IgE_prob_CXCR5down < 0.) {
         // upregulate CXCR4 instead of CXCR5
         if (drandom() < -1.0 * ctmp.IgE_prob_CXCR5down) {
            ctmp.CXCR5failure = 2;
            ctmp.responsive2signal[CXCL13] = false;
            ctmp.responsive2signal[CXCL12] = true;
         }
      } else if (ctmp.IgE_prob_CXCR5down > 0.) {
         // just downregulate CXCR5
         if (drandom() < ctmp.IgE_prob_CXCR5down) {
            ctmp.CXCR5failure = 1;
            ctmp.responsive2signal[CXCL13] = false;
         }
      }
   }
   // +++++DEC+++++
   if (inject_antiDEC205OVA) {
     if (ctmp.DEC205 && 
	 (time >= inject_antiDEC205OVA_t0) && 
	 (time < DEC205_ova_activity)) {
       ctmp.DEC205_ova = true;
     }
     // If this line above is active, CB should not be activated upon injection right away.
     // The thought is that CB get only dec205 activated when they differentiate to CC
     // Activity must then also be stopped here if the time window has passed:
     if (time > DEC205_ova_activity) {
       ctmp.DEC205_ova = false;
     }
     /* Thus, those who came from a TB interaction now get inactivated 
      * if not within the time window. */
   }
   // +++++DEC+++++
   // New CC is in state unselected
   ctmp.state = unselected;
   // If ag is not retained in CB, collection of ag has to be started from zero
   if (retain_ag == false) {
      /* In this case, CB did not transport the Ag collected in the previous round
       * and here, the Ag memory (nFDCcontacts and collected_ag_portions) is reset
       * to zero in the new CC. */
      ctmp.delete_antigen();
   } else {
      /* Here, CB retained the Ag collected in the LZ and eventually divided it
       * asymmetrically. Now there might still be some left. Depending on the option
       * to keep the remaining Ag to the next round of LZ selection, Ag memory
       * nFDCcontacts and collected_ag_portions is reset to zero or not. */
      ctmp.process_antigen();
   }
   // cout<<"nFDCcontacts in diff CB2CC="<<ctmp.nFDCcontacts<<"\n";
   // Delete CB
   if (l.cellknot[i].cell != CB) {
      cout << "Fehler in differ2CC_CB: no CB! ";
      exit(1);
   }
   del_CB(i, li, l, shape);
   // Erzeuge CC
   if (put_CC(i, ctmp, l, shape) != 0) {
      exit(1);
   }
   return 0;
}
void cellman::write_log_out_aff(double bestaff, long ASpos, AffinitySpace &shape) {
   log_out_aff << time << "  " << bestaff;
   int n_ags = shape.get_n_Antigen();
   for (int ag = 0; ag < n_ags; ag++) {
      double aff = shape.affinity_norm(ASpos, shape.get_Antigen(ag)); /// Philippe : in replacement
                                                                      // of get_affinity2ag() that
                                                                      // was just calling
                                                                      // getAffinity_norm
      log_out_aff << "  " << aff;
   }
   log_out_aff << "\n";
}
void cellman::write_log_bc_aff_this(short celltyp, long &ASpos, AffinitySpace &AS) {
   double aff = AS.best_affinity_norm(ASpos);
   log_bc_aff << time << "  " << celltyp << "  " << aff << "  ";
   int n_ags = AS.get_n_Antigen();
   for (int ag = 0; ag < n_ags; ag++) {
      aff = AS.affinity_norm(ASpos, AS.get_Antigen(ag)); /// Philippe : in replacement of
                                                         // get_affinity2ag(
      log_bc_aff << "  " << aff;
   }
   log_bc_aff << "\n";
}
void cellman::write_log_bc_aff(AffinitySpace &AS) {
   long ASpos;
   long ncell = CB_list.benutzt();
   for (int i = 0; i < ncell; i++) {
      ASpos = CB_list[i].pos_ss;
      write_log_bc_aff_this(0, ASpos, AS);
   }
   ncell = CC_list.benutzt();
   for (int i = 0; i < ncell; i++) {
      ASpos = CC_list[i].pos_ss;
      write_log_bc_aff_this(1, ASpos, AS);
   }
   ncell = OUT_list.benutzt();
   for (int i = 0; i < ncell; i++) {
      ASpos = OUT_list[i].pos_ss;
      write_log_bc_aff_this(2, ASpos, AS);
   }
}
short int cellman::CB_differ2OUT(long i, long li, space &l, AffinitySpace &shape) {

    cerr<<"LEDA NOT CODED FOR TFR\n";exit(1);

   // adjust global counters
   ++n_outs;
   if (CB_list[li].n_recycling > 0) {
      ++n_outs_with_recycling;
   }
   n_muts += CB_list[li].n_mutation;
   n_recmuts += CB_list[li].n_recandmute;
   if (CB_list[li].n_mutation < max_mutation_bin) {
      ++mutation_frequency[CB_list[li].n_mutation];
   } else {
      ++mutation_frequency[max_mutation_bin - 1];
   }

   // Save old cell states
   cellOUT newOUT;
   newOUT = CB_list[li];
   newOUT.index = i;
   newOUT.get_radius(l.dim);
   newOUT.responsive2signal[CXCL13] = false;
   newOUT.responsive2signal[CXCL12] = false;   // shall we have CXCL12 activity here ? ###
   for (short d = 0; d < l.dim; d++) {
      newOUT.barycenter[d] = CB_list[li].barycenter[d];
   }


   // Delete CB
   if (l.cellknot[i].cell != CB) {
      cout << "Fehler in differ2CC_CB: no CB! ";
      exit(1);
   }
   del_CB(i, li, l, shape);

   // Create new OUT-cell:
   l.set_knot(i, out, OUT_list.benutzt());
   // Write OUT-cell on cell_list und save list-index on the lattice
   OUT_list.add(newOUT);
   // correct shape space:
   shape.add_cell(sout, newOUT.pos_ss);
   // add to shape space if DEC205 positive
   if (newOUT.DEC205) {
      shape.add_cell(soutdec, newOUT.pos_ss);
   }
   //MSchips
   if (newOUT.selfMutation) {
       shape.add_cell(soutSelf,newOUT.pos_ss);
   } else {shape.add_cell(soutNonSelf, newOUT.pos_ss);}
   // save this output event in the Ig-class output counter
   ++integrate_out[newOUT.IgX.Ig_class];
   ++integrate_out[nIg_classes];

   double aff = shape.best_affinity_norm(newOUT.pos_ss);
   // correct average affinity of output cells
   shape.correct_average_affinity(sout, newOUT.pos_ss, cellOUT::average_affinity);
   if (aff > cellOUT::max_affinity) {
      cellOUT::max_affinity = aff;
   }

   // Write to log_out_aff
   write_log_out_aff(aff, newOUT.pos_ss, shape);

   newOUT.deliberate_memory();

   return 0;
}
short cellman::BCinflux(space &l, AffinitySpace &shape, sigs &s) {
   if (p_BCinflux == 0) {
      return 0;    // don't do anything, no error
   }
   short err = 1;
   // set the probability of BCinflux
   double p_in = p_BCinflux;
   if (do_smooth_stopBCinflux) {
      p_in = p_BCinflux / (1.0 + exp((time - t_stopBCinflux) / smooth_stopBCinflux));
   } else if (time > t_stopBCinflux) {
      p_in = 0;
   }
   // eventually add another new BC to the GC reaction:
   if (drandom() < p_in) {
      // generate a new BC in CB state
      cellCB newCB;
      // determine an available target place:
      long getpos = l.get_random_empty_position();
      if (getpos < 0) {
	err = 1;
	++pp[3];
	// write this failure to prolog.out
	prolog << time << "  " 
	       << pp[0] << "  " << pp[1] << "  " << pp[2] << "  " << pp[3] << "  "
	       << -1 << "\n"; // -1 allows to distinguish the entry from division.
      } else {
	// use random seeder cells within the list of seeder cells
	newCB.pos_ss = shape.get_Seeder();    // this uses a random choice
	/*
        if (shape.best_affinity_norm(newCB.pos_ss)>0.5) { 
	  newCB.DEC205 = true; 
	  cout << "WARNING! Code specific for Oliver Bannard to tag SWHEL cells\n"
	       << "         and know their fraction in the course of the reaction.\n"
	       << "         This option is set in cellman::BCinflux and MUST be\n"
	       << "         inactivated when DEC205-experiments are performed.\n"
	       << "WARNING END!\n\n\n";
	}
	*/
	// cout<<"ss-position="<<newCB.pos_ss<<"\n";
	newCB.make_CB_new();    // sets the state of the cell to CB and resets its internal clock
	// by default n_divisions2do=0
	if (cellCB::fixed_number_of_divisions()) {
	  // set the number of divisions
	  newCB.n_divisions2do = get_founder_times_of_division();
	}
	if (s.signal_use[CXCL12]) {
	  newCB.responsive2signal[CXCL12] = true;
	}
	if (def_DEC205 && def_DEC205_t0 < 0) {
	  newCB.attribute_DEC205(p_DEC205);
	}
	// The following was used to mimick low and high affinity founder cells
	// for the vanishing cells of Thomas Winkler:
	//if (newCB.DEC205) { newCB.pos_ss = shape.get_Seeder(0); }
	//else { newCB.pos_ss = shape.get_Seeder(1); }
	long newindex = put_CB(getpos, newCB, l, shape);
	// Adapt the average seeder cell affinity
	cellCB::average_seeder_affinity
	  = (n_founder * cellCB::average_seeder_affinity
	     + shape.best_affinity_norm(newCB.pos_ss)) / (n_founder + 1);
	// Set cellOUT-parameters if needed:
	if ((cellOUT::use_threshold > 0) && (cellOUT::initial_ab_affinity < 0)) {
	  /* These values of cellOUT::average_ and max_affinity are used to define thresholds
	   * of binding antigen (in cellman::calc_CC()).
	   * If cellOUT::initial_ab_affinity>=0 this value is used as initial threshold
	   * and was set in cellOUT::set_statics() and should not be changed here.
	   * The new seeder cell has to be included to define the new threshold,
	   * but only before start of output production: */
	  if (n_outs == 0) {
            cellOUT::average_affinity = cellCB::average_seeder_affinity;
            double aff = shape.best_affinity_norm(newCB.pos_ss);
            if (aff > cellOUT::max_affinity) {
	      cellOUT::max_affinity = aff;
            }
	  }
	}
	// cout<<"cellCB::average_seeder_affinity="<<cellCB::average_seeder_affinity<<"\n";
	// increase the cummulative number of founder cells
	++n_founder;
	// write the new cell as founder cell to the brainbow class
	if (track_mutations) {
	  CB_list[newindex].brainbow_index
            = trackmutations.write(time, true, true, -1, CB_list[newindex].pos_ss);
	}
	/* Eventually, preload the CB with antigen:
	 * Normally this would be used for the in vitro BC asymmetric division of antigen
	 * experiment by Batista. BC influx for in vitro doesn't really make sense.
	 */
	CB_list[newindex].preload_with_ag();
        ///This is the experiment in Benet et al. Fr.Imm2018
        /// Half of the cells are CCL3KO, i.e. they interact at lower freq with Tfr due to lack of CCL3 secretion
        /// Here, selfMutation are actually representing the WT
        /// Tfr can then inteeract only with SelfMutation-Bcells
        if (drandom()<CCL3KO) {
            CB_list[newindex].ccl3KO=1;
        }
	// in case of multi-fragment cells, temporary arrays may be removed
	newCB.deliberate_memory();
	err = 0;
	long ag_pos = shape.get_Antigen(shape.get_nearest_Antigen(CB_list[newindex].pos_ss));
	cout << "... new founder at " << CB_list[newindex].pos_ss << " with distance "
	     << sqrt(shape.Abstandquad(CB_list[newindex].pos_ss, ag_pos)) << "; " << n_founder
	     << " founder cells incorporated\n";
	// New founder cells are not monitored in the analysis file.
	// If track_mutations==true, they are monitored in the brainbow file.
      } // end if (getpos >=0)
   }   // end if (drandom()<p_in)
   return err;
}
short cellman::proliferate_CB(long i, long li, space &l, AffinitySpace &shape) {
   short err = 1;
   long newli;
   long j = CB_list[li].ask_mitosis(pp, l);
   // note :  in case of one fragments, j is the found position.
   //        in case of multi-fragment cell, then j = 0/1 depending on if can divide or not.
   if (cellCB::target_volume == 1) {
      if (j >= 0) {
         // if (li==50) cout<<"n_divisions2do="<<CB_list[li].n_divisions2do<<",
         // state="<<CB_list[li].state;
         // write in log-file for proliferation
         if (outputfiles == 0) {
            prolog << time << "  " << pp[0] << "  " << pp[1] << "  " << pp[2] << "  "
                   << pp[3] << "  "
                   << l.Abstand(i, j) << "\n";
         }
         err = 0;
         // new cell has the same state as the dividing one (apart of state)
         cellCB newCB = CB_list[li];
         // cout<<"in proliferate_CB-1: CB["<<li<<"].p_mutation="<<CB_list[li].p_mutation<<"; "
         // <<"newCB.p_mutation="<<newCB.p_mutation<<"\n";
         // cerr<<"BEFORE: newCB.state="
	 // <<newCB.state<<" n_divisions2do="<<newCB.n_divisions2do<<",";
         // cerr<<"CBlist.state="<<CB_list[li].state<<"
         // n_divisions2do="<<CB_list[li].n_divisions2do<<"; ";
         newCB.make_CB_new();     // newCB.state=cb_normal;
         CB_list[li].make_CB_new();
         CB_list[li].set_remaining_divisions();
         newCB.set_remaining_divisions();
         //      if (li==50) cout<<" ... --> old.n_div="<<CB_list[li].n_divisions2do
         //	      <<", new.n_div="<<newCB.n_divisions2do<<"\n";
         // cerr<<"AFTER: newCB.state="<<newCB.state<<" n_divisions2do="<<newCB.n_divisions2do<<",
         // ";
         // cerr<<"CBlist.state="<<CB_list[li].state<<"
         // n_divisions2do="<<CB_list[li].n_divisions2do<<"\n";
         // ##### The following might be shifted to make_CB_new()! also for fragmented cell below!
         // +++++DEC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         // The new daughter does not inherit the DEC205_ova state
         if (time > DEC205_ova_activity) {
            newCB.DEC205_ova = false;
         }
         // In addition one might ask to induce differentiation upon division events in both
         // daughters:
         if (inject_antiDEC205OVA && 
	     DEC205_induce_differentiation && 
	     (time >= inject_antiDEC205OVA_t0) && 
	     (time < DEC205_ova_activity)) {
	   if (newCB.DEC205) {
	     newCB.state = newCB.set2differentiation();
	   }
	   if (CB_list[li].DEC205) {
	     CB_list[li].state = CB_list[li].set2differentiation();
	   }
         }
         // +++++DEC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         // newCB.responsive2signal[CXCL13]=false;
         newCB.trackit = false;
         // ###### shouldn't the daughter inherit fluorescence??? (see below for DEC205+
         // photoactivation)
         newCB.trackno = -1;
         newli = put_CB(j, newCB, l, shape);
         newCB.deliberate_memory(); // ### is this necessary in the one node one cell case?
      }
   } else if (j == 0) {
      // case of more fragment object and do proliferate
      cellCB newCB = CB_list[li];
      // reset state, do not track, and preserve all other variables
      newCB.make_CB_new();    // newCB.state=cb_normal;
      CB_list[li].make_CB_new();
      CB_list[li].set_remaining_divisions();
      newCB.set_remaining_divisions();
      // +++++DEC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      if (inject_antiDEC205OVA) {
	// The new daughter does not inherit the DEC205_ova state
	if (time > DEC205_ova_activity) {
	  newCB.DEC205_ova = false;
	}
	// In addition one might ask to induce differentiation upon division in both daughters:
	if (DEC205_induce_differentiation && 
	    (time >= inject_antiDEC205OVA_t0) && 
	    (time < DEC205_ova_activity)) {
	  if (newCB.DEC205) {
	    newCB.state = newCB.set2differentiation();
	  }
	  if (CB_list[li].DEC205) {
	    CB_list[li].state = CB_list[li].set2differentiation();
	  }
	}
      }
      // +++++DEC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      // newCB.responsive2signal[CXCL13]=false;
      newCB.trackit = false;
      newCB.trackno = -1;
      // add a new element on the list
      newli = CB_list.add(newCB);
      // Hier kann man CB_list[newli].ini(...) nicht verwenden, da newCB=CB_list[li] wichtig!

      err = CB_list[li].mitosis(i, li, newli, CB_list[newli], l);
      if (err == 0) {
	// Actualize shape space statistics
	/// Philippe: 
	///I don't understand why there is not the same thing in the first case (volume = 1)
	shape.add_cell(sCB, CB_list[newli].pos_ss);
        shape.add_cell(total, CB_list[newli].pos_ss);
	// check_connection of the result
	if (checkit == 2) {
	  CB_list[li].check_connection(CB, li, l);
	  CB_list[newli].check_connection(CB, newli, l);
	}
	// end of preliminary check
	++pp[0];
	if (outputfiles == 0) {
	  prolog << time << "  " << pp[0] << "\n";
	}
      } else {
	exit(1);
      }
      newCB.deliberate_memory();
   }
   if (err == 0) {
      // if division successful new cells in li and newli:
      // Divide BrdU symmetrically on the daughters:
      CB_list[li].BrdU /= 2.;
      CB_list[newli].BrdU /= 2.;
      // Option how to set the mutaton probability in the daughter cells:
      CB_list[newli].p_mutation = CB_list[li].p_mutation;
      /* Note that p_mutation in [newli] is set to the value in [li].
       * This implies that both daughters will have the same mutation rate
       * (which is transferred with the =operater in the cell-class)
       * and that an increased or reduced mutation rate will be inherited
       * in the spirit of a symmetric division.
       * If the above line was deactivated, thus, p_mutation was not touched,
       * then p_mutation would be set to the standard value p_mut because
       * cellCB::ini() is called via cellman::put_CB(). Then, one daughter
       * is always reset to the standard value irrespective of the number
       * of selections of the mother. This is true in a more node CB,
       * where no setting implies a symmetric daughter.
       */
      // if antigen is retained, decide whether division will be symmetric or asymmetric
      bool asymmetric_division = false;
      if (retain_ag) {
         if (CB_list[li].DEC205_ova || (drandom() < divide_ag_asymmetric)) {
            asymmetric_division = true;
         }
      }
      /*
       * if (CB_list[li].DEC205_ova) {
       * if (asymmetric_division )
       *  cout<<"asymmetric ... retained="<<CB_list[li].retained_ag<<" ... ";
       * else cout<<" symmetric ... retained="<<CB_list[li].retained_ag<<" ... ";
       * cout<<"\n";
       * }
       */
      // Now switch Ig-class of both daughters if set:
      if (immunoglobulin_class::do_switch_classes == 2) {
         CB_list[li].IgX.class_switch();
         CB_list[newli].IgX.class_switch();
      }
      //MSc
      CB_list[newli].nTFRcontacts=CB_list[li].nTFRcontacts;
      CB_list[newli].ccl3KO = CB_list[li].ccl3KO;

      // Now mutate both daughters if appropriate
      if (time >= mutation_start_time) {
          //MSchips

         // cout<<"in proliferate_CB-2: CB["<<li<<"].p_mutation="<<CB_list[li].p_mutation<<"; "
         // <<"CB["<<newli<<"].p_mutation="<<CB_list[newli].p_mutation<<"\n";
         // suppress mutation if antigen was retained and the flag to stop mutation is set
         if (not (asymmetric_division)
             || not (retain_ag && (CB_list[li].retained_ag > 0.)
                     && cellCB::ag_loaded_CB_stop_mutation)) {
            // NOTE: if cellCB::asymmetric_polarity_index<1. every little retained antigen
            //       will induce a stop of mutations if cellCB::ag_loaded_CB_stop_mutation is set.
            //       Added a warning in setparam.C
            // mutate the old and new cell
            long oldpss_mother=CB_list[li].pos_ss;
            bool isSelf=CB_list[li].selfMutation, isRedeemed=CB_list[li].redeemed;
            CB_list[newli].redeemed=isRedeemed;
            double pred=p_redeem;

            CB_list[li].mutate(shape);
            CB_list[newli].mutate(shape);

            bool cell1_mutated=0, cell2_mutated=0;
            if (CB_list[li].pos_ss!=oldpss_mother) { //daughter1 did mutate
                    cell1_mutated=1;
            }
            if (CB_list[newli].pos_ss!=oldpss_mother) { //daughter2 did mutate
                    cell2_mutated=1;
            }

            if (SelfReact_mode==0) { // 0: With probability p_selfMut at mutation
                if (isSelf) { //-3 --> no inheritance
                    if (p_redeem>0) {
                        if (cell1_mutated &&
                                drandom()<pred) {CB_list[li].selfMutation=0; CB_list[li].redeemed=1;++newly_generated_redeemed;/*redeemed*/}
                        else {CB_list[li].selfMutation=1;}
                        if (cell2_mutated &&
                                drandom()<pred) {CB_list[newli].selfMutation=0; CB_list[newli].redeemed=1;++newly_generated_redeemed;/*redeemed*/}
                        else {CB_list[newli].selfMutation=1;}
                    } else {CB_list[li].selfMutation=1; CB_list[newli].selfMutation=1;}
                } else {
                    if (cell1_mutated && !isRedeemed) { //mutated and not already redeemed cell
                        if (drandom()<p_selfMut) {CB_list[li].selfMutation=1;++newly_generated_sCBs;}
                        else {CB_list[li].selfMutation=isSelf;}
                    } else {CB_list[li].selfMutation=isSelf;}
                    if (cell2_mutated && !isRedeemed) {
                        if (drandom()<p_selfMut) {CB_list[newli].selfMutation=1;++newly_generated_sCBs;}
                        else {CB_list[newli].selfMutation=isSelf;}
                    } else {CB_list[newli].selfMutation=isSelf;}
                }
            }
            // cout<<"mutate ... ";
         }
         // OPTION ++++++++++++++++++++
         // Do symmetric divisions also lead to suppression of mutation if antigen was retained?
         // If the first line in the if condition is active (not(asymmetric_division) ||
         //  mutation may only be suppressed for asymmetric division
         // If the first line is inactive
         //  mutation is also suppressed if antigen is retained but cells divided symmetrically!
         // end OPTION
      }


      // Control the handling of collected/retained antigen
      if (retain_ag) {
         // otherwise nothing is to do and retained_ag is set zero
         // By default both daughters have the full amount of nFDCcontacts
         // divide_ag_asymmetric is a % which allows to make a fraction of divisions symmetric:
         // cout<<"total-ag="<<CB_list[li].retained_ag;
         if (asymmetric_division) {
            // CB_list[newli].retained_ag=0; // and keep the other daughter as is
            // cout<<CB_list[li].retained_ag<<"; ";
            CB_list[newli].iamhighag = false;
            double all_ag = CB_list[li].retained_ag;
            double pitmp = cellCB::asymmetric_polarity_index;
            // +++++++++++ OPTION +++++++++++++++++++++++++++++++++++++++++++++
            // vary the PI value around its mean
            double shift = 0.;
            double gaussf;
            double width = 0.;
            // Define the width of a smooth distribution around the PI set in the parameter file:
            if (pitmp < 1. - 4.0 * cellCB::smooth_PI) {
               width = cellCB::smooth_PI * pitmp;
            } else {
               width = cellCB::smooth_PI * ((1.0 - pitmp)) * pitmp;       // linear switch from
                                                                          // smooth_PI to 0%
            }
            // ... this generates a fairly symmetric distribution around the set PI value (pitmp).
            // However, it becomes asymmetric when smooth_PI is too large compared to pitmp.
            // Note that the asymmetric distribution is generated because
            // large values are cut off below in the while loop.
            double tmp = 3.0;
	    // ##### make a proper sampling from a Gaussian here.
	    // ##### introduce the polarity index of 0.88 (Batista) into the simulations
            while (tmp<0. || tmp> 1.) {
               gaussf = drandom();
               if ((gaussf == 1.) || (gaussf == 0.)) {
                  shift = 0;
               } else {
                  shift = width * log((1. - gaussf) / gaussf);
               }
               tmp = (pitmp + shift);
               // Note that with cellCB::smooth_PI=0., everything still works and results in
               // tmp=pitmp
               // here.
               // cout<<shift<<", "<<tmp<<"; ";
            }
            // ### Can't we just use the sampling function I defined somewhere ? (Michael)
            // cout<<"\n";
            pitmp = tmp;
            // +++++++ end OPTION +++++++++++++++++++++++++++++++++++++++++++++
            // pitmp contains the fraction of ag which is kept in CB_list[li]
            CB_list[li].retained_ag = pitmp * all_ag;
            CB_list[li].adapt_specific_ag(pitmp);
            CB_list[newli].retained_ag = all_ag - CB_list[li].retained_ag;
            CB_list[newli].adapt_specific_ag(1 - pitmp);
         } else {
            // divide symmetric:
            // OPTION +++++++++++++++++++
            // Either distribute retained antigen on both daughters as symmetric as possible:
            double all_ag = CB_list[li].retained_ag;
            CB_list[newli].retained_ag = all_ag / 2.;
            CB_list[newli].adapt_specific_ag(0.5);
            CB_list[li].retained_ag = all_ag - CB_list[newli].retained_ag;
            CB_list[li].adapt_specific_ag(0.5);
            CB_list[newli].iamhighag = false;
            CB_list[li].iamhighag = false;
            // Note that this prohibits differentiation to output after symmetric
            // divisions when retained_ag is used and ag_loaded_CB_diff2output=true
            // (see cellthis::ask\_differentiate(..)).
            // As this also doesn't make sense, the combination of symmetric division
            // and ag_loaded_CB_diff2output=true is prohibited in setparam.C.
            // end OPTION ++++++++++++++++++++
         }
         // cout<<"; CB_list[li].retained_ag="<<CB_list[li].retained_ag<<";
         // CB_list[newli].retained_ag="<<CB_list[newli].retained_ag<<"\n";
      }

      // Here the two daughter cells are written as events to the brainbow class
      if (track_mutations) {
         // The two daughters are CB_list[li] and CB_list[newli]
         // The mother cell was CB_list[li] before division
         // Thus, CB_list[li].brainbow_index contains the event index of the birth of the mother
         // cell
         long mother_index = CB_list[li].brainbow_index;
         /* Call for both new cells the brainbow writing routine 
	  * and save the returned event index:
	  * parameters: 
	  * (double time, bool founder, bool birth, long mother_index, long ss_position)
	  */
         CB_list[li].brainbow_index = trackmutations.write(time,
                                                           false,
                                                           true,
                                                           mother_index,
                                                           CB_list[li].pos_ss);
         CB_list[newli].brainbow_index
            = trackmutations.write(time, false, true, mother_index, CB_list[newli].pos_ss);
      }

      // The following is for the case that DEC205+ cells are tracked:
      if ((CB_list[newli].DEC205)    // || CB_list[li].trackit)
          &&
          photoactivation && (time > photoactivation_t0) && (time < TRACK::TRACKUNTIL)) {
         // track it
         cout << "TRACKit ................................!!!!!!!!!\n";
         double r[l.dim];
         l.get_koord(CB_list[newli].index, r);
         CB_list[newli].trackit = true;
         double l2saxis = CB_list[newli].get_long2short_axis(newli, CB, l);
         double elongat = CB_list[newli].get_elongation(newli, CB, l);
         CB_list[newli].trackno = trackdata.add_new_track(CB,
                                                          time,
                                                          r,
                                                          CB_list[newli].polarity,
                                                          elongat,
                                                          l2saxis);
      }
   }
   return err;
}
void cellman::calc_CB(long * m,
                      long mlang,
                      short int ss_save,
                      space &l,
                      sigs &s,
                      AffinitySpace &shape,
                      dynarray<long> &redo) {
   /* m die random Reihenfolge der Indizes der CB auf dem Gitter,
    * also von CB_list[r].index-Werten */
   long li, n;
   short int errd, errp;
   int chose;
   cellCB::set_differentiation(time);

   // +++ OPTION: diminish proliferation rate with time; rewrite with cellCB::CB_list.p_pro-=...
   // CB_list[0].p_pro-=CB_list[0].delta_p_pro;
   // +++ OPTION end

   // for (n=0; n<mlang; n++) if (l.cellknot[m[n]].listi==-1)
   // cout<<"at n="<<n<<",m[n]="<<m[n]<<",listi=-1;;";
   for (n = 0; n < mlang; n++) {
      /* Hier wird random gewaehlt, welche Aktion zuerst probiert wird.
       * Die uebrigen dann nur probieren, wenn keine Aktion durchgefuehrt
       * wurde. Diffusion geht immer solange die Zelle noch CB ist.
       */

      // cerr<<"CB"<<m[n]<<"(max"<<mlang<<")->";

      /*
       * cout<<"n="<<n<<",m[n]="<<m[n]<<",listi="<<l.cellknot[m[n]].listi<<": ";
       * for (int nn=n; nn<mlang; nn++) if (l.cellknot[m[nn]].listi==-1)
       * cout<<"at nn="<<nn<<",m[nn]="<<m[nn]<<",listi=-1;;";
       * cout<<" -- ";
       */

      // Get index on the cell list
      li = l.cellknot[m[n]].listi;
      // cerr<<"n="<<n<<",m[n]="<<m[n]<<",li="<<li
      //	 <<"(x,y)=("<<l.knot[m[n]].x[0]<<","<<l.knot[m[n]].x[1]<<");;\n";

      // Set clock
      CB_list[li].set_clock();
      // if (CB_list[li].state>cb_G2)
      // cerr<<"s="<<CB_list[li].state<<" a="<<CB_list[li].time_of_cycle_state_switch<<";;; ";

      // Set adhesion
      /* Later on, when the expression of adhesion molecules becomes dynamical,
       * the calculation will be initiated by
       * "CB_list[li].set_adhesion();"
       * For the moment this is simply waste of time.
       ###
       */

      // do a lot of checks (delete after a while ###)
      if (checkit == 2) {
         if (li == -1) {
            cout << "call for listindex = -1 at lattice p= " << m[n] << "!\n";
            exit(1);
         }
         if (l.cellknot[m[n]].cell != CB) {
            cout << "Fehler bei lattice point " << m[n] << " and cell#" << li
                 << " in calc_CB: Zelltyp = " << l.cellknot[m[n]].cell << "\n";
            exit(1);
         }
         // Barycenter ueberpruefen
         CB_list[li].check_barycenter(l);
      }

      /*
       * // Output zu cbqual.out
       * if (ss_save==1) cbquality<<time<<"  "<<m[n]<<"  "<<CB_list[li].pos_ss<<"  "
       * <<shape.best_affinity(CB_list[li].pos_ss)<<"\n";
       */

      // =====================================================
      // ============== Start of Actions: ==============================
      // =====================================================
      // Zustand bestimmen (cycling_time erhöhen, Differenzierungs-Potenzial bestimmen):
      // cerr<<"vor .get_new_state() ...\n";
      CB_list[li].get_new_state(m[n], dt, l, s);
      // if (CB_list[li].state>cb_G2)
      // cerr<<"after get_new_state: s="<<CB_list[li].state<<"
      // a="<<CB_list[li].time_of_cycle_state_switch<<";;; ";
      // if (li==20) cout<<"CT="<<CB_list[li].cycling_time<<"; ";
      if ((s.signal_use[glucose] == 1) && (s.signal_use[oxygen] == 1)) {
         CB_list[li].use_nutrient(s, l.knot[CB_list[li].index].guest_grid_pointer);
      }

      if (CB_list[li].status == necrotic) {
         macrophagocyte_CB(m[n], li, l, shape);
      } else {
         // Mache Proliferation, Wachstum, Differenzierung und Diffusion
         CB_list[li].set_CXCR4expression();
         CB_list[li].resensitise4CXCL12(s);

         // cerr<<"vor proliferate and differentiate ...\n";
         if (cellCB::fixed_number_of_divisions()) {
            errd = 1;
            if (CB_list[li].state == cb_divide) {
               // cerr<<"will proliferate ... ";
               errp = proliferate_CB(m[n], li, l, shape);
               // cerr<<"errp="<<errp<<" s="<<CB_list[li].state<<"
               // a="<<CB_list[li].time_of_cycle_state_switch<<";;; ";
            } else if (CB_list[li].state == cb_stop_dividing) {
               errd = CB_list[li].ask_differentiate();       // errd==0 if, yes, differentiate!
               if (errd == 0) {
                  // cerr<<"will differentiate ... ";
                  // if (CB_list[li].DEC205_ova && CB_list[li].diff2output)
                  // cout<<CB_list[li].retained_ag<<": XXXXXXXXXXX\n";
                  if (CB_list[li].diff2output) {
                     CB_differ2OUT(m[n], li, l, shape);
                  } else {
                     differ2CC_CB(m[n], li, l, shape);
                  }
               }
            }
         } else {
            errd = 1;
            errp = 1;
            cerr<<"'here1'"<<endl;exit(1);
            chose = irandom(2);
            if (chose == 0) {
               errd = CB_list[li].ask_differentiate();       // errd==0 if, yes, differentiate!
               if (errd == 0) {
                  if (CB_list[li].diff2output) {
                     CB_differ2OUT(m[n], li, l, shape);
                  } else {
                     differ2CC_CB(m[n], li, l, shape);
                  }
               }
            }
            if (chose == 1) {
               errp = proliferate_CB(m[n], li, l, shape);
            }
            // Differentiate wurde probiert. Falls misserfolg, try proliferate
            if ((chose == 0) && (errd != 0)) {
               errp = proliferate_CB(m[n], li, l, shape);
            }
            // Proliferate wurde probiert.   Falls misserfolg, try differentiate
            if ((chose == 1) && (errp != 0)) {
               errd = CB_list[li].ask_differentiate();
               if (errd == 0) {
                  if (CB_list[li].diff2output) {
                     CB_differ2OUT(m[n], li, l, shape);
                  } else {
                     differ2CC_CB(m[n], li, l, shape);
                  }
               }
            }
         }     // end else
               // errd==0 if CB differentiated to another cell subset.

         // Wenn die Zelle noch CB ist, kann sie diffundieren oder wachsen:
         if ((errd != 0) && (CB_list[li].state != cb_M)) {
            // call movement:
            double moved = CB_list[li].move(li, l, s, trackdata, time);

            // Treat the case of initiated move but suppression by lack of space:
            if (CB_list[li].contact_inhibited) {
               // this also implies volume==1 and moved==0.
               // reset contact_inhibited in order to ensure that it remains false in general
               CB_list[li].contact_inhibited = false;
               // save this cell for later processing
               redo.add(CB_list[li].index);
               // moved=try2exchange_cells(m[n],li,CB,CB_list[li].polarity,l);
               // cout<<"  back in calc_CB if contact_inhibited ...\n";
            }
            // write to movement file ...
            if (CB_list[li].writethis2track != trackini) {
               // cerr<<"Write this 2 track ...\n";
               if (int (v_resolution * moved + 0.5) <= v_resolution) {
                  ++velocity[int (v_resolution * moved + 0.5)];
               }
               // calculate total movements:
               long tmp[l.dim];
               l.get_koord(CB_list[li].born_index, tmp);
               double pathlength = l.get_2norm(CB_list[li].barycenter, tmp);

               // not used for trackdata!
               double tmpi[l.dim];
               l.get_koord(CB_list[li].index, tmpi);

               // do not use index here!
               // Index is the lattice point nearest of the barycenter.
               // It is actualised after every call of frag_cell::get_barycenter(...)
               // in frag_cell::fragdiffuse(...), thus, automatically.

               double l2saxis = CB_list[li].get_long2short_axis(li, CB, l);
               double elongat = CB_list[li].get_elongation(li, CB, l);
               trackdata.Write_movement(CB_list[li].trackno,
                                        CB,
                                        time,
                                        CB_list[li].barycenter,
                                        CB_list[li].polarity,
                                        elongat,
                                        l2saxis,       // movement,
                                        CB_list[li].writethis2track);
               CB_list[li].writethis2track = trackini;

               movements << time * 60. << "   " << sqrt(time * 60.) << "   "  // in min
                         << CB_list[li].volume << "   "                       // in fragments
                         << moved << "   "   // distance in this time step in lattice constants
                         << moved * l.dx / (60. * dt * double (CB_list[li].n_immobile))
                         << "   "           // v [microns/min]
                         << pathlength * l.dx << "   " // total walk from born position [microns]
                         << pathlength * pathlength * l.dx * l.dx
                  / (l.dim2 * pow(sqrt(60. * time) - sqrt(60. * CB_list[li].born_time), 2.))
                         << "    "       // associated diffusion constant
                         << (tmpi[0] - tmp[0]) * l.dx << "   " << (tmpi[1] - tmp[1]) * l.dx;
               // distance to born position in each dimension in microns
               if (l.dim > 2) {
                  movements << "   " << (tmpi[2] - tmp[2]) * l.dx;
               }
               movements << "\n";
               CB_list[li].n_immobile = 1;
               // cout<<"wrote to movement\n";
            }

            if ((outputfiles == 0)
                && (60. * time / CB_list[0].deltat_v
                    - double (long (60. * time / CB_list[0].deltat_v))
                    < 60. * dt / CB_list[0].deltat_v)) {
               // cerr<<"in if ((outputfiles==0) ...) { ...\n";
               long tmp[l.dim];
               l.get_koord(CB_list[li].born_index, tmp);
               double pathlength = l.get_2norm(CB_list[li].barycenter, tmp);
               long tmpi[l.dim];
               l.get_koord(CB_list[li].index, tmpi);
               // calculate velocity in microns/min
               double v = CB_list[li].add_moves * l.dx / (CB_list[0].deltat_v);
               // calculate also the velocity on the basis of the displacement since last
               // measurement
               double v2 = l.get_2norm(CB_list[li].barycenter, CB_list[li].last_position);
               v2 *= (l.dx / CB_list[0].deltat_v);
               // write to file
               move_vdt << time * 60. << "   " << sqrt(time * 60.) << "   "
                        << CB_list[li].volume << "   "
                        << CB_list[li].add_moves << "   " << v << "   " << v2 << "   "
                        << pathlength * l.dx
                        << "   "
                        << pathlength * pathlength * l.dx * l.dx
                  / (l.dim2
                     * pow(sqrt(60. * time) - sqrt(60. * CB_list[li].born_time),
                           2.)) << "    "
                        << (tmpi[0] - tmp[0]) * l.dx << "   " << (tmpi[1] - tmp[1]) * l.dx;
               if (l.dim > 2) {
                  move_vdt << "   " << (tmpi[2] - tmp[2]) * l.dx;
               }
               move_vdt << "\n";
               // write to movement array (v/delta_v is the index of the velocity-array):
               if (int (v / delta_v + 0.5) <= v_resolution) {
                  ++delta_velocity[int (v / delta_v + 0.5)];
               }
               if (int (v2 / delta_v + 0.5) <= v_resolution) {
                  ++delta_velocity_euklid[int (v2 / delta_v + 0.5)];
               }
               // Note: no round, done for bar diagrams: 0..delta_v = 0; delta_v..2delta_v=1
               // This implies, that the v-value in the output file is to be set to
               // n*delta_v+delta_v/2

               // Reset add_moves of this cell:
               CB_list[li].add_moves = 0.;
               for (short a = 0; a < l.dim; a++) {
                  CB_list[li].last_position[a] = CB_list[li].barycenter[a];
               }

               // get long2short axis parameter ...
               double longaxis, shortaxis;
               double l2s = CB_list[li].get_long2short_axis(li,
                                                            l.cellknot[m[n]].cell,
                                                            l,
                                                            longaxis,
                                                            shortaxis);
               // and write into the axis-file
               axis << time * 60. << "   " << longaxis * l.dx << "   " << shortaxis * l.dx
                    << "   " << l2s << "\n";
            }
            // end of write to movement files

            if (CB_list[li].changed_polarity == 1) {
               CB_list[li].changed_polarity = 0;
               if (outputfiles == 0) {
                  polarities << li << "   " << 60. * time << "   0\n";
               }
            }
         }

         if ((errd != 0) && (CB_list[li].state != cb_G0) && (CB_list[li].state != cb_M)
             && (CB_list[li].state != cb_S)
             && (CB_list[li].state != cb_stop_dividing)) {
            // call cell growth:
            CB_list[li].grow(li, l);
         }
      }    // not necrotic
   }
   // now that all actions with pre-existing CB are finished, try having new ones
   BCinflux(l, shape, s);
}
// =======================================================
// Centrocytes ===========================================
// =======================================================

short int cellman::put_CC(long i, cellCC &newCC, space &l, AffinitySpace &shape) {
   // Test if point i empty
   if (l.cellknot[i].cell != nocell) {
      return 1;
   }
   // save index on CC_list to be used:
   long li = CC_list.benutzt();
   // Write .cell .pos_ss .ccstate
   l.set_knot(i, CC, li);
   // l.knot[i].cell=CC;
   // Write on cell_list und save list-index on the lattice
   // l.knot[i].listi=CC_list.benutzt();
   newCC.index = i;
   newCC.selectable = 1;
   newCC.mobile = 1;
   newCC.clock = 0;
   newCC.affinity = shape.best_affinity_norm(newCC.pos_ss);
   CC_list.add(newCC);
   // shape space:
   shape.add_cell(sCCunselected, newCC.pos_ss);
   shape.add_cell(sCC, newCC.pos_ss);
   shape.add_cell(total, newCC.pos_ss);
   //MS
//   if (newCC.selfMutation) {
//       shape.add_cell(sCCself,newCC.pos_ss);
//   }
   CC_list[li].trackfate_initialize();
   CC_list[li].trackfate(time,true);
   //CC_list[li].trackfate_show();
   // extra treatment if new CC is supposed to directly interact with TFH:
   if (cellCC::ag_loaded_CC_directly2TFH && (CC_list[li].nFDCcontacts > 0)) {
      CC_list[li].go2TCselection(shape);
      /* i.e. when pre-loaded antigen is present from a previous round of selection
       * and still there after (asymmetric) division in the CB state,
       * these cells (when back to the CC state) are assumed to directly
       * search for TFH and not to search for further antigen on FDC
       * if the flag ag_loaded_CC_directly2TFH is set true
       */
      // cout<<"Pre-loaded CB with "<<CC_list[li].nFDCcontacts<<" ag-portions goes for TFH!\n";
   }
   // Fuer die Statistik, die Zahl der erzeugten CCs zaehlen:
   ++CC_total;

   return 0;
}
/*
 * long int cellman::find_CC(long i) {
 * long int found=-1;
 * long int n=0;
 * long int bis=CC_list.benutzt();
 * while (n<bis && found==-1) {
 *    if (CC_list[n].index==i) { found=n; }
 *    else { n++; }
 * }
 * if (found==-1) { cout<<"Error: Index "<<i<<" not found in find_CC!"; exit(1); }
 * return found;
 * }
 */

short int cellman::del_CC(long i, long li, space &l) {
  CC_list[li].trackfate(time,true);
  CC_list[li].writetrackfate();
  // remove on lattice
  l.clear_knot(i);
  // l.knot[i].cell=empty;
  // l.knot[i].listi=-1;
  // Delete cell on position li on the list
  CC_list.erase_jump(li);
  // Correct on the lattice for the shift of last cell on the list to position li
  // Don't if the deleted cell is the last cell on the list
  if (li < CC_list.benutzt()) {
    l.cellknot[CC_list[li].index].listi = li;
  }
  return 0;
}
short cellman::CC_differ2out(long i, long li, space &l, AffinitySpace &shape) {
   // adjust global counters
   ++n_outs;
   if (CC_list[li].n_recycling > 0) {
      ++n_outs_with_recycling;
   }
   n_muts += CC_list[li].n_mutation;
   n_recmuts += CC_list[li].n_recandmute;
   if (CC_list[li].n_mutation < max_mutation_bin) {
      ++mutation_frequency[CC_list[li].n_mutation];
   } else {
      ++mutation_frequency[max_mutation_bin - 1];
   }

   // Save old cell states
   cellOUT newOUT;
   newOUT = CC_list[li];
   // use the position of CC as barycenter of the new OUT cell.
   for (short a = 0; a < l.dim; a++) {
      newOUT.barycenter[a] = l.knot[i].x[a];
   }
   newOUT.index = i;
   newOUT.volume = 1;
   newOUT.get_radius(l.dim);
   newOUT.responsive2signal[CXCL13] = false;
   newOUT.responsive2signal[CXCL12] = false;

   // Delete CC from shape space, lattice and CC_list
   shape.rem_cell(sCCselected, CC_list[li].pos_ss);
   shape.rem_cell(sCC, CC_list[li].pos_ss);
   del_CC(i, li, l);

   // Create new OUT-cell:
   l.set_knot(i, out, OUT_list.benutzt());
   // Write OUT-cell on cell_list und save list-index on the lattice
   OUT_list.add(newOUT);
   // correct shape space:
   shape.add_cell(sout, newOUT.pos_ss);
   if (newOUT.DEC205) {
      shape.add_cell(soutdec, newOUT.pos_ss);
   }
   //MSchips
   if (newOUT.selfMutation) {
       shape.add_cell(soutSelf,newOUT.pos_ss);
   } else {shape.add_cell(soutNonSelf, newOUT.pos_ss);}

   // save this output event in the Ig-class output counter
   ++integrate_out[newOUT.IgX.Ig_class];
   ++integrate_out[nIg_classes];

   double aff = shape.best_affinity_norm(newOUT.pos_ss);
   // correct average affinity of output cells
   shape.correct_average_affinity(sout, newOUT.pos_ss, cellOUT::average_affinity);
   if (aff > cellOUT::max_affinity) {
      cellOUT::max_affinity = aff;
   }

   // Write to log_out_aff
   write_log_out_aff(aff, newOUT.pos_ss, shape);

   newOUT.deliberate_memory();

   return 0;
}
short cellman::CC_differ2CB(long i, long li, space &l, AffinitySpace &shape) {
   ++n_recycling_events;
   // alte Werte in neuen CB uebernehmen
   cellCB newCB;
   // Werte vom CC uebernehmen (volume, fragment ... definieren)
   newCB = CC_list[li];

   // ###### some of the following might be shifted to the operator=
   newCB.responsive2signal[CXCL13] = false;
   newCB.responsive2signal[CXCL12] = true;

   // Speichere, dass recycled wurde:
   ++newCB.n_recycling;
   // if ag is not retained, reset this variable (by default it is transmitted from nFDCcontacts)
   if (retain_ag == false) {
      newCB.retained_ag = 0.;
   } else if (newCB.retained_ag > 0.) {
      newCB.iamhighag = true;
   }
   // this last condition defines newly differentiated CB by default as highag
   // unless they were selected without having collected ag (which is not possible?)
   // cout<<"retained antigen in diff CC to CB = "<<newCB.retained_ag<<"\n";
   // Set the mutation probability to the one after TC contact
   newCB.set_mutation_after_TC(shape);
   // cerr<<"CC selected: CB["<<CB_list.benutzt()<<"].p_mutation="<<newCB.p_mutation<<"\n";
   // cerr<<"CC selected: CB["<<CB_list.benutzt()<<"].state="<<newCB.state<<"\n";
   for (short a = 0; a < l.dim; a++) {
      newCB.barycenter[a] = l.knot[i].x[a];
   }
   newCB.get_radius(l.dim);
   // force differentiation to output after division if chosen:
   if (newCB.DEC205_ova) {
      if ((DEC205_forces_output > 0.) && (drandom() < DEC205_forces_output)) {
         newCB.diff2output = true;
      }
      if (retain_DEC205_ag) {
         // note that this also implies retain_ag=true (see setparam.C)
         double factor = (dec205_max_antigen + newCB.retained_ag) / newCB.retained_ag;
         newCB.retained_ag += dec205_max_antigen;     // force selection
         newCB.adapt_specific_ag(factor);
         newCB.iamhighag = true;                      // and save that this cell state
         // for eventual final differentiation, see cellthis::ask_differentiate(...)
      }
   } else {
      // (if not(newCB.DEC205_ova), so in the normal case
      if ((p_CB2OUT > 0.) && (drandom() < p_CB2OUT)) {
         newCB.diff2output = true;
      }
   }
   // if differentiation was delayed and option is set start later in the cell cycle:
   if (((CC_list[li].DEC205_ova == false) && (CC_list[li].dif_delay > 0.))
       || ((CC_list[li].DEC205_ova == true) && (CC_list[li].dif_delay_DEC > 0.))) {
      newCB.transmit_CCdelay2cycle(CC_list[li].selected_clock);
   }

   del_CC(i, li, l);

   // CC wird wieder zu CB:
   l.set_knot(i, CB, CB_list.benutzt());
   // l.knot[i].cell=CB;
   // Write on cell_list und save list-index on the lattice
   // l.knot[i].listi=CB_list.benutzt();
   // cout<<"differ2CB: newCB.n_divisions2do="<<newCB.n_divisions2do<<"\n";
   CB_list.add(newCB);
   // CB_list[l.knot[i].listi].get_barycenter(CB,l.knot[i].listi,l);
   // CB_list[l.cellknot[i].listi].get_barycenter(l);
   // !!! is the barycenter not necessarily the position of the single fragment of the CC???
   // Yes it is!

   // shape space:
   shape.add_cell(sCB, newCB.pos_ss);
   shape.rem_cell(sCCselected, newCB.pos_ss);
   shape.rem_cell(sCC, newCB.pos_ss);

   newCB.deliberate_memory();

   return 0;
}

short int cellman::get_into_FDC_contact(long li, AffinitySpace &shape, AntibodyDyn &Ab,
					space &l, TRACK &td, double &time) {
  short int err = 1;
  long fdc_i;
  if (CC_list[li].selectable == 1) {
    fdc_i = CC_list[li].contact2FDC(l);
  } else {
    fdc_i = -1;
  }
  //MSc
  ///if TfrModel==3 S-CCs that bound Tfr should not be allowed ag-consumption
  bool tfr_allow_bind=1;
  if ((TFR_mode==3||TFR_mode==16) && (CC_list[li].selfMutation && CC_list[li].nTFRcontacts>0)) {
      tfr_allow_bind=0;
  }
  if (fdc_i != -1 && tfr_allow_bind) {
    // Get the index of the FDC in question on the list
    long int fdc_li = l.cellknot[fdc_i].FDClisti;
    // Determine the index of the FDC-fragment on the FDC fragment list
    int FDC_frag_index = FDC_list[fdc_li].get_fragment_index(fdc_i);
    /* Determine the antigen on the FDC fragment for interaction with the BCR,
     * where FDC_ag_index is the index on the antigen list in AffinitySpace */
    int FDC_ag_index
      = FDC_list[fdc_li].local_interaction_with_ag(FDC_frag_index,
						   CC_list[li].pos_ss,
						   shape);
    if (FDC_ag_index == -1) { // FDC_ag_index == -1 if no antigen was there.
      //cerr << "local_interaction_with_ag on FDC failed!\n"; 
      ++TestedFDCsiteVoidOfAg;
    }
    // Determine a threshold for binding
    double threshold = 0.;
    if (cellOUT::use_threshold > 0) {
      if (cellOUT::use_threshold == 1) {
	threshold = cellOUT::average_affinity;
      } else if (cellOUT::use_threshold == 2) {
	threshold = cellOUT::max_affinity;
      } else if (cellOUT::use_threshold == 3) {
	threshold = Ab.average_ab_affinity(FDC_ag_index);
	// returns threshold=1 (=no binding in bind_antigen(...)) if
	// FDC_ag_index==-1
      }
    }
    //cerr << cellOUT::use_threshold << "->" << threshold << "; ";
    // Try to bind_antigen with these parameters
    err = CC_list[li].bind_antigen(FDC_list[fdc_li],
				   FDC_frag_index,
				   FDC_ag_index,
				   shape,
				   threshold);
    /* err==0 if binding was successful, then suppress movement, and write incontact to TRACK;
     * err==1 if binding was unsuccessful, allow movement but suppress interactions
     */
    if (err == 0) {
      CC_list[li].trackfate(time,true);      
      if (CC_list[li].trackit) {
	double tmpi[l.dim];
	l.get_koord(CC_list[li].index,tmpi);
	td.Write_movement(CC_list[li].trackno, CC, time, tmpi,
			  CC_list[li].polarity, 1, 1,
			  CC_list[li].fdc_clock, incontact);
      }
    }
    /* If CXCR5 down-regulation is not regulated by a rate
     * do it by hand upon first FDC encounter: */
    // if (cellCC::p_CXCR5down<=0.) CC_list[li].responsive2signal[CXCL13]=false;
  }
  return err;
}



void cellman::calc_CC(long * m,
                      long mlang,
                      space &l,
                      sigs &s,
                      AntibodyDyn &Ab,
                      AffinitySpace &shape,
                      dynarray<long> &redo) {
    long li, n;
    // Randomize the sequence of points that are actualized
    short int err;
    bool CCexists;
    cellCC::set_differentiation(time);   // to update time-dependent probabilities
    for (n = 0; n < mlang; n++) {
        // Get index on the cell list
        li = l.cellknot[m[n]].listi;
        if (TFR_mode!=CC_list[li].TFR_CC_interaction_mode) {cerr<<"wronge mode\n";exit(1);}
        CCexists = true;
        CC_list[li].aging(dt);
        CC_list[li].progress_signalling(); // mTORC1 and FoxO
        CC_list[li].addTFHsignalAUC(dt);
        CC_list[li].decay_TFHsignal();
        CC_list[li].decay_cMyc();
        CC_list[li].trackfate(time,false);

        //MS newout for correlation-like plots
        int status=-1;
        bool is_self=0;
        if (time>double(first_record+time_window-dt) && time<=double(first_record+time_window)) {
            //here no accumulation is needed --> just write at the end of the time window
            double tfhsig=CC_list[li].tc_signal,
                    pmhc=CC_list[li].nFDCcontacts,
                    aff=shape.best_affinity_norm(CC_list[li].pos_ss);
            //note for the status:
            // 0=cc; 1=apo; 2=selected
            if (CC_list[li].state==apoptosis) {
                status=1;
            } else if (CC_list[li].state==selected) {
                status=2;
            } else {status=0;}
            if (CC_list[li].selfMutation) {is_self=1;}
            write_info(time, tfhsig, pmhc, aff, status, is_self);
        }

        if (CC_list[li].trackit) {
            /* Here the contact times to FDCs are controlled.
             * When a BC has no physical contact to FDCs it is considered as detached from FDCs.
             * This really monitors physical contacts to a FDC
             * irrespective of actually being bound to the FDC.
             * This enters the read-out file fdccontact.out generated in the TRACK class.
             *
             * One might consider to monitor the times bound to FDC as well. TOBEDONE.
             * The information is in the TRACK class as well and need to be based on the
             * time between the events <incontact> and <offcontact>.
             */
            long fdc_in_contact = CC_list[li].contact2FDC(l);
            if (fdc_in_contact == -1) {
                // if FDC contact clock is running: save and reset
                if (CC_list[li].fdc_clock > 0.) {
                    double tmpi[l.dim];
                    l.get_koord(CC_list[li].index, tmpi);
                    trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                             CC_list[li].polarity, 1, 1,
                                             CC_list[li].fdc_clock, fdcdetachment);
                    CC_list[li].fdc_clock = 0.;
                }
                // ++++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++
                // Activate this if CC should be stoped to be tracked after positive FDC selection
                /*
                 * *if (CC_list[li].state==FDCselected
                 * *   || CC_list[li].state==TCcontact
                 * *   || CC_list[li].state==selected) {
                 * * double tmpi[l.dim];
                 * * l.get_koord(CC_list[li].index,tmpi);
                 * * trackdata.Stop_movement(CC_list[li].trackno,nocell,time,
                 * *                         tmpi,CC_list[li].polarity);
                 * * CC_list[li].trackit=false;
                 * *}
                 * */
                // ++++++++++++++++++++++ end OPTION ++++++++++++++++++++++++
            } else {
                // FDC contact
                /* +++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++
                 * * This version is used to analyse CC-FDC contact time excluding refractory times:
                 * * if ((CC_list[li].state==unselected && CC_list[li].selectable==1) ||
                 * * CC_list[li].state==contact)
                 * * This version is used to analyse CC-FDC contact time including refractory times:
                 * * if ((CC_list[li].state == unselected) || (CC_list[li].state == contact)) {
                 * * From 2018-01-02 both versions have to be complemented by the state FDCselected
                 * in order to account for simultaneous Tfh and FDC search.
                 * +++++++++++++++++++ end OPTION ++++++++++++++++++++++++++++++++
                 */
                 if ( (CC_list[li].state == unselected) || (CC_list[li].state == contact)
                      || (CC_list[li].state == FDCselected && cellCC::simultaneousTfhFDC) ) {
                     /* Increase clock of FDC contact
                      * * but only if the cell has the potential to interact with FDC */
                     CC_list[li].fdc_clock += dt;
                 }
            }
        }

        // Centrocytes do in dependence of their state
    switch (CC_list[li].state) {
    case unselected: {
        // cerr<<"unselected ... ";
        if ( ( cellCC::CC_FDC_selection == 0 ) || ( CC_list[li].DEC205_ova ) ) {
            // if (cellCC::p_CXCR5down<=0.) CC_list[li].responsive2signal[CXCL13]=false;
            CC_list[li].go2TCselection(shape);
            /* i.e. jump directly to TC selection! This has to be changed, if the binding
             * * process of soluble antigen and cell receptors is treated dynamically.
             * Then a binding routine has to be called here. ###
             *
             * Note that DEC205_ova positive cells are thought to acquire antigen
             * independent of FDC and BCR, thus, no FDC search in this case.
             * */
        } else {
            // Case that selection of CC has to go through FDC contact:
            // Possible switch to state apoptosis.
            if (not (cellCC::collectFDCsignals)) {
                err = CC_list[li].apoptose(time, shape);
            }
            // If the cell lives: err=1
            else {
                err = 1;  // in case collectFDCsignals is true.
            }
            if (err == 1) {
                // Continue here if cell still alive.
                /* If in contact to antigen then switch state to antigen-contact with
                 * affinity weight probability (depending on AffinitySpace position).
                 * Otherwise move the CC.
                 * */
                // cout<<" no apoptosis! ";
                CC_list[li].set_selectable();
                if (CC_list[li].CXCR5failure == 0) {
                    CC_list[li].set_CXCR5expression(); // time desensitisation
                    CC_list[li].resensitise4CXCL13(s); // undercritical CXCL13 resensitisation
                } else if (CC_list[li].CXCR5failure == 2) {
                    CC_list[li].set_CXCR4expression(); // time desensitisation
                    CC_list[li].resensitise4CXCL12(s); // undercritical CXCL12 resensitisation
                } // else if CC_list[li].CXCR5failure==1 do random walk
                // Ask whether serial collection of FDC signals shall be stopped:
                err = CC_list[li].stop_collecting_FDCsignals(shape, time, dt);
                // err==0: nothing happened; err==1: cell either dead or switched state
                if (err == 1) { CC_list[li].trackfate(time,true); }
                if (err == 0) { // continue collecting FDC signals:
                    err = get_into_FDC_contact(li, shape, Ab, l, trackdata, time);
                    // err==0 if FDC found and bound; err==1 otherwise.
                }
                else { err = 0; } // cell is dead or progressed in selection state --> no movement
                if ((err != 0) && (CC_list[li].mobile == 1)) {
                    CC_list[li].move(li, l, s, trackdata, time);
                }
            }
        }
        break; // if state has changed only work on this state in the next time step.
    }// end of case unselected.

    case contact: {
        // this means contact to FDC.
        // cout<<"contact ... ";
        CC_list[li].selected_clock += dt;
        // Assumption: No apoptosis if in contact to FDC.
        // The following controls the duration of the interaction with the FDC
        err = CC_list[li].FDCdetach(shape);
        /* Case simultaneousTfhFDC==true:
         * Currently, the first round of antigen uptake is done as before,
         * i.e. only interactions with FDCs are possible for the time of
         * antigen uptake as set in the parameter file.
         * Alternatively, one might switch to the state FDCselected upon first antigen uptake,
         * which is not yet programmed.
         */
        // if err == 0, the cell has detached for whatever reason, monitor in trackdata:
        if (err == 0) {
            CC_list[li].trackfate(time, true);
            if (CC_list[li].trackit) {
                double tmpi[l.dim];
                l.get_koord(CC_list[li].index,tmpi);
                trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                         CC_list[li].polarity, 1, 1,
                                         CC_list[li].fdc_clock, offcontact);
            }
        }
        // If not already done, downregulate CXCR5 expression here at the latest:
        CC_list[li].responsive2signal[CXCL13] = false;
        // No movement possible as in affin contact to FDC.
        break;
    }// end of case contact (to FDC).

    case FDCselected: {
        CC_list[li].time_before_tfr+=dt;
        // Case of search for Tfh help with terminated or simultaneous search for antigen from FDC.
        // Test for selection and progress FDCselected_clock
        short wasselected = CC_list[li].try_selection(time, dt, shape);
        if (wasselected == 1) { // either selected or apoptosis
            CC_list[li].trackfate(time,true);
        }
        if (wasselected == 0) {
            // Only continue when the CC was neither selected nor killed by apoptosis:
            if (CC_list[li].apoptose(time, shape) == 0/* && !exp_stop_apo*/) {//MS: returns 1 if Bcl2 exp is on
                // if apoptosis was done monitor this:
                CC_list[li].trackfate(time,true);
            } else { // i.e. the CC survived
                CC_list[li].set_selectable();
                if (CC_list[li].CXCR5failure == 0) {
                    CC_list[li].set_CXCR5expression();     // time desensitisation
                    CC_list[li].resensitise4CXCL13(s);     // undercritical CXCL13 resensitisation
                } else if (CC_list[li].CXCR5failure == 2) {
                    CC_list[li].set_CXCR4expression();     // time desensitisation
                    CC_list[li].resensitise4CXCL12(s);     // undercritical CXCL12 resensitisation
                } // else if CC_list[li].CXCR5failure==1 do random walk
                short int err = 2;
                if (cellCC::simultaneousTfhFDC) {
                    // Do a random decision whether FDC is tested first:
                    if ( drandom() < CC_list[li].get_prob2bindFDCfirst() ) {
                        // Test whether there is an FDC present and bind it
                        err = get_into_FDC_contact(li, shape, Ab, l, trackdata, time);
                        /*
                         * cout << "FDCtest1[" << li << "]: selectable=" << CC_list[li].selectable
                         * << "; err=" << err
                         * << "; CXCL13responsive=" << CC_list[li].responsive2signal[CXCL13]
                         * << "\n";
                         * */
                    }
                    /* err == 0 if FDC found and bound;
                     * err == 1 if failure to bind;
                     * err == 2 if FDC were not tested.
                     */
                }
                if (err > 0) {
                    //MSchips
                    //check the temporary tfr
                    long tfr_i=-1;
                    long tc_i = -1;
                    bool canbindtfr=0, canbindtfh=1, isSelf=CB_list[li].selfMutation;
                    short tobind=0; //=1->tfr =2->tfh
                    //check whether there is a temporary binding of tfr and cc
                    if (CC_list[li].canbindTFR(time)) {
                        tfr_i = CC_list[li].get_contact(TFR,l);
                    }
                    //CC scans neighb to find TFH
                    if (CC_list[li].try2findTFH()) {
                        tc_i = CC_list[li].get_contact(TC, l);
                    }
                    if (((tfr_i!=-1)
                            && !(CC_list[li].same_TFR_as_before(tfr_i,TFR_list[l.cellknot[tfr_i].listi].id))
                            && TFR_list[l.cellknot[tfr_i].listi].state == TFRnormal)) {
                        canbindtfr=1; //CC can bind the found tfr
                    }
                    if (tc_i == -1 || CC_list[li].same_TC_as_before(tc_i, TC_list[l.cellknot[tc_i].listi].id) ) {
                        canbindtfh=0; //Either no tfh or cc cannot bind the tfh
                    }
                    if (isSelf) {
                        if (canbindtfr) {tobind=1;} //self cell preferentially polise to tfr
                        else if (canbindtfh) {tobind=2;} //self binds tfh if no tfr
                    } else {//non self preferentially binds to tfh and binds tfr if no tfh available
                        if (canbindtfh) {tobind=2;}
                        else if (canbindtfr) {tobind=1;}
                    }
                    if (tobind==1) {
                        //free TFR was found
                        //                        if (CC_list[li].ccInd_TFRbound==-1) {
                        CC_list[li].bind_TFR(time, TFR_list[l.cellknot[tfr_i].listi], l, shape);
                        //                            CC_list[li].tfr_fate(TFR_list[l.cellknot[tfr_i].listi]);

                        CC_list[li].trackfate(time, true);
                        if (CC_list[li].trackit && err == 0) {
                            double tmpi[l.dim];
                            l.get_koord(CC_list[li].index,tmpi);
                            trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                                     CC_list[li].polarity, 1, 1,
                                                     CC_list[li].fdc_clock, inTFRcontact);
                        }
                    } else if (tobind==2) { //bind to tfh
                        if (tc_i!=TC_list[l.cellknot[tc_i].listi].index) {cerr<<"???wringtcind";exit(1);}

                        CC_list[li].bind_TC(time, TC_list[l.cellknot[tc_i].listi], l, shape);
                        /* This is always done when get_contact (to TC) was successful,
                         * * thus trackdata has to be set to incontact and fate is tracked: */
                        CC_list[li].trackfate(time, true);
                        if (CC_list[li].trackit && err == 0) {
                            double tmpi[l.dim];
                            l.get_koord(CC_list[li].index,tmpi);
                            trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                                     CC_list[li].polarity, 1, 1,
                                                     CC_list[li].fdc_clock, incontact);
                        }
                    } else { //no tfr -- notfh --> cell has to be tested for fdc if parameters require and move!
                        // Failure to find Tfh:
                        if (cellCC::simultaneousTfhFDC && err == 2) {
                            // Test whether there is an FDC present if not already tested before.
                            err = get_into_FDC_contact(li, shape, Ab, l, trackdata, time);
                            // err == 0 if FDC found and bound
                            /*
                             * cout << "FDCtest2[" << li << "]: selectable=" << CC_list[li].selectable
                             * << "; err=" << err
                             * << "; CXCL13responsive=" << CC_list[li].responsive2signal[CXCL13]
                             * << "\n";
                             */
                        }
                        if (err > 0) {
                            CC_list[li].move(li, l, s, trackdata, time);
                        }
                    }
                }
            } // end of if no apoptosis was done
        } // end of if (wasselected == 0)
        break;
    }// end of case FDCselected.

    case TCcontact: {
        if (CC_list[li].index!=CC_list[li].ccInd_TFHbound) {
            cerr<<" in TCcont index now: "<<CC_list[li].index
               <<" index at binding "<<CC_list[li].ccInd_TFHbound<<endl;
            exit(1);
        }
        // cout<<"TCcontact ... ";
        CC_list[li].responsive2signal[CXCL13] = false;
        err = CC_list[li].got_tc_signals(time, dt,
                                         TC_list[l.cellknot[CC_list[li].tc_index].listi],
                l, shape);
        if (err == 0 || err==2) {
            // the cell got deliberated from the TC, save in trackdata:
            CC_list[li].trackfate(time, true);
            if (CC_list[li].trackit) {
                double tmpi[l.dim];
                l.get_koord(CC_list[li].index,tmpi);
                trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                         CC_list[li].polarity, 1, 1,
                                         CC_list[li].fdc_clock, offcontact);
            }
        }
        break;
    } // end of case TCcontact.

    //MSchips
    case TFRcontact: {
    //        cerr<<"TFRcontact ... ";
        CC_list[li].responsive2signal[CXCL13] = false;
        if (CC_list[li].index!=CC_list[li].ccInd_TFRbound) {
            cerr<<"in case TFRcont index now: "<<CC_list[li].index
               <<" index at binding "<<CC_list[li].ccInd_TFRbound<<endl;
            exit(1);
        }
        short err_tfr=1;

        err_tfr = CC_list[li].got_tfr_interaction(time, dt,
                                              TFR_list[l.cellknot[CC_list[li].tfr_index].listi],l, shape);
        if (err_tfr == 0) {
            // the cell got deliberated from the TFR, save in trackdata:
            CC_list[li].trackfate(time, true);
            if (CC_list[li].trackit) {
                double tmpi[l.dim];
                l.get_koord(CC_list[li].index,tmpi);
                trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                         CC_list[li].polarity, 1, 1,
                                         CC_list[li].fdc_clock, offTFRcontact);
            }
        }
        break;
    } // end of case TFRcontact.

    case selected: {
        if  ((CC_list[li].tfr_tmp_index!=-1)) {cerr<<"not properlyupd\n";exit(1);}
        bool tfrbound=0;
      // cout<<"selected ... ";
      // +++++++++++++++++ OPTION ++++++++++++++++++++++++++++++
      // either keep CXCR5 de- and re-sensitisation mechanism ...
      // CC_list[li].set_CXCR5expression();
      // CC_list[li].resensitise4CXCL13(s);
      // ... or just switch chemotaxis to CXCL13 off:
      CC_list[li].responsive2signal[CXCL13] = false; 
      // +++++++++++++ end OPTION ++++++++++++++++++++++++++++++
      // Probability to differentiate to output cell and probability to recycle to CB.
      // If already preselected for differentiation to output:
      if (CC_list[li].selected4output) {
	if (CC_list[li].final_differentiation()) {
	  CC_differ2out(m[n], li, l, shape);
	  CCexists = false;
	} else {
	  // nomove=
	  CC_list[li].move(li, l, s, trackdata, time);
	}
      } else {
          if (((TFR_mode==13)
                  || (TFR_mode==19 && CC_list[li].CD138))  && CC_list[li].nTFRcontacts==0) {
              long tfr_i=-1;
//              bool canbindtfr=0, isSelf=CB_list[li].selfMutation;
//              short tobind=0; //=1->tfr =2->tfh
//              //CC scans neighb to find TFR
              if (CC_list[li].canbindTFR(time)) {
                  tfr_i = CC_list[li].get_contact(TFR,l);
              }
              if (((tfr_i!=-1)
                      && !(CC_list[li].same_TFR_as_before(tfr_i,TFR_list[l.cellknot[tfr_i].listi].id))
                      && TFR_list[l.cellknot[tfr_i].listi].state == TFRnormal)) {
                  if (tfr_i!=TFR_list[l.cellknot[tfr_i].listi].index) {cerr<<"???wrongtfrind";exit(1);}
                  //free TFR was found
                  //                        if (CC_list[li].ccInd_TFRbound==-1) {
                  CC_list[li].bind_TFR(time, TFR_list[l.cellknot[tfr_i].listi], l, shape);
                  //                            CC_list[li].tfr_fate(TFR_list[l.cellknot[tfr_i].listi]);
                  tfrbound=1;
                  CC_list[li].trackfate(time, true);
                  if (CC_list[li].trackit) {
                      double tmpi[l.dim];
                      l.get_koord(CC_list[li].index,tmpi);
                      trackdata.Write_movement(CC_list[li].trackno, CC, time, tmpi,
                                               CC_list[li].polarity, 1, 1,
                                               CC_list[li].fdc_clock, inTFRcontact);
                  }
              }
          }
        // If not yet differentiation type selected and no tfr to bind
          if (!tfrbound){
              err = CC_list[li].dif2OUTorCB(dt);
              if (err == 0) {
                  // nomove=
                  CC_list[li].move(li, l, s, trackdata, time);
              } else {
                  if (err == 2) {
                      // case differentiate to CB
                      CC_differ2CB(m[n], li, l, shape);
                      CCexists = false;
                  } else {
                      if (err == 1) {
                          // case differentiate to output without further delay
                          CC_differ2out(m[n], li, l, shape);
                          CCexists = false;
                      }
                      // if (err==3) do nothing now. Delay of differentiation is active.
                  }
              }
          }
      }
      // cout<<" selected: err="<<err<<", nomove="<<nomove;
      break;
    } // end of case selected.

    case apoptosis: {
      // cout<<" apoptosis: ";
      // Probability to leave the system through macrophage phagocytosis
      if (CC_list[li].macrophagocyte(shape) == 1) {
	// (switch to cell=empty).
	if (CC_list[li].trackit == true) {
	  double tmpi[l.dim];
	  l.get_koord(CC_list[li].index, tmpi);
	  trackdata.Stop_movement(CC_list[li].trackno, nocell, time, tmpi, CC_list[li].polarity);
	}
	if (track_mutations) {
	  // Write the deletion event to the brainbow class. Parameter are: 
	  // (double time, bool founder, bool birth, long mother_index, long ss_position)
	  trackmutations.write(time, false, false, 
			       CC_list[li].brainbow_index, CC_list[li].pos_ss);
	}
	// cout<<" deleted CC! contact_inhibited="<<CC_list[li].contact_inhibited<<"\n";
	del_CC(m[n], li, l);
	CCexists = false;
      } else {
	// If not yet eaten:
	// CC determined for apoptosis may continue to move
	bool keep_motility = CC_list[li].set_apoptotic_motility(s);
	if (keep_motility) {
	  // nomove=
	  CC_list[li].move(li, l, s, trackdata, time);
	}
      }
      break;
    } // end of case apoptosis.

    } // end switch
    // cout<<" CCexists="<<CCexists;
    // If cell was contact inhibited in its movement save for exchange:
    if (CCexists && CC_list[li].contact_inhibited) {
      // Note that deleted cells are excluded to fulfill the criterion
      // cout<<" contact inhibited ";
      // Reset contact_inhibited in order to ensure that it remains false in general
      CC_list[li].contact_inhibited = false;
      // this also implies volume==1 and moved==0.
      // Save this cell for later processing
      redo.add(CC_list[li].index);
    }

    // cout<<", nomove="<<nomove;
    // if the cell moved track it
    if (CCexists && (CC_list[li].writethis2track != trackini)) {
      // cout<<" in trackit ";
      double tmpi[l.dim];
      l.get_koord(CC_list[li].index, tmpi);
      trackdata.Write_movement(CC_list[li].trackno,
			       CC,
			       time,
			       tmpi,
			       CC_list[li].polarity,
			       CC_list[li].writethis2track);
      CC_list[li].writethis2track = trackini;
    }

    // cout<<" end this CC.\n\n";
  } // end for loop for each CC.
}


void cellman::free_tfr(long cc_ind, long tfr_ind) {
    CC_list[cc_ind].tfr_tmp_index=-1;
    TFR_list[tfr_ind].tmp_blocked=0;
}
// =======================================================
// TCs ===================================================
// =======================================================

long cellman::put_TC(long i, cellTC &newTC, space &l, AffinitySpace &shape) {
   /* This routine is to be understood as construction of a
    * first TC consisting of one fragment only, i.e. which is possibly
    * smaller than its asymptotic final state.
    * Returns the TC_list-index of the new cell.
    */
   if (i < 0) {
      cout << "Negative Index in put_TC!!!\n";
      exit(1);
   }
   if (l.cellknot[i].cell != nocell) {
      return -1;
   }
   // Add new cell on the list
   long li = TC_list.add(newTC);
   // Initialize it as new cell:
   TC_list[li].ini(i, li, time, l, shape);
   // That's it.
   return li;
}
void cellman::addTC(long nTC, space &l, double restrict, 
		    AffinitySpace &shape, ofstream &anafile) {
   anafile << "Put " << nTC << " T cells ...";
   long n = 0;
   long k[l.dim];
   while (n < nTC) {
      cellTC seederTC;   // Variable for T cells
      for (short i = 0; i < l.dim; i++) {
         k[i] = irandom(l.prodimvec[i]);
      }
      k[l.dim - 1] = irandom(int ((l.prodimvec[l.dim - 1]) / restrict));
      long getpos = l.Index(k);
      // cout<<"iTC="<<getpos<<"; ";
      // Assume TC to be specific for the Antigen!
      seederTC.pos_ss = shape.get_Antigen();
      if (put_TC(getpos, seederTC, l, shape) != -1) {
         anafile << " " << getpos << " ";
         ++n;
      }
   }
   anafile << " done.\n";
}
short int cellman::del_TC(long i, long li, space &l) {
   if (l.cellknot[i].cell != TC) {
      cout << "Error in del_TC: Zelltyp ist = " << l.cellknot[i].cell << " !!!\n";
      exit(1);
   }
   // remove on lattice
   l.clear_knot(i);
   // l.knot[i].cell=empty;
   // l.knot[i].listi=-1;
   // Delete cell on position li on the list
   TC_list.erase_jump(li);
   // Correct on the lattice for the shift of last cell on the list to position li
   // Don't if the deleted cell is the last cell on the list
   if (li < TC_list.benutzt()) {
      l.cellknot[TC_list[li].index].listi = li;
   }
   return 0;
}
void cellman::divide_TC(long i, long li, long j, space &l, AffinitySpace &shape) {
   /* Division of TC is only possible when the TC switched from the state TCnormal to TCdivide.
    * This implies that we can assume here, that this TC is bound to no BC.
    * In particular, the dividing cell has n_CC_nn==0 and this is copied to the daughter.
    * Also it can be assumed that space is available as the target position <j> was determined
    * before.
    */
   long newli;
   // write in log-file for proliferation
   if (outputfiles == 0) {
      prolog << time << "  " << pp[0] << "  " << pp[1] << "  " << pp[2] << "  " << pp[3]
             << "  " << l.Abstand(i, j)
             << "\n";
   }

   // Reset the times (this will be copied to the new cell as well)
   TC_list[li].reset_cycle_times();

   // Declare and define a new TC:
   cellTC newTC;
   // note that many parameters are set in put_TC() --> cellTC::ini() below, so define later

   // ##### Note that the old TC is not copied to the new one, 
   // which is needed when they become specific!

   // ##### first write then activate these
   // newTC.set_remaining_divisions();
   // TC_list[li].set_remaining_divisions();

   // newCB.responsive2signal[CXCL13]=false; // ### ??? no chemotaxis concept right now
   newTC.trackit = false;
   newTC.trackno = -1;

   // put the new TC on the lattice --> this calls cellTC::ini() which sets pos_ss and others
   newli = put_TC(j, newTC, l, shape);
   
   if (newli != -1) {

     //newTC.deliberate_memory(); // activate if TC are treated as frag_cells

     // copy the affinityspace position from the dividing cell
     TC_list[newli].pos_ss = TC_list[li].pos_ss;
     
     // Divide BrdU symmetrically on the daughters:
     TC_list[li].BrdU /= 2.;
     TC_list[newli].BrdU /= 2.;

     // Added by Marta Schips:
     TC_list[li].state = TCnormal;
     TC_list[newli].state = TCnormal; // when a copy-operator is activated this is necessary.
     
     /* Further options to handle the cells can be found in cellman::proliferate_CB(...),
      * e.g. mutation, cell tracking, etc.
      */
   }

   /*
   cerr << "benutzt=" << TC_list.benutzt() << "\n";
   cerr << "{l["<<i<<"]="<<l.cellknot[i].listi<<"="<<li<<" "<<l.cellknot[i].cell
	<< " / l["<<j<<"]="<<l.cellknot[j].listi<<"="<<newli<<" "<<l.cellknot[i].cell
	<< "}\n";
   char c;
   cin >> c;
   */
}
void cellman::calc_TC(long * m,
                      long mlang,
                      space &l,
                      sigs &s,
                      AffinitySpace &shape,
                      dynarray<long> &redo) {
  long li, n;
  // Randomize the sequence of points that are actualized
  for (n = 0; n < mlang; n++) {
    // Get index on the cell list
    li = l.cellknot[m[n]].listi;
    //cerr << "(" << n << "," << li << "):\n";

    // TC act in dependence on their state
    switch (TC_list[li].state) {
    case TCnormal:
        //MSchips
        ///WARNING: I think in each routine, move should be the last action
        /// otherwise 'li' is lost
//      TC_list[li].move(li, l, s, trackdata, time);
      if (cellTC::do_division) {
          // check whether a division-related action is necessary
          TC_list[li].ask_enter_cell_cycle();  // resets state and times if answer is yes
      }
      if (cell::north_weight==-5) {
          TC_list[li].set_CXCR5expression(); // time desensitisation
          TC_list[li].resensitise4CXCL13(s); // undercritical CXCL13 resensitisation
      }

      TC_list[li].move(li, l, s, trackdata, time);
      break;
      
    case TC_CCcontact:
      // cerr << "contact" << li << "; ";
      // No division or motility possible, for in contact to CC.
      // Evolve time variables:
      TC_list[li].evolve_polarisation_state(dt);
      // Set polarity for the direction of rescue signals.
      TC_list[li].set_polarity(l);
      break;
      
    case TCdivide:
      // cerr << "divide" << li << "; ";
      // move or divide
      bool divide_now = TC_list[li].progress_cell_cycle(dt);
      if (divide_now) {
	//cerr << "dividenow" << li << "(";
	/* in ask_mitosis, division is forced by the state==TCdivide 
	 * (<cellTC::proliferation> is ignored)
	 */
	long j = TC_list[li].ask_mitosis(pp, l);
	//cerr<<"j="<<j<<",m[n]="<<m[n]<<"); ";
	if (j >= 0) {
	  divide_TC(m[n], li, j, l, shape);
	}
	/* if this fails, in the next time step <divide_now==true> will happen again as
	 * clock is unchanged */
      } else {
	TC_list[li].move(li, l, s, trackdata, time);
      }
      break;
    }    // end switch
    
    //cerr<<"end of "<<li<<" switch, ";
    // if cell was contact inhibited in its movement try exchange:
    if (TC_list[li].contact_inhibited) {
      // reset contact_inhibited in order to ensure that it remains false in general
      TC_list[li].contact_inhibited = false;
      // this also implies volume==1 and moved==0.
      // save this cell for later processing
      redo.add(TC_list[li].index);
      // short moved=try2exchange_cells(m[n],li,TC,TC_list[li].polarity,l);
      // if (moved==1) nomove=0;
    }
    //cerr<<"end of "<<li<<" contactinhibited, ";
    
    if (TC_list[li].writethis2track != trackini) {
      double tmpi[l.dim];
      l.get_koord(TC_list[li].index, tmpi);
      trackdata.Write_movement(TC_list[li].trackno,
			       TC,
			       time,
			       tmpi,
			       TC_list[li].polarity,
			       TC_list[li].writethis2track);
      TC_list[li].writethis2track = trackini;
    }
    //cerr<<"end of "<<li<<" track.\n";
  }   // end for
  //cerr << "End of calc_TC!\n";
}

//MSchips
// =======================================================
// TFRcells===============================================
// =======================================================

long cellman::put_TFR(long i, cellTFR &newTFR, space &l, AffinitySpace &shape) {
   /* Returns the TFR_list-index of the new cell.
    */
   if (i < 0) {
      cout << "Negative Index in put_TFR!!!\n";
      exit(1);
   }
   if (l.cellknot[i].cell != nocell) {
      return -1;
   }

   // Add new cell on the list
   long li = TFR_list.add(newTFR);
   // Initialize it as new cell:
   TFR_list[li].ini(i, li, time, l, shape);
   // That's it.
   return li;
}

void cellman::addTFR(short tfr_mode, long nTFR, space &l, double restrict,
                     AffinitySpace &shape, ofstream &anafile) {
  anafile << "Put " << nTFR << " Tfr cells ...";
  bool inipos=0;
   long n = 0;
   long k[l.dim];
   while (n < nTFR) {
      cellTFR seederTFR;   // Variable for T cells
      for (short i = 0; i < l.dim; i++) {
         k[i] = irandom(l.prodimvec[i]);
      }
      if (cellTFR::TFR_north_weight==-1) {
          restrict=1.0;
          inipos=1;
      }
      int maxcoord=int ( (l.prodimvec[l.dim - 1]) / restrict);
      k[l.dim - 1] = irandom(maxcoord);
//      cerr<<k[l.dim-1]<<endl;
      long getpos = l.Index(k);
      // cout<<"iTC="<<getpos<<"; ";
      // Assume TC to be specific for the Antigen!
      seederTFR.pos_ss = shape.get_Antigen();
      if (!inipos || (l.no_border(getpos)==0) ) {
          if (put_TFR(getpos, seederTFR, l, shape) != -1) {
              anafile << " " << getpos << " ";
              ++n;
          }
      }
   }
   anafile << " done.\n";
}
short int cellman::del_TFR(long i, long li, space &l) {
   if (l.cellknot[i].cell != TFR) {
      cout << "Error in del_TFR: celltype is = " << l.cellknot[i].cell << " !!!\n";
      exit(1);
   }
   cerr<<"tfr deleted in state "<<TFR_list[li].state<<endl;
   // remove on lattice
   l.clear_knot(i);

   // Delete cell on position li on the list
   TFR_list.erase_jump(li);
   // Correct on the lattice for the shift of last cell on the list to position li
   // Don't if the deleted cell is the last cell on the list
   if (li < TFR_list.benutzt()) {
      l.cellknot[TFR_list[li].index].listi = li;
   }
   return 0;
}
//not used for now
void cellman::divide_TFR(long i, long li, long j, space &l, AffinitySpace &shape) {
   /* Division of TFR is only possible when the TFR switched from the state TFRnormal to TFRdivide.
    * This implies that we can assume here, that this TFR is bound to no BC.
    * In particular, the dividing cell has n_CC_nn==0 and this is copied to the daughter.
    * Also it can be assumed that space is available as the target position <j> was determined
    * before.
    */
   long newli;
   // write in log-file for proliferation
   if (outputfiles == 0) {
      prolog << time << "  " << pp[0] << "  " << pp[1] << "  " << pp[2] << "  " << pp[3]
             << "  " << l.Abstand(i, j)
             << "\n";
   }

   // Reset the times (this will be copied to the new cell as well)
   TFR_list[li].reset_cycle_times();

   // Declare and define a new TFR:
   cellTFR newTFR;
   // note that many parameters are set in put_TFR() --> cellTFR::ini() below, so define later

   // ##### Note that the old TFR is not copied to the new one,
   // which is needed when they become specific!

   // ##### first write then activate these
   // newTFR.set_remaining_divisions();
   // TFR_list[li].set_remaining_divisions();

   newTFR.trackit = false;
   newTFR.trackno = -1;

   // put the new TFR on the lattice --> this calls cellTFR::ini() which sets pos_ss and others
   newli = put_TFR(j, newTFR, l, shape);

   if (newli != -1) {

     //newTFR.deliberate_memory(); // activate if TFR are treated as frag_cells

     // copy the affinityspace position from the dividing cell
     TFR_list[newli].pos_ss = TFR_list[li].pos_ss;

     // Divide BrdU symmetrically on the daughters:
     TFR_list[li].BrdU /= 2.;
     TFR_list[newli].BrdU /= 2.;

     TFR_list[li].state = TFRnormal;
     TFR_list[newli].state = TFRnormal;

     /* Further options to handle the cells can be found in cellman::proliferate_CB(...),
      * e.g. mutation, cell tracking, etc.
      */
   }

   /*
   cerr << "benutzt=" << TFR_list.benutzt() << "\n";
   cerr << "{l["<<i<<"]="<<l.cellknot[i].listi<<"="<<li<<" "<<l.cellknot[i].cell
        << " / l["<<j<<"]="<<l.cellknot[j].listi<<"="<<newli<<" "<<l.cellknot[i].cell
        << "}\n";
   char c;
   cin >> c;
   */
}
//ms
void cellman::calc_TFR(long * m,
                      long mlang,
                      space &l,
                      sigs &s,
                      AffinitySpace &shape,
                      dynarray<long> &redo) {
//    cerr << "In calc_TFR...\n";
  long li, n;
  // Randomize the sequence of points that are actualized
  for (n = 0; n < mlang; n++) {
    // Get index on the cell list
    li = l.cellknot[m[n]].listi;

    // TFR act in dependence on their state
    switch (TFR_list[li].state) {
    case TFRnormal: {
        if (TFR_list[li].tmp_blocked) {
            cerr<<"blocked.not properly tested\n"; exit(1);
        }
        /** @brief: DTR experiment --
         * 3 dtr injections starting from tfr_depletion_time
         * a fraction of Tfr are depleted at each injection
         **/
        if ((tfr_depletion_time>0 && ((time>=tfr_depletion_time && time<(tfr_depletion_time+dt))
              || (time>=(tfr_depletion_time+24) && time<(tfr_depletion_time+24+dt))
              || (time>=(tfr_depletion_time+48) && time<(tfr_depletion_time+48+dt)))
             && (drandom()<(0.4)))
                || (TFR_list[li].to_delete)) {
            cerr<<time;
            del_TFR(m[n],li,l);
        } else {
            //division is not active
            if (cellTFR::do_division) {
                // check whether a division-related action is necessary
                TFR_list[li].ask_enter_cell_cycle();  // resets state and times if answer is yes
            }

            //        cerr<< "THIS IS NORTHTFR" <<cell::TFR_north_weight<<endl;
            if (cell::TFR_north_weight==-5) {
                TFR_list[li].set_CXCR5expression(); // time desensitisation
                TFR_list[li].resensitise4CXCL13(s); // undercritical CXCL13 resensitisation
            } else if (cell::TFR_north_weight==-4) {
                TFR_list[li].set_CXCR4expression(); // time desensitisation
                TFR_list[li].resensitise4CXCL12(s); // undercritical CXCL12 resensitisation
            }
            //50: tfr actively scans for CC
            // --> can bind free CCs or disrupt CC-Tfh bound for CC not polarised
            long cc_i=-1;
            if (TFR_mode>=50 && TFR_mode<100) {
                cc_i = TFR_list[li].get_contact(CC,l);
                cerr<<"not properly tested\n"; exit(1);
            }
            if (!TFR_list[li].tmp_blocked) {
                TFR_list[li].move(li, l, s, trackdata, time);
            }
        }
        break;
    }

    case TFR_CCcontact: {
        if ((tfr_depletion_time>0 && ((time>=tfr_depletion_time && time<(tfr_depletion_time+dt))
                                      || (time>=(tfr_depletion_time+24) && time<(tfr_depletion_time+24+dt))
                                      || (time>=(tfr_depletion_time+48) && time<(tfr_depletion_time+48+dt)))/*(time>tfr_depletion_time && time<(tfr_depletion_time+48))*/
             && (drandom()<(0.4)))) {
            TFR_list[li].to_delete=1;
        }
        short nbound = TFR_list[li].get_n_boundCC();
        if (nbound>1) {
            cerr<<"more than one cc bound at index "
               <<TFR_list[li].index
              <<" and id "<<TFR_list[li].id <<endl;
            exit(1);
        }
        break;
    }


    case TFRdivide: {
        // cerr << "divide" << li << "; ";
        // move or divide
        bool divide_now = TFR_list[li].progress_cell_cycle(dt);
        if (divide_now) {
            //cerr << "dividenow" << li << "(";
            /* in ask_mitosis, division is forced by the state==TFRdivide
         * (<cellTFR::proliferation> is ignored)
         */
            long j = TFR_list[li].ask_mitosis(pp, l);
            //cerr<<"j="<<j<<",m[n]="<<m[n]<<"); ";
            if (j >= 0) {
                divide_TFR(m[n], li, j, l, shape);
            }
            /* if this fails, in the next time step <divide_now==true> will happen again as
         * clock is unchanged */
        } else {
            TFR_list[li].move(li, l, s, trackdata, time);
        }
    }
        break;
    }    // end switch

//    cerr<<"end of "<<li<<" switch\n ";
    // if cell was contact inhibited in its movement try exchange:
    if (TFR_list[li].contact_inhibited) {
      // reset contact_inhibited in order to ensure that it remains false in general
      TFR_list[li].contact_inhibited = false;
      // this also implies volume==1 and moved==0.
      // save this cell for later processing
      redo.add(TFR_list[li].index);
      // short moved=try2exchange_cells(m[n],li,TFR,TFR_list[li].polarity,l);
      // if (moved==1) nomove=0;
    }
    //cerr<<"end of "<<li<<" contactinhibited, ";
    ///Warning: track is not happening for tfr (TODO!) // MMH2Marta: yes

    if (TFR_list[li].writethis2track != trackini) {
      double tmpi[l.dim];
      l.get_koord(TFR_list[li].index, tmpi);
      trackdata.Write_movement(TFR_list[li].trackno,
                               TFR,
                               time,
                               tmpi,
                               TFR_list[li].polarity,
                               TFR_list[li].writethis2track);
      TFR_list[li].writethis2track = trackini;
    }
    //cerr<<"end of "<<li<<" track.\n";
  }   // end for
//  cerr << "End of calc_TFR!\n";
}


// =======================================================
// FDCs ==================================================
// =======================================================

short int cellman::put_FDC(long i, space &l, AffinitySpace &shape, short transparence) {
   int r, n;
   long ind;
   long k[l.dim];
   long ktmp[l.dim];
   l.get_koord(i, k);
   for (n = 0; n < l.dim; n++) {
      ktmp[n] = k[n];
   }
   if ((l.cellknot[i].cell == nocell) && 
       ((l.cellknot[i].FDClisti == -1) || (transparence == 1))) {
      // erster Teil der Bedingung ueberfuessig ##
      // cout<<"start...\n";
      // save list-index on the lattice
      l.set_knot(i, FDC, FDC_list.benutzt());
      // cout<<"did set_knot\n";
      // l.knot[i].cell=FDC;
      l.cellknot[i].FDClisti = FDC_list.benutzt();
      // cout<<"saved FDClisti\n";
      // l.knot[i].listi=l.knot[i].FDClisti;
      cellFDC tmp;
      tmp.index = i;
      tmp.born_index = i;
      tmp.born_time = time;
      // tmp.center=i;
      tmp.state = soma;
      /*multiAg:
       * With multiple Ags per FDC and per FDC site, this following doesn't make sense.
       * Thus, monitoring sFDC on AffinitySpace is useful only in the 1 Ag case.
       * Take a random Antigen anyway in the 1 Ag case.
       * It will define pos_ss on all fragments of this
       * FDC.
       * But it will not be used further on, as everything is shifted to cellFDC::antigen_amount**
       */
      if (shape.get_n_Antigen() == 1) {
         tmp.pos_ss = shape.get_Antigen();
         // Add on shape space
         shape.add_cell(sFDC, tmp.pos_ss);
         // cout<<"added on shape\n";
      } else {
         tmp.pos_ss = -1;
         /*###multiAg: For multiple Ag, pos_ss should not be used anywhere!
          * Setting it to -1 is a test and will lead to program abortions in case
          * pos_ss would be used accidentally anywhere. Remove later.*/
      }
      // Write on cell_list
      FDC_list.add(tmp);
      // add this as first fragment:
      FDC_list[FDC_list.benutzt() - 1].fragments[FDC_list[FDC_list.benutzt() - 1].volume] = i;
      ++FDC_list[FDC_list.benutzt() - 1].volume;
      // total wird nicht erhoeht
      for (n = 0; n < l.dim; n++) {
         for (r = 1; r <= cellFDC::DendriteLength; r++) {
            ktmp[n] = k[n] + r;
            ind = l.Index(ktmp);
            if (ind != -1) {
               if (l.knot[ind].status != external) {
                  // ### Note: eventually already existing other cells are overwritten here!
                  // i.e. no check for CB CC or other cells is done!
                  // This is fine only as long as FDC are introduced first.

                  // if no FDC transparence save cell list index in listi also:
                  if (transparence == 0) {
                     if (l.cellknot[ind].cell != FDC) {
                        l.set_knot(ind, FDC, l.cellknot[i].FDClisti);
                        // l.cellknot[ind].listi=l.cellknot[i].FDClisti;
                        // (write only if no FDC was there before!)
                        // l.cellknot[ind].cell=FDC;
                     }
                  }
                  // nothing to do in the case of transparence==1
                  // else
                  // l.cellknot[ind].cell=nocell; // isn't that redundant? It should be nocell
                  // anyway
                  // It is even wrong, as a previously produced FDC soma at this point would
                  // be deleted by this. -> delete it!
                  // in any case of transparence save it in FDClisti:
                  if (l.cellknot[ind].FDClisti == -1) {
                     // condition: no FDC was there before
                     l.cellknot[ind].FDClisti = l.cellknot[i].FDClisti;
                  }
                  /*###multiAg -- I see a problem here when both FDCs are on the same grid-node
                   * but the grid-node only refers to one of the FDCs. Implications? One FDC
                   *fragment
                   * will never be used in terms of antigen? Yes, I think the new one which was not
                   *saved
                   * on the grid. No worries, just a bit of antigen will not be used. Other
                   * inconsistencies?
                   */
                  // and add new lattice index to fragments of the FDC
                  FDC_list[FDC_list.benutzt()
                           - 1].fragments[FDC_list[FDC_list.benutzt() - 1].volume] = ind;
                  ++FDC_list[FDC_list.benutzt() - 1].volume;
               }
            }
            ktmp[n] = k[n] - r;
            ind = l.Index(ktmp);
            if (ind != -1) {
               if (l.knot[ind].status != external) {
                  // comments as above
                  if (transparence == 0) {
                     if (l.cellknot[ind].cell != FDC) {
                        l.set_knot(ind, FDC, l.cellknot[i].FDClisti);
                        // l.cellknot[ind].listi=l.cellknot[i].FDClisti;
                        // l.cellknot[ind].cell=FDC;
                     }
                  }
                  // else l.cellknot[ind].cell=nocell; // as above
                  if (l.cellknot[ind].FDClisti == -1) {
                     // condition: no FDC was there before
                     l.cellknot[ind].FDClisti = l.cellknot[i].FDClisti;
                  }
                  FDC_list[FDC_list.benutzt()
                           - 1].fragments[FDC_list[FDC_list.benutzt() - 1].volume] = ind;
                  ++FDC_list[FDC_list.benutzt() - 1].volume;
               }
            }
            ktmp[n] = k[n];
         }
      }
      // Distribute the Antigen on the fragments:
      if (cellFDC::use_antigen == 1) {
         FDC_list[FDC_list.benutzt() - 1].add_antigen(shape.get_n_Antigen());
      }

      /*
       * cout<<"Ag=(";
       * for (int a=0; a<FDC_list[FDC_list.benutzt()-1].volume; a++)
       * cout<<FDC_list[FDC_list.benutzt()-1].antigen_amount[a]<<" ";
       * cout<<")\n";
       */
      return 0;
   } else {
      return 1;
   }
}
void cellman::calc_FDC(space &l, sigs &s, AntibodyDyn &Ab, AffinitySpace &AS) {
   long n;
   // Randomize the sequence of points that are actualized
   long max = FDC_list.benutzt();
   long m[max];
   random2_sequence(m, max);
   for (n = 0; n < max; n++) {
      m[n] = FDC_list[int (m[n])].index;
   }

   for (n = 0; n < max; n++) {
      // Production of diff2CC signal molecules in the center of FDC
      // This corresponds to the picture of a signal production at the FDC soma
      FDC_list[l.cellknot[m[n]].FDClisti].signal_production(l.knot[m[n]].guest_grid_pointer, s);
      // rate equation for production of immune complex:
      if (s.signal_use[antibody] == 1) {
         if (AS.get_n_Antigen() == 1) {
            FDC_list[l.cellknot[m[n]].FDClisti].mk_immune_complex(dt, s);
         } else {
            cerr << "Calling cellFDC::mk_immune_complex(const double&, sigs&):\n"
                 << "This version of mk_immune_complex only works with single Ags.\n"
                 << "Abort.\n";
            exit(1);
         }
      } else if (use_antibody_bins) {
         FDC_list[l.cellknot[m[n]].FDClisti].mk_immune_complex(dt, Ab, AS);
      }
      // Do not use the following as it is not yet functional!!!
      // FDC_list[l.cellknot[m[n]].FDClisti].mk_immune_complex(s);
   }
}
// =======================================================
// Stroma cells ==========================================
// =======================================================

void cellman::calc_stroma(sigs &s) {
   if (s.signal_use[CXCL12] == 1) {
      for (int n = 0; n < STROMA_list.benutzt(); n++) {
         s.signal_put(STROMA_list[n], CXCL12, p_mkCXCL12);
      }
   }
}
// =======================================================
// Output cells ==========================================
// =======================================================

/*
 * long int cellman::find_OUT(long i) {
 * long int found=-1;
 * long int n=0;
 * long int bis=OUT_list.benutzt();
 * while (n<bis && found==-1) {
 *    if (OUT_list[n].index==i) { found=n; }
 *    else { ++n; }
 * }
 * if (found==-1) { cout<<"Error: Index "<<i<<" not found in find_OUT!"; exit(1); }
 * return found;
 * }
 */

short int cellman::del_OUT(long i, long li, space &l) {
   // Erase from tracking
   if (OUT_list[li].trackit == true) {
      double tmpi[l.dim];
      l.get_koord(OUT_list[li].index, tmpi);
      trackdata.Stop_movement(OUT_list[li].trackno, out, time, tmpi, OUT_list[li].polarity);
   }
   if (track_mutations) {
      // write the deletion event to the brainbow class
      // parameters: (double time, bool founder, bool birth, long mother_index, long ss_position)
      trackmutations.write(time, false, false, OUT_list[li].brainbow_index, OUT_list[li].pos_ss);
   }
   // Write .cell .pos_ss .center .ccstate
   // remove on lattice
   l.clear_knot(i);
   //   l.knot[i].cell=empty;
   //   l.knot[i].listi=-1;
   // remove from list
   // Delete cell on position li on the list
   OUT_list.erase_jump(li);
   // Correct on the lattice for the shift of last cell on the list to position li
   // Don't if the deleted cell is the last cell on the list
   if (li < OUT_list.benutzt()) {
      l.cellknot[OUT_list[li].index].listi = li;
   }
   return 0;
}
void cellman::calc_OUT(space &l, sigs &s, AffinitySpace &shape, dynarray<long> &redo) {
   long li, n;
   // Randomize the sequence of points that are actualized
   long max = OUT_list.benutzt();
   // cout<<"OUT#="<<max<<"\n";
   long m[max];
   random2_sequence(m, max);
   // cout<<"back from randomize.\n";
   for (n = 0; n < max; n++) {
      m[n] = OUT_list[int (m[n])].index;
   }

   for (n = 0; n < max; n++) {
      if (l.cellknot[m[n]].cell != out) {
         cout << "In calc_OUT: l.cellknot[" << m[n] << "].cell=" << l.cellknot[m[n]].cell
              << "\n";
         cout << "Falscher Zelltyp!\n";
         exit(1);
      }
      li = l.cellknot[m[n]].listi; /// Philippe : this should be int (m[n]), no ?

      switch (OUT_list[li].state) {
      case Outfree:
          if (OUT_list[li].try_eject(m[n], l, shape) == 0) {
              long pZ=l.knot[OUT_list[li].index].x[2],
                      pX=l.knot[OUT_list[li].index].x[0],
                      pY=l.knot[OUT_list[li].index].x[1];
              ofstream exitpoint_xyz;
              exitpoint_xyz.open("exitpoint_xyz.out",ofstream::app);
              for (int i=0; i<l.prodimvec[0] ;i++) {
                  exitpoint_xyz << time
                       << " " << pX << " " << pY
                          << " " << pZ << "\n";
              }
              exitpoint_xyz.close();
              if (OUT_list[li].out_to_die) {
                  cerr<<"self out "<<OUT_list[li].selfMutation<<" with tfrcon "<<OUT_list[li].nTFRcontacts
                     <<" erased from ss\n";
                  //selfmut with tfrcont leaving the GC is dying --> will not be part of total produced output
                  //should be removed by the total number of cells
                  shape.rem_cell(sout, OUT_list[li].pos_ss);
                  shape.rem_cell(soutSelf, OUT_list[li].pos_ss);
              }
              del_OUT(m[n], li, l);
          } else {
              //MSchips
              //check the temporary tfr
              long tfr_i=-1;
              if ((TFR_mode==17) && OUT_list[li].canbindTFR(time)) {
                  tfr_i = OUT_list[li].get_contact(TFR,l);
              }
              if (tfr_i!=-1 && TFR_list[l.cellknot[tfr_i].listi].state == TFRnormal) {
                  OUT_list[li].bind_TFR(time, TFR_list[l.cellknot[tfr_i].listi], l, shape);
                  ++asc_got_contact;
              } else {
                  OUT_list[li].move(li, l, s, trackdata, time);  // result not saved as not needed in the
                                                                 // following
                  /// Philippe : suggestion : make move return a boolean, and then the field
                  // contact_inhibited is not required
                  // if cell was contact inhibited in its movement try exchange:
                  if (OUT_list[li].contact_inhibited) {
                     // reset contact_inhibited in order to ensure that it remains false in general
                     OUT_list[li].contact_inhibited = false;
                     // this also implies volume==1 and moved==0.
                     // save this cell for later processing
                     redo.add(OUT_list[li].index);
                     // short moved=try2exchange_cells(m[n],li,out,OUT_list[li].polarity,l);
                     // if (moved==1) nomove=0;
                  }
                  if (OUT_list[li].writethis2track != trackini) {
                     double tmpi[l.dim];
                     l.get_koord(OUT_list[li].index, tmpi);
                     trackdata.Write_movement(
                        OUT_list[li].trackno,
                        out,
                        time,
                        tmpi,
                        OUT_list[li].polarity,
                        OUT_list[li].writethis2track);
                     OUT_list[li].writethis2track = trackini;
                  }
                  OUT_list[li].signal_production(l.knot[m[n]].guest_grid_pointer, s);
              }
          }

          break;
      case OutTFRcontact:
          cerr<<"in cont "<<OUT_list[li].index<< " with "<<OUT_list[li].tfr_index<<endl;
          short err_tfr=0;
          err_tfr = OUT_list[li].got_tfr_interaction(time,dt,
                                                     TFR_list[l.cellknot[OUT_list[li].tfr_index].listi],l, shape);
          break;
      }
   }
   // cout<<"end calc_OUT.\n";
}
// =======================================================
// Beta-cells   ==========================================
// =======================================================

long cellman::put_BETA(long i, cellbeta &newBETA, space &l) {
   /* This routine is to be understood as construction of a
    * first betacell which is smaller
    * compared to its asymptotic final state. put_BETA cannot
    * be used to restore large betacells after their deletion.
    * This is done by the shift of single fragments following
    * Pott's ideas.
    * ==> no generalization to betacells with large radii is necessary!
    * Returns the BETA_list-index of the new cell.
    */
   if (i < 0) {
      cout << "Negative Index in put_BETA!!!\n";
      exit(1);
   }
   if (l.cellknot[i].cell != nocell) {
      return -1;
   }
   // Add new cell on the list
   newBETA.get_radius(l.dim);
   // cout<<"Add betacell at position i="<<i<<"; ";
   long li = BETA_list.add(newBETA);
   // cout<<"list-index is li="<<li<<"; ";
   // Initialize it as new cell:
   BETA_list[li].ini(i, li, time, l);
   // cout<<"initialized the new cell!\n";
   // That's it.
   return li;
}
short int cellman::del_BETA(long int i, long li, space &l) {
   if (l.cellknot[i].cell != BETA) {
      cout << "Error in del_BETA: Zelltyp ist = " << l.cellknot[i].cell << " !!!\n";
      exit(1);
   }
   // Delete all fragments from the lattice
   for (int a = 0; a < BETA_list[li].volume; a++) {
      l.clear_knot(BETA_list[li].fragments[a]);
      // l.cellknot[BETA_list[li].fragments[a]].cell=nocell;
      // l.cellknot[BETA_list[li].fragments[a]].listi=-1;
   }
   // Delete cell on position li on the list
   BETA_list.erase_jump(li);
   // Correct on the lattice for the shift of last cell on the list to position li
   // Don't if the deleted cell is the last cell on the list
   if (li < BETA_list.benutzt()) {
      for (int a = 0; a < BETA_list[li].volume; a++) {
         l.cellknot[BETA_list[li].fragments[a]].listi = li;
      }
   }
   return 0;
}
short cellman::macrophagocyte_BETA(long i, long li, space &l) {
   if (drandom() < p_macrophage_BETA) {
      del_BETA(i, li, l);
      return 1;
   }
   return 0;
}
short cellman::proliferate_BETA(long i, long li, space &l) {
   short err = 1;
   long newli;
   long j = BETA_list[li].ask_mitosis(pp, l);
   if (cellbeta::target_volume == 1) {
      if (j >= 0) {
         // write in log-file for proliferation
         if (outputfiles == 0) {
            prolog << time << "  " << pp[0] << "  " << pp[1] << "  " << pp[2] << "  "
                   << pp[3] << "  "
                   << l.Abstand(i, j) << "\n";
         }
         err = 0;
         // new cell has the same state as the dividing one (apart of state)
         cellbeta newBETA = BETA_list[li];
         // newBETA.state=cb_normal;
         newBETA.trackit = false;
         newBETA.trackno = -1;
         newli = put_BETA(j, newBETA, l);
         newBETA.deliberate_memory();
      }
   } else if (j == 0) {
      // case of more fragment object and do proliferate
      cellbeta newBETA = BETA_list[li];
      // reset state, do not track, and preserve all other variables
      // newBETA.state=cb_normal;
      newBETA.trackit = false;
      newBETA.trackno = -1;
      // add a new element on the list
      newli = BETA_list.add(newBETA);
      // Hier kann man BETA_list[newli].ini(...) nicht verwenden, da newBETA=BETA_list[li] wichtig!

      err = BETA_list[li].mitosis(i, li, newli, BETA_list[newli], l);
      if (err == 0) {
         // check_connection of the result
         if (checkit == 2) {
            BETA_list[li].check_connection(BETA, li, l);
            BETA_list[newli].check_connection(BETA, newli, l);
         }
         // end of preliminary check
         ++pp[0];
         if (outputfiles == 0) {
            prolog << time << "  " << pp[0] << "\n";
         }
      } else {
         exit(1);
      }
      newBETA.deliberate_memory();
   }
   return err;
}
void cellman::show_BETA() {
   // cerr<<"output at t="<<time<<" hours.\n";
   if (cellbeta::LOCAL_FILES) {
      for (long n = 0; n < BETA_list.benutzt(); n++) {
         BETA_list[n].show_all(3600. * time);
      }
   } else {
      for (long n = 0; n < BETA_list.benutzt(); n++) {
         BETA_list[n].show_all(3600. * time, beta_a, beta_b, beta_i, beta_r, beta_n, beta_g);
      }
   }
}
long cellman::get_global_index(long &li) {
   // li is the index of the node in the betacell list
   return li * cellbeta::N_equations;
   // returns the index in the global equation array y of the first equation to cell li
}
void cellman::get_gap_junction_static(double * y,
                                      long startat,
                                      const long &celli,
                                      dynarray<cellbeta> &list,
                                      long xypos,
                                      space &l) {
   /* This version of get_gap_junction is used in a global Runge-Kutta setting.
    * The quantities of all cells are passed in a single array y.
    * The specific array starts at position startat.
    * The cell under consideration has position xypos on the spatial lattice.
    * Note that the calling array is modified at gap_K,Na,Ca only.
    * These quantities are not used in other cells such that the Runge-Kutta
    * philosophy is guaranteed.
    * Note also that no control 1-cell-1-node is done as in the local non-static routine.
    */
   // cerr<<"in get_gap ... ";
   // initialise the variables for all ions with the own potential
   y[startat + cellbeta::gap_K] = 0.;    // K
   y[startat + cellbeta::gap_Na] = 0.;   // Na
   y[startat + cellbeta::gap_Ca] = 0.;   // Ca
   // go through all neighbours and sum the contributions:
   for (short i = 0; i < l.dim2; i++) {
      long j = l.knot[xypos].near_n[i];
      // cerr<<"index="<<index<<"  j="<<j<<"  ";
      if ((j != -1) && (l.cellknot[j].cell == BETA)) {
         long cellj = l.cellknot[j].listi;
         long startnn = get_global_index(cellj);
         // find the minimum of both gap-junction densities and use it for the current
         double rho_gap = list[celli].rho[betaWerte::gap];
         if (list[cellj].rho[betaWerte::gap] < rho_gap) {
            rho_gap = list[cellj].rho[betaWerte::gap];
         }
         // cout<<rho_gap<<"\n";
         // cerr<<"cellj="<<cellj<<"  ";
         y[startat + cellbeta::gap_K]
            += rho_gap * (y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                          - cellbeta::get_Nernst(y[startat + cellbeta::K],
                                                 y[startnn + cellbeta::K], cellbeta::z_K));
         y[startat + cellbeta::gap_Na]
            += rho_gap * (y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                          - cellbeta::get_Nernst(y[startat + cellbeta::Na],
                                                 y[startnn + cellbeta::Na], cellbeta::z_Na));
         y[startat + cellbeta::gap_Ca]
            += rho_gap * (y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                          - cellbeta::get_Nernst(y[startat + cellbeta::Ca],
                                                 y[startnn + cellbeta::Ca], cellbeta::z_Ca));
      }
   }
   // multiply by the factor which is identical for all ions:
   // double factor=p.gbar_gap*p.rho[betaWerte::gap]/l.dim2;
   //  double factor=cellbeta::p.gbar_gap*cellbeta::p.rho[betaWerte::gap]/6.;
   double factor = cellbeta::p.gbar_gap / 6.;
   y[startat + cellbeta::gap_K] *= factor;
   y[startat + cellbeta::gap_Na] *= factor;
   y[startat + cellbeta::gap_Ca] *= factor;
   /*  if (startat==258)
    * cerr<<"old: "<<y[startat+cellbeta::gap_K]<<"  "<<y[startat+cellbeta::gap_Na]<<"  "
    * <<y[startat+cellbeta::gap_Ca]<<"\n";*/
}
void cellman::get_gap_junction_dynamic(double * y,
                                       long startat,
                                       const long &celli,
                                       dynarray<cellbeta> &list,
                                       long xypos,
                                       space &l) {
   /* This version of get_gap_junction is used in a global Runge-Kutta setting
    * WITH DYNAMIC GAP-JUNCTION CONDUCTANCE.
    * The quantities of all cells are passed in a single array y.
    * The specific array starts at position startat.
    * Note that the calling array is modified at gap_K,Na,Ca only.
    * These quantities are not used in other cells such that the Runge-Kutta
    * philosophy is guaranteed.
    * Note also that no control 1-cell-1-node is done as in the local non-static routine.
    */
   // cerr<<"in get_gap ... ";
   // initialise the variables for all ions with the own potential
   y[startat + cellbeta::gap_K] = 0.;    // K
   y[startat + cellbeta::gap_Na] = 0.;   // Na
   y[startat + cellbeta::gap_Ca] = 0.;   // Ca
   // go through all neighbours and sum the contributions:
   for (short i = 0; i < l.dim2; i++) {
      long j = l.knot[xypos].near_n[i];
      // cerr<<"index="<<index<<"  j="<<j<<"  ";
      if ((j != -1) && (l.cellknot[j].cell == BETA)) {
         long cellj = l.cellknot[j].listi;
         // find the minimum of both gap-junction densities and use it for the current
         double rho_gap = list[celli].rho[betaWerte::gap];
         if (list[cellj].rho[betaWerte::gap] < rho_gap) {
            rho_gap = list[cellj].rho[betaWerte::gap];
         }
         // cout<<rho_gap<<"\n";
         // cout<<"rho_gap XXX= "<<rho_gap<<"\n";
         // cerr<<"li="<<li<<"  ";
         y[startat + cellbeta::gap_K] += (rho_gap * y[startat + cellbeta::gap_K + i + 1]);
         y[startat + cellbeta::gap_Na] += (rho_gap * y[startat + cellbeta::gap_Na + i + 1]);
         y[startat + cellbeta::gap_Ca] += (rho_gap * y[startat + cellbeta::gap_Ca + i + 1]);
      }
   }
   // multiply by the factor which is identical for all ions:
   // double factor=p.gbar_gap*p.rho[betaWerte::gap]/l.dim2;
   // double factor=cellbeta::p.gbar_gap*cellbeta::p.rho[betaWerte::gap]/6.;
   double factor = cellbeta::p.gbar_gap / 6.;
   y[startat + cellbeta::gap_K] *= factor;
   y[startat + cellbeta::gap_Na] *= factor;
   y[startat + cellbeta::gap_Ca] *= factor;
   /*
    * if (startat==430)
    * cerr<<"new: "<<y[startat+cellbeta::gap_K]<<"  "<<y[startat+cellbeta::gap_Na]<<"  "
    * <<y[startat+cellbeta::gap_Ca]<<"\n";
    */
}
void cellman::get_gap_rhs(double * y,
                          double * derivative,
                          const long &startat,
                          const long &celli,
                          dynarray<cellbeta> &list,
                          space &l) {
   // cerr<<" ... in cellman::get_gap_rhs(...) ... \n";
   double Vbar_gap;
   derivative[startat + cellbeta::gap_K] = 0.;
   derivative[startat + cellbeta::gap_Na] = 0.;
   derivative[startat + cellbeta::gap_Ca] = 0.;
   for (short i = 0; i < l.dim2; i++) {
      long j = l.knot[list[celli].index].near_n[i];
      // cerr<<"index="<<index<<"  j="<<j<<"  ";
      if ((j != -1) && (l.cellknot[j].cell == BETA)) {
         long li = l.cellknot[j].listi;
         long startnn = get_global_index(li);
         Vbar_gap = y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                    - cellbeta::get_Nernst(y[startat + cellbeta::K],
                                           y[startnn + cellbeta::K],
                                           cellbeta::z_K);
         /*
          * if (startat==430)// && y[startat+cellbeta::V]>-0.0699)
          * cerr<<"j="<<j<<", li="<<li<<", startat="<<startat<<", V_i="
          *    <<y[startat+cellbeta::V]<<", V_j="
          *    <<y[startnn+cellbeta::V]<<", K_i="
          *    <<y[startat+cellbeta::K]<<", K_j="
          *    <<y[startnn+cellbeta::K]<<", V_Nernst="
          *    <<cellbeta::get_Nernst(y[startat+cellbeta::K],y[startnn+cellbeta::K],cellbeta::z_K)
          *    <<", tau="<<cellbeta::p.tau_gap
          *    <<", y[430+..]="<<y[startat+cellbeta::gap_K+i+1]
          *    <<", Vbar_gap="<<Vbar_gap
          *    <<", dV_gap["<<i<<"]="
          *    <<(Vbar_gap-y[startat+cellbeta::gap_K+i+1])/cellbeta::p.tau_gap
          *    <<"\n";
          */
         derivative[startat + cellbeta::gap_K + i + 1]
            = (Vbar_gap - y[startat + cellbeta::gap_K + i + 1]) / cellbeta::p.tau_gap;
         // +1 is needed in the index because of gap_K pointing to the summed gap-junctions
         Vbar_gap = y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                    - cellbeta::get_Nernst(y[startat + cellbeta::Na],
                                           y[startnn + cellbeta::Na],
                                           cellbeta::z_Na);
         derivative[startat + cellbeta::gap_Na + i + 1]
            = (Vbar_gap - y[startat + cellbeta::gap_Na + i + 1]) / cellbeta::p.tau_gap;
         Vbar_gap = y[startat + cellbeta::V] - y[startnn + cellbeta::V]
                    - cellbeta::get_Nernst(y[startat + cellbeta::Ca],
                                           y[startnn + cellbeta::Ca],
                                           cellbeta::z_Ca);
         derivative[startat + cellbeta::gap_Ca + i + 1]
            = (Vbar_gap - y[startat + cellbeta::gap_Ca + i + 1]) / cellbeta::p.tau_gap;
      } else {
         derivative[startat + cellbeta::gap_K + i + 1] = 0.;
         derivative[startat + cellbeta::gap_Na + i + 1] = 0.;
         derivative[startat + cellbeta::gap_Ca + i + 1] = 0.;
      }
   }
}
void cellman::rhs(double t, double * y, double * derivative, dynarray<cellbeta> &list, space &l) {
   /* The y contains all cells in the list in the sequence of their position
    * in the list. rhs does not change the value of y but for the gap-junction
    * values. However, these are never used again but are calculated in every
    * call from the beginning.
    */
   // run through all cells and calculate cellbeta::rhs()
   double ycell[cellbeta::N_equations];
   double dydtcell[cellbeta::N_equations];
   for (long i = 0; i < list.benutzt(); i++) {
      // cerr<<"in rhs(t="<<t<<",...); work on cell li="<<i<<"of "<<list.benutzt()<<"; ";
      // get the position of the first equation to the cell under consideration:
      long startat = get_global_index(i);
      // cerr<<"startat="<<startat<<" y[glu]="<<y[startat+cellbeta::glu]<<".\n";;
      /* This could be used to get in between steps of glucose changes
       * // get the right glucose level
       * long k[l.dim];
       * l.get_koord(list[i].index,k);
       * y[startat+cellbeta::glu]=sigs::get_const_signal_value(t,cellbeta::p.glu_0,k);
       */
      // Calculate the gap-junction currents to cell i

      for (long j = 0; j < cellbeta::N_equations; j++) {
         derivative[startat + j] = 0.;
      }
      if (cellbeta::p.gap_dynamic) {
         get_gap_junction_dynamic(y, startat, i, list, list[i].index, l);
         /*
          * cerr<<"startat="<<startat
          * <<", y[V]="<<y[startat+cellbeta::V]
          * <<", kold="<<ktmp<<" vs knew="<<y[startat+cellbeta::gap_K]<<"; "
          * <<"naold="<<natmp<<" vs nanew="<<y[startat+cellbeta::gap_Na]<<"; "
          * <<"caold="<<catmp<<" vs canew="<<y[startat+cellbeta::gap_Ca]<<"\n";*/
         // get the rhs for the dynamic
         get_gap_rhs(y, derivative, startat, i, list, l);
      } else {
         get_gap_junction_static(y, startat, i, list, list[i].index, l);
      }
      /*      if (startat==258 && t>3.-1.e-10 && t<3+1.e-10) {
       *  cerr<<"new:\n";
       *  cerr<<"  y[V]="<<y[startat+cellbeta::V]
       *      <<", y[K]="<<y[startat+cellbeta::K]<<"\n  ";
       *  for (short ik=0; ik<l.dim2; ik++) {
       *    long j=l.knot[list[i].index].near_n[ik];
       *    //cerr<<"index="<<index<<"  j="<<j<<"  ";
       *    if (j!=-1 && l.cellknot[j].cell==BETA) {
       *      long li=l.cellknot[j].listi;
       *      long startnn=get_global_index(li);
       *      cerr<<"  ynn"<<ik<<"[V]="<<y[startnn+cellbeta::V]
       *          <<", ynn"<<ik<<"[K]="<<y[startnn+cellbeta::K]
       *          <<", ynn"<<ik<<"[gap_K]="<<y[startnn+cellbeta::gap_K+ik+1]
       *          <<", nernst="
       *
       *        <<cellbeta::get_Nernst(y[startat+cellbeta::K],y[startnn+cellbeta::K],cellbeta::z_K)
       *          <<"\n  ";
       *    }
       *  }
       *  cerr<<"  factor="<<cellbeta::p.gbar_gap*cellbeta::p.rho[betaWerte::gap]/6.
       *      <<"  result: "<<y[startat+cellbeta::gap_K]<<"  "<<y[startat+cellbeta::gap_Na]<<"  "
       *      <<y[startat+cellbeta::gap_Ca]<<"\n";
       * }
       */

      // save the values of the quantities for the cell-specific call
      for (int j = 0; j < cellbeta::N_equations; j++) {
         ycell[j] = y[startat + j];
      }
      // load the cell-specific protein densities to the cellbeta::p parameter variable
      cellbeta::p.rho = list[i].rho;
      // call the cell-specific ode-routine
      cellbeta::rhs(t, ycell, dydtcell);
      // save the cell-specific results in derivative for all cells
      // on the global array (but only up to index gap_K-1 -- the other ones
      // are calculated before in get_gap_rhs(...)!
      for (int j = 0; j < cellbeta::gap_K; j++) {
         derivative[startat + j] = dydtcell[j];
      }
      for (int j = cellbeta::g_fNa_V; j < cellbeta::N_equations; j++) {
         derivative[startat + j] = dydtcell[j];
      }
      // derivative[startat+cellbeta::g_fNa_V]=dydtcell[cellbeta::g_fNa_V];
      /*
       * if (startat==430)// && y[startat+cellbeta::V]>-0.0699)
       * cerr<<"derivative[K]="<<derivative[startat+cellbeta::gap_K_2]<<"\n";
       */
   }
}
void cellman::call_solver(dynarray<cellbeta> &list, space &l) {
   long Nglobal = cellbeta::N_equations * list.benutzt();
   // cout<<"benutzt="<<list.benutzt()<<"; N_equations="<<cellbeta::N_equations<<";
   // Nglobal="<<Nglobal<<"\n";
   double y[Nglobal];
   double y1[Nglobal];
   // save all betacells in a single array
   for (long i = 0; i < list.benutzt(); i++) {
      long itmp = get_global_index(i);
      for (int j = 0; j < cellbeta::N_equations; j++) {
         y[itmp + j] = list[i].y_n[j];
      }
   }
   double tsend = 3600. * time;
   double dtsend = 3600. * dt;
   // now the external time step is also used for betacell-electrophysiology!
   solver.ode_step(y, y1, list, l, Nglobal, tsend, dtsend, cellbeta::p.dy, prhs, method);
   // if (y[258+cellbeta::V]>-0.0699) cerr<<"new y0006[gap_K_2]="<<y[258+cellbeta::gap_K_2]<<"\n";
   // save result in the list-variables
   for (long i = 0; i < list.benutzt(); i++) {
      long itmp = get_global_index(i);
      for (int j = 0; j < cellbeta::N_equations; j++) {
         list[i].y_n[j] = y[itmp + j];
      }
   }
}
void cellman::calc_BETA(long * m, long mlang, space &l, sigs &s, dynarray<long> &redo) {
   /* m die random Reihenfolge der Indizes der betacells auf dem Gitter,
    * also von BETA_list[r].index-Werten */
   long li, n;

   // for (n=0; n<mlang; n++) if (l.cellknot[m[n]].listi==-1)
   // cout<<"at n="<<n<<",m[n]="<<m[n]<<",listi=-1;;";

   if (cellbeta::FULL_RUNGE == false) {
      // use the global solver for consistent Runge-Kutta
      //     call_solver(BETA_list,l);
      //   else // save the current values as old for use of explicit solution of gap-junctions
      for (n = 0; n < mlang; n++) {
         BETA_list[l.cellknot[m[n]].listi].synchronise();
      }
   }

   for (n = 0; n < mlang; n++) {
      // cerr<<"BETA"<<m[n]<<"(max"<<mlang<<")->";

      /*
       * cout<<"n="<<n<<",m[n]="<<m[n]<<",listi="<<l.cellknot[m[n]].listi<<": ";
       * for (int nn=n; nn<mlang; nn++) if (l.cellknot[m[nn]].listi==-1)
       * cout<<"at nn="<<nn<<",m[nn]="<<m[nn]<<",listi=-1;;";
       * cout<<" -- ";
       */

      // Get index on the cell list
      li = l.cellknot[m[n]].listi;
      // cout<<"n="<<n<<",m[n]="<<m[n]<<",li="<<li<<";;";

      // Set clock
      BETA_list[li].set_clock();

      // Set adhesion
      /* Later on, when the expression of adhesion molecules becomes dynamical,
       * the calculation will be initiated by
       * "BETA_list[li].set_adhesion();"
       * For the moment this is simply waste of time.
       ###
       */

      // do a lot of checks (delete after a while ###)
      if (checkit == 2) {
         if (li == -1) {
            cout << "call for listindex = -1 at lattice p= " << m[n] << "!\n";
            exit(1);
         }
         if (l.cellknot[m[n]].cell != BETA) {
            cout << "Fehler bei lattice point " << m[n] << " and cell#" << li
                 << " in calc_BETA: Zelltyp = " << l.cellknot[m[n]].cell << "\n";
            exit(1);
         }
         // Barycenter ueberpruefen
         BETA_list[li].check_barycenter(l);
      }

      // =====================================================
      // ============== Start of Actions: ==============================
      // =====================================================
      // Zustand bestimmen:
      BETA_list[li].get_new_state(m[n], s);
      /* // Activate the following if betacells shall consume glucose!
       * if (s.signal_use[glucose]==1)
       * BETA_list[li].use_nutrient(s,l.knot[BETA_list[li].index].guest_grid_pointer);
       */

      if (BETA_list[li].status == necrotic) {
         macrophagocyte_BETA(m[n], li, l);
      } else {
         if (cellbeta::FULL_RUNGE == false) {
            BETA_list[li].electrophysiology(time, dt, s, BETA_list, l);
         }

         // Mache Proliferation, Wachstum und Diffusion
         proliferate_BETA(m[n], li, l);

         // call movement:
         BETA_list[li].move(li, l, s, trackdata, time);

         // Treat the case of initiated move but suppression by lack of space:
         if (BETA_list[li].contact_inhibited) {
            // this also implies volume==1 and moved==0.
            // reset contact_inhibited in order to ensure that it remains false in general
            BETA_list[li].contact_inhibited = false;
            // save this cell for later processing
            redo.add(BETA_list[li].index);
            // moved=try2exchange_cells(m[n],li,BETA,BETA_list[li].polarity,l);
            // cout<<"  back in calc_BETA if contact_inhibited ...\n";
         }

         // write to movement file ...
         if (BETA_list[li].writethis2track != trackini) {
            // cout<<"Write this 2 track ...\n";

            // do not use index here!
            // Index is the lattice point nearest of the barycenter.
            // It is actualised after every call of frag_cell::get_barycenter(...)
            // in frag_cell::fragdiffuse(...), thus, automatically.

            double l2saxis = BETA_list[li].get_long2short_axis(li, BETA, l);
            double elongat = BETA_list[li].get_elongation(li, BETA, l);
            trackdata.Write_movement(BETA_list[li].trackno,
                                     BETA,
                                     time,
                                     BETA_list[li].barycenter,
                                     BETA_list[li].polarity,
                                     elongat,
                                     l2saxis,      // movement,
                                     BETA_list[li].writethis2track);
            BETA_list[li].writethis2track = trackini;
            BETA_list[li].n_immobile = 1;
            // cout<<"wrote to movement\n";
         }

         // call cell growth:
         BETA_list[li].grow(li, l);
      }    // not necrotic
   }
   if (cellbeta::FULL_RUNGE) {
      // use the global solver for consistent Runge-Kutta
      call_solver(BETA_list, l);
   }
}
// =======================================================
// Marker Ki67 ===========================================
// =======================================================

void cellman::inject_Ki67() {
   // Delete all Ki67 markers in every cell:
   for (long n = 0; n < CC_list.benutzt(); n++) {
      CC_list[n].Ki67 = 0;
   }
   for (long n = 0; n < OUT_list.benutzt(); n++) {
      OUT_list[n].Ki67 = 0;
   }
   // Mark all CBs now
   for (long n = 0; n < CB_list.benutzt(); n++) {
      CB_list[n].Ki67 = 1;
   }
   for (long n = 0; n < BETA_list.benutzt(); n++) {
      BETA_list[n].Ki67 = 1;
   }
   cout << "Injection of Ki67\n";
}
// =======================================================
// Data analysis to files ================================
// =======================================================

void cellman::zone_add(long i, long &n_CB, long &n_CB_nr, long * n_CC, long &n_all, space &l) {
   // Es werden nur CBs und CCs in all gezaehlt!
   long li = l.cellknot[i].listi;
   /*if (l.knot[i].cell!=ext && l.knot[i].cell!=FDC
    * && l.knot[i].cell!=empty ) { */
   if (((l.cellknot[i].cell == CB) && (CB_list[l.cellknot[i].listi].index == i))
       || ((l.cellknot[i].cell == CC) && (CC_list[l.cellknot[i].listi].index == i)
           && ((ignore_apoptotic_CC == false) || (CC_list[li].state != apoptosis)))) {
      // This condition restricts the count to the centers of the cells
      if (li == -1) {
         cout << "list-index=-1, but lattice point not empty.\n";
         exit(1);
      }
      ++n_all;
      if (l.cellknot[i].cell == CB) {
         ++n_CB;
         if (CB_list[li].n_recycling == 0) {
            ++n_CB_nr;
         }
      } else if (l.cellknot[i].cell == CC) {
         ++n_CC[CC_list[li].state];
         ++n_CC[apoptosis + 1];
      }
   }
}
void cellman::zone_put(double t, long n_CB, long n_CB_nr, long * n_CC, long n_all,
                       ofstream &datfile) {
   // Schreibt die Werte nach datfile
   double all = double (n_all);
   double n_allCC = double (n_CC[apoptosis + 1]);
   if (n_all > 0) {
      datfile << t << "    ";
      datfile << double (n_CB) / all << "   ";
      datfile << double (n_CB_nr) / all << "   ";
      datfile << double (n_allCC) / all << "   ";
      if (n_allCC > 0) {
         datfile << double (n_CB) / double (n_allCC) << "   " << double (n_CB_nr)
            / double (n_allCC) << "   ";
      } else {
         if (n_CB > 0) {
            datfile << "20000.0   ";      // d.h. infty
         } else {
            datfile << "1.0   ";
         }
         if (n_CB_nr > 0) {
            datfile << "20000.0   ";      // d.h. infty
         } else {
            datfile << "1.0   ";
         }
      }
   } else {
      datfile << t << "   0   0   0   ";
      if (n_allCC > 0) {
         datfile << double (n_CB) / double (n_allCC) << "   " << double (n_CB_nr)
            / double (n_allCC) << "   ";
      } else {
         if (n_CB > 0) {
            datfile << "20000.0   ";      // d.h. infty
         } else {
            datfile << "0.0   ";
         }
         if (n_CB_nr > 0) {
            datfile << "20000.0   ";      // d.h. infty
         } else {
            datfile << "0.0   ";
         }
      }
   }
   datfile << "  " << n_all << "   " << n_CB << "   " << n_CC[apoptosis + 1] << "   " << n_CB_nr;
   datfile << "\n";
}
void cellman::zone_files(suffix tnr, space &l) {
   // +++ OPTION: Make groups of lines/planes in order to reduce noise in particular at the border
   int group = 5;
   // end OPTION
   // Spatial l.Distribution of CBs and CCs relative to the FDCs
   long k[l.dim];
   long n_CB, n_CB_nr, n_all;
   long n_CC[apoptosis + 2];
   // Eroeffne das Zielfile:
   char datname[30] = "";
   strcat(datname, "zone");
   strcat(datname, tnr);
   strcat(datname, ".out");
   cout << "Write to " << datname << "\n";
   // Stream oeffnen
   ofstream ff(datname);
   ff.setf(ios::scientific, ios::floatfield);
   // Kopfzeile
   ff << "! time=" << time << "h\n";
   ff << "! z #CB/#all #norecyclingCB/#all #CC/#all #CB/#CC #nrCB/#CC\n";
   // Durchlaufe die Ebenen/Linen senkrecht zur Ausbildung von zones
   int n_group = 0;
   for (int n = 0; n < l.prodimvec[l.dim - 1]; n++) {
      k[l.dim - 1] = n;
      // start in the corner:
      for (short d = 0; d < l.dim - 1; d++) {
         k[d] = 0;
      }
      // Gesamtzahl der Punkte
      // long max=long(pow(double(l.prodim),l.dim-1));
      long max = 1;
      for (short d = 0; d < l.dim - 1; d++) {
         max *= l.prodimvec[d];
      }
      if (n_group == 0) {
         // Initialisiere die Zaehl-Variablen
         n_CB = 0;        // Total number of centroblasts
         n_CB_nr = 0;     // Number of not recycled centroblasts
         for (short d = unselected; d < apoptosis + 2; d++) {
            n_CC[d] = 0;      // number of centrocytes
         }
         n_all = 0;           // Total number of other in-GC points
      }
      ++n_group;
      // Durchlaufe alle Punkte der Linie/Ebene
      for (long m = 0; m < max; m++) {
         long i = l.Index(k);
         zone_add(i, n_CB, n_CB_nr, n_CC, n_all, l);
         // Erhoehe die entsprechende Koordinate
         if (k[0] == l.prodimvec[0] - 1) {
            k[0] = 0;
            if (l.dim == 3) {
               ++k[1];
            }
         } else {
            ++k[0];
         }
      }
      if (n_group == group) {
         // Write the values into the file
         zone_put((double (k[l.dim - 1]) - (double (group) - 1.) / 2.) * l.dx,
                  n_CB,
                  n_CB_nr,
                  n_CC,
                  n_all,
                  ff);
         n_group = 0;
      }
   }   // next line/plane
       // close file:
   ff.close();
}
short cellman::get_bin(double value) {
   // assumes the lowest value to be 0 and the maximum value of 1
   // rescale the value to the number of bins:
   // short bin=short(double(affinity_resolution)*value/max_value);
   short bin = short (double (affinity_resolution) * value);
   // the max_value is outside the range and is added to the highest bin:
   if (bin >= affinity_resolution) {
      bin = affinity_resolution - 1;
   }
   // this returns all values in each interval to the lower bin index
   // i.e. with max_value=1 and resolution=10 it is
   // [0,0.1[ -> 0; ... ; [0.9,1.0] -> 9
   return bin;
}
/*******************************************************************************************/
void cellman::save_foxo(ofstream& foxofile, short& foxores, int totalCC, int* nfoxo) {
  double foxofrac = 0.;
  foxofile << totalCC << "  ";
  for (short ifoxo = 0; ifoxo < foxores; ifoxo++) { 
    if (totalCC > 0) { foxofrac = double(nfoxo[ifoxo])/double(totalCC); }
    foxofile << foxofrac << "  "; 
  }
  foxofile << "  ";
}
void cellman::show_foxo(double time, long* nCC) {
  // +++++++++++ OPTION +++++++++++++++++++++++++++++
  short foxores = 5;
  // ++++++++++++++++++++++++++++++++++++++++++++++++
  int foxoCC[foxores], foxoCCnoapo[foxores], foxoCCunsel[foxores], foxoCCFDCsel[foxores];
  for (short ifoxo = 0; ifoxo < foxores; ++ifoxo) {
    foxoCC[ifoxo] = 0;
    foxoCCnoapo[ifoxo] = 0;
    foxoCCunsel[ifoxo] = 0;
    foxoCCFDCsel[ifoxo] = 0;
  }
  short foxores1 = foxores - 1;
  short foxobin = 0;
  for (long icc = 0; icc < CC_list.benutzt(); icc++) {
    foxobin = CC_list[icc].get_FoxO() * double(foxores1);
    if (foxobin >= foxores) { foxobin = foxores1; }
    ++foxoCC[foxobin];
    if (CC_list[icc].state != apoptosis) { ++foxoCCnoapo[foxobin]; }
    if (CC_list[icc].state == unselected || CC_list[icc].state == contact) {
      ++foxoCCunsel[foxobin];
    }
    if (CC_list[icc].state == FDCselected || CC_list[icc].state == TCcontact) {
      ++foxoCCFDCsel[foxobin];
    }    
  }
  ofstream foxofile;
  foxofile.open("foxoCC.out", ofstream::app);
  foxofile << time << "  "
	   << time/24. << "   ";
  save_foxo(foxofile, foxores, nCC[apoptosis + 1], foxoCC);
  save_foxo(foxofile, foxores, nCC[apoptosis + 1] - nCC[apoptosis], foxoCCnoapo);
  save_foxo(foxofile, foxores, nCC[unselected] + nCC[contact], foxoCCunsel);
  save_foxo(foxofile, foxores, nCC[FDCselected] + nCC[TCcontact], foxoCCFDCsel);
  foxofile << "\n";
  foxofile.close();
}
/*******************************************************************************************/
void cellman::show_ag_collected() {
   const int total_nCC = CC_list.benutzt();
   int * collected_ag_fdc = new int[total_nCC + 1];
   int * collected_ag_tc = new int[total_nCC + 1];
   int * collected_ag_tc_noapo = new int[total_nCC + 1];
   for (int aa = 0; aa < total_nCC; aa++) {
      collected_ag_fdc[aa] = 0;
      collected_ag_tc[aa] = 0;
      collected_ag_tc_noapo[aa] = 0;
   }
   int nCC_betweenFDCandTC = 0, nCC_include_afterTC = 0, nCC_include_afterTC_noapo = 0;

   // now load all nFDCcontacts from these cells (sorted according to state) into these new arrays
   short CCprogress;
   for (int li = 0; li < total_nCC; li++) {
      CCprogress = 0;
      if ((CC_list[li].state == FDCselected) || (CC_list[li].state == TCcontact)) {
         CCprogress = 1;
      }
      if (CC_list[li].state == selected) {
         CCprogress = 2;
      }
      if (CC_list[li].state == apoptosis) {
         CCprogress = 3;
      }
      if (CCprogress > 0) {
         // FDCselected, TCcontact, selected, apoptosis
         // if CCprogress is 1, 2 or 3
         collected_ag_tc[nCC_include_afterTC] = CC_list[li].nFDCcontacts;
         ++nCC_include_afterTC;
         if (CCprogress < 3) {
            // FDCselected, TCcontact, selected
            collected_ag_tc_noapo[nCC_include_afterTC_noapo] = CC_list[li].nFDCcontacts;
            ++nCC_include_afterTC_noapo;
         }
         if (CCprogress == 1) {
            // FDCselected, TCcontact
            collected_ag_fdc[nCC_betweenFDCandTC] = CC_list[li].nFDCcontacts;
            ++nCC_betweenFDCandTC;
         }
      }
   }

   // Find the average values and sd in these new arrays:
   double averageFDC = 0, sdFDC = 0;
   get_average_and_sd(collected_ag_fdc, nCC_betweenFDCandTC, averageFDC, sdFDC);
   double averageTC = 0, sdTC = 0;
   get_average_and_sd(collected_ag_tc, nCC_include_afterTC, averageTC, sdTC);
   double averageTCnoapo = 0, sdTCnoapo = 0;
   get_average_and_sd(collected_ag_tc_noapo, nCC_include_afterTC_noapo, 
		      averageTCnoapo, sdTCnoapo);

   // write this to files:
   ag_collected << time << "    " 
		<< nCC_betweenFDCandTC << "  " << averageFDC << "  " << sdFDC
                << "     "
                << nCC_include_afterTC_noapo << "  " << averageTCnoapo << "  " << sdTCnoapo
                << "     "
                << nCC_include_afterTC << "  " << averageTC << "  " << sdTC << "\n";
}
/*******************************************************************************************/
long cellman::getNofVoidFDCsites() {
  double voidsites = 0;
  for (int i = 0; i < FDC_list.benutzt(); i++) {
    voidsites += FDC_list[i].get_voidsites();
  }
  return voidsites;
}
long cellman::getFDCnetwork_volume() {
  long no = 0;
  for (int i = 0; i < FDC_list.benutzt(); i++) {
    no += FDC_list[i].volume;
  }
  return no;
}
void cellman::showFDCfailure(double& t) {
  ofstream failure;
  failure.open("FDCagfailure.out", ofstream::app);
  failure << t << "  " << t/24. << "    "
	  << TestedFDCsiteVoidOfAg << "  ";
  double ratio = 0.;
  if (CC_list.benutzt() > 0) {
    ratio = double(TestedFDCsiteVoidOfAg)/double(CC_list.benutzt());
  }
  failure << CC_list.benutzt() << "  " << ratio << "    ";
  // Reset counter for next time period:
  TestedFDCsiteVoidOfAg = 0;
  long totalsites = getFDCnetwork_volume();
  long emptysites = getNofVoidFDCsites();
  failure << emptysites << "  " << totalsites << "   " 
	  << double(emptysites)/double(totalsites) << endl;
  failure.close();
}

/*******************************************************************************************/
void cellman::show_K_signal_intensity_pMHC() {
  double mean = 0.;
  double sd = 0.; 
  double v = 0.;
  long n = TC_list.benutzt();
  for (long i = 0; i < n; i++) { mean += TC_list[i].K_signal_intensity_pMHC; }
  if (n > 0) { mean /= double(n); }
  for (long i = 0; i < n; i++) { 
    v = TC_list[i].K_signal_intensity_pMHC - mean;
    sd += v * v;
  }
  if ( n > 1 ) { 
    sd /= double(n - 1); 
    sd = sqrt(sd);
  } else { sd = 0.; }
  ofstream K("TFHsignalintensityK.out", ofstream::app);
  K << time << "   " << mean << "   " << sd << "   " << n << "\n";
  K.close();
}
/*******************************************************************************************/
void cellman::show_cycle_phases(long zone_separator) {
  /* writes into a file that collects the time course of BC numbers in different
   * phases of the cell cycle. 
   * This is also differentiated to the spatially defined zones.
   * zone_separator is the first index attributed to the DZ, used to attribute to zones.
   */
  long phases[cb_statenumber];
  long phasesDZ[cb_statenumber];
  long phasesLZ[cb_statenumber];
  long totalBC = 0, totalBCDZ = 0, totalBCLZ = 0;
  for (short j = 0; j < cb_statenumber; j++) 
    { phases[j] = 0; phasesDZ[j] = 0; phasesLZ[j] = 0; }
  for (long i = 0; i < CB_list.benutzt(); i++ ) {
    ++phases[CB_list[i].state];
    ++totalBC;
    if (CB_list[i].index < zone_separator) 
      { ++phasesLZ[CB_list[i].state]; ++totalBCLZ; }
    else
      { ++phasesDZ[CB_list[i].state]; ++totalBCDZ; }
  }
  // Write this to the output file ...
  ofstream cycle_phases("cellcycle_phases.out", ofstream::app);
  cycle_phases << time << "    "
	       << totalBC << "    ";
  //short first = int(cb_G1), last = int(cb_M) + 1; for (short c = first; c < last; c++) {
  for (short c = 0; c < short(cb_statenumber); c++) {
    cycle_phases << phases[c] << "  ";
  }
  /* Now add from CC those who are in waiting mode.
   * The thought is that these cells actually are already in cell cycle
   * but are still sensitive to CXCL13, thus not passing to the DZ yet.
   * In order to reproduce the number of cells in different cell cycle
   * phases, ignoring CC which are already dividing would shift the
   * numbers towards later phases, which is a mistake. Therefore,
   * selected CC will be attributed to a cell cycle phase depending
   * on the time they are already in the selected state.
   * TODO: It might be better to let selected CC differentiate to
   * CBs right away, which then could enter cell cycle normally,
   * and to delay only the up-regulation of CXCR4. Then, this artificial
   * correction could be deleted again.
   */
  if ( cellCB::shiftCCdelay2CBcycle() ) {
    for (long i = 0; i < CC_list.benutzt(); i++ ) {
      if ( CC_list[i].state == selected ) {
	centroblasts phase = cellCB::get_virtual_cell_cycle_phase(CC_list[i].selected_clock);
	++phases[phase];
	++totalBC;
	if (CC_list[i].index < zone_separator) 
	  { ++phasesLZ[phase]; ++totalBCLZ; }
	else
	  { ++phasesDZ[phase]; ++totalBCDZ; }
      }
    }
  }
  cycle_phases << "  " << totalBC << "    ";
  //short first = int(cb_G1), last = int(cb_M) + 1; for (short c = first; c < last; c++) {
  for (short c = 0; c < short(cb_statenumber); c++) {
    cycle_phases << phases[c] << "  ";
  }
  cycle_phases << "\n";
  cycle_phases.close();
  /* Write the data specific for the zone to output file. Here, only
   * those eventually combined with the CCs in a started cell cycle are shown.
   */
  ofstream cycle_phases_zones("cellcycle_phases_zones.out", ofstream::app);
  cycle_phases_zones << time << "  " << totalBCDZ << "    ";
  for (short c = 0; c < short(cb_statenumber); c++) {
    cycle_phases_zones << phasesDZ[c] << "  ";
  }
  cycle_phases_zones << "  " << totalBCLZ << "    ";
  for (short c = 0; c < short(cb_statenumber); c++) {
    cycle_phases_zones << phasesLZ[c] << "  ";
  }
  cycle_phases_zones << "  ";
  for (short c = 0; c < short(cb_statenumber); c++) {
    double fraction = phasesDZ[c];
    if (fraction > 0) { fraction /= (phasesDZ[c]+phasesLZ[c]); }
    cycle_phases_zones << fraction << "  ";
  }
  cycle_phases_zones << "\n";
  cycle_phases_zones.close();
}
/*********************************************************************************************/
void write_MHCdeficient(double& time, 
			long& mhcdef_CB, long& mhcdef_CC, 
			long& dec205_CB, long& dec205_CC) {
  long dec205_BC = dec205_CB + dec205_CC;
  long mhcdef_BC = mhcdef_CB + mhcdef_CC;
  ofstream MHCdeficient_out;
  MHCdeficient_out.open("MHCdeficient.out", ofstream::app);
  MHCdeficient_out << time << "   "
		   << mhcdef_BC << "   "
		   << mhcdef_CB << "   "
		   << mhcdef_CC << "     ";
  if (dec205_BC > 0 ) { MHCdeficient_out << double(mhcdef_BC)/double(dec205_BC) << "   "; }
  else { MHCdeficient_out << "0    "; }
  if (dec205_CB > 0 ) { MHCdeficient_out << double(mhcdef_CB)/double(dec205_CB) << "   "; }
  else { MHCdeficient_out << "0    "; }
  if (dec205_CC > 0 ) { MHCdeficient_out << double(mhcdef_CC)/double(dec205_CC) << "   "; }
  else { MHCdeficient_out << "0    "; }
  MHCdeficient_out << "\n";
  MHCdeficient_out.close();
}
/*********************************************************************************************/
void cellman::mk_cell_sum(space &l, sigs &s, AntibodyDyn &Ab, AffinitySpace &shape) {
   // kann nicht nach lattice, sollte in uebergeordneter Zellklasse bleiben ?
    //MSchips
    ///Warning: check if this function needs to be modified
    /// for the compatibility with TFR model (TODO!)
  long li;
  long nd_CB, nd_CB_nr, nd_all;
  long nd_CB_nr_ha, nd_CB_r_ha;
  long nl_CB, nl_CB_nr, nl_all;
  long nl_CB_nr_ha, nl_CB_r_ha;
  long nd_out, nl_out;
  long dec_pos_dz_CB = 0, dec_pos_lz_CB = 0, dec_neg_dz_CB = 0, dec_neg_lz_CB = 0;
  long dec_pos_dz_CC = 0, dec_pos_lz_CC = 0, dec_neg_dz_CC = 0, dec_neg_lz_CC = 0;
  long dec_pos_dz_CC_sel = 0, dec_pos_lz_CC_sel = 0, 
    dec_neg_dz_CC_sel = 0, dec_neg_lz_CC_sel = 0;
  long dec_pos_dz_CC_apo = 0, dec_pos_lz_CC_apo = 0, 
    dec_neg_dz_CC_apo = 0, dec_neg_lz_CC_apo = 0;
  long dec_pos_dz_out = 0, dec_pos_lz_out = 0, dec_neg_dz_out = 0, dec_neg_lz_out = 0;
  long mhcdef_CB = 0, mhcdef_CC = 0;
  long n_diff_CB = 0;
  long notempty = 0;
  long nouts = 0;
  long nl_CC[apoptosis + 2];
  long nd_CC[apoptosis + 2];
  long cb_aff_lz[affinity_resolution], 
    cc_aff_lz[affinity_resolution],
    out_aff_lz[affinity_resolution],
    cbdec_aff_lz[affinity_resolution], 
    ccdec_aff_lz[affinity_resolution],
    outdec_aff_lz[affinity_resolution],
    cb_aff_dz[affinity_resolution], 
    cc_aff_dz[affinity_resolution],
    out_aff_dz[affinity_resolution],
    cbdec_aff_dz[affinity_resolution], 
    ccdec_aff_dz[affinity_resolution],
    outdec_aff_dz[affinity_resolution];
  for (short i = 0; i < affinity_resolution; i++) {
    cb_aff_lz[i] = 0;
    cc_aff_lz[i] = 0;
    out_aff_lz[i] = 0;
    cb_aff_dz[i] = 0;
    cc_aff_dz[i] = 0;
    out_aff_dz[i] = 0;
    cbdec_aff_lz[i] = 0;
    ccdec_aff_lz[i] = 0;
    outdec_aff_lz[i] = 0;
    cbdec_aff_dz[i] = 0;
    ccdec_aff_dz[i] = 0;
    outdec_aff_dz[i] = 0;
  }
  long cb_ig_dz[nIg_classes], 
    cc_ig_dz[nIg_classes], 
    out_ig_dz[nIg_classes],
    cb_ig_lz[nIg_classes],
    cc_ig_lz[nIg_classes], 
    out_ig_lz[nIg_classes], 
    cc_ig_apo_dz[nIg_classes],
    cc_ig_apo_lz[nIg_classes];
  double cb_ig_aff[nIg_classes], cc_ig_aff[nIg_classes], out_ig_aff[nIg_classes];
  for (int i = 0; i < nIg_classes; i++) {
      cb_ig_dz[i] = 0;
      cc_ig_dz[i] = 0;
      out_ig_dz[i] = 0;
      cb_ig_lz[i] = 0;
      cc_ig_lz[i] = 0;
      out_ig_lz[i] = 0;
      cb_ig_aff[i] = 0.0;
      cc_ig_aff[i] = 0.0;
      out_ig_aff[i] = 0.0;
      cc_ig_apo_dz[i] = 0;
      cc_ig_apo_lz[i] = 0;
  }

   // Light zone:
   nl_CB = 0;         // Total number of centroblasts
   nl_CB_nr = 0;      // Number of not recycled centroblasts
   nl_CB_r_ha = 0;    // Number of high affinity recycled centroblasts
   nl_CB_nr_ha = 0;   // Number of high affinity non-recycled centroblasts
   nl_out = 0;
   for (short n = unselected; n < apoptosis + 2; n++) {
      nl_CC[n] = 0;    // number of centrocytes
   }
   nl_all = 0;         // Total number of CB or CC in-GC points
   // Determine the index that separates the DZ from the LZ (assumed 50:50 DZ and LZ):
   long max = l.zone_separator;
   for (long n = 0; n < max; n++) {
      li = l.cellknot[n].listi;
      zone_add(n, nl_CB, nl_CB_nr, nl_CC, nl_all, l);
      if ((l.cellknot[n].cell != nocell) && (l.knot[n].status != external)) {
         ++notempty;
      }
      if (l.cellknot[n].cell == out) {
         ++nouts;
         ++nl_out;
         ++out_aff_lz[get_bin(shape.best_affinity_norm(OUT_list[li].pos_ss))];
         ++out_ig_lz[OUT_list[li].IgX.Ig_class];
         out_ig_aff[OUT_list[li].IgX.Ig_class] += shape.best_affinity_norm(OUT_list[li].pos_ss);
         if (OUT_list[li].DEC205) {
            ++dec_pos_lz_out;
            ++outdec_aff_lz[get_bin(shape.best_affinity_norm(OUT_list[li].pos_ss))];
         } else {
            ++dec_neg_lz_out;
         }
      }
      if ((l.cellknot[n].cell == CB) && (CB_list[li].index == n)) {
         // count only if center of cell
         if ((CB_list[li].state == cb_differentiate)
             || (CB_list[li].state == cb_stop_dividing)) {
            ++n_diff_CB;
         }
         if (shape.best_affinity_norm(CB_list[li].pos_ss) > 0.3) {
            if (CB_list[li].n_recycling == 0) {
               ++nl_CB_nr_ha;
            } else {
               ++nl_CB_r_ha;
            }
         }
         ++cb_aff_lz[get_bin(shape.best_affinity_norm(CB_list[li].pos_ss))];
         ++cb_ig_lz[CB_list[li].IgX.Ig_class];
         cb_ig_aff[CB_list[li].IgX.Ig_class] += shape.best_affinity_norm(CB_list[li].pos_ss);
         if (CB_list[li].DEC205) {
            ++dec_pos_lz_CB;
            ++cbdec_aff_lz[get_bin(shape.best_affinity_norm(CB_list[li].pos_ss))];
         } else {
            ++dec_neg_lz_CB;
         }
	 if (CB_list[li].MHCdeficient) { ++mhcdef_CB; }
      }
      if ((l.cellknot[n].cell == CC)
          && ((ignore_apoptotic_CC == false) || (CC_list[li].state != apoptosis))) {
         ++cc_aff_lz[get_bin(shape.best_affinity_norm(CC_list[li].pos_ss))];
         ++cc_ig_lz[CC_list[li].IgX.Ig_class];
         cc_ig_aff[CC_list[li].IgX.Ig_class] += shape.best_affinity_norm(CC_list[li].pos_ss);
         if (CC_list[li].state == apoptosis) {
            ++cc_ig_apo_lz[CC_list[li].IgX.Ig_class];
         }
         if (CC_list[li].DEC205) {
            ++dec_pos_lz_CC;
            if (CC_list[li].state == selected) {
               ++dec_pos_lz_CC_sel;
            }
            if (CC_list[li].state == apoptosis) {
               ++dec_pos_lz_CC_apo;
            }
            ++ccdec_aff_lz[get_bin(shape.best_affinity_norm(CC_list[li].pos_ss))];
         } else {
            ++dec_neg_lz_CC;
            if (CC_list[li].state == selected) {
               ++dec_neg_lz_CC_sel;
            }
            if (CC_list[li].state == apoptosis) {
               ++dec_neg_lz_CC_apo;
            }
         }
	 if (CC_list[li].MHCdeficient) { ++mhcdef_CC; }
      }
   }
   if (show_mode != islet) {
      zone_put(time, nl_CB, nl_CB_nr, nl_CC, nl_all, light_zone);
   }

   // Dark zone
   nd_CB = 0;         // Total number of centroblasts
   nd_CB_nr = 0;      // Number of not recycled centroblasts
   nd_CB_r_ha = 0;    // Number of high affinity recycled centroblasts
   nd_CB_nr_ha = 0;   // Number of high affinity non-recycled centroblasts
   nd_out = 0;
   for (short n = unselected; n < apoptosis + 2; n++) {
      nd_CC[n] = 0;    // number of centrocytes
   }
   nd_all = 0;         // Total number of CB or CC in-GC points
   // der Rest gibt die dark zone:
   for (long n = max; n < l.pointnumber; n++) {
      li = l.cellknot[n].listi;
      zone_add(n, nd_CB, nd_CB_nr, nd_CC, nd_all, l);
      if ((l.cellknot[n].cell != nocell) && (l.knot[n].status != external)) {
         ++notempty;
      }
      if (l.cellknot[n].cell == out) {
         ++nouts;
         ++nd_out;
         ++out_aff_dz[get_bin(shape.best_affinity_norm(OUT_list[li].pos_ss))];
         ++out_ig_dz[OUT_list[li].IgX.Ig_class];
         out_ig_aff[OUT_list[li].IgX.Ig_class] += shape.best_affinity_norm(OUT_list[li].pos_ss);
         if (OUT_list[li].DEC205) {
            ++dec_pos_dz_out;
            ++outdec_aff_dz[get_bin(shape.best_affinity_norm(OUT_list[li].pos_ss))];
         } else {
            ++dec_neg_dz_out;
         }
      }
      if ((l.cellknot[n].cell == CB) && (CB_list[li].index == n)) {
         if ((CB_list[li].state == cb_differentiate)
             || (CB_list[li].state == cb_stop_dividing)) {
            ++n_diff_CB;
         }
         if (shape.best_affinity_norm(CB_list[li].pos_ss) > 0.3) {
            if (CB_list[li].n_recycling == 0) {
               ++nd_CB_nr_ha;
            } else {
               ++nd_CB_r_ha;
            }
         }
         ++cb_aff_dz[get_bin(shape.best_affinity_norm(CB_list[li].pos_ss))];
         ++cb_ig_dz[CB_list[li].IgX.Ig_class];
         cb_ig_aff[CB_list[li].IgX.Ig_class] += shape.best_affinity_norm(CB_list[li].pos_ss);
         if (CB_list[li].DEC205) {
            ++dec_pos_dz_CB;
            ++cbdec_aff_dz[get_bin(shape.best_affinity_norm(CB_list[li].pos_ss))];
         } else {
            ++dec_neg_dz_CB;
         }
	 if (CB_list[li].MHCdeficient) { ++mhcdef_CB; }
      }
      if ((l.cellknot[n].cell == CC)
          && ((ignore_apoptotic_CC == false) || (CC_list[li].state != apoptosis))) {
         ++cc_aff_dz[get_bin(shape.best_affinity_norm(CC_list[li].pos_ss))];
         ++cc_ig_dz[CC_list[li].IgX.Ig_class];
         cc_ig_aff[CC_list[li].IgX.Ig_class] += shape.best_affinity_norm(CC_list[li].pos_ss);
         if (CC_list[li].state == apoptosis) {
            ++cc_ig_apo_dz[CC_list[li].IgX.Ig_class];
         }
         if (CC_list[li].DEC205) {
            ++dec_pos_dz_CC;
            if (CC_list[li].state == selected) {
               ++dec_pos_dz_CC_sel;
            }
            if (CC_list[li].state == apoptosis) {
               ++dec_pos_dz_CC_apo;
            }
            ++ccdec_aff_dz[get_bin(shape.best_affinity_norm(CC_list[li].pos_ss))];
         } else {
            ++dec_neg_dz_CC;
            if (CC_list[li].state == selected) {
               ++dec_neg_dz_CC_sel;
            }
            if (CC_list[li].state == apoptosis) {
               ++dec_neg_dz_CC_apo;
            }
         }
	 if (CC_list[li].MHCdeficient) { ++mhcdef_CC; }
      }
   }

   if (show_mode != islet) {
      zone_put(time, nd_CB, nd_CB_nr, nd_CC, nd_all, dark_zone);
   }

   // Merke die Zeit des Endes der dark zone
   if ((time < t_dark_end) && (nd_CB < nd_CC[apoptosis + 1])) {
      t_dark_end = time;
   }

   // write the dec205-file and the affinity files
   // !!! NOTE: out-data are read above but not part of the read-out into stream xdec205
   // !!! add them below if output cells in the GC (not external) are to be included
   if (show_mode != islet) {
     // MHC-deficiency (mhcdef_{CB,CC} were counted above)
     long tmp1 = dec_pos_dz_CB + dec_pos_lz_CB; // all DEC205+ CBs
     long tmp2 = dec_pos_dz_CC + dec_pos_lz_CC; // all DEC205+ CCs
     write_MHCdeficient(time, mhcdef_CB, mhcdef_CC, tmp1, tmp2);
      // dec205-file:
      tmp1 = dec_pos_dz_CB + dec_pos_lz_CB + dec_pos_dz_CC + dec_pos_lz_CC;
      tmp2 = dec_neg_dz_CB + dec_neg_lz_CB + dec_neg_dz_CC + dec_neg_lz_CC;
      xdec205 << time << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      tmp1 = dec_pos_dz_CB + dec_pos_dz_CC;
      tmp2 = dec_neg_dz_CB + dec_neg_dz_CC;
      xdec205 << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      tmp1 = dec_pos_lz_CB + dec_pos_lz_CC;
      tmp2 = dec_neg_lz_CB + dec_neg_lz_CC;
      xdec205 << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      tmp1 = dec_pos_dz_CB + dec_pos_lz_CB;
      tmp2 = dec_neg_dz_CB + dec_neg_lz_CB;
      xdec205 << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      tmp1 = dec_pos_dz_CC + dec_pos_lz_CC;
      tmp2 = dec_neg_dz_CC + dec_neg_lz_CC;
      xdec205 << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      tmp1 = dec_pos_dz_CC_sel + dec_pos_lz_CC_sel;
      tmp2 = dec_neg_dz_CC_sel + dec_neg_lz_CC_sel;
      xdec205 << "   " << tmp1 << "   " << tmp2;
      if (tmp1 + tmp2 > 0) {
         xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
      } else {
         xdec205 << "   0";
      }

      if (ignore_apoptotic_CC) {
         xdec205 << "   0   0   0   0   0";
      } else {
         tmp1 = dec_pos_dz_CC_apo + dec_pos_lz_CC_apo;
         tmp2 = dec_neg_dz_CC_apo + dec_neg_lz_CC_apo;
         xdec205 << "   " << tmp1 << "   " << tmp2;
         if (tmp1 + tmp2 > 0) {
            xdec205 << "   " << double (tmp1) / double (tmp1 + tmp2);
         } else {
            xdec205 << "   0";
         }
         xdec205 << "   " << dec_pos_dz_CC_apo << "    " << dec_neg_dz_CC_apo;
      }

      xdec205 << "\n";

      // save all the LZ affinity files:
      aff_cb_lz << time;
      aff_cc_lz << time;
      aff_out_lz << time;
      aff_cb_dz << time;
      aff_cc_dz << time;
      aff_out_dz << time;
      for (short i = 0; i < affinity_resolution; i++) {
         aff_cb_lz << "   " << cb_aff_lz[i] << "  " << cbdec_aff_lz[i];
         aff_cc_lz << "   " << cc_aff_lz[i] << "  " << ccdec_aff_lz[i];
         aff_out_lz << "   " << out_aff_lz[i];
         aff_cb_dz << "   " << cb_aff_dz[i] << "  " << cbdec_aff_dz[i];
         aff_cc_dz << "   " << cc_aff_dz[i] << "  " << ccdec_aff_dz[i];
         aff_out_dz << "   " << out_aff_dz[i];
      }
      aff_cb_lz << "\n";
      aff_cc_lz << "\n";
      aff_out_lz << "\n";
      aff_cb_dz << "\n";
      aff_cc_dz << "\n";
      aff_out_dz << "\n";

      // ===============================
      // Save immunoglobulin class data:
      // ===============================
      // xsumBCig<<"! time : CB+CC : DZ : LZ : DZ/LZ: IgM : %IgM : IgM DZ : IgM LZ : IgM DZ/LZ : "
      //        <<"IgG... : IgE... : IgA...\n";
      // ... first the time:
      xsumBCig << time << "    ";
      xsumBCOig << time << "    ";
      xapoig << time << "    ";
      xapoigdz << time << "    ";
      xapoiglz << time << "    ";
      xaffig << time << "    ";
      // ... get the total populations:
      double sum_cb_ig_dz = 0, sum_cb_ig_lz = 0, sum_cb_ig_aff = 0, sum_cc_ig_dz = 0,
             sum_cc_ig_lz = 0,
             sum_cc_ig_aff = 0, sum_out_ig_dz = 0, sum_out_ig_lz = 0, sum_out_ig_aff = 0,
             sum_cc_ig_apo_dz = 0,
             sum_cc_ig_apo_lz = 0, sum_cc_ig_apo_tot = 0;
      for (int k = 0; k < nIg_classes; k++) {
         sum_cb_ig_dz += cb_ig_dz[k];
         sum_cb_ig_lz += cb_ig_lz[k];
         sum_cb_ig_aff += cb_ig_aff[k];
         sum_cc_ig_dz += cc_ig_dz[k];
         sum_cc_ig_lz += cc_ig_lz[k];
         sum_cc_ig_aff += cc_ig_aff[k];
         sum_out_ig_dz += out_ig_dz[k];
         sum_out_ig_lz += out_ig_lz[k];
         sum_out_ig_aff += out_ig_aff[k];
         sum_cc_ig_apo_dz += cc_ig_apo_dz[k];
         sum_cc_ig_apo_lz += cc_ig_apo_lz[k];
         sum_cc_ig_apo_tot += cc_ig_apo_dz[k] + cc_ig_apo_lz[k];
      }
      // ... write these to file:
      xsumBCig << sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz << "  "
               << sum_cb_ig_dz + sum_cc_ig_dz
               << "  " << sum_cb_ig_lz + sum_cc_ig_lz << "  ";
      xapoigdz << sum_cc_ig_apo_dz << "  ";
      xapoiglz << sum_cc_ig_apo_lz << "  ";
      xapoig << sum_cc_ig_apo_tot << "  ";
      if (sum_cb_ig_lz + sum_cc_ig_lz > 0) {
         xsumBCig << (sum_cb_ig_dz + sum_cc_ig_dz) / (sum_cb_ig_lz + sum_cc_ig_lz) << "    ";
      } else {
         xsumBCig << "0    ";
      }
      xsumBCOig << sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
         + sum_out_ig_lz << "  "
                << sum_cb_ig_dz + sum_cc_ig_dz + sum_out_ig_dz << "  " << sum_cb_ig_lz
         + sum_cc_ig_lz + sum_out_ig_lz
                << "  ";
      if (sum_cb_ig_lz + sum_cc_ig_lz + sum_out_ig_lz > 0) {
         xsumBCOig << (sum_cb_ig_dz + sum_cc_ig_dz + sum_out_ig_dz)
            / (sum_cb_ig_lz + sum_cc_ig_lz + sum_out_ig_lz)
                   << "    ";
      } else {
         xsumBCOig << "0    ";
      }
      // ... apoptosis:
      if (sum_cb_ig_dz + sum_cc_ig_dz + sum_out_ig_dz > 0) {
         xapoigdz << double (sum_cc_ig_apo_dz)
            / double (sum_cb_ig_dz + sum_cc_ig_dz + sum_out_ig_dz);
      } else {
         xapoigdz << "0";
      }
      xapoigdz << "    ";
      if (sum_cb_ig_lz + sum_cc_ig_lz + sum_out_ig_lz > 0) {
         xapoiglz << double (sum_cc_ig_apo_lz)
            / double (sum_cb_ig_lz + sum_cc_ig_lz + sum_out_ig_lz);
      } else {
         xapoiglz << "0";
      }
      xapoiglz << "    ";
      if (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
          + sum_out_ig_lz > 0) {
         xapoig << double (sum_cc_ig_apo_tot)
            / double (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz
                      + sum_out_ig_dz + sum_out_ig_lz);
      } else {
         xapoig << "0";
      }
      xapoig << "    ";
      // ... affinity
      if (sum_cb_ig_dz + sum_cb_ig_lz > 0) {
         xaffig << double (sum_cb_ig_aff) / double (sum_cb_ig_dz + sum_cb_ig_lz);
      } else {
         xaffig << "0";
      }
      xaffig << "  ";
      if (sum_cc_ig_dz + sum_cc_ig_lz > 0) {
         xaffig << double (sum_cc_ig_aff) / double (sum_cc_ig_dz + sum_cc_ig_lz);
      } else {
         xaffig << "0";
      }
      xaffig << "  ";
      if (sum_out_ig_dz + sum_out_ig_lz > 0) {
         xaffig << double (sum_out_ig_aff) / double (sum_out_ig_dz + sum_out_ig_lz);
      } else {
         xaffig << "0";
      }
      xaffig << "  ";
      if (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz > 0) {
         xaffig << double (sum_cb_ig_aff + sum_cc_ig_aff)
            / double (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz);
      } else {
         xaffig << "0";
      }
      xaffig << "  ";
      if (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
          + sum_out_ig_lz > 0) {
         xaffig << double (sum_cb_ig_aff + sum_cc_ig_aff + sum_out_ig_aff)
            / double (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
                      + sum_out_ig_lz);
      } else {
         xaffig << "0";
      }
      xaffig << "    ";

      // ... write total, % of total, DZ, LZ, DZ/LZ for each class
      for (int k = 0; k < nIg_classes; k++) {
         // ... first for CB+CC
         xsumBCig << cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] << "  ";
         if (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz > 0) {
            xsumBCig << double (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k])
               / (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz);
         } else {
            xsumBCig << "0";
         }
         xsumBCig << "  ";
         xsumBCig << cb_ig_dz[k] + cc_ig_dz[k] << "  " << cb_ig_lz[k] + cc_ig_lz[k] << "  ";
         if (cb_ig_lz[k] + cc_ig_lz[k] > 0) {
            xsumBCig << double (cb_ig_dz[k] + cc_ig_dz[k])
               / double (cb_ig_lz[k] + cc_ig_lz[k]) << "    ";
         } else {
            xsumBCig << "0    ";
         }
         // ... then for CB+CC+OUT
         xsumBCOig << cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] + out_ig_dz[k]
            + out_ig_lz[k] << "  ";
         if (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
             + sum_out_ig_lz > 0) {
            xsumBCOig
               << double (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] + out_ig_dz[k]
                          + out_ig_lz[k])
               / (sum_cb_ig_dz + sum_cb_ig_lz + sum_cc_ig_dz + sum_cc_ig_lz + sum_out_ig_dz
                  + sum_out_ig_lz);
         } else {
            xsumBCOig << "0";
         }
         xsumBCOig << "  ";
         xsumBCOig << cb_ig_dz[k] + cc_ig_dz[k] + out_ig_dz[k] << "  " << cb_ig_lz[k]
            + cc_ig_lz[k] + out_ig_lz[k]
                   << "  ";
         if (cb_ig_lz[k] + cc_ig_lz[k] + out_ig_lz[k] > 0) {
            xsumBCOig << double (cb_ig_dz[k] + cc_ig_dz[k] + out_ig_dz[k])
               / double (cb_ig_lz[k] + cc_ig_lz[k] + out_ig_lz[k]) << "    ";
         } else {
            xsumBCOig << "0    ";
         }

         // ... write the number of apoptotic BC on the lattice of this class for each zone and
         // total:
         xapoigdz << cc_ig_apo_dz[k] << "  ";
         // ... and the frequency with respect to the total population of this class
         if (cb_ig_dz[k] + cc_ig_dz[k] + out_ig_dz[k] > 0) {
            xapoigdz << double (cc_ig_apo_dz[k])
               / double (cb_ig_dz[k] + cc_ig_dz[k] + out_ig_dz[k]);
         } else {
            xapoigdz << "0";
         }
         xapoigdz << "    ";
         xapoiglz << cc_ig_apo_lz[k] << "  ";
         // ... and the frequency with respect to the total population of this class
         if (cb_ig_lz[k] + cc_ig_lz[k] + out_ig_lz[k] > 0) {
            xapoiglz << double (cc_ig_apo_lz[k])
               / double (cb_ig_lz[k] + cc_ig_lz[k] + out_ig_lz[k]);
         } else {
            xapoiglz << "0";
         }
         xapoiglz << "    ";
         xapoig << cc_ig_apo_dz[k] + cc_ig_apo_lz[k] << "  ";
         // ... and the frequency with respect to the total population of this class
         if (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] + out_ig_dz[k]
             + out_ig_lz[k] > 0) {
            xapoig << double (cc_ig_apo_dz[k] + cc_ig_apo_lz[k])
               / double (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k]
                         + cc_ig_lz[k]
                         + out_ig_dz[k] + out_ig_lz[k]);
         } else {
            xapoig << "0";
         }
         xapoig << "    ";

         // ... and the affinities
         if (cb_ig_dz[k] + cb_ig_lz[k] > 0) {
            xaffig << double (cb_ig_aff[k]) / double (cb_ig_dz[k] + cb_ig_lz[k]);
         } else {
            xaffig << "0";
         }
         xaffig << "  ";
         if (cc_ig_dz[k] + cc_ig_lz[k] > 0) {
            xaffig << double (cc_ig_aff[k]) / double (cc_ig_dz[k] + cc_ig_lz[k]);
         } else {
            xaffig << "0";
         }
         xaffig << "  ";
         if (out_ig_dz[k] + out_ig_lz[k] > 0) {
            xaffig << double (out_ig_aff[k]) / double (out_ig_dz[k] + out_ig_lz[k]);
         } else {
            xaffig << "0";
         }
         xaffig << "  ";
         if (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] > 0) {
            xaffig << double (cb_ig_aff[k] + cc_ig_aff[k])
               / double (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k]);
         } else {
            xaffig << "0";
         }
         xaffig << "  ";
         if (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] + out_ig_dz[k]
             + out_ig_lz[k] > 0) {
            xaffig << double (cb_ig_aff[k] + cc_ig_aff[k] + out_ig_aff[k])
               / double (cb_ig_dz[k] + cb_ig_lz[k] + cc_ig_dz[k] + cc_ig_lz[k] + out_ig_dz[k]
                         + out_ig_lz[k]);
         } else {
            xaffig << "0";
         }
         xaffig << "    ";
      }
      // ... finally end line:
      xsumBCig << "\n";
      xsumBCOig << "\n";
      xapoig << "\n";
      xapoiglz << "\n";
      xapoigdz << "\n";
      xaffig << "\n";
   }


   //MSc
   ofstream newly_gen;
   newly_gen.open("newly_gen.out",ofstream::app);
   newly_gen << time << " " << newly_generated_sCBs << " " << newly_generated_redeemed << "\n";
   newly_gen.close();
   newly_generated_redeemed=0;
   newly_generated_sCBs=0;

   // Zone independent read-outs
   nd_CB += nl_CB;
   nd_CB_nr += nl_CB_nr;
   nd_CB_nr_ha += nl_CB_nr_ha;
   nd_CB_r_ha += nl_CB_r_ha;
   for (short n = unselected; n < apoptosis + 2; n++) {
      nd_CC[n] += nl_CC[n];
   }

   // Calculate the high affinity fractions:
   double CBha_nr;
   if (nd_CB_nr > 0) {
      CBha_nr = double (nd_CB_nr_ha) / double (nd_CB_nr);
   } else {
      CBha_nr = 0.;
   }
   double CBha_r;
   if (nd_CB - nd_CB_nr > 0) {
      CBha_r = double (nd_CB_r_ha) / double (nd_CB - nd_CB_nr);
   } else {
      CBha_r = 0.;
   }
   double CBha;
   if (nd_CB > 0) {
      CBha = double (nd_CB_r_ha + nd_CB_nr_ha) / double (nd_CB);
   } else {
      CBha = 0.;
   }
   cbhighaff << time << "   " << CBha << "   " << CBha_nr << "   " << CBha_r << "\n";
   // Docuementation during run:
   // cout<<time<<"  CBha="<<CBha<<"   "<<CBha_nr<<"   "<<CBha_r<<"\n";

   double fraction_of_recycled_all = 1.;
   long tmp
      = long (n_recycling_events + shape.get_sum_cell(sallapoptosis) + shape.get_sum_cell(sout));
   if (tmp > 0.) {
      fraction_of_recycled_all = double (n_recycling_events) / double (tmp);
   }
   double fraction_of_recycled = 1.;
   long tmpold
      = long (n_recycling_events_last + shape.get_oldsum_cell(sallapoptosis)
              + shape.get_oldsum_cell(sout));
   if (tmp - tmpold > 0) {
      fraction_of_recycled = double (n_recycling_events - n_recycling_events_last)
                             / double (tmp - tmpold);
   }

   xsums[CB] << time << "   " 
	     << nd_CB << "   " 
	     << nd_CB_nr << "   " 
	     << n_recycling_events << "   "
             << n_recycling_events - n_recycling_events_last << "   "
             << fraction_of_recycled_all << "   "
             << fraction_of_recycled << "\n";
   n_recycling_events_last = n_recycling_events;

   xsums[CC] << time << "   " 
	     << nd_CC[apoptosis + 1] << "   " // this is the total CC
	     << nd_CC[unselected] << "   "
             << nd_CC[contact] << "   "
             << nd_CC[FDCselected] << "   " 
	     << nd_CC[TCcontact] << "   " 
	     << nd_CC[selected] << "   "
             << nd_CC[apoptosis] << "\n";
   
   show_foxo(time, nd_CC);

   long total_out = nd_out + nl_out;
   xsums[out] << time << "   " << total_out << "   " << nd_out << "   " << nl_out << "    ";
   // note that nd_CB/CC contain total numbers
   if (nd_CB + nd_CC[apoptosis + 1] + total_out > 0) {
      xsums[out] << double (total_out) / double (nd_CB + nd_CC[apoptosis + 1] + total_out);
   } else {
      xsums[out] << "0";
   }
   xsums[out] << "\n";

   CB_end = nd_CB;
   // Mittelung der CBs ueber die letzten zwei Tage:
   if (time > last2days + 0.5) {
      CB_average += nd_CB;
      CB_variance += (nd_CB * nd_CB);
   }

   // Write the number of occupied lattice points,
   //       the number ignoring FDCs
   //       the occupied fraction of the total GC volume,
   //       the occupied fraction ignoring FDCs
   //       the occupied fraction ignoring FDCs and OUTs
   BC_integral += (nd_CB + nd_CC[apoptosis + 1]);
   cout << "t=" << time << " --> #BCs=" << nd_CB + nd_CC[apoptosis + 1] << "\n";
   xvolume << time << "   " << notempty << "   " << notempty - FDCpointnumber << "   " << nd_CB
      + nd_CC[apoptosis + 1]
           << "   " << double (notempty) / double (GCpointnumber) << "   "
           << double (notempty - FDCpointnumber) / (double (GCpointnumber - FDCpointnumber))
           << "   "
           << double (notempty - FDCpointnumber - nouts)
      / (double (GCpointnumber - FDCpointnumber)) << "   "
           << nd_CB + nd_CC[apoptosis + 1] + total_out << "   " << BC_integral << "\n";
   cum_vol << time << "   " << nd_CB + nd_CC[apoptosis + 1] << "\n";

   // Zone-independent read-outs, read off the full lattice
   // TC:
   long ntc = 0;
   const unsigned short len = 7;   // the number of columns is always l.dim2+1 for dim=3
   long tc_nbc[len];
   for (short i = 0; i < len; i++) {
      tc_nbc[i] = 0;
   }
   for (long n = 0; n < l.pointnumber; n++) {
      li = l.cellknot[n].listi;
      if ((l.knot[n].status != external) && (l.cellknot[n].cell == TC)) {
         // At first add this to the number of TC
         ++ntc;
         // Look through the neighbours and count all those that are CC in TC contact:
         short nbc_contact = 0;
         // Pruefe Zelltyp der Nachbarn
         for (short i = 0; i < l.dim2; i++) {
            long j = l.knot[n].near_n[i];
            if ((j != -1) && (l.cellknot[j].cell == CC)
                && (CC_list[l.cellknot[j].listi].state == TCcontact)) {
               ++nbc_contact;
            }
         }
         // add this one to the array of TC with contact to a number of BC:
         ++tc_nbc[nbc_contact];
      }
   }
   // write the result to the xtc-file:
   xtc << time << "   " << ntc;
   for (short i = 0; i < len; i++) {
      xtc << "   " << tc_nbc[i];
   }
   xtc << "\n";

   if (checkit == 2) {
      long cbsum = CB_list.benutzt();
      if (cbsum != nd_CB) {
         cout << "Error: total number of CBs inconsistent in mk_cell_sum!\n";
         exit(1);
      }
      long cbnrsum = 0;
      long cbrsum = 0;
      for (int j = 0; j < cbsum; j++) {
         if (CB_list[j].n_recycling > 0) {
            ++cbrsum;
         } else {
            ++cbnrsum;
         }
      }
      if (cbnrsum != nd_CB_nr) {
         cout << "Error: number of non-recycled CBs inconsistent in mk_cell_sum!\n";
         exit(1);
      }
   }

   // Mutations-Statistik am Output:
   mutation_out << time << "   " << n_outs << "   " << n_muts << "   " << n_recmuts;
   if (n_outs > 0) {
      mutation_out << "   " << double (n_muts) / double (n_outs) << "   " 
		   << double (n_recmuts) / double (n_outs);
   } else {
      mutation_out << "   0   0";
   }
   mutation_out << "   " << n_outs_with_recycling;
   if (n_outs_with_recycling > 0) {
      mutation_out << "   " << double (n_recmuts) / double (n_outs_with_recycling);
   } else {
      mutation_out << "    0";
   }
   mutation_out << "\n";
   // I did this strange sequence in order to keep compatibility with previous versions

   // Signal production
   s.write_siglog(time);

   if (cellFDC::use_antigen == 1) {
     // Antigen presentation:
     double free_ag = 0.;
     double ics = 0.;
     int n_ags = shape.get_n_Antigen();
     double free_ags[n_ags];
     double ics_ags[n_ags];
     for (int a = 0; a < n_ags; a++) {
       free_ags[a] = 0;
       ics_ags[a] = 0;
     }
     for (int f = 0; f < FDC_list.benutzt(); f++) {
       // total free antigen of any kind:
       free_ag += FDC_list[f].get_total_antigen();
       //free_ag += FDC_list[f].get_total_antigen(n_ags);
       // free antigen for each antigen separately
       for (int a = 0; a < n_ags; a++) {
	 free_ags[a] += FDC_list[f].get_total_antigen(a);
       }
       // total amount of ICs for all Ags together:
       ics += FDC_list[f].get_total_immune_complex();
       // and for each antigen separately:
       for (int a = 0; a < n_ags; a++) {
	 ics_ags[a] += FDC_list[f].get_total_immune_complex(a);
       }
     }
     double total_ab = 0.;
     if (s.signal_use[antibody] == 1) {
       total_ab = s.get_signal_total(antibody);
     } else if (use_antibody_bins) {
       total_ab = Ab.get_total_antibodies();
     }
     ag_presentation << time << "   " << free_ag << "   " << ics / cellFDC::ag_threshold << "   "
		     << total_ab / cellFDC::ag_threshold << "   " << free_ag
       * cellFDC::ag_threshold << "   " << ics
		     << "   " << total_ab
       // <<"   "<<ics*cellFDC::ag_threshold
       // <<"   "<<total_ab*cellFDC::ag_threshold
		     << "\n";
     // now write out ag and ic for each antigen separately as well as the fraction of both:
     ofstream multiag_file("antigens.out", ofstream::app);
     for (int a = 0; a < n_ags; a++) {
       multiag_file << free_ags[a] << "   ";
     }
     multiag_file << "\n";
     multiag_file.close();
     ofstream multiic_file("ics.out", ofstream::app);
     for (int a = 0; a < n_ags; a++) {
       multiic_file << ics_ags[a] << "   ";
     }
     multiic_file << "\n";
     multiic_file.close();
     ofstream multiicfrag_file("ics_frac.out", ofstream::app);
     for (int a = 0; a < n_ags; a++) {
       multiicfrag_file << ics_ags[a] / (ics_ags[a] + free_ags[a] * cellFDC::ag_threshold)
			<< "   ";
     }
     // note that ics_ags is in Mol while free_ags is number of antigen portions of siez
     // ag_threshold
     multiicfrag_file << "\n";
     multiicfrag_file.close();
     
     if (use_antibody_bins) {
       total_ab = Ab.get_total_antibodies();
       Ab.show_antibodies(time);
       Ab.write_antibodies(time);
     }
     
     cum_ag << time << "   " << free_ag << "   " << ics << "   " << total_ab << "   "
	    << free_ag * cellFDC::ag_threshold << "   " << ics * cellFDC::ag_threshold << "   "
	    << total_ab * cellFDC::ag_threshold << "\n";
   }
   
   // Documentation during run
   /* cout<<" #CB_aktu="<<nd_CB
    *  <<" #CB_dif_state="<<n_diff_CB
    *  <<" #produced_CC="<<CC_total
    *  <<" sig_total="<<signal_total[sig_differ2CC]
    *  <<" sig_aktu="<<signal_aktuell[sig_differ2CC]
    *  <<"\n"; */

   // Sicherheitscheck:
   if (checkit == 1) {
      if (s.signal_use[sig_differ2CC] == 1) {
         double sig_differ2CC_aktuell = s.get_signal_total(sig_differ2CC);
         // #erzeugterCC+#CBimStatusDifferentiate=#erzeugterSignale-#aktuellerSignale
         double sigs;
         if (s.diffusion_mode == 2) {
            sigs = CC_total + n_diff_CB - s.signal_used[sig_differ2CC];
         } else {
            // cases l.diffusion_mode in [0,1]:
            sigs = s.signal_total[sig_differ2CC] - sig_differ2CC_aktuell;
            //                -l.signal_diffout[sig_differ2CC];
            // diffout ist bereits in _total enthalten! deswegen nicht abziehen.
            if (cellCB::receptor_use == 0) {
               sigs -= (CC_total + n_diff_CB);
            } else {
               sigs -= s.signal_used[sig_differ2CC];
            }
         }

         if (sigs < 0) {
            sigs *= -1.0;
         }
         if (sigs > 1.e-03) {
            cout << "FEHLER!\n"
                 << "CC_total=" << CC_total << " n_diff_CB=" << n_diff_CB
                 << " sig_total=" << s.signal_total[sig_differ2CC] << " sig_aktuell="
                 << sig_differ2CC_aktuell
                 << " sig_diffout=" << s.signal_diffout[sig_differ2CC]
                 << " sig_used=" << s.signal_used[sig_differ2CC] << "\n"
                 << "Zahl erzeugter Signale mit Zahl der Differenzierungsbefehle "
                 << "inkompatibel!\n";
            exit(1);
         }
      }
   }
}
void cellman::count_dec205_ova_positive() {
   long noCB = CB_list.benutzt();
   long noCC = CC_list.benutzt();
   long b_dec205p = 0, b_dec205n = 0, b_dec205ovap = 0, b_dec205ovan = 0;
   long c_dec205p = 0, c_dec205n = 0, c_dec205ovap = 0, c_dec205ovan = 0;

   // go through all centroblasts
   for (long i = 0; i < noCB; ++i) {
      if (CB_list[i].DEC205) {
         ++b_dec205p;
         if (CB_list[i].DEC205_ova) {
            ++b_dec205ovap;
         } else {
            ++b_dec205ovan;
         }
      } else {
         ++b_dec205n;
      }
   }
   // go through all centroblasts
   for (long i = 0; i < noCC; ++i) {
      if (CC_list[i].DEC205) {
         ++c_dec205p;
         if (CC_list[i].DEC205_ova) {
            ++c_dec205ovap;
         } else {
            ++c_dec205ovan;
         }
      } else {
         ++c_dec205n;
      }
   }
   // ++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++
   /*
    * cout<<b_dec205n<<" DEC-CB, "<<b_dec205p<<" DEC+CB of which "
    *  <<b_dec205ovap<<" OVA-activated and "<<b_dec205ovan<<" OVA-inactivated.\n";
    * cout<<c_dec205n<<" DEC-CC, "<<c_dec205p<<" DEC+CC of which "
    *  <<c_dec205ovap<<" OVA-activated and "<<c_dec205ovan<<" OVA-inactivated.\n";
    * cout<<b_dec205n+c_dec205n<<" DEC-BC, "<<b_dec205p+c_dec205p<<" DEC+BC of which "
    *  <<b_dec205ovap+c_dec205ovap<<" OVA-activated and "
    *  <<b_dec205ovan+c_dec205ovan<<" OVA-inactivated.\n";
    */
   // ++++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++++++++++++
}
void cellman::set_pov_sphere(int x,
                             int y,
                             int z,
                             double r,
                             double cr,
                             double cg,
                             double cb,
                             ofstream &fff) {
   fff << "sphere { <" << double (x) + drandom(0.4) - 0.2 << "," << double (y) + drandom(0.4)
      - 0.2 << ","
       << double (z) + drandom(0.4) - 0.2 << ">, " << r << " texture { pigment { color "
       << "red " << cr << " green " << cg << " blue " << cb << " } } }\n";
}
void cellman::xfiles(suffix tnr, space &l) {
   bool opendx = false;

   // cell files  ==============================
   char datname[30] = "";
   strcat(datname, "xy");
   strcat(datname, tnr);
   strcat(datname, ".ppm");
   cout << "Write to " << datname << "\n";
   char opendxname[30] = "";
   if (opendx) {
      strcat(opendxname, "xy");
      strcat(opendxname, tnr);
      strcat(opendxname, ".dat");
      cout << "Write to " << opendxname << "\n";
   }
   char povname[30] = "";
   strcat(povname, "xy");
   strcat(povname, tnr);
   strcat(povname, ".pov");
   cout << "Write to " << povname << "\n";

   // File movie ergaenzen:
   // movie<<datname<<" ";
   /* Troubling: The use of this call leads to an empty movie-file!!??
    * Putting the same call in hyphasma.C the movie-file is o.k.
    * But only if the movie.close() command is deleted in the close_files routine.
    * I guess I didn't get an important feature of the ofstreams.
    */

   // Stream oeffnen
   ofstream ff(datname);
   ff.setf(ios::scientific, ios::floatfield);
   ofstream ffdx;
   if (opendx) {
      ffdx.open(opendxname);
      ffdx.setf(ios::floatfield);
   }
   ofstream ffpov(povname);
   ffpov.setf(ios::floatfield);

   // Zeige Ebenen (je nach Dimension unterschiedlich)
   long k[l.dim];
   for (short d = 0; d < l.dim; d++) {
      k[d] = 0;
   }
   long d1 = 0;
   long d2 = 1;
   if (l.dim == 3) {
      d1 = 1;
      d2 = 2;
      k[0] = l.prodimvec[0] / 2;
   }

   // Kopf:
   ff << "P3\n";
   ff << l.prodimvec[d1] << "\n";
   ff << l.prodimvec[d2] << "\n";
   if (l.system == 0) {
      ff << "255\n";
   }
   ffpov << "#include \"colors.inc\"\n"
         << "#include \"textures.inc\"\n"
         << "#include \"shapes.inc\"\n\n"
      // ###<<"background { color Black }\n\n"
         << "background { color red 0.8 green 0.8 blue 0.8 }\n\n"
         << "camera { location <"
      // for r=160 and dx=5 micron in 3D
         << l.prodimvec[d1] / 2 << ", " << l.prodimvec[d2] / 2
         << ", -40> up 1.8*y right 1.8*x look_at <"
      // for r=220 and dx=5 micron in 2D
      //	<<l.prodimvec[d1]/2<<", "<<l.prodimvec[d2]/2<<", -40> up 2.3*y right 2.3*x look_at <"
         << l.prodimvec[d1] / 2 << ", " << l.prodimvec[d2] / 2 << ", 0> }\n\n"
      //	<<"light_source { <-20, -20, -80> color White }\n\n"
         << "light_source { <-30, -30, -100> color White }\n\n"
         << "global_settings { ambient_light rgb <3, 3, 3> }\n\n";

   // CBs Farben zuordnen:
   long ncb = CB_list.benutzt();
   short colcb[ncb];
   // +++ OPTION: color variance within CBs
   for (long a = 0; a < CB_list.benutzt(); a++) {
      colcb[a] = irandom(70);
   }

   // colcb[a]=50;
   // +++ OPTION end
   // +++ OPTION: chose true if two subsequent layer are to be visualised
   int layer_min = 0;
   int layer_max = 0;
   if (l.dim == 3) {
      // layer_min=-1;
      // layer_max=+1;
      layer_min = -5;
      layer_max = +5;
      if (photoactivation) {
         layer_min = int (-1. * double (l.prodimvec[0]) / 2.);
         layer_max = int (double (l.prodimvec[0]) / 2.);
      }
   }
   // +++ OPTION end

   // Gitterpunkte durchlaufen
   long i, li;
   double x, y, z;
   bool found;
   x = 0.;
   y = 0.;
   z = 0.;
   for (long i1 = 0; i1 < l.prodimvec[d2]; i1++) {
      for (long i2 = 0; i2 < l.prodimvec[d1]; i2++) {
         found = false;
         for (long i0 = layer_min; i0 <= layer_max; i0++) {
            k[d1] = i2;
            k[d2] = i1;
            if (l.dim == 3) {
               k[0] = int (l.prodimvec[0] / 2) + i0;
            }
            x = (i2 + 0.5) * l.dx;
            y = (i1 + 0.5) * l.dx;
            z = (i0 + 0.5) * l.dx;
            i = l.Index(k);
            li = l.cellknot[i].listi;

            if (l.knot[i].status != external) {
               switch (l.cellknot[i].cell) {
                  case nocell:
                     if (l.dim == 2) {
                        ff << "0 0 0\n";          // schwarz
                     }
                     break;

                  case CB:
                     if (show_mode == tumour) {
                        if (CB_list[li].status == necrotic) {
                           if (found == false) {
                              ff << "200 0 0\n";            // hellrot
                           }
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 7 0.9 0 0\n";
                           }
                        } else {
                           if (CB_list[li].get_contact(nocell, l)
                               == -1) {
                              // case no free neighbour -> quiescent
                              if (found == false) {
                                 ff << "150 150 150\n";             // hellrot
                              }
                              if (opendx) {
                                 ffdx << x << " " << y << " " << z << " 7 0.9 0.9 0.9\n";
                              }
                           } else {
                              if (found == false) {
                                 ff << "0 200 0\n";             // hellrot
                              }
                              if (opendx) {
                                 ffdx << x << " " << y << " " << z << " 7 0 0.9 0\n";
                              }
                           }
                        }
                     } else {
                        if (photoactivation) {
                           if (CB_list[li].trackit) {
                              if (found == false) {
                                 ff << "255 0 255"
                                    << "\n";             // hellpink
                              }
                              set_pov_sphere(
                                 k[d1],
                                 l.prodimvec[d2] - k[d2],
                                 0,
                                 1.23,
                                 1.0,
                                 0,
                                 1.0,
                                 ffpov);                                                                   //
                                                                                                           //
                                                                                                           //
                                                                                                           //
                                                                                                           // hellpink
                              // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                              if (opendx) {
                                 ffdx << x << " " << y << " " << z << " 7 1 0 0\n";
                              }
                              found = true;
                           } else if (l.dim == 2) {
                              ff << "0 0 0\n";            // schwarz
                           }
                        } else if (def_DEC205) {
                           if (CB_list[li].DEC205) {
                              if (found == false) {
                                 ff << "255 0 255"
                                    << "\n";             // hellpink
                              }
                              set_pov_sphere(
                                 k[d1],
                                 l.prodimvec[d2] - k[d2],
                                 0,
                                 1.23,
                                 1.0,
                                 0,
                                 1.0,
                                 ffpov);                                                                   //
                                                                                                           //
                                                                                                           //
                                                                                                           //
                                                                                                           // hellpink
                              if (opendx) {
                                 ffdx << x << " " << y << " " << z << " 7 1 0 0\n";
                              }
                              found = true;
                           } else if (l.dim == 2) {
                              ff << "0 0 0\n";            // schwarz
                           }
                        } else if (CB_list[li].n_recycling == 0) {
                           if ((CB_list[li].state != cb_differentiate)
                               && (CB_list[li].state != cb_stop_dividing)) {
                              if (show_Ki67 == 0) {
                                 // if (found==false) ff<<"0 0 "<<150+colcb[li]<<"\n"; //
                                 // dunkelblau-toene
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,0,0.7,ffpov);
                                 // //
                                 // dunkelblau
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.8,0,0,0.7,ffpov);
                                 // //
                                 // dunkelblau
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                1.23,
                                                0,
                                                0,
                                                0.7,
                                                ffpov);             // dunkelblau
                                 // if (found==false) ff<<150+colcb[li]<<" 0 0\n"; // dunkelrot
                                 if (found == false) {
                                    switch (CB_list[li].state) {
                                       case cb_normal: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_differentiate: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_G1: {
                                          ff << "100 0 255\n";
                                          break;
                                       }

                                       case cb_S: {
                                          ff << "150 0 255\n";
                                          break;
                                       }

                                       case cb_G2: {
                                          ff << "200 0 255\n";
                                          break;
                                       }

                                       case cb_M: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_G0: {
                                          ff << "50 0 255\n";
                                          break;
                                       }

                                       case cb_divide: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_stop_dividing: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_statenumber:
                                          break;
                                    }
                                 }
                                 // if (found==false) ff<<"255 255 255\n";
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z
                                         << " 7 0.25 0 0\n";
                                 }
                              } else {
                                 if (CB_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    // if (found==false) ff<<"0 0 255\n";
                                    // if (found==false) ff<<"255 255 255\n";
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.25 0 0\n";
                                    }
                                 }
                              }
                              // sig<<"150 0 0\n";
                           } else {
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "255 0 " << 150 + colcb[li] << "\n";              // dunkelblau-toene
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0,0.7,ffpov);
                                 // //
                                 // dunkelblau
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.8,0,0,0.7,ffpov);
                                 // //
                                 // dunkelblau
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                1.23,
                                                0,
                                                0,
                                                0.7,
                                                ffpov);             // dunkelblau
                                 // if (found==false) ff<<150+colcb[li]<<" 0 0\n"; // dunkelrot
                                 // if (found==false) ff<<"255 0 0\n";
                                 // if (found==false) ff<<"255 255 255\n";
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z
                                         << " 7 0.25 0 0\n";
                                 }
                              } else {
                                 if (CB_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "0 0 200\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.25 0 0\n";
                                    }
                                 }
                              }
                              // sig<<"100 100 0\n"; // dunkelrot gruen
                           }
                        } else {
                           if (CB_list[li].n_recycling == 1) {
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    switch (CB_list[li].state) {
                                       case cb_normal: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_differentiate: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_G1: {
                                          ff << "100 0 255\n";
                                          break;
                                       }

                                       case cb_S: {
                                          ff << "150 0 255\n";
                                          break;
                                       }

                                       case cb_G2: {
                                          ff << "200 0 255\n";
                                          break;
                                       }

                                       case cb_M: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_G0: {
                                          ff << "50 0 255\n";
                                          break;
                                       }

                                       case cb_divide: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_stop_dividing: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_statenumber:
                                          break;
                                    }
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,0,1.0,ffpov);
                                 // //
                                 // hellblau
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.8,0,0,1.0,ffpov);
                                 // //
                                 // hellblau
                                 set_pov_sphere(
                                    k[d1],
                                    l.prodimvec[d2] - k[d2],
                                    0,
                                    1.23,
                                    0,
                                    0,
                                    1.0,
                                    ffpov);                                                                  //
                                                                                                             //
                                                                                                             //
                                                                                                             //
                                                                                                             // hellblau
                                 // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 7 1 0 0\n";
                                 }
                              } else {
                                 if (CB_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "0 0 200\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.25 0 0\n";
                                    }
                                 }
                              }
                              // sig<<"255 0 0\n"; // hellrot
                           } else {
                              // mehrfaches recycling
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    switch (CB_list[li].state) {
                                       case cb_normal: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_differentiate: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_G1: {
                                          ff << "100 0 255\n";
                                          break;
                                       }

                                       case cb_S: {
                                          ff << "150 0 255\n";
                                          break;
                                       }

                                       case cb_G2: {
                                          ff << "200 0 255\n";
                                          break;
                                       }

                                       case cb_M: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_G0: {
                                          ff << "50 0 255\n";
                                          break;
                                       }

                                       case cb_divide: {
                                          ff << "255 0 255\n";
                                          break;
                                       }

                                       case cb_stop_dividing: {
                                          ff << "0 255 255\n";
                                          break;
                                       }

                                       case cb_statenumber:
                                          break;
                                    }
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,0,1.0,ffpov);
                                 // //
                                 // hellblau
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.8,0,0,1.0,ffpov);
                                 // //
                                 // hellblau
                                 set_pov_sphere(
                                    k[d1],
                                    l.prodimvec[d2] - k[d2],
                                    0,
                                    1.23,
                                    0,
                                    0,
                                    1.0,
                                    ffpov);                                                                  //
                                                                                                             //
                                                                                                             //
                                                                                                             //
                                                                                                             // hellblau
                                 // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 7 1 0 0\n";
                                 }
                              } else {
                                 if (CB_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "0 0 200\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 7 0.25 0 0\n";
                                    }
                                 }
                              }
                              // sig<<"255 0 0\n"; // hellrot
                           }
                        }
                     }
                     break;

                  case CC:
                     if (photoactivation) {
                        if (CC_list[li].trackit) {
                           if (found == false) {
                              ff << "0 255 0"
                                 << "\n";            // hellgruen
                           }
                           set_pov_sphere(
                              k[d1], l.prodimvec[d2] - k[d2], 0, 0.93, 0, 1.0, 0, ffpov);            //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     // hellgruen
                           // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                           }
                           found = true;
                        } else if (l.dim == 2) {
                           ff << "0 0 0\n";           // schwarz
                        }
                     } else if (def_DEC205) {
                        if (CC_list[li].DEC205) {
                           if (found == false) {
                              ff << "0 255 0"
                                 << "\n";            // hellgruen
                           }
                           set_pov_sphere(
                              k[d1], l.prodimvec[d2] - k[d2], 0, 0.93, 0, 1.0, 0, ffpov);            //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     // hellgruen
                           // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                           }
                           found = true;
                        } else if (l.dim == 2) {
                           ff << "0 0 0\n";           // schwarz
                        }
                     } else {
                        switch (CC_list[li].state) {
                           case unselected:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 100 0\n";              // dunkelgruen
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0.6,0,ffpov);
                                 // //
                                 // dundelgruen
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,0.6,0,ffpov);
                                 // //
                                 // dundelgruen
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                0.93,
                                                0,
                                                0.6,
                                                0,
                                                ffpov);             // dunkelgruen
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                                 }
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.25 0 0\n";
                                    }
                                 }
                              }
                              break;

                           case contact:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 200 0\n";              // hellgruen
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0.6,0,ffpov);
                                 // //
                                 // dundelgruen
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,0.6,0,ffpov);
                                 // //
                                 // dundelgruen
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                0.93,
                                                0,
                                                1.0,
                                                0,
                                                ffpov);             // hellgruen
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                                 }
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.25 0 0\n";
                                    }
                                 }
                              }
                              break;

                           case TCcontact:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 200 0\n";              // hellgruen
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0.8,0,ffpov);
                                 // //
                                 // hellgruen
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,1.0,0,ffpov);
                                 // //
                                 // hellgruen
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                0.93,
                                                0,
                                                1.0,
                                                0,
                                                ffpov);             // hellgruen
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                                 }
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.25 0 0\n";
                                    }
                                 }
                              }
                              break;

                           case FDCselected:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 100 0\n";              // dunkelgruen
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0.8,0,ffpov);
                                 // //
                                 // hellgruen
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,1.0,0,ffpov);
                                 // //
                                 // hellgruen
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                0.93,
                                                0,
                                                0.6,
                                                0,
                                                ffpov);             // dunkelgruen
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                                 }
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.25 0 0\n";
                                    }
                                 }
                              }
                              break;

                           case selected:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 200 200\n";              // cyan
                                 }
                                 // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,1.0,0,ffpov);
                                 // //
                                 // hellstesgruen
                                 // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,1.0,0,ffpov);
                                 // //
                                 // hellstesgruen
                                 set_pov_sphere(
                                    k[d1],
                                    l.prodimvec[d2] - k[d2],
                                    0,
                                    0.93,
                                    0,
                                    1.0,
                                    1.0,
                                    ffpov);                                                                    //
                                                                                                               //
                                                                                                               //
                                                                                                               //
                                                                                                               // cyan
                                 if (opendx) {
                                    ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                                 }
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.9 0.9 0.9\n";
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                    if (opendx) {
                                       ffdx << x << " " << y << " " << z
                                            << " 3 0.25 0 0\n";
                                    }
                                 }
                              }
                              break;

                           case apoptosis:
                              if (show_Ki67 == 0) {
                                 if (found == false) {
                                    ff << "0 30 0\n";              // ddunkelgruen
                                 }
                                 set_pov_sphere(k[d1],
                                                l.prodimvec[d2] - k[d2],
                                                0,
                                                0.69,
                                                0,
                                                0.3,
                                                0,
                                                ffpov);             // ddunkelgruen
                              } else {
                                 if (CC_list[li].Ki67 == 0) {
                                    if (found == false) {
                                       ff << "150 150 150\n";               // hellgrau
                                    }
                                 } else {
                                    if (found == false) {
                                       ff << "200 0 0\n";               // dunkelrot fuer
                                                                        // Ki67-markiert
                                    }
                                 }
                              }
                              break;
                              // MMH2Marta: we have to show TFR on the lattice picture
                                            ///Marta2MMH done BUT colors need to be differently set (TODO!)
            case TFRcontact:
                            if (show_Ki67 == 0) {
                               if (found == false) {
                                  ff << "0 200 0\n";              // hellgruen
                               }
                               // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.5,0,0.8,0,ffpov);
                               // //
                               // hellgruen
                               // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0,1.0,0,ffpov);
                               // //
                               // hellgruen
                               set_pov_sphere(k[d1],
                                              l.prodimvec[d2] - k[d2],
                                              0,
                                              0.93,
                                              0,
                                              1.0,
                                              0,
                                              ffpov);             // hellgruen
                               if (opendx) {
                                  ffdx << x << " " << y << " " << z << " 3 0 1 1\n";
                               }
                            } else {
                               if (CC_list[li].Ki67 == 0) {
                                  if (found == false) {
                                     ff << "150 150 150\n";               // hellgrau
                                  }
                                  if (opendx) {
                                     ffdx << x << " " << y << " " << z
                                          << " 3 0.9 0.9 0.9\n";
                                  }
                               } else {
                                  if (found == false) {
                                     ff << "200 0 0\n";               // dunkelrot fuer
                                                                      // Ki67-markiert
                                  }
                                  if (opendx) {
                                     ffdx << x << " " << y << " " << z
                                          << " 3 0.25 0 0\n";
                                  }
                               }
                            }
                            break;
                        }
                     }
                     break;

                  case TC:
                     if ((photoactivation == false) && (def_DEC205 == false)) {
                        if (found == false) {
                           ff << "255 0 0\n";           // rot
                        }
                        // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0.7,0,0,ffpov); //
                        // hellrot
                        set_pov_sphere(k[d1],
                                       l.prodimvec[d2] - k[d2],
                                       0,
                                       0.93,
                                       1.0,
                                       0,
                                       0,
                                       ffpov);                                                              //
                                                                                                            //
                                                                                                            //
                                                                                                            //
                                                                                                            // hellrot
                        // if (found==false) ff<<"0 0 255\n"; // hellblau
                        if (opendx) {
                           ffdx << x << " " << y << " " << z << " 3 0 0 1\n";
                        }
                     } else if (l.dim == 2) {
                        ff << "0 0 0\n";          // schwarz
                     }
                     break;

                     //MSchips
                     case TFR:
                        if ((photoactivation == false) && (def_DEC205 == false)) {
                           if (found == false) {
                              ff << "128 128 0\n";           //128 128 0 olive green //255	0	255fuchsia
                           }
                           // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0.7,0,0,ffpov); //
                           // hellrot
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          1.0,
                                          0,
                                          0,
                                          ffpov);                                                              //
                                                                                                               //
                                                                                                               //
                                                                                                               //
                                                                                                               // hellrot
                           // if (found==false) ff<<"0 0 255\n"; // hellblau
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0 0 1\n";
                           }
                        } else if (l.dim == 2) {
                           ff << "0 0 0\n";          // schwarz
                        }
                        break;

                  case BETA:
                     if ((photoactivation == false) && (def_DEC205 == false)) {
                        //	 cout<<"in beta!;  ";
                        if (BETA_list[li].y_n[cellbeta::V] < -0.06) {
                           if (found == false) {
                              ff << "0 0 255\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0,
                                          0,
                                          1.0,
                                          ffpov);
                        } else if (BETA_list[li].y_n[cellbeta::V] < -0.055) {
                           if (found == false) {
                              ff << "125 0 130\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0.5,
                                          0,
                                          0.5,
                                          ffpov);
                        } else if (BETA_list[li].y_n[cellbeta::V] < -0.05) {
                           if (found == false) {
                              ff << "150 0 110\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0.6,
                                          0,
                                          0.4,
                                          ffpov);
                        } else if (BETA_list[li].y_n[cellbeta::V] < -0.045) {
                           if (found == false) {
                              ff << "170 0 80\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0.7,
                                          0,
                                          0.3,
                                          ffpov);
                        } else if (BETA_list[li].y_n[cellbeta::V] < -0.04) {
                           if (found == false) {
                              ff << "200 0 60\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0.8,
                                          0,
                                          0.2,
                                          ffpov);
                        } else if (BETA_list[li].y_n[cellbeta::V] < -0.03) {
                           if (found == false) {
                              ff << "230 0 30\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          0.9,
                                          0,
                                          0.1,
                                          ffpov);
                        } else {
                           if (found == false) {
                              ff << "255 0 0\n";
                           }
                           set_pov_sphere(k[d1],
                                          l.prodimvec[d2] - k[d2],
                                          0,
                                          0.93,
                                          1.0,
                                          0,
                                          0,
                                          ffpov);
                        }
                        if (opendx) {
                           ffdx << x << " " << y << " " << z << " 3 0 0 1\n";
                        }
                     } else if (l.dim == 2) {
                        // cout<<"in black!;  ";
                        ff << "0 0 0\n";          // schwarz
                     }
                     break;

                  case FDC:
                     if ((photoactivation == false) && (def_DEC205 == false)) {
                        if (found == false) {
                           ff << "255 255 0\n";           // gelb
                        }
                        // ###set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.8,0.8,0.8,0.,ffpov); //
                        // gelb
                        // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,1.2,0.7,0.7,0.,ffpov); //
                        // gelb
                        set_pov_sphere(k[d1],
                                       l.prodimvec[d2] - k[d2],
                                       0,
                                       1.83,
                                       0.7,
                                       0.7,
                                       0.,
                                       ffpov);                                                                 //
                                                                                                               //
                                                                                                               //
                                                                                                               //
                                                                                                               // gelb
                        // if (found==false) ff<<"191 191 0\n"; // gelb
                        if (li == -1) {
                           cout << "li for FDC ==-1\n;";
                           cout << "FDClisti=" << l.cellknot[i].FDClisti << "\n";
                           exit(1);
                        }
                        if (FDC_list[li].state == soma) {
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 7 1 1 0\n";
                              for (short xxx = 0; xxx < 15; xxx++) {
                                 ffdx << x + l.dx / 2. + 2. * xxx << " " << y << " "
                                      << z << " 2 1 1 0\n";
                                 ffdx << x - l.dx / 2. - 2. * xxx << " " << y << " "
                                      << z << " 2 1 1 0\n";
                                 ffdx << x << " " << y + l.dx / 2. + 2. * xxx << " "
                                      << z << " 2 1 1 0\n";
                                 ffdx << x << " " << y - l.dx / 2. - 2. * xxx << " "
                                      << z << " 2 1 1 0\n";
                              }
                           }
                        }
                        // sig<<"191 191 0\n"; // gelb
                     } else if (l.dim == 2) {
                        ff << "0 0 0\n";          // schwarz
                     }
                     break;

                  case out:
                     if (def_DEC205) {
                        if (l.dim == 2) {
                           ff << "0 0 0\n";           // schwarz
                        }
                     } else if (photoactivation) {
                        if (OUT_list[li].trackit) {
                           if (found == false) {
                              ff << "255 255 255"
                                 << "\n";            // weiss
                           }
                           set_pov_sphere(
                              k[d1],
                              l.prodimvec[d2] - k[d2],
                              0,
                              1.23,
                              0.8,
                              0.8,
                              0.8,
                              ffpov);                                                                    //
                                                                                                         //
                                                                                                         //
                                                                                                         //
                                                                                                         // grau
                           // if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0 1 0\n";
                           }
                           found = true;
                        } else if (l.dim == 2) {
                           ff << "0 0 0\n";           // schwarz
                        }
                     } else
                     /*
                      * if (def_DEC205) {
                      * if (OUT_list[li].DEC205) {
                      * if (found==false) ff<<"255 255 255"<<"\n"; // weiss
                      * set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,1.23,0.8,0.8,0.8,ffpov); //
                      *grau
                      * //if (found==false) ff<<"255 0 "<<50+colcb[li]<<"\n"; // hellrot
                      * if (opendx) ffdx<<x<<" "<<y<<" "<<z<<" 3 0 1 0\n";
                      * found=true;
                      * } else if (l.dim==2) ff<<"0 0 0\n"; // schwarz
                      * }
                      * else
                      */
                     if (show_Ki67 == 0) {
                        if (found == false) {
                           ff << "191 191 191\n";           // gelb
                        }
                        // set_pov_sphere(k[d1],l.prodimvec[d2]-k[d2],0,0.6,0.8,0.8,0.8,ffpov); //
                        // hellgruen
                        set_pov_sphere(
                           k[d1], l.prodimvec[d2] - k[d2], 0, 1.23, 0.8, 0.8, 0.8, ffpov);           //
                                                                                                     //
                                                                                                     //
                                                                                                     //
                                                                                                     // hellgrau
                        if (opendx) {
                           ffdx << x << " " << y << " " << z << " 3 0 1 0\n";
                        }
                     } else {
                        if (OUT_list[li].Ki67 == 0) {
                           if (found == false) {
                              ff << "150 150 150\n";            // hellgrau
                           }
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0.9 0.9 0.9\n";
                           }
                        } else {
                           if (found == false) {
                              ff << "200 0 0\n";            // dunkelrot fuer Ki67-markiert
                           }
                           if (opendx) {
                              ffdx << x << " " << y << " " << z << " 3 0.25 0 0\n";
                           }
                        }
                     }
                     break;

                  /*
                   * case ext:
                   * if (found==false) ff<<"63 63 63\n"; // dunkelgrau
                   * break;
                   */
                  case blast1:
                     break;

                  case blast2:
                     if ((photoactivation == false) && (found == false)) {
                        ff << "63 63 63\n";          // dunkelgrau
                     }
                     break;

                  case N_cells:
                     break;
               }
            } else     // i.e. if (l.knot[i].status==external)
            if (l.dim == 2) {
               ff << "63 63 63\n";       // dunkelgrau
            }
            if ((l.cellknot[i].FDClisti > -1)
                && (FDC_list[l.cellknot[i].FDClisti].state == dendrite)) {
               if (opendx) {
                  ffdx << x << " " << y << " " << z << " 2 1 1 0\n";
               }
            }
            // if photoactivation is active don't adjust found, as it is done above
            if ((photoactivation == false) && (def_DEC205 == false)
                && (l.knot[i].status != external)
                && (l.cellknot[i].cell != nocell)) {
               found = true;
            }
         }
         // if no cell was written write something in ppm:
         if ((found == false) && (l.dim == 3)) {
            // Recalculate the index of the central layer
            k[d1] = i2;
            k[d2] = i1;
            k[0] = int (l.prodimvec[0] / 2);      // used i0=0!
            x = (i2 + 0.5) * l.dx;
            y = (i1 + 0.5) * l.dx;
            z = (0.5) * l.dx;      // used i0=0!
            i = l.Index(k);
            if (l.knot[i].status == external) {
               ff << "63 63 63\n";       // dunkelgrau
            } else {
               // if l.cellknot[i].cell==nocell
               ff << "0 0 0\n";          // black
            }
         }
      }
   }

   // Stream schliessen
   ff.close();
   if (opendx) {
      ffdx.close();
   }
   ffpov.close();
   //cerr << "Finished writing to xy....pov.\n";
}
void cellman::show_BrdU() {
   int BrdU_bins = 50;
   double BrdU_max = 100.;
   double BrdU_min = 0.09765625;
   double BrdU_level[BrdU_bins];
   int BrdU_array[BrdU_bins];
   // initialise x values on basis of log_2(BrdU)
   double dBrdU = ((log(BrdU_max) / log(2.)) - (log(BrdU_min) / log(2.))) / double (BrdU_bins);
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_max) / log(2.)) - ((double (j) + 0.5) * dBrdU)));
      // note that the +0.5 makes that BrdU_level contains the lower bound of the bin-interval
      BrdU_array[j] = 0;
   }
   // fill the bins with the cells on CB_list[]
   int bin;
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      bin = 0;
      while ((bin < BrdU_bins) && (CB_list[i].BrdU < BrdU_level[bin])) {
         ++bin;
      }
      if (bin == BrdU_bins) {
         --bin;
      }
      ++BrdU_array[bin];
   }
   // shift levels back to the centre of the intervals
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_level[j]) / log(2.)) + (0.5 * dBrdU)));
   }
   ofstream brdu_out;
   brdu_out.open("brdu.out");
   brdu_out << "! BrdU bin (centre) : # of cells\n";
   for (int j = 0; j < BrdU_bins; ++j) {
      brdu_out << BrdU_level[j] << "   " << BrdU_array[j] << "\n";
   }
   brdu_out.close();
}
void cellman::show_cell_in_BrdU() {
   // searches in CB_list for BrdU bins defined above
   // and attribute high and low antigen loaded cells
   // ++++++++++++++++ OPTION +++++++++++++++++++++++
   int BrdU_bins = 10;
   double BrdU_max = 100.;
   double BrdU_min = 0.09765625;
   double agx_threshold1 = 20.0;
   double agx_threshold2 = 5.0;
   // ++++++++++++ end OPqTION +++++++++++++++++++++++
   double BrdU_level[BrdU_bins];
   enum cell_ag_class {
      all, ag_low, ag_middle, ag_high, n_ag_classes
   };
   int BrdU_array[BrdU_bins][n_ag_classes];
   // initialise x values on basis of log_2(BrdU)
   double dBrdU = ((log(BrdU_max) / log(2.)) - (log(BrdU_min) / log(2.))) / double (BrdU_bins);
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_max) / log(2.)) - ((double (j) + 0.5) * dBrdU)));
      // note that the +0.5 makes that BrdU_level contains the lower bound of the bin-interval
      for (short k = 0; k < short (n_ag_classes); ++k) {
         BrdU_array[j][k] = 0;
      }
   }
   // fill the bins with the cells on CB_list[]
   int bin;
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      bin = 0;
      while ((bin < BrdU_bins) && (CB_list[i].BrdU < BrdU_level[bin])) {
         ++bin;
      }
      if (bin == BrdU_bins) {
         --bin;
      }
      ++BrdU_array[bin][all];
      if (CB_list[i].retained_ag > agx_threshold1) {
         ++BrdU_array[bin][ag_high];
      } else if (CB_list[i].retained_ag > agx_threshold2) {
         ++BrdU_array[bin][ag_middle];
      } else {
         ++BrdU_array[bin][ag_low];
      }
   }
   // shift levels back to the centre of the intervals
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_level[j]) / log(2.)) + (0.5 * dBrdU)));
   }
   ofstream brdu_cell_out;
   brdu_cell_out.open("brdu_cell.out");
   brdu_cell_out
      << "! Generation : BrdU bin (centre) : total # of cells : # of high ag "
      << ": # of middle ag : # of low ag\n";
   for (int j = 0; j < BrdU_bins; ++j) {
      brdu_cell_out << j << "    " << BrdU_level[j];
      for (short k = 0; k < short (n_ag_classes); ++k) {
         brdu_cell_out << "   " << BrdU_array[j][k];
      }
      brdu_cell_out << "\n";
   }
   brdu_cell_out.close();
}
void cellman::show_BrdU_at_time() {
  /* Goes through all BC and attributes the cellular BrdU level to BrdU bins */
   // and attribute high and low antigen loaded cells
   // ++++++++++++++++ OPTION +++++++++++++++++++++++
   int BrdU_bins = 10;
   double BrdU_max = 100.;
   double BrdU_min = 0.09765625;
   // ++++++++++++ end OPTION +++++++++++++++++++++++
   double BrdU_level[BrdU_bins];
   int BrdU_array[BrdU_bins];
   // initialise x values on basis of log_2(BrdU)
   double dBrdU = ((log(BrdU_max) / log(2.)) - (log(BrdU_min) / log(2.))) / double (BrdU_bins);
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_max) / log(2.)) - ((double (j) + 0.5) * dBrdU)));
      // note that the +0.5 makes that BrdU_level contains the lower bound of the bin-interval
      BrdU_array[j] = 0;
   }
   // fill the bins with the cells on CB_list[] and CC_list[]
   int bin;
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      bin = 0;
      while ((bin < BrdU_bins) && (CB_list[i].BrdU < BrdU_level[bin])) { ++bin; }
      if (bin == BrdU_bins) { --bin; }
      ++BrdU_array[bin];
   }
   for (long i = 0; i < CC_list.benutzt(); ++i) {
      bin = 0;
      while ((bin < BrdU_bins) && (CC_list[i].BrdU < BrdU_level[bin])) { ++bin; }
      if (bin == BrdU_bins) { --bin; }
      ++BrdU_array[bin];
   }
   // shift levels back to the centre of the intervals for read-out
   for (int j = 0; j < BrdU_bins; ++j) {
      BrdU_level[j] = pow(2.0, ((log(BrdU_level[j]) / log(2.)) + (0.5 * dBrdU)));
   }
   double tt = time - tfirst_inject_BrdU;
   if (tt < 0) { tt = 0; }
   ofstream brdu_cell_t;
   brdu_cell_t.open("brdu_t.out", ofstream::app);
   for (int j = 0; j < BrdU_bins; ++j) {
     brdu_cell_t << tt << "  " << time << "   " 
		 << j << "    " << BrdU_level[j]
		 << "    " << BrdU_array[j] << "\n";
   }
   brdu_cell_t.close();
}
void cellman::show_BCstate_in_BrdUpositive() {
  /* Shows the BC states in the subset of BrdU positive cells 
   * above the threshold BrdU_detection_threshold.
   * The distinguished states are G1,G0,S,G2,M,CCnotyetselected.
   * Note that this is only called when cell cycle phases are explicitly followed.
   */
  long phases[cb_statenumber];
  long totalBC = 0, CCnotyetselected = 0;
  for (short j = 0; j < cb_statenumber; j++) 
    { phases[j] = 0; }
  for (long i = 0; i < CB_list.benutzt(); i++ ) {
    if (CB_list[i].BrdU > BrdU_detection_threshold) {
      ++phases[CB_list[i].state];
      ++totalBC;
    }
  }
  /* Now add from CC those who are in waiting mode.
   * The thought is that these cells actually are already in cell cycle
   * but are still sensitive to CXCL13, thus not passing to the DZ yet.
   * In order to reproduce the number of cells in different cell cycle
   * phases, ignoring CC which are already dividing would shift the
   * numbers towards later phases, which is a mistake. Therefore,
   * selected CC will be attributed to a cell cycle phase depending
   * on the time they are already in the selected state.
   * TODO: It might be better to let selected CC differentiate to
   * CBs right away, which then could enter cell cycle normally,
   * and to delay only the up-regulation of CXCR4. Then, this artificial
   * correction could be deleted again.
   */
  if ( cellCB::shiftCCdelay2CBcycle() ) {
    for (long i = 0; i < CC_list.benutzt(); i++ ) {
      if ( CC_list[i].state == selected ) {
	centroblasts phase = cellCB::get_virtual_cell_cycle_phase(CC_list[i].selected_clock);
	if (CC_list[i].BrdU > BrdU_detection_threshold) {
	  ++phases[phase];
	  ++totalBC;
	}
      }
      else {
	if (CC_list[i].BrdU > BrdU_detection_threshold) {
	  ++CCnotyetselected;
	  ++totalBC;
	}
      }
    }
  }
  ofstream brduphases;
  brduphases.open("cellcycle_phases_brdu.out", ofstream::app);
  double tt = time - tfirst_inject_BrdU;
  if (tt < 0) { tt = 0; }
  brduphases << tt << "  " << time << "  " << totalBC << "    ";
  short first = int(cb_G1), last = int(cb_M) + 1; 
  for (short c = first; c < last; c++) {
    brduphases << phases[c] << "  ";
  }
  brduphases << CCnotyetselected << "     ";
  double frac = 0.;
  for (short c = first; c < last; c++) {
    frac = 0.; 
    if (totalBC > 0) { frac = phases[c]/double(totalBC); }
    brduphases << frac << "  ";
  }
  frac = 0.;
  if (totalBC > 0) { frac = CCnotyetselected/double(totalBC); }
  brduphases << frac << "     ";
  brduphases << "\n";
  brduphases.close();
}
void cellman::show_ag_loaded() {
   // ++++++++++++++++ OPTION +++++++++++++++++++++++
   bool use_log = true;
   int ag_bins = 39;
   double ag_max = 400.;
   double ag_min = 0.09765625 / 2;
   // 0.0030517578125; //100*2^{-15}
   // 100*2^(-10)=0.09765625;
   // ++++++++++++ end OPTION +++++++++++++++++++++++
   double ag_level[ag_bins];
   int ag_array[ag_bins];
   double dag;
   if (use_log) {
      // initialise x values on basis of log_2(ag)
      dag = ((log(ag_max) / log(2.)) - (log(ag_min) / log(2.))) / double (ag_bins);
      for (int j = 0; j < ag_bins; ++j) {
         ag_level[j] = pow(2.0, ((log(ag_max) / log(2.)) - ((double (j) + 0.5) * dag)));
         // note that the +0.5 makes that BrdU_level contains the lower bound of the bin-interval
         ag_array[j] = 0;
      }
   } else {
      // use linear bins
      dag = (ag_max - ag_min) / double (ag_bins);
      for (int j = 0; j < ag_bins; ++j) {
         ag_level[j] = ag_max - ((double (j) + 0.5) * dag);
         ag_array[j] = 0;
      }
   }
   // fill the bins with the cells on CB_list[]
   int bin;
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      bin = 0;
      while ((bin < ag_bins) && (CB_list[i].retained_ag < ag_level[bin])) {
         ++bin;
      }
      if (bin == ag_bins) {
         --bin;
      }
      ++ag_array[bin];
   }
   // shift levels back to the centre of the intervals
   for (int j = 0; j < ag_bins; ++j) {
      if (use_log) {
         ag_level[j] = pow(2.0, ((log(ag_level[j]) / log(2.)) + (0.5 * dag)));
      } else {
         ag_level[j] = ag_level[j] + (0.5 * dag);
      }
   }
   ofstream ag_out;
   ag_out.open("ag_loaded.out");
   ag_out << "! Shows retained antigen in dividing BCs, i.e. CBs\n"
          << "! antigen bin (centre) : # of cells\n";
   for (int j = 0; j < ag_bins; ++j) {
      ag_out << ag_level[j] << "   " << ag_array[j] << "\n";
   }
   ag_out.close();
}
void cellman::show_ag_BrdU(bool ini) {
   ofstream ag_brdu;
   if (ini) {
      ag_brdu.open("ag_brdu_ini.out");
   } else {
      ag_brdu.open("ag_brdu.out");
   }
   ag_brdu << "! brdu : ag\n";
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      // the following makes the scale linear for small values and log_2 for large ones
      double tmp = CB_list[i].retained_ag;
      double hill_k = 1.0;
      double hill = tmp * tmp / (hill_k * hill_k + tmp * tmp);
      if (tmp > 0.) {
         tmp = (1.0 - hill) * tmp + hill * log(tmp) / log(2.);
      } // see supplement of Cell Reports 2012 for more details
      ag_brdu << CB_list[i].BrdU << "   " << CB_list[i].retained_ag << "    " << tmp << "\n";
   }
   // add CC in virtual cell cycle here?
   ag_brdu.close();
}
//MS
void cellman::get_average_mutation_prob_SELF(double &muta_avS, double &muta_sdS) {
   muta_avS = 0;
   muta_sdS = 0;
   long muta_nn = CB_list.benutzt();
   long count_self=0;
   for (long i = 0; i < muta_nn; ++i) {
       if (CB_list[i].selfMutation) {
           muta_avS += CB_list[i].p_mutation;
           ++count_self;
       }
   }
   if (count_self > 0) {
      muta_avS /= double (count_self);
   }
   for (long i = 0; i < muta_nn; ++i) {
       if (CB_list[i].selfMutation) {
           muta_sdS += (CB_list[i].p_mutation - muta_avS) * (CB_list[i].p_mutation - muta_avS);
       }
   }
   if (count_self > 1) {
      muta_sdS /= double (count_self - 1);
      muta_sdS = sqrt(muta_sdS);
   } else {
      muta_sdS = 0;
   }
}
void cellman::get_average_mutation_prob(double &muta_av, double &muta_sd) {
   muta_av = 0;
   muta_sd = 0;
   long muta_nn = CB_list.benutzt();
   for (long i = 0; i < muta_nn; ++i) {
      muta_av += CB_list[i].p_mutation;
   }
   if (muta_nn > 0) {
      muta_av /= double (muta_nn);
   }
   for (long i = 0; i < muta_nn; ++i) {
      muta_sd += (CB_list[i].p_mutation - muta_av) * (CB_list[i].p_mutation - muta_av);
   }
   if (muta_nn > 1) {
      muta_sd /= double (muta_nn - 1);
      muta_sd = sqrt(muta_sd);
   } else {
      muta_sd = 0;
   }
}
void cellman::write_mutations(double time) {
  //cerr << "Start write_mutations(double) ... \n";
   // +++++++++++++++ OPTION ++++++++++++++++
   int max_mutation = 20;
   // +++++++++++ end OPTION ++++++++++++++++
   // define a set of counters for each mutation number
   int nbc[max_mutation + 1];
   for (int n = 0; n <= max_mutation; n++) {
      nbc[n] = 0;
   }
   int mutationhere = 0;
   // Go through all GC-BC = CB + CC and add to the counters:
   for (long i = 0; i < CB_list.benutzt(); ++i) {
      // get the total number of GC-BCs at this point
      mutationhere = CB_list[i].n_mutation;
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      if (mutationhere > max_mutation) {
         mutationhere = max_mutation;
      }
      ++nbc[mutationhere];
   }
   for (long i = 0; i < CC_list.benutzt(); ++i) {
      // get the total number of GC-BCs at this point
      mutationhere = CC_list[i].n_mutation;
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      if (mutationhere > max_mutation) {
         mutationhere = max_mutation;
      }
      ++nbc[mutationhere];
   }
   // open the output file (append)
   ofstream mutation_histo;
   mutation_histo.open("mutation_histo.out", ofstream::app);
   for (int n = 0; n <= max_mutation; n++) {
      mutation_histo << time / 24. << "   " << time << "   "
                     << n << "   "
                     << nbc[n]
                     << "\n";
   }
   mutation_histo.close();
   //cerr << "Finished write_mutations.\n";
}
void cellman::write_Ig_class_output() {
  ofstream igoutput;
  igoutput.open("xsumout_ig.out", ofstream::app);
  igoutput << time << "   " << time / 24. << "     " 
	   << integrate_out[nIg_classes] << "  "
	   << integrate_out[nIg_classes] - integrate_out_old[nIg_classes] << "     ";
  for (int n = 0; n < nIg_classes; n++) {
    igoutput << integrate_out[n] << "  "
	     << integrate_out[n] - integrate_out_old[n] << "     ";
  }
  igoutput << "\n";
  igoutput.close();
  for (int n = 0; n <= nIg_classes; n++) {
    integrate_out_old[n] = integrate_out[n];
  }
}
void cellman::show_ag_distribution_on_FDC() {
  /* Histogram of remaining antigen in % of initial per FDC site */
  int resolution = 25;
  int free_ag_site[resolution+1] = {0};
  int tot_ag_site[resolution+1] = {0};
  for (int i = 0; i < FDC_list.benutzt(); i++) {
    FDC_list[i].add2histo_ag_per_site(free_ag_site, tot_ag_site, resolution);
  }
  // write ag_site to file
  ofstream aghisto;
  aghisto.open("ag_FDCsite_histo.out");
  aghisto.setf(ios::scientific, ios::floatfield);
  aghisto << "! bin : remaining frac of initial ag t=end : freq of total : free on FDC-sites \n";
  for (int i = 0; i <= resolution; i++) {
    aghisto << i << "   " << double(i) / double(resolution) 
	    << "   " << tot_ag_site[i] 
	    << "   " << free_ag_site[i] 
	    << "\n";
  }
  aghisto.close();
}
void cellman::write_final(const double &deltax) {
  /* write some data to output files after the last time step */
  // write into velocity files:
  //  if (outputfiles==0) {
  double average_v = 0;
  ofstream vmean;
  vmean.open("vmean.out");
  vmean.setf(ios::scientific, ios::floatfield);
  vmean << "! mean v (single moves) : #moves : mean v (within intervalls) : #moves : 0"
	<< " : mean v (intervalls+euklid) : # moves\n";
  
  long all_moves = 0;
  ofstream vdist;
  vdist.open("vdist.out");
  vdist.setf(ios::scientific, ios::floatfield);
  vdist << "! walked distance in dx : v in microns/min : "
	<< "# of occurrences : same per # of time steps\n";
  for (int i = 0; i <= v_resolution; i++) {
    all_moves += velocity[i];
    average_v += double (velocity[i]) * double (i) * deltax
      / (60. * dt * double (v_resolution));
  }
  average_v /= all_moves;
  vmean << average_v << "   " << all_moves << "   ";
  for (int i = 0; i <= v_resolution; i++) {
    vdist << i << "   " << double (i) * deltax / (60. * dt * double (v_resolution)) << "   "
	  << velocity[i] << "   "
	  << double (velocity[i]) / double (all_moves) << "\n";
  }
  vdist.close();
  
  vdist.open("vdistd.out");
  vdist.setf(ios::scientific, ios::floatfield);
  vdist << "! array-index : v in microns/min : # of occurrences : probability distribution"
	<< " : # of occs euklid : prob-dist euklid\n";
  
  all_moves = 0;
  long all_moves_euklid = 0;
  average_v = 0.;
  double average_v_euklid = 0.;
  for (int i = 0; i <= v_resolution; i++) {
    all_moves += delta_velocity[i];
    all_moves_euklid += delta_velocity_euklid[i];
    //    average_v+=double(delta_velocity[i])*(double(i)*delta_v+delta_v/2.);
    average_v += double (delta_velocity[i]) * double (i) * delta_v;
    if (i == 0) {
      average_v += double (delta_velocity[i]) * 0.25 * delta_v;
    }
    average_v_euklid += double (delta_velocity_euklid[i]) * double (i) * delta_v;
    if (i == 0) {
      average_v_euklid += double (delta_velocity_euklid[i]) * 0.25 * delta_v;
    }
  }
  /* Note, that if delta_v/2. is added here, this has to be done in the output below
   * as well as in the counting routine.
   */
  average_v /= all_moves;
  average_v_euklid /= all_moves_euklid;
  vmean << average_v << "   " << all_moves << "   " << 0.0 << "   " << average_v_euklid
	<< "   " << all_moves_euklid
	<< "   "
	<< "\n";
  for (int i = 0; i <= v_resolution; i++) {
    vdist << i << "   "
      //	 <<double(i)*delta_v+delta_v/2.<<"   "
	  << double (i) * delta_v << "   " << delta_velocity[i] << "   "
	  << double (delta_velocity[i]) / double (all_moves) << "   "
	  << delta_velocity_euklid[i] << "   "
	  << double (delta_velocity_euklid[i]) / double (all_moves_euklid) << "\n";
  }
  vdist.close();
  
  vmean.close();

  // Write fdc-encounters
  ofstream encount;
  encount.open("fdcencount.out");
  encount << "! Histogram for the number of FDC encounters in dying CC\n";
  encount << "! number of encounters & frequence of occurrence\n";
  for (int a = 0; a < cellCC::fdc_max_encounters; a++) {
    encount << a << "   " << cellCC::fdc_encounters[a] << "\n";
  }
  encount.close();

  ofstream muthisto("mutation_out_histo.out");
  for (int i = 0; i < max_mutation_bin; i++) {
    muthisto << i << "   " << mutation_frequency[i] << endl;
  }
  muthisto.close();
  
  cellCB::show_cummulative_number_of_divisions();
  cellCC::show_cummulative_ag_collection();
  cellCB::show_cummulative_mutation_prob();
  cellCB::show_cell_cycle_phase_duration();
  show_ag_distribution_on_FDC();
  
   //  }
}
void cellman::mk_beta1file(const long &ind, const char filename[30], const int columns) {
   bool target_open = false;
   char finput[30];
   double dat[columns + 1];
   long dati = -9;
   ofstream target;
   strcpy(finput, filename);
   strcat(finput, ".out");
   // cerr<<finput<<":\n";
   ifstream fff(finput);
   while (fff.eof() == 0) {
      // get one data line
      fff >> dati;
      // cerr<<"ind="<<ind<<", dati="<<dati<<"\n";
      for (int j = 1; j < columns; j++) {
         fff >> dat[j];
      }
      if (dati == ind) {
         if (target_open == false) {
            char foutput[30];
            strcpy(foutput, filename);
            strcat(foutput, cellbeta::beta_index);
            strcat(foutput, ".out");
            target.open(foutput);
            target_open = true;
         }
         for (int j = 1; j < columns; j++) {
            target << dat[j] << "  ";
         }
         target << "\n";
         dati = -9;
      }
   }
   fff.close();
   if (target_open) {
      target.close();
   }
}
void cellman::mk_single_beta_files(space &l) {
   cerr << "Write " << BETA_list.benutzt() << " betacells in single files ";
   /*
    * for (long i=0; i<BETA_list.benutzt(); i++) {
    * cerr<<".";
    * // position on the lattice:
    * long& ind=BETA_list[i].index;
    */

   for (long ind = 0; ind < l.pointnumber; ind++) {
      if (l.cellknot[ind].cell == BETA) {
         // long ind=l.cellknot[i].listi;
         cerr << ".";

         char inputfile[30];
         strcpy(inputfile, "beta_a");
         mk_beta1file(ind, inputfile, 36);
         strcpy(inputfile, "beta_b");
         mk_beta1file(ind, inputfile, 5);
         strcpy(inputfile, "beta_i");
         mk_beta1file(ind, inputfile, 22);
         strcpy(inputfile, "beta_r");
         mk_beta1file(ind, inputfile, 22);
         strcpy(inputfile, "beta_n");
         mk_beta1file(ind, inputfile, 7);
         strcpy(inputfile, "beta_g");
         mk_beta1file(ind, inputfile, 23);

         addchar(cellbeta::beta_index);
      }
   }
   cout << " done.\n";
}
void cellman::write_cell_specific() {
  cellCC::write_histograms();
}
// =======================================================
// Main routines =========================================
// =======================================================

void cellman::set_pars(Parameter &par, space &l) {
   cellCB::set_statics(time, par);
   cellCC::set_statics(time, par);
   cellFDC::set_statics(time, par);
}
short cellman::check_all(space &l, AffinitySpace &shape) {
   long n;

   // Testrun: consistency of l.knot[i].listi with cell_list[listi].index=i
   cout << "Run consistency checks of cell-list and lattice ... ";
   for (n = 0; n < CB_list.benutzt(); n++) {
      if (l.cellknot[CB_list[n].index].listi != n) {
         cout << "\nCB_list index=" << n << " while l.cellknot[CB_list[" << n
              << "].listi=" << l.cellknot[CB_list[n].index].listi << "\n";
         exit(1);
      }
      CB_list[n].check_connection(CB, n, l);
   }
   for (n = 0; n < TC_list.benutzt(); n++) {
      if (l.cellknot[TC_list[n].index].listi != n) {
         cout << "\nTC_list index=" << n << " while l.cellknot[TC_list[" << n
              << "].listi=" << l.cellknot[TC_list[n].index].listi << "\n";
         exit(1);
      }
   }
   for (n = 0; n < CC_list.benutzt(); n++) {
      if (l.cellknot[CC_list[n].index].listi != n) {
         cout << "\nCC_list index=" << n << " while l.cellknot[CC_list[" << n
              << "].listi=" << l.cellknot[CC_list[n].index].listi << "\n";
         exit(1);
      }
   }
   for (n = 0; n < TFR_list.benutzt(); n++) {
      if (l.cellknot[TFR_list[n].index].listi != n) {
         cout << "\nTFR_list index=" << n << " while l.cellknot[TFR_list[" << n
              << "].listi=" << l.cellknot[TFR_list[n].index].listi << "\n";
         exit(1);
      }
   }
   /*
    * Note, that the test of FDC is deleted because of a change the 4.8.2004:
    * I prefered to have connected FDCs on the FDC fragment-lists. This implies
    * different FDCs to occupy the same lattice point. On the lattice only one of them
    * is pointed to, which leads to an error detection in the following routine.
    * for (n=0; n<FDC_list.benutzt(); n++) {
    * if (l.knot[FDC_list[n].index].listi!=n) {
    *  cout<<"\nFDC_list index="<<n<<" while l.knot[FDC_list["<<n
    *      <<"].listi="<<l.knot[FDC_list[n].index].listi<<"\n";
    *  exit(1);
    * }
    * for (int m=0; m<FDC_list[n].volume; m++)
    *  if ((FDCtransparent==0 && l.knot[FDC_list[n].fragments[m]].FDClisti!=n) ||
    *      (FDCtransparent==1 && l.knot[FDC_list[n].fragments[m]].FDClisti==-1)
    *      ) {
    *  cout<<"\nFDC_list index="<<n<<" while l.knot[FDC_list["<<n
    *      <<"].fragments["<<m<<"]].FDClisti="
    *      <<l.knot[FDC_list[n].fragments[m]].FDClisti<<"\n";
    *  exit(1);
    * }
    * }
    */

   /*
    * long fdcn=0;
    * long fdcc=0;
    * for (long x=0; x<l.pointnumber; x++) {
    * if (l.knot[x].cell==FDC) ++fdcc;
    * if (l.knot[x].FDClisti!=-1) ++fdcn;
    * }
    * cout<<"FDC number: "<<fdcc<<" cell=FDC and "<<fdcn<<" FDClisti not -1!\n";
    */

   for (n = 0; n < OUT_list.benutzt(); n++) {
      if (l.cellknot[OUT_list[n].index].listi != n) {
         cout << "\nOUT_list index=" << n << " while l.cellknot[OUT_list[" << n
              << "].listi=" << l.cellknot[OUT_list[n].index].listi << "\n";
         exit(1);
      }
   }
   cout << "done.\n";
   fragment_consistency(l);
   check_listi(l);
   // check_TC_neighbours(l);
   shape.sum_check();
   for (n = 0; n < l.pointnumber; n++) {
      if ((l.knot[n].status == external) && (l.cellknot[n].cell != N_cells)) {
         cout << "At grid point " << n << " status is external but cell="
              << l.cellknot[n].cell << ".\n";
         exit(1);
      }
   }
   // cout<<"done.\n";

   return 0;
}
long cellman::get_GCvolume() {
   /* ### For the moment differences in cell volume are ignored.
    * It may be important for cells with more than one subunit to
    * include the "real" volume. However, now the total volume
    * occupied by cells is approximated by the total number of
    * cells.
    */
   long cb_plus_cc = CB_list.benutzt() + CC_list.benutzt();
   return cb_plus_cc;
}
/*
 * int CChowoften2619(long* x, long l) {
 * int n=0;
 * for (long i=0; i<l; i++) if (x[i]==2619) ++n;
 * if (n>1) {
 *  cout<<"Values in mCC are: ";
 *  for (long i=0; i<l; i++) cout<<i<<"->"<<x[i]<<";  ";
 *  cout<<"\n";
 * }
 * return n;
 * }
 */
//MSchips
/* #time : CB_list : CC_list : OUT_list : apop (total:self:normal) : out_in_GC*/
void cellman::BcAnalysis(double time, space &l , AffinitySpace &shape) {
    ofstream Bc_anal;
    int cb=0, cc=0, apo=0;
    int cb_self=0, cc_self=0, apo_self=0;
    int outingc_self=0, outingc=0;
    Bc_anal.open("BcAnal.out",ofstream::app);
    for (int b = 0; b < CB_list.benutzt(); b++) {
        if (CB_list[b].selfMutation) {++cb_self;}
        else ++cb;
    }
    for (int c = 0; c < CC_list.benutzt(); c++) {
        if (CC_list[c].state == apoptosis) {
            if (CC_list[c].selfMutation) {++apo_self;}
            else ++apo;
        }
        else {
            if (CC_list[c].selfMutation) {++cc_self;}
            else ++cc;
        }
    }

    for (int o = 0; o < OUT_list.benutzt(); o++) {
        if (OUT_list[o].selfMutation) { ++outingc_self; }
        else { ++outingc; }
    }

    double out_self = shape.get_sum_cell(soutSelf);
     double out = (shape.get_sum_cell(sout) - out_self);
//    Bc_anal << time << "  " << CB_list.benutzt() << "  "
//            << (CC_list.benutzt() - apo) << "  " <<
//               "\n";
    Bc_anal << time << " " << CB_list.benutzt() << " "
            << CC_list.benutzt() -(apo_self+apo) << " "
            << (out+out_self) << " "
            << (apo_self+apo) << " "
            << cb_self << " "
            << cc_self << " "
            << out_self << " "
            << apo_self << " "
            << cb << " "
            << cc << " "
            << out << " "
            << apo << " "
            << OUT_list.benutzt() << " "
            << outingc_self << " "
            << outingc << "\n";

    Bc_anal.close();
}
void cellman::BcAnalysis2(double time, space &l , AffinitySpace &shape) {
    int live_Bcells=0, CB_DZ=0, CB_LZ=0, liveCC_DZ=0, liveCC_LZ=0, apopt=0, apopt_DZ=0, apopt_LZ=0,
            self_CB=0, self_liveCC=0, self_apopt=0, CB_redeemed=0;
    double live_cell_ratio=0, cell_ratio=0, live_zone_ratio=0, zone_ratio=0;
    int out_inGC=0, self_out_inGC=0;
    double prop_selfOut_inGC=0, prop_selfOut=0, propGC_selfOut=0, cd138=0, cd138fr=0, TfrS_p_TfhS=0;
    double total_out_inGC=0;

    double ccl3ko_cells=0, ccl3ko_fraction=0;

    double select_frac=0;
    double Aselect_frac_self=0,Aselect_frac_ns=0, Bselect_self=0,Bselect_ns=0;

    double inddiffdel_tocb=0, inddiffdel_toout=0;

    double death_by_tfh=0,death_by_fdc=0;


    long LZ_sep = l.zone_separator;
    for (int b = 0; b < CB_list.benutzt(); b++) {
        ++live_Bcells;
        if (CB_list[b].selfMutation) {++self_CB; if (CB_list[b].redeemed) {cerr<<"mistakeinredemp";exit(1);}}
        else if (CB_list[b].redeemed) {++CB_redeemed;}
        if (CB_list[b].index<LZ_sep) {++CB_LZ;}
        else {++CB_DZ;}
        if (CB_list[b].ccl3KO) {++ccl3ko_cells;}
    }
    for (int c = 0; c < CC_list.benutzt(); c++) {
        if (CC_list[c].state!=apoptosis) {
            ++live_Bcells;
            if (CC_list[c].CD138) {
                ++cd138; inddiffdel_toout+=CC_list[c].individual_dif_delay;
            } else {
                if (CC_list[c].state==selected){
                    inddiffdel_tocb+=CC_list[c].individual_dif_delay;
                }
            }
            if (CC_list[c].selfMutation) {++self_liveCC;}
            if (CC_list[c].index<LZ_sep) {++liveCC_LZ;}
            else {++liveCC_DZ;}
            if (CC_list[c].ccl3KO) {
                ++ccl3ko_cells;
                //check
                if (exp_CCL3KO&& CC_list[c].nTFRcontacts>0) {cerr<<"should have no cont with tfr\n";exit(1);}
            }
            if (CC_list[c].state==selected) {
                ++select_frac;
                if (CC_list[c].selfMutation) {++Aselect_frac_self;++Bselect_self;} else {++Aselect_frac_ns;++Bselect_ns;}
            }
        } else {
            ++apopt;
            if (CC_list[c].selfMutation) {++self_apopt;}
            if (CC_list[c].index<LZ_sep) {++apopt_LZ;}
            else {++apopt_DZ;}
            if (CC_list[c].tfr_signal>=CC_list[c].tc_signal) {++TfrS_p_TfhS;}
            if (CC_list[c].nFDCcontacts==0) {++death_by_fdc;}else{++death_by_tfh;}
        }
    }
    TfrS_p_TfhS/=apopt;
    inddiffdel_tocb/=(select_frac-cd138);
    inddiffdel_toout/=(cd138);

    //ns_selBCs/selBC
    Aselect_frac_ns/=select_frac;
    //s_selBC/selBC
    Aselect_frac_self/=select_frac;

    select_frac/=(liveCC_DZ+liveCC_LZ);

    ofstream INDDIFFDEL;
    INDDIFFDEL.open("INDDIFFDELFILE.out",ofstream::app);
    INDDIFFDEL << time << " " << inddiffdel_tocb
                << " " << inddiffdel_toout << "\n";
    INDDIFFDEL.close();
//    cerr<< time << " " << inddiffdel_tocb
//        << " " << inddiffdel_toout << "\n";


    for (int o = 0; o < OUT_list.benutzt(); o++) {
        ++out_inGC;
        if (OUT_list[o].selfMutation) {
            ++self_out_inGC;
        }
    }

    double self_out = (shape.get_sum_cell(soutSelf));
    double out = (shape.get_sum_cell(sout));

    if (CC_list.benutzt()>0) {
        live_cell_ratio=double(CB_list.benutzt())/(liveCC_DZ+liveCC_LZ);
        cell_ratio=double(CB_list.benutzt()/CC_list.benutzt());
        live_zone_ratio=double(CB_DZ+liveCC_DZ)/double(CB_LZ+liveCC_LZ);
        zone_ratio=double(CB_DZ+liveCC_DZ+apopt_DZ)/double(CB_LZ+liveCC_LZ+apopt_LZ);
    } else {
        if (apopt>0) {
            cell_ratio=double(CB_list.benutzt())/double(apopt);
            zone_ratio=double(CB_DZ+liveCC_DZ+apopt_DZ)/double(CB_LZ+liveCC_LZ+apopt_LZ);
        }
    }
    if (out>0) {
        if (out_inGC>0) {prop_selfOut_inGC=double(self_out_inGC)/double(out_inGC);}
        prop_selfOut=self_out/out;
        if (out_inGC>0) {propGC_selfOut=double(out_inGC)/double(live_Bcells+out_inGC);} else {propGC_selfOut=0;}
    }
    if (live_Bcells>0) {
        cd138fr=cd138/double(live_Bcells);
        total_out_inGC=(cd138+out_inGC)/(live_Bcells+out_inGC);
    } else {cd138fr=0; total_out_inGC=0;}

    bool err=0;
    if ((CB_LZ+CB_DZ)!=shape.get_sum_cell(sCB)) {err=1;cerr<<"CB_list="<<CB_LZ+CB_DZ<<" sCB="<<shape.get_sum_cell(sCB)<<endl;}
    if ((liveCC_LZ+liveCC_DZ+apopt)!=shape.get_sum_cell(sCC)) {err=1;cerr<<"CC_list="<<liveCC_LZ+liveCC_DZ<<" sCC="<<shape.get_sum_cell(sCC)<<endl;}
    if ((apopt)!=shape.get_sum_cell(sCCapoptosis)) {err=1;cerr<<"apo_list="<<apopt<<" sapo="<<shape.get_sum_cell(sCCapoptosis)<<endl;}
//    if (self_liveCC!=shape.get_sum_cell(sCCself)) {err=1; cerr<<"counted self="<<self_liveCC<<" sccself="<<shape.get_sum_cell(sCCself)<<endl;}
    if (err) {exit(1);}

//    for (long n = 0; n < max; n++) {//LZ
//       li = l.cellknot[n].listi;
//       if ((l.cellknot[n].cell == CB) && (CB_list[li].index == n)) {

//       }
//       if ((l.cellknot[n].cell == CC) && (CC_list[li].index == n)) {}
//    }
//    for (long n = max; n < l.pointnumber; n++) {//DZ
//       li = l.cellknot[n].listi;
//    }

    ofstream BC_analysis_TFR;
//    time(1): live_Bcells(2): CB_DZ(3): liveCC_DZ(4): apopt_DZ(5):
//      CB_LZ(6): liveCC_LZ(7): apopt_LZ(8):
//      (self_CB+self_liveCC)(9): self_apopt(10):
//      out_inGC(11): self_out_inGC(12): out(13): self_out(14):
//      live_cell_ratio(15): live_zone_ratio(16): cell_ratio(17): zone_ratio(18):
//      propGC_selfOut(19): prop_selfOut_inGC(20): prop_selfOut(21): CB_redeemed(22):
//      cd138/liveBC(23) : TfrSplTfh(/apop) (24) : total_out_inGC(25)

    BC_analysis_TFR.open("BC_analysis_TFR.out",ofstream::app);
    BC_analysis_TFR << time << " " << live_Bcells << " "
                    << CB_DZ << " " << liveCC_DZ << " " << apopt_DZ << " "
                    << CB_LZ << " " << liveCC_LZ << " " << apopt_LZ << " "
                    << (self_CB+self_liveCC) << " " << self_apopt << " "
                    << out_inGC << " " << self_out_inGC << " " << out << " " << self_out << " "
                    << live_cell_ratio << " " << live_zone_ratio << " "
                    << cell_ratio << " " << zone_ratio << " "
                    << propGC_selfOut << " " << prop_selfOut_inGC << " " << prop_selfOut << " "
                    << CB_redeemed << " " << cd138fr << " " << TfrS_p_TfhS << " " << total_out_inGC << " "
                    <<cd138<<" "<< select_frac <<"\n";
    BC_analysis_TFR.close();


    ofstream asc_contact;
    double freq=0;
    if (OUT_list.benutzt()>0) {freq=double(asc_got_contact/OUT_list.benutzt());} else {freq=0;}
    asc_contact.open("asc_contact.out",ofstream::app);
    asc_contact << time << " " << asc_got_contact
                << " " << freq << "\n";
    asc_contact.close();

    ccl3ko_fraction = double((ccl3ko_cells/live_Bcells)*100);
    ofstream ccl3_file;
    ccl3_file.open("ccl3_file.out",ofstream::app);
    ccl3_file << time << " " << ccl3ko_cells
                << " " << ccl3ko_fraction << "\n";
    ccl3_file.close();

    ofstream FRAC_SEL;
    FRAC_SEL.open("FRAC_SEL.out",ofstream::app);
    FRAC_SEL << time << " " << Aselect_frac_ns
                << " " << Aselect_frac_self << " " << Bselect_ns
                << " " << Bselect_self << "\n";
    FRAC_SEL.close();

    ofstream CAUSE_OF_DEATH;
    CAUSE_OF_DEATH.open("CAUSE_OF_DEATH.out",ofstream::app);
    CAUSE_OF_DEATH << time << " " << death_by_fdc
                << " " << (death_by_fdc/apopt) << " " << death_by_tfh
                << " " << (death_by_tfh/apopt) << "\n";
    CAUSE_OF_DEATH.close();
}
void cellman::writeZPositions(double time, space &l){
    long double zTfh[l.prodimvec[0]], zTfr[l.prodimvec[0]],
            zCB[l.prodimvec[0]], zCC[l.prodimvec[0]],
            zCD138[l.prodimvec[0]], zASC[l.prodimvec[0]],
            xTfh[l.prodimvec[0]], xTfr[l.prodimvec[0]],
            xCB[l.prodimvec[0]], xCC[l.prodimvec[0]],
            xCD138[l.prodimvec[0]], xASC[l.prodimvec[0]],
            yTfh[l.prodimvec[0]], yTfr[l.prodimvec[0]],
            yCB[l.prodimvec[0]], yCC[l.prodimvec[0]],
            yCD138[l.prodimvec[0]], yASC[l.prodimvec[0]];

    int cd138_n=0;
    for (int i=0; i<l.prodimvec[0]; i++) {
        zTfh[i]=0; zTfr[i]=0;
        zCB[i]=0; zCC[i]=0;
        zCD138[i]=0; zASC[i]=0;
        xTfh[i]=0; xTfr[i]=0;
        xCB[i]=0; xCC[i]=0;
        xCD138[i]=0; xASC[i]=0;
        yTfr[i]=0; yTfh[i]=0;
        yCB[i]=0; yCC[i]=0;
        yCD138[i]=0; yASC[i]=0;
    }
    for (int b = 0; b < CB_list.benutzt(); b++) {
        long pZ=l.knot[CB_list[b].index].x[2], pX=l.knot[CB_list[b].index].x[0], pY=l.knot[CB_list[b].index].x[1];
        ++zCB[pZ]; ++xCB[pX]; ++yCB[pY];
    }

    for (int b = 0; b < CC_list.benutzt(); b++) {
        long pZ=l.knot[CC_list[b].index].x[2], pX=l.knot[CC_list[b].index].x[0], pY=l.knot[CC_list[b].index].x[1];
        ++zCC[pZ]; ++xCC[pX]; ++yCC[pY];
        if (CC_list[b].CD138) {
            ++cd138_n; ++zCD138[pZ]; ++xCD138[pX]; ++yCD138[pY];
        }
        ofstream cd138_xyz;
        cd138_xyz.open("cd138_xyz.out",ofstream::app);
        for (int i=0; i<l.prodimvec[0] ;i++) {
            cd138_xyz << time
                 << " " << pX << " " << pY
                    << " " << pZ << "\n";
        }
        cd138_xyz.close();
    }

    for (int b = 0; b < OUT_list.benutzt(); b++) {
        long pZ=l.knot[OUT_list[b].index].x[2], pX=l.knot[OUT_list[b].index].x[0], pY=l.knot[OUT_list[b].index].x[1];
        ++zASC[pZ]; ++xASC[pX]; ++yASC[pY];
        ofstream ASC_xyz;
        ASC_xyz.open("ASC_xyz.out",ofstream::app);
        for (int i=0; i<l.prodimvec[0] ;i++) {
            ASC_xyz << time
                 << " " << pX << " " << pY
                    << " " << pZ << "\n";
        }
        ASC_xyz.close();
    }

    ofstream tfh_xyz;
    tfh_xyz.open("tfh_xyz.out",ofstream::app);
    for (int b = 0; b < TC_list.benutzt(); b++) {
        long pZ=l.knot[TC_list[b].index].x[2], pX=l.knot[TC_list[b].index].x[0], pY=l.knot[TC_list[b].index].x[1];
        ++zTfh[pZ]; ++xTfh[pX]; ++yTfh[pY];
        tfh_xyz<<time<<" "<<pX<<" "<<pY<<" "<<pZ<<"\n";
    }
    tfh_xyz.close();

    ofstream tfr_xyz;
    tfr_xyz.open("tfr_xyz.out",ofstream::app);
    for (int b = 0; b < TFR_list.benutzt(); b++) {
        long pZ=l.knot[TFR_list[b].index].x[2], pX=l.knot[TFR_list[b].index].x[0], pY=l.knot[TFR_list[b].index].x[1];
        ++zTfr[pZ]; ++xTfr[pX]; ++yTfr[pY];
        tfr_xyz<<time<<" "<<pX<<" "<<pY<<" "<<pZ<<"\n";
    }
    tfr_xyz.close();

    ofstream Zpos;
    Zpos.open("Zpos.out",ofstream::app);
    for (int i=0; i<l.prodimvec[0] ;i++) {
        Zpos << time << " " << i
             << " " << zCB[i]/CB_list.benutzt() << " " << zCC[i]/CC_list.benutzt()
                << " " << zTfh[i]/TC_list.benutzt() << " " << zTfr[i]/TFR_list.benutzt()
                << " " << zCD138[i]/cd138_n << " " << zASC[i]/OUT_list.benutzt() << "\n";
    }
    Zpos.close();
//    ofstream PTfr;
//    PTfr.open("PTfr.out",ofstream::app);
//    for (int i=0; i<l.prodimvec[0] ;i++) {
//        PTfr << time << " " << i
//             << " " << xTfr[i]/TFR_list.benutzt() << " " << yTfr[i]/TFR_list.benutzt()
//                << " " << zTfr[i]/TFR_list.benutzt() << "\n";
//    }
//    PTfr.close();
}

void cellman::writeAffinity2(double time, space &l, AffinitySpace &shape) {
    short size=4;
    double aff_mean[size], aff_sd[size];
    for (int i=0; i<size; i++) {
        aff_mean[i]=0;aff_sd[i]=0;
    }
    shape.meanSD_out_affinity(aff_mean,aff_sd);

    //0 : 1 : 2 --> sout (self : total : nonSelf)
    //3 --> total bcs
//    time(1): sout--Self[m:sd](2:3): total[m:sd](4:5): nonself[m:sd](6:7)
//      liveBCs[m:sd](8:9)
    ofstream BCaff_TFR;
    BCaff_TFR.open("BCaff_TFR.out",ofstream::app);
    BCaff_TFR << time << " ";
    for (int i=0; i<size; i++) {
        BCaff_TFR << aff_mean[i] << " " << aff_sd[i] << " ";
    }
    BCaff_TFR << "\n";
    BCaff_TFR.close();
}
void cellman::write_info(double time, double tfhsig, double pmhc, double affinity, short status, bool self) {
    ofstream outCCInfo;
    outCCInfo.open("outCCInfo.out", ofstream::app);
    outCCInfo << time << " " <<
                 tfhsig << " " <<
                 pmhc << " " <<
                 affinity << " " <<
                 status << " " << self << endl;
    outCCInfo.close();
}
void cellman::writeAffinity(double time, space &l, AffinitySpace &shape) {
    int CC=0, CB=0, apo=0, Out=0, Self=0, SelfOut=0;
    int tot = CB_list.benutzt()+CC_list.benutzt()+OUT_list.benutzt();
    int totout=OUT_list.benutzt(), totBC=CB_list.benutzt()+CC_list.benutzt();

    // affinity variables
    int aff_bins = 10;
    double aff_max = 1;
    double aff_min = 0;
    double aff_level[aff_bins];
    int aff_array[aff_bins];
    double stepAff = (aff_max-aff_min)/aff_bins;
    for (int j = 0; j < aff_bins; ++j) {
        aff_level[j] = aff_max-(j*stepAff);
        aff_array[j] = 0;
    }
    int bin=0;

    double avg_tot=0, sd_tot=0,
            avg_out=0, sd_out,
            avg_self=0, sd_self,
            avg_selfOut=0, sd_selfOut,
            avg_bc=0, sd_bc;
    double affsTOT[tot], affsOUT[totout], affsBC[totBC],
            affsSelf[totout], affsSelfOut[totout];
    for (int j = 0; j < tot; j++) {
        affsTOT[j] = 0;
        if (j<totout) {
            affsOUT[j] = 0;
            affsSelf[j] = 0;
            affsSelfOut[j] = 0;
        }
        if (j<totBC) affsBC[j] = 0;
    }

    int all_index = -1, indexN=-1, bc_index=-1, selfIndex=-1, selfOutIndex=-1;

    for (int c = 0; c < CB_list.benutzt(); c++) {
        all_index+=1; //0
        indexN+=1; //0
        bc_index+=1; //0
        ++CB; //1
        bin = 0;
        while ((bin < aff_bins) && (shape.best_affinity_norm(CB_list[c].pos_ss) < aff_level[bin])) {
            ++bin;
        }
        if (bin == aff_bins) {
            --bin;
        }
        if (CB_list[c].selfMutation) {
            selfIndex+=1;
            affsSelf[selfIndex]=shape.best_affinity_norm(CB_list[c].pos_ss);
            ++Self;
        }
        ++aff_array[bin];
        affsTOT[all_index] = shape.best_affinity_norm(CB_list[c].pos_ss);
        affsBC[bc_index] = shape.best_affinity_norm(CB_list[c].pos_ss);
    }
    for (int c = 0; c < CC_list.benutzt(); c++) {
            if (CC_list[c].state == apoptosis) {
                ++apo;
            } else {
                indexN+=1;
                all_index+=1;
                bc_index+=1;
                ++CC;
                bin = 0;
                while ((bin < aff_bins) && (shape.best_affinity_norm(CC_list[c].pos_ss) < aff_level[bin])) {
                    ++bin;
                }
                if (bin == aff_bins) {
                    --bin;
                }
                if (CC_list[c].selfMutation) {
                    selfIndex+=1;
                    affsSelf[selfIndex]=shape.best_affinity_norm(CC_list[c].pos_ss);
                    ++Self;
                }
                ++aff_array[bin];
                affsTOT[all_index] = shape.best_affinity_norm(CC_list[c].pos_ss);
                affsBC[bc_index] = shape.best_affinity_norm(CC_list[c].pos_ss);
             }
        }
    for (int o = 0; o < OUT_list.benutzt(); o++) {
        indexN+=1;
        all_index+=1;
        ++Out;
        bin = 0;
        while ((bin < aff_bins) && (shape.best_affinity_norm(OUT_list[o].pos_ss) < aff_level[bin])) {
           ++bin;
        }
        if (bin == aff_bins) {
            --bin;
        }
        if (OUT_list[o].selfMutation) {
            selfIndex+=1;
            selfOutIndex+=1;
            affsSelf[selfIndex]=shape.best_affinity_norm(OUT_list[o].pos_ss);
            affsSelfOut[selfOutIndex]=shape.best_affinity_norm(OUT_list[o].pos_ss);
            ++Self; ++SelfOut;
        }
        ++aff_array[bin];
        affsTOT[all_index] = shape.best_affinity_norm(OUT_list[o].pos_ss);
        affsOUT[o] = shape.best_affinity_norm(OUT_list[o].pos_ss);
    }

    int NS=tot-apo, NSN=(CB+CC);
    double sum1=0, sum2=0, sum3=0, //aff of BC // aff of OUT cells // aff total
            sum4=0, sum5=0; //aff of Mem //aff of PC

    for (int j = 0; j < NS; j++) { //CB+CC+OUT-apo > OUT /// < / <= ???
        avg_tot += affsTOT[j]; //sum of the affinities of all BCs
        if (j<totout) avg_out+=affsOUT[j]; //sum of the affinities of all BCs in OUTlist
        if (j<selfIndex) avg_self+=affsSelf[j]; //sum of the affinities of all self GCBC cells
        if (j<selfOutIndex) avg_selfOut+=affsSelfOut[j]; //sum of the affinities of all self Out cells
        if (j<NSN) avg_bc+=affsBC[j]; //sum of the affinities of all BCs in CBs and CCs
    }

    /// --> the 'else' statements should not even be necessary
    /// if any cell number is 0 => the affs* will not be updated and will just stay as initialised, i.e.=0
    /// => avg_* value will also stay as initialised, i.e.=0
    if (NS>0) avg_tot/=NS; else avg_tot=0;
    if (NSN>0) avg_bc/=NSN; else avg_bc=0;
    if (Self>0) {avg_self/=Self;} else avg_self=0;
    if (totout>0) {
        avg_out/=totout;
        if (SelfOut>0) {avg_selfOut/=SelfOut;} else avg_selfOut=0;
    } else {
        avg_out=0;
    }

    for (int j = 0; j < NS; j++) {
        if (j<NSN) sum1+=( (affsBC[j]-avg_bc) * (affsBC[j]-avg_bc) );
        if (j<totout) sum2+=( (affsOUT[j]-avg_out) * (affsOUT[j]-avg_out) );
        sum3+=((affsTOT[j]-avg_tot) * (affsTOT[j]-avg_tot));
        if (j<Self) sum4+=( (affsSelf[j]-avg_self) * (affsSelf[j]+avg_self) );
        if (j<SelfOut) sum5+=( (affsSelfOut[j]-avg_selfOut) * (affsSelfOut[j]-avg_selfOut) );
    }
    if (NS>0) sd_tot=sqrt(sum3/NS); else sd_tot=0;
    if (NSN>0) sd_bc=sqrt(sum1/NSN); else sd_bc=0;
    if (totout>0) {
        sd_out=sqrt(sum2/totout);
        if (Self>0) {sd_self=sqrt(sum4/Self);} else sd_self=0;
        if (SelfOut>0) {sd_selfOut=sqrt(sum5/SelfOut);} else sd_selfOut=0;
    } else {
        sd_out=0;
    }

     ofstream writeAff;
     writeAff.open("writeAffinity.out",ofstream::app);

     // time(1) : CB(2) : CC(3) : apo(4) :
     // OUT(5) : Memory(6) : Plasma(7) :
     // Affinities (avg : sd) :
     // CB+CC+OUT(8:9) : CB+CC(10:11) : OUT(12:13) :
     // memory (14:15) : plasma (16:17)

     writeAff << time << " "
//                         << CB << " " << CC << " "
//                         << apo << " " << Out << " "
//                         << Memory << " " << Plasma << " "
                         << avg_tot << " " << sd_tot << " "
                         << avg_bc << " " << sd_bc << " "
                         << avg_out << " " << sd_out << " "
                         << avg_self << " " << sd_self << " "
                         << avg_selfOut << " " << sd_selfOut << "\n";
     writeAff.close();
}

void cellman::contactscheck(double time, space &l) {
    // # time : #TFR : #freeTFR : #TFR in nonselfCC : #TFR in selfCC

//    ofstream tfrposition;
//    tfrposition.open("tfrposition.out",ofstream::app);

    int tfrcc=0, cctfr=0, cctfr_self=0, cctfr_nonself=0;
    for (int tt=0; tt<TFR_list.benutzt(); tt++) {
        if (TFR_list[tt].state==TFR_CCcontact) {++tfrcc;}
//        tfrposition << time << " " << TFR_list[tt].state << " "
//             << l.knot[TFR_list[tt].index].x[0] << " "
//             << l.knot[TFR_list[tt].index].x[1] << " "
//             << l.knot[TFR_list[tt].index].x[2] << "\n";
    }
//    tfrposition.close();
    int tfh=0, freetfh=0;
    std::array<long, 6> numTFHcc;
    numTFHcc = std::array<long, 6> { 0,0,0,0,0,0 };

    for (int th=0; th<TC_list.benutzt(); th++) {
        ++tfh;
        if (TC_list[th].state==TC_CCcontact) {
            int NN=(TC_list[th].get_n_boundCC() -1);
            ++numTFHcc[NN];
        } else { ++freetfh; }
    }

    ofstream tfhcontact;
    tfhcontact.open("tfh.out",ofstream::app);
    tfhcontact << time << " " << tfh << " " << freetfh << " ";
    for (int ii=0; ii<6; ii++) {
        tfhcontact << numTFHcc[ii] << " ";
    }
    tfhcontact << "\n";
    tfhcontact.close();


    for (int cc=0; cc<CC_list.benutzt(); cc++) {
//        if (CC_list[cc].ccInd_TFRbound!=-1) {
        if (CC_list[cc].state==TFRcontact) {
            ++cctfr;
            if (CC_list[cc].selfMutation) ++cctfr_self;
            else ++cctfr_nonself;
        }
    }
    int freetfr=TFR_list.benutzt()-cctfr;


    ofstream tfrcontcheck;
    tfrcontcheck.open("tfr.out",ofstream::app);

    int bound_to_asc=0;
    bound_to_asc=tfrcc-cctfr;
    if (bound_to_asc<0) {bound_to_asc=-bound_to_asc;}

    tfrcontcheck << time
                 << " " <<TFR_list.benutzt()
                 << " " << freetfr
                 << " " << cctfr_nonself
                 << " " << cctfr_self
                 << " " << bound_to_asc
                 << "\n";

    if (tfrcc!=cctfr) {
        //WARNING THIS IS NOT VALID ANYMORE:
        //tfr HAVE ONLY ONE 'BINDING' STATE FOR CC AND ASC!!!
//        cerr<<" tfrcc "<<tfrcc <<" cctfr "<<cctfr<<endl;
//        exit(1);
    }

    tfrcontcheck.close();

}


void cellman::time_step(short int ss_save, short int day_save,
                        space &l,
                        sigs &s,
                        AntibodyDyn &Ab,
                        AffinitySpace &shape,
                        ofstream &ana) {
  // cerr<<"Start timestep ... (time="<<time<<")\n";
  long n;
  
  /* Restrict global counting of attributed numbers of division to
   * events after injection of anti-DEC205-OVA when injection is active: */
  if (inject_antiDEC205OVA && time < inject_antiDEC205OVA_t0) {
    for (int k=0; k <= cellCB::max_n_of_divisions; k++) 
      { cellCB::cummulative_attributed_n_of_divisions[k] = 0; }
  }
  
  // if (time>=TRACK::TRACKFROM+0.5*dt && time<TRACK::TRACKFROM+1.5*dt)
  if ((time >= TRACK::TRACKFROM - 0.5 * dt) && (time < TRACK::TRACKFROM + 0.5 * dt)) {
    make_fluorescent(l);
  }
  if (tmx.get_tamoxifen_action(time) && tmx.MHC) { recombine_cells(time); }

   /* ### If new Ag are supposed to be generated on the run, this has to be initiated here:
    *  Add the new Ag with AffinitySpace::add_new_Antigen(..).
    *  Distribute the new Ag onto the FDCs (unclear how to combine with pre-existing ones),
    *  where there is no routine prepared for this.
    *  The FDC-class is prepared to handle additional Ags that appear later.
    *  The number of additional Ags is limited to a factor (static int n_Antigen_dim_factor)
    *  times the initial Ags, which is set in the FDC-class in cellthis.cpp.
    *  All this should be handled by an extra routine cellman::add_antigen() (not yet there).
    */
   Ab.check4new_antigen(shape, ana);        // eventually, open new set of ab-affinity-bins
   Ab.ask2inject_antibody(time, dt, ana);   // eventually, inject antibodies
   // DEC205 attribution is done at the time creation of new BCs if def_DEC205_t0 is negative
   if (def_DEC205 && (time >= def_DEC205_t0 - 0.5 * dt) && (time < def_DEC205_t0 + 0.5 * dt)) {
      attribute_DEC205();
   }
   if (inject_antiDEC205OVA && 
       (time >= inject_antiDEC205OVA_t0 - 0.5 * dt) && 
       (time < inject_antiDEC205OVA_t0 + 0.5 * dt)) 
     { do_inject_antiDEC205OVA(l, shape, ana); }
   if (photoactivation && (time >= photoactivation_t0 - 0.5 * dt)
       && (time < photoactivation_t0 + 0.5 * dt)) 
     { make_photoactivation(l); }

   // Set up an array of grid positions that are to be recalculated at the
   // end because cell movements of single node cells were suppressed by
   // contact inhibition.
   dynarray<long> redo(100, 1, 0);

   // Put the signals if they are set by hand
   // ++++++++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++
   // This option for precalculation of glucose is switched on and off in signals.h
   // Don't touch this call here!
   // Put glucose for the betacells -> define the function in signals.C mk_const_signal(...)
   if ((show_mode == islet) && (s.signal_use[glucose] == 1)) 
     { s.mk_const_signal(glucose, time * 3600., cellbeta::glucose_rest); }
   // ++++++++++++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++++

   // Berechne schon mal die Zelltypen, die nicht zu anderen Zelltypen werden koennen
   // where the sequence is randomized in calc_xxx(...)

   // cerr<<"vor OUT ... \n";
   calc_OUT(l, s, shape, redo);

   // differentiate external plasma cells to antibody producing PM and make them produce:
   shape.PM_differentiate(soutext, soutextproduce, dt);
   // ##### note that this might be done always -- not only for usage of antibody bins.
   if (use_antibody_bins) 
     { Ab.produce_antibodies_outside(shape); }
   /// Philippe : I don't understand how antibodies are produces in the other case (not bins)

   // cerr<<"vor FDC ... \n";
   calc_FDC(l, s, Ab, shape);
   // cerr<<"nach FDC.\n";
   calc_stroma(s);

   // Jetzt die Typen, die sich ineinander umwandeln koennen
   // Randomize the sequence of cells of all those types beforehand in order to
   // avoid double calculation of the same cell.
   long mCBlang = CB_list.benutzt();
   long mCB[mCBlang];
   random2_sequence(mCB, mCBlang);
   // cerr<<"nach random2_sequence(CB).\n";
   for (n = 0; n < mCBlang; n++) {
      // if (l.knot[CB_list[int(mCB[n])].index].listi==-1) {
      // cout<<"cell="<<int(mCB[n])<<" total="<<CB_list.benutzt(); exit(1); }
      mCB[n] = CB_list[int (mCB[n])].index;
      // cerr<<mCB[n]<<"; ";
      // for (int mm=0; mm<CB_list[l.knot[mCB[n]].listi].volume; mm++)
      // cout<<CB_list[l.knot[mCB[n]].listi].fragments[mm]<<", ";
      // cout<<"\n";
   }

   /// Philippe : putting that in the Calc_XX function ?
   long mCClang = CC_list.benutzt();
   long mCC[mCClang];
   random2_sequence(mCC, mCClang);
   for (n = 0; n < mCClang; n++) {
      mCC[n] = CC_list[int (mCC[n])].index;
   }
   // cerr<<"nach random CC.\n";

   long mTClang = TC_list.benutzt();
   long mTC[mTClang];
   random2_sequence(mTC, mTClang);
   for (n = 0; n < mTClang; n++) {
      mTC[n] = TC_list[int (mTC[n])].index;
   }
   // if (time>380) cout<<"nach random TC.\n";

   //MSchips
   long mTFRlang = TFR_list.benutzt();
   long mTFR[mTFRlang];
   random2_sequence(mTFR, mTFRlang);
   for (n = 0; n < mTFRlang; n++) {
      mTFR[n] = TFR_list[int (mTFR[n])].index;
   }

   // Berechne diese
   // cerr<<"vor CC ... "<<mCClang<<"\n";
   calc_CC(mCC, mCClang, l, s, Ab, shape, redo);

   // cerr<<"vor CB ... "<<mCBlang<<"\n";
   calc_CB(mCB, mCBlang, ss_save, l, s, shape, redo);
   // cerr<<"nach CB. \n";

   // cerr<<"vor TC ... "<<mTClang<<"\n";
   calc_TC(mTC, mTClang, l, s, shape, redo);
   // if (time>380) cout<<"nach TC. \n";

   calc_TFR(mTFR, mTFRlang, l, s, shape, redo);
//    cerr<<"endof tfr\n";


   long mBETAlang = BETA_list.benutzt();
   long mBETA[mBETAlang];
   if (cellbeta::RANDOMISE_SEQUENCE) {
      random2_sequence(mBETA, mBETAlang);
   } else {
      for (n = 0; n < mBETAlang; n++) {
         mBETA[n] = n;
      }
   }
   for (n = 0; n < mBETAlang; n++) {
      mBETA[n] = BETA_list[int (mBETA[n])].index;
   }
   // cerr<<"nach random BETA.\n";
   // cerr<<"vor BETA ... "<<mBETAlang<<"\n";
   calc_BETA(mBETA, mBETAlang, l, s, redo);
   // cerr<<"nach BETA. \n";

   if (allow_exchange) {
      // Now process the contact inhibited movements as saved in redo[]:
      // Note that no further random sequence is needed as the cell-positions where already
      // written in a random order into redo[].
      // cerr<<"in allow_exchange ...\n";
      while (redo.benutzt() > 0) {
         // process the first cell
         long partner = retry_movement(redo, l);
         // cerr<<"partner="<<partner<<"\n";
         // partner contains -1 if no movement was performed
         // or the grid-index of the partner if the cells were exchanged.
         // Delete the processed cell and eventually its partner from the list
         redo.erase_jump(0);
         if (partner > -1) {
            long where_partner = redo.find(partner);
            if (where_partner >= 0) {
               redo.erase_jump(where_partner);
            }
            // Note that processing the partner anyway is not an option, because before the
            // position in redo[] would have to be adopted to the novel position.
         }
      }
      // cerr<<"After redo!\n";
   }

   // Testroutine zur Deposition von Signalen in einer Zelle
   /*
    * if (time>=0.0 && time<0.6) {
    * l.knot[1860].signal_new[sig_differ2CC]=1000;
    * l.signal_total[sig_differ2CC]+=1000;
    * }
    */

   // Diffusion der Signalmolekuele:
   // cerr<<"Vor signals ... \n";
   s.signal_diffuse(l);
   // cerr<<"nach signals.\n";

   /* BrdU staining is done after the time step calculations to make the number
    * of cells stained compatible with the number of S-phase BCs in the cellcycle_phases.out
    * file, which is written right below.
    */
   if (do_inject_BrdU && 
       (time >= tnext_inject_BrdU - 0.5 * dt) && 
       (time < tnext_inject_BrdU + 0.5 * dt)) { 
     label_BrdU(false); 
     //     show_BrdU_at_time();
     //     show_BCstate_in_BrdUpositive();
   }
   //MSchips
//   BcAnalysis(time,l,shape);
   BcAnalysis2(time,l,shape);
//   writeAffinity(time, l, shape);
   writeAffinity2(time, l, shape);
   contactscheck(time, l);
   if ((time >= 24 && time<(24+dt))
           || (time >= 48 && time<(48+dt))
           || (time >= 72 && time<(72+dt))
           || (time >= 120 && time<(120+dt))
           || (time >= 144 && time<(144+dt))
           || (time >= 168 && time<(168+dt))
           || (time >= 216 && time<(216+dt))
           || (time >= 240 && time<(240+dt))
           || (time >= 312 && time<(312+dt))
           || (time >= 360 && time<(360+dt)) ) {
       writeZPositions(time,l);
   }

   if ((show_mode != islet) && (outputfiles == 0)) {
      // corrs<<time<<"   "<<2.-CB_list[0].performed2aimed_move_now<<"\n";
      /*
       * The previous variant makes no sense if more cells are around! The plotted variable is
       * defined for every cell separately (as it should be), and therefore an
       * average value over all cells would have to be calculated here (the new
       * version below works equally for one or more cells:
       */
      double p2a = 0.;
      if (CB_list.benutzt() > 0) {
         for (int i = 0; i < CB_list.benutzt(); i++) {
            p2a += 2. - CB_list[i].performed2aimed_move_now;
         }
         p2a /= double (CB_list.benutzt());
      } else {
         p2a = 999.;
      }
      corrs << time << "   " << p2a << "\n";
      /*
       ### The other point: Is it more reasonable to write every hour or
       ##################something only? Then the corrs-writing has to be shifted into
       ##################the if (ss_save==1) block below.
       */
   }
   // cout<<"t="<<time*60.<<" vol="<<CB_list[0].volume<<"\n";
   // if (time>=144. && time<144.+dt) FDC_list[0].set_antigen_amount(time,dt);

   if ((show_mode != islet) && (ss_save == 1)) {
      mk_cell_sum(l, s, Ab, shape);    // makes no changes to shape.
      shape.to_multiag_files(time);
      shape.to_ssfiles(time);          // changes shape.oldsum_cell[xx] !!
      write_Ig_class_output();  // writes integrate_out[] to xsumout_ig.out
      count_dec205_ova_positive();
      show_ag_collected();
      if (do_inject_BrdU && time >= tfirst_inject_BrdU - 0.5 * dt) {
          show_BrdU_at_time();
          show_BCstate_in_BrdUpositive();
      }
      showFDCfailure(time);
      show_cycle_phases(l.zone_separator);
      show_K_signal_intensity_pMHC();
      cellCB::show_number_of_divisions(time, ndivtime);
      cellCC::monitor_TFHintensity.write_mean2file(time);
      double muta_av = 0, muta_sd = 0;
      get_average_mutation_prob(muta_av, muta_sd);
      cellCB::show_mutation_prob(time, muta_av, muta_sd, mutation_time);
      //MSchips
      cellCB::selfBC_number_of_divisions(time, selfBC_ndivtime);
      double muta_avS = 0, muta_sdS = 0;
      get_average_mutation_prob_SELF(muta_avS, muta_sdS);
      cellCB::show_mutation_prob_SELF(time, muta_avS, muta_sdS, Selfmutation_time);
      if (TFR_mode!=-1) {
          cellCC::show_freq_of_TfrCont(time, freq_nonself, freq_self);
      }

   }
   if ((show_mode != islet) && (day_save == 1)) {
     cellCC::update_histograms();
   }
   if ((time >= 288.) && (time < 288. + dt)) {
      CC2CBratio = double (CC_list.benutzt()) / double (CB_list.benutzt());
   }
   if ((time > 72.) && (t_1st_under100 < 0.) && (CB_list.benutzt() + CC_list.benutzt() < 100)) {
      t_1st_under100 = time;
   }
   if (checkit == 2) {
      check_all(l, shape);
      //  shape.sum_check();
      Ab.check4consistency(shape, ana);
   }

   // compare to volume experiments
   if (show_mode != islet) {
      kinetics.check_GCvolume(get_GCvolume(), time);
   }

   // if (time>=TRACK::TRACKUNTIL-0.5*dt && time<TRACK::TRACKUNTIL+0.5*dt)
   if ((time > TRACK::TRACKUNTIL - 0.5 * dt) && (time <= TRACK::TRACKUNTIL + 0.5 * dt)) {
      stop_fluorescent();
   }

   //MS newout for correlation-like plots
   double record_gap=1*24; //every 'record_gap' hours
   if (time>double(first_record+time_window)
              && time<=double(first_record+time_window+dt)) {
       first_record+=record_gap; //if first record is a global then i can use that to jump to th next recording time
   }

   // cout<<"glucose(t="<<time<<")="<<l.get_signal_total(glucose)<<"\n";
   // cout<<"oxygen(t="<<time<<")="<<l.get_signal_total(oxygen)<<"\n";

   // if (time>394) cout<<"end of timestep.\n\n";
}
