// #include <math.h>
#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include "cellman.h"

/*
 * ====================================================================================
 * ================= output files =====================================================
 * ====================================================================================
 * value || function
 * ------------------------------------------------------------------------------------
 * 0        write all standard files
 * 1        write less standard files
 * 2        write a minimal set of standard files (no signal files)
 * 3        write as 2 plus files for analysis of two-photon data about cell motility
 * ====================================================================================
 */

int run_t(Parameter &p,
          cellman &c,
          space &l,
          sigs &g,
          AntibodyDyn &Ab,
          AffinitySpace &s,
          ofstream &ana) {
   // Die Zeiten:
   // +++ OPTION: time gap between global writings (standard 24 hours)
   double t_gap = 24.;
   // end OPTION
   // Zeitschritt-Koordinate
   long int n = 0;
   // Gesamtzahl der Zeitschritte
   long int nmax = long ((p.Value.tmax - p.Value.tmin) / p.Value.deltat + 0.5);

   // Sonderfall bei tmin<0:
   long switch_toGCphase = -1;
   if (p.Value.tmin < 0) {
      switch_toGCphase = long ((0. - p.Value.tmin) / p.Value.deltat + 0.5) + 1;
   }
   // mutation is now switched on in the CB-proliferation routine in cellman.C: delete here
   //   long switch_toMutation;
   // switch_toMutation=long((p.Value.Start_Mutation-p.Value.tmin)/p.Value.deltat+0.5)+1;

   /// Philippe : I have the impression that these variables are not required
   /// because set_differentiation is called at each time step
   long switch_toDifferentiation;
   if (cellCB::SMOOTH_DIFFERENTIATION) {
      switch_toDifferentiation = 1;
   } else {
      switch_toDifferentiation
         = long ((p.Value.Start_Differentiation - p.Value.tmin) / p.Value.deltat + 0.5) + 1;
   }
   long switch_toOutput;
   if (cellCC::SMOOTH_DIFFERENTIATION) {
      switch_toOutput = 1;
   } else {
      switch_toOutput = long ((p.Value.StartOutput - p.Value.tmin) / p.Value.deltat + 0.5) + 1;
   }
   // Initialize data output:
   // shape space and lattice writing depending on p.Value.ToFileStep
   short int write_x = 0;
   // Initialisierung des Dateinamens fuer den Output:
   suffix tnr = "0000";
   if (p.Value.tmin > 0) {
      for (int push = 0; push < int (double (p.Value.tmin) / t_gap + 0.5); push++) {
         addchar(tnr);
      }
   }

   // various (continuous in time) characteritic data (automatic intervalls)
   short write_t = 0, write_day = 0;
   // write each hour (number of time steps for one hour):
   long write_int = long (nmax / (p.Value.tmax - p.Value.tmin));
   if (write_int == 0) {
      write_int = 1;
   }

   // For betacell output:
   if (p.Value.show_mode == islet) {
      c.show_BETA();
   }
   long write_beta = long (p.betaValue.dt_output / (3600. * p.Value.deltat) + 0.5);
   long beta_count = 0;
   if ((p.Value.show_mode == islet) && (cellbeta::LOCAL_FILES == false) && (write_beta == 0)) {
      write_beta = 1;
      cout << "WARNING: beta-output-dt is too small for the external timestep!\n";
   }

   // number of time steps for 12 and 24 hours
   long write12_int = long (12 * nmax / (p.Value.tmax - p.Value.tmin));
   long write24_int = long (24 * nmax / (p.Value.tmax - p.Value.tmin));
   short int inject = 0;

   // write zone-files each day
   short int write_z = 0;
   suffix tnz = "0000";
   if (p.Value.tmin > 0) {
      for (int push = 0; push < int (double (p.Value.tmin) / t_gap + 0.5); push++) {
         addchar(tnz);
      }
   }

   // for signal tests:
   if (sigs::TEST_MODE) {
      g.write_TEST(c.time);
   }

   // cerr<<"Before c.set_pars ...\n";
   // G.tofile(tnr,"b",logit);
   // G.tofile(tnr,"a",logit);
   // G.toTfiles();
   c.set_pars(p, l);

   // cerr<<"Before for (n=1; ...) { ...\n";
   // Dynamik: Durchlauf der Zeitschritte
   for (n = 1; n <= nmax; n++) {
      ++beta_count;
      // should I write to files? if so: value =1 (=0 otherwise)
      write_x = (n == (n / p.Value.ToFileStep) * p.Value.ToFileStep);
      write_t = ((n == (n / write_int) * write_int) || (n == nmax));
      inject = (n == (n / write12_int) * write12_int);
      write_day = (n == (n / write24_int) * write24_int);
      write_z = ((n == (n / (int (t_gap) * write_int)) * int (t_gap) * write_int) || (n == nmax));

      c.time += p.Value.deltat;
      // Documentation during run:
      // cerr<<"Timestep "<<n<<" von "<<nmax<<" = "<<c.time<<" h :\n";

      if ((n == switch_toGCphase)
          ||   //	  n==switch_toMutation ||
          (n == switch_toDifferentiation)
          || (n == switch_toOutput)) {
         cout << "Adjust parameters ... ";
         c.set_pars(p, l);  /// Philippe : Why is this function not called at each time point ???
         cout << "\n";
      }

      c.time_step(write_t, write_day, l, g, Ab, s, ana);

      if (write_x == 1) {
         addchar(tnr);
         c.movie << "$1\"" << tnr << ".gif\" ";
         // c.mkgifs<<"ppmtogif "<<"xy"<<tnr<<".ppm > xy"<<tnr<<".gif\n";
         c.xfiles(tnr, l);
         c.write_mutations(c.time);
         s.write_gcbc_hamming(c.time);
         s.write_gcbc_affinity(c.time);
         g.write_files(tnr, false);
         if (sigs::TEST_MODE) {
            g.write_TEST(c.time);
         }

         if (p.Value.safety_checks == 1) {
            c.check_all(l, s);
            Ab.check4consistency(s, ana);
         }
      }
      if ((c.outputfiles < 2) && (write_z == 1) && (c.show_mode != islet)) {
         addchar(tnz);
         c.zone_files(tnz, l);
      }
      if (inject == 1) {
         // every 12 hours
         if (n >= switch_toDifferentiation - 1) {
            c.inject_Ki67();
         }
         c.write_log_bc_aff(s);
      }
      if ((cellbeta::LOCAL_FILES == false) && (beta_count == write_beta)) {
         c.show_BETA();
         beta_count = 0;
      }
   }
   c.write_final(l.dx);
   c.write_cell_specific();
   cout << "Close movie-files ... \n";
   c.movie << "> $1\".gif\"\n";
   c.movie.close();
   // c.mkgifs.close();
   // cout<<"done.\n";
   c.check_all(l, s);
   Ab.check4consistency(s, ana);

   return 1;
}
int main(int argn, char * * argument) {
  // cout << sequence::testeAffinityFunctions(50, 8, 25, seqAffNorm);

  //  cout<<"argn="<<argn<<" "<<argument[0]<<" "<<argument[1]<<" "<<argument[2]<<"\n";
  char cross[2] = "#";
  if (argn == 1) {
    argument[1] = cross;
  }
  //  cout<<"argn="<<argn<<" "<<argument[0]<<" "<<argument[1]<<" "<<argument[2]<<"\n";

  // Definiere die Parameter-Variablen
  Parameter par;
  int done = 0; 
  if (argn <= 2) {
    done = par.wahl(argument[1],true,true);
  } else if (argn == 3) {
    par.wahl(argument[1],false,true);
    done = par.wahl(argument[2],true,false);
  } else {
    cout << "ERROR: more than 2 command line parameters not supported --> exit\n";
    exit(1);
  }

   // Falls Parameter eingelesen wird gestartet
   if (done == 5) {
      // Zufallsgenerator initialisieren
      if (par.Value.ini_random == -1) {
         srand(time(NULL));
      } else { srand(par.Value.ini_random); }
      for (long n = 0; n < par.Value.ini_random; n++) {
         int x = irandom(100);
         double y = drandom();
         y = double (x) + y;
      }

      // Eroeffnung der Analyse Dateien:
      ofstream analyze("ana_ini.out");
      // +++ OPTION: Use an additional analysis-file with special results
      // The file is designed for several runs stored in the same file
      ofstream sel_analyze("sel_analyze.out", std::ios_base::app);
      // The stored values are explained in a separate text-file:
      ofstream sel_ana_text("sel_readme.out");
      if (par.Value.FDClength >= par.Value.dx) {
         sel_analyze << 12 * int (par.Value.FDClength / par.Value.dx + 0.5)
            * par.Value.FDCnumber;
      } else {
         sel_analyze << 4 * par.Value.FDCnumber;
      }
      sel_analyze << "   " << par.Value.FDCnumber << "   " << par.Value.FDClength << "   ";
      if (par.Value.ag_per_FDC < 0) {
         sel_analyze << "-1   ";
      } else {
         sel_analyze << par.Value.ag_per_FDC << "   ";
      }
      sel_analyze << par.Value.ag_saturation_FDC << "   " << par.Value.totalTC << "   ";
      if (par.Value.TC_CC_selection == 1) {
         sel_analyze << par.Value.TC_time << "   " << par.Value.TC_rescue_time << "   ";
      } else {
         sel_analyze << "0   -1   ";
      }
      if (par.Value.mk_ab < 0.) {
         sel_analyze << "0   ";
      } else {
         sel_analyze << par.Value.mk_ab << "   ";
      }
      sel_analyze << par.Value.ic_k_on << "   " << par.Value.ic_k_off << "   "
                  << par.Value.ag_threshold << "   "
                  << log(2.) / par.Value.tolight << "   " << par.Value.mksignal << "   ";

      sel_ana_text
         << "Definition of the values stored in the file sel_analyze.out:\n"
         << "------------------------------------------------------------\n"
         << "\n"
         << "The values shown in a single line are the following\n"
         << "from left to right:\n\n"
         << "column 01: Number of FDC sites\n"
         << "column 02: Number of FDC\n"
         << "column 03: Length of FDC arms in microns\n"
         << "column 04: Antigen presentation per FDC (in Binding-Quanta)\n"
         <<
         "column 05: Antigen amount per FDC fragment below which binding probability is linearly reduced\n"
         << "column 06: Number of TC\n"
         << "column 07: Duration of CC-TC contact (hours)\n"
         << "column 08: Minimum polarisation duration of TC towards CC for selection\n"
         << "column 09: Antibody production of OUTPUT in Mol per hour and cell\n"
         << "column 10: IC=Ab-Ag dynamics: k_on\n"
         << "column 11: IC=Ab-Ag dynamics: k_off\n"
         << "column 12: Threshold of Ag concentration for CC binding in Mol\n"
         << "column 13: CB2CC differentiation time in hours\n"
         << "column 14: Differentiation signal production in Quanta per FDC and hour\n"
         // and now results:
         << "column 15: First time of population below 100 cells (hours)\n"
         << "column 16: Number of cells at the end\n"
         << "column 17: Fraction of high affinity cells at the end (average CB+CC)\n"
         << "column 18: Number of produced OUTPUT cells\n"
         << "column 19: Average binding probability of produced OUTPUT cells\n"
         << "column 20: End of dark zone (hours)\n"
         << "column 21: Ratio of CC2CB at t=288 hours\n"
         << "column 22: Ratio of output cells at day 12 to day 6\n"
         << "column 23: Standard deviation from measured GC kinetics\n";
      sel_ana_text.close();
      // end OPTION

      // lattices definieren
      space l(par, analyze);
      sigs s(l, par, analyze);

      // Before starting, a few settings in the cell-class need to be done,
      // which are needed in derived classes:
      immunoglobulin_class::load_matrix(par, analyze);
      cell::set_statics(par,analyze);
      cellCB::set_statics(par, l, analyze);
      cellCC::set_statics(par, l, analyze);
      cellTC::set_statics(par, l, analyze);
      cellOUT::set_statics(par, l, analyze);
      cellFDC::set_statics(par, l, analyze);
      cellbeta::set_statics(par, l, analyze);
      AntibodyDyn::set_statics(par, analyze);
      arupProtein::set_statics(par, analyze);
      //MSchips
      cellTFR::set_statics(par, l, analyze);

      // Shapespace definieren
      AffinitySpace * newShape = NULL;
      if (par.Value.use_sequence_space == 1) {
         newShape = (AffinitySpace*) new sequenceSpace(par,
                                                       analyze);
      } else if (par.Value.use_sequence_space == 0) {
         newShape = (AffinitySpace*) new SS(par,analyze);
      } else if (par.Value.use_arup_space == 1) {
         newShape = (AffinitySpace*) new arupSpace(par,analyze);
      }
      if ((newShape == NULL) || (par.Value.use_sequence_space * par.Value.use_arup_space > 0)) {
         cerr << "ERR: conflicting options for use_sequence_space and use_arup_space" << endl;
      }

      // Antibody definition
      AntibodyDyn Ab(par, *newShape, analyze);

      // define cells
      cellman c(par, l, s, *newShape, analyze);

      // show initital BrdU and pre-loaded antigen FACS plot (in vitro BC for Cell Reports 2012)
      c.show_ag_BrdU(true);

      analyze << "\nGerminal Center initialisiert!\n";
      analyze << "Start calculation ...\n";
      cout << "Start calculation ...\n";

      // ==================================================================
      // ==================================================================
      // Zeitschleife durchlaufen:
      run_t(par, c, l, s, Ab, *newShape, analyze);
      // ==================================================================
      // ==================================================================

      cout << " ... end of calculation.\n\n";
      analyze << " ... end of calculation.\n\n";
      analyze << "Data analysis:\n";

      if (c.outputfiles == 0) {
         // if at least one cell is left do this analysis
         if (c.CB_list.benutzt() > 0) {
            // deformation parameter:
            long nmoves = 0;
            c.CB_list[0].frac_average = 0.;
            for (long n = 0; n < c.CB_list.benutzt(); ++n) {
               nmoves += c.CB_list[n].n_moves;
               c.CB_list[0].frac_average += (c.CB_list[n].n_moves * c.CB_list[n].alpha_mean);
            }
            c.CB_list[0].frac_average /= double (nmoves);
            analyze << "Average fraction of considered free neighbors for movement = "
                    << c.CB_list[0].frac_average
                    << "\n";
         }

         // If no CB is left add one artificially to allow access to the data
         if (c.CB_list.benutzt() == 0) {
            cellCB tmpcell;
            c.CB_list.add(tmpcell);
         }

         /* OPTION: MOVE_ANALYSIS
          * analyze <<"Fraction of successful movement tries = "
          * <<double(c.CB_list[0].n_move_done)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * analyze <<"Fraction of movements on itself = "
          * <<double(c.CB_list[0].n_move_self)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * analyze <<"Fraction of movements on itself of back fragments = "
          * <<double(c.CB_list[0].n_move_self_back)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * analyze <<"Fraction of movement tries removed by previous movement = "
          * <<double(c.CB_list[0].n_move_removed)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * analyze <<"Fraction of movement tries that are forbidden = "
          * <<double(c.CB_list[0].n_move_forbidden)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * analyze <<"Fraction of movement tries that are forbidden for back fragments= "
          * <<double(c.CB_list[0].n_move_forbidden_back)/double(c.CB_list[0].n_try_move)
          * <<"\n";
          * double factor=-1.0;
          * if (c.CB_list[0].n_try_move!=0
          * && c.CB_list[0].n_move_self!=0
          * && c.CB_list[0].n_move_forbidden!=0) {
          * double denom=double(c.CB_list[0].n_move_done)/double(c.CB_list[0].n_try_move)
          * +(1.0-double(c.CB_list[0].n_move_self_back)/double(c.CB_list[0].n_move_self))
          * double(c.CB_list[0].n_move_self)/double(c.CB_list[0].n_try_move)
          ******+(1.0-double(c.CB_list[0].n_move_forbidden_back)/double(c.CB_list[0].n_move_forbidden))
          * double(c.CB_list[0].n_move_forbidden)/double(c.CB_list[0].n_try_move);
          * if (denom>0.)
          * factor=(1.0-double(c.CB_list[0].n_move_removed)/double(c.CB_list[0].n_try_move))/denom;
          * }
          * analyze<<"Correction factor for diffusion at the end = ";
          * if (factor==-1.0) analyze<<"error\n";
          * else analyze<<factor<<"\n";
          */// end of OPTION: MOVE_ANALYSIS
         analyze << "Ratio of performed to aimed move per time step (average over "
                 << c.CB_list[0].n_directed_moves
                 << " moves) = " << c.CB_list[0].performed2aimed_move << "\n";
      }

      analyze << "Number of CBs at the end = " << c.CB_end << "\n";
      c.CB_average /= 48.;
      analyze << "Average number of CBs during the last 2 days = " << c.CB_average << " +- "
              << sqrt(c.CB_variance / 48. - c.CB_average * c.CB_average) << "\n";
      double _OUT_haffinity, _OUT_steepness, _CB_haffinity, _CC_haffinity;
      newShape->getStatistics(_OUT_haffinity, _OUT_steepness, _CB_haffinity, _CC_haffinity);
      analyze << "High affinity of CBs at the end = " << 100. * _CB_haffinity << "%\n";
      analyze << "High affinity of CCs at the end = " << 100. * _CC_haffinity << "%\n";
      analyze << "High affinity of OUTs at the end = " << 100. * _OUT_haffinity << "%\n";
      analyze << "Total number of OUTs at the end = " << newShape->get_sum_cell(sout) << "\n";
      analyze << "Steepness of OUT at times 288h/144h = " << _OUT_steepness << "\n";
      analyze << "Effectivity index out/mid_CB = " << double (newShape->get_sum_cell(sout))
         / c.CB_average << "\n";
      analyze << "End of dark zone at " << c.t_dark_end << " h\n";
      if (cellFDC::ic_calculations > 0) {
         double ic_errs
            = double (cellFDC::ab_sign_errors + cellFDC::ag_sign_errors
                      + cellFDC::ic_sign_errors)
              / double (cellFDC::ic_calculations);
         analyze << "Fraction of errors in mk_immune_complex(..): " << ic_errs << "\n";
         analyze << "  ab: " << double (cellFDC::ab_sign_errors)
            / double (cellFDC::ic_calculations)
                 << ", ag: " << double (cellFDC::ag_sign_errors)
            / double (cellFDC::ic_calculations)
                 << ", ic: " << double (cellFDC::ic_sign_errors)
            / double (cellFDC::ic_calculations) << ".\n";
         if (ic_errs > 0.05) {
            cout << "Note large number of errors in cellFDC::mk_immune_complex(...)!\n";
         }
      }
      analyze << "... end of data analysis.\n";
      analyze.close();

      // +++ OPTION: results of optional analysis
      sel_analyze << c.t_1st_under100 << "   " << newShape->get_sum_cell(sCB)
         + newShape->get_sum_cell(sCC) << "   ";
      sel_analyze << (double (newShape->get_sum_cell(sCB)) * _CB_haffinity
                      + double (newShape->get_sum_cell(sCC)) * _CC_haffinity)
         / double (newShape->get_sum_cell(sCB) + newShape->get_sum_cell(sCC)) << "   "
                  << newShape->get_sum_cell(sout) << "   " << _OUT_haffinity << "   "
                  << c.t_dark_end << "   "
                  << c.CC2CBratio << "   " << _OUT_steepness << "   ";
      // if (c.n_chi_values>1) sel_analyze<<c.chi_sum/(double(c.n_chi_values));
      double sigma2 = c.kinetics.get_sigma2_GCvolume();
      if (sigma2 >= 0.) {
         sel_analyze << sigma2;
      } else {
         sel_analyze << "*";
      }
      sel_analyze << "\n";
      sel_analyze.close();
      // end OPTION.

      // +++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++++
      // Write all signals to a file
      bool writeallsignals = false;
      if (writeallsignals) {
         cout
            <<
            "Write signals to signal.out that can be used as input if copied to signal.sig.\n";
         s.write_all_signals(cell::chemo_max, cell::chemo_steep, cell::chemo_half);
      } else {
         cout << "Do not write signals to signal.out that could be used as input signal.sig.\n";
      }
      // +++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++++

      c.close_files();
      newShape->close_files();

      c.trackdata.Write_files();
      c.trackmutations.save_to_file();
      // activate this for brainbow analysis
      //c.trackmutations.analysis();
      // The following was used for the analysis in Cell Reports 2012 Fig. S1
      c.show_BrdU();
      c.show_ag_loaded();
      c.show_cell_in_BrdU();
      c.show_ag_BrdU(false);  // FACS: ag versus BrdU

      if ((c.show_mode == islet) && (cellbeta::LOCAL_FILES == false)) {
         c.mk_single_beta_files(l);
      }
   }
   cerr << "HYPHASMA terminated.\n";
   return 0;
}
