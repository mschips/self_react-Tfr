#include "ss.h"
#include <math.h>
#include "random.h"

SSpoint::SSpoint() { }
void SSpoint::init(Parameter &par) {
   external_cell = false;
   ab_producer = false;
   k = new long[par.Value.DimShapeSpace];
   // Es gibt doppelt soviele Nachbarn wie die Dimension des shapespace
   nn = new long[2 * par.Value.DimShapeSpace];
   index = 0;
   short int i;
   for (i = 0; i < par.Value.DimShapeSpace; i++) {
      k[i] = 0;
   }
   for (i = 0; i < 2 * par.Value.DimShapeSpace; i++) {
      nn[i] = 0;
   }
   for (i = 0; i < number; i++) {
      n_cell[i] = 0;
   }
}
SSpoint::~SSpoint() {
   delete[] k;
   delete[] nn;
}
char operator ==(const SSpoint &a, const SSpoint &b) {
   char back = 1;
   if (a.index != b.index) {
      back = 0;
   }
   return back;
}
char operator !=(const SSpoint &a, const SSpoint &b) { return !((a == b) == 1); }

// =============================================================

SS::SS() { }
double SS::get_sum_cell(cells celltype) { return sum_cell[celltype]; }
double SS::get_oldsum_cell(cells celltype) { return oldsum_cell[celltype]; }
long SS::getMutation(long from_position) {
   long int nnn = -1;
   while (nnn == -1) {
      nnn = this->ssp[from_position].nn[irandom(2 * this->dim)];
   }
   return nnn;
}
SS::SS(Parameter &par, ofstream &ana) {
   long int n;
   ana << "Initialize Shapespace fields ... \n";

   // +++++++++++++++++++++++++++++++++++++++++++++++
   // Put in parameter file
   pm_differentiation_rate = log(2.) / par.Value.pm_differentiation_time;
   //   pm_differentiation_rate=log(2.)/24.;
   // +++++++++++++++++++++++++++++++++++++++++++++++

   // Range of space
   dim = par.Value.DimShapeSpace;
   PointsPerDim = par.Value.SSRangePerDim;
   PointsTotal = par.Value.SSStates;
   measure = par.Value.metrik;
   width2 = par.Value.GammaGauss * par.Value.GammaGauss;
   amplitude = par.Value.amplitudeGauss;
   ana << "  Gauss affinity weight: amplitude=" << amplitude << "; width^2=" << width2 << "\n";
   /*cout<<"\n  dim="<<dim<<"  PointsPerDim="<<PointsPerDim
    *  <<"  PointsTotal="<<PointsTotal<<"  measure="<<measure
    *  <<"  width="<<width<<"\n";*/
   // Analysis vars
   OUT_haffinity = 0.;
   OUT_steepness = 0.;
   CB_haffinity = 0.;
   CC_haffinity = 0.;
   // Zunaechst alle Felder initialisieren
   ssp = new SSpoint[par.Value.SSStates];
   for (int nn = 0; nn < par.Value.SSStates; nn++) {
      ssp[nn].init(par);
   }
   SSpoint p;
   p.init(par);
   for (n = 0; n < PointsTotal; n++) {
      //       p.index=n;
      ssp[n].index = n;
      get_koord(n, p.k);
      // Dann die Nachbarn bestimmen
      // fix_neighbors(n,p.nn);
      for (short i = 0; i < dim; i++) {
         ssp[n].k[i] = p.k[i];
      }
   }
   // Dann die Nachbarn bestimmen
   for (n = 0; n < PointsTotal; n++) {
      fix_neighbors(n);
   }
   ana << " done.\n";

   // Antigen-pool definieren
   ana << "Initialize Antigen pool ...\n ";
   n_Antigen = par.Value.APeakNumber;
   Antigen = new long[n_Antigen];
   for (n = 0; n < n_Antigen; n++) {
      //   for (n=0; n<Antigen.lang(); n++) {
      Antigen[n] = par.Value.takeA[n];
      if (Antigen[n] == -1) {
         // short int schonda=0;
         long x[dim];
         // ana << "PointsPerDim="<<int(PointsPerDim)<<"\n";
         ana << "  Random coordinates=(";
         for (int j = 0; j < dim; j++) {
            x[j] = irandom(int (PointsPerDim));
            ana << x[j] << "  ";
         }
         ana << "): ";

         Antigen[n] = Index(x);
         // Comment:
         ana << "Antigen-Peak: Antigen[" << n << "]= " << Antigen[n] << "\n";
         // schonda=0;
         // for (j=0; j<n; j++) { if (AIndex[n]==AIndex[j]) schonda=1; }
         // if (schonda==1) --i;
      } else {
         ana << "  Fixed Antigen-Peak: Antigen[" << n << "]= " << Antigen[n] << "\n";
      }
   }
   ana << "... done.\n";

   // initialise ab_producers to none:
   external_cells.reserve(5000);
   external_cells.clear();
   ab_producers.reserve(5000);
   ab_producers.clear();

   // define pool of seeder cells
   ana << "Initialize Seeder cell pool ... \n";
   n_Seeder = par.Value.totalBss;
   Seeder = new long[n_Seeder];
   int n_fixed_seeders = 0;
   while (par.Value.takeB[n_fixed_seeders] != -1) { ++n_fixed_seeders; }
   if (n_fixed_seeders == 0 && par.Value.fixedfraction >= 0) {
     cerr << "WARNING!!! fraction of clones from the list is " 
	  << par.Value.fixedfraction 
	  << " but no elements are provided in the list.\n"
	  << "           Reset the fraction to -1 !\n\n";
     par.Value.fixedfraction = -1;
   }
   for (n = 0; n < n_Seeder; n++) {
     if (par.Value.fixedfraction < 0) {
       // use each CB value from the parameter file once
       Seeder[n] = par.Value.takeB[n]; 
     } else { 
       if ( double(n)/double(n_Seeder) < par.Value.fixedfraction ) {
	 // use the CB values from the par file to account for <fixedfraction> seeder-positions
	 Seeder[n] = par.Value.takeB[irandom(n_fixed_seeders)];
       }
       else { Seeder[n] = -1; } // which means attribute random SS-position below
     }
     if (Seeder[n] == -1) {
       // short int schonda=0;
       ana << "  Random coordinates --> ";
       Seeder[n] = get_random_SS_position(par.Value.min_seeder_dist, par.Value.max_seeder_dist);
       ana << "Seeder cell: Seeder[" << n << "]=  " << Seeder[n] << "\n";
       // schonda=0;
       // for (j=0; j<n; j++) { if (AIndex[n]==AIndex[j]) schonda=1; }
       // if (schonda==1) --i;
     } else {
       ana << "  Fixed seeder cell : Seeder[" << n << "]= " << Seeder[n] << "\n";
     }
     //cout << "Seeder[" << n << "] = " << Seeder[n] << endl;
   }
   ana << "... done.\n";

   // Output-Dateien eroeffnen:
   ana << "Define shape space output files ...";
   for (short j = 0; j < SSlogs; j++) {
      sum_cell[j] = 0.;
      oldsum_cell[j] = 0.;
   }
   logdata[sCB].open("ssumcb.out");
   logdata[sCB] << "! Summe aller Centroblasten im SS\n";
   logdata[sCC].open("ssumcc.out");
   logdata[sCC] << "! Summe aller Centrocyten im SS\n";
   logdata[sCCunselected].open("ssumccun.out");
   logdata[sCCunselected] << "! Summe aller unselektierten Centrocyten im SS\n";
   logdata[sCCcontact].open("ssumccco.out");
   logdata[sCCcontact] << "! Summe aller Antigen-selektierten Centrocyten im SS\n";
   logdata[sCCselected].open("ssumccse.out");
   logdata[sCCselected] << "! Summe aller positiv selektierten Centrocyten im SS\n";
   logdata[sCCapoptosis].open("ssumapo1.out");
   logdata[sCCapoptosis] << "! Summe aller gestorbenen Zellen auf dem Gitter im SS\n";
   logdata[sallapoptosis].open("ssumapo2.out");
   logdata[sallapoptosis] << "! Summe aller jemals gestorbenen Zellen im SS : increment\n";
   // Lasse FDC und Tcell weg!
   logdata[sout].open("ssumout.out");
   logdata[sout] << "! Summe aller erzeugter Output Zellen im SS "
                 << ": DEC205+ : fraction of DEC205+\n";
   logdata[soutext].open("ssumoute.out");
   logdata[soutext] << "! Summe aller ausgeworfener Output Zellen im SS\n";
   logdata[soutextproduce].open("ssumoutep.out");
   logdata[soutextproduce] << "! Summe aller ausgeworfener Output Zellen im SS, "
                           << "die Antikoerper produzieren\n";
   logdata[total].open("ssum.out");
   logdata[total] << "! Summe aller Zellen im SS\n";

   logmeanaff.open("saffin.out");
   logmeanaff << "! Mean affinity of CB+CC : CB : CC : out\n";
   loghighaff.open("shaffin.out");
   loghighaff << "! Fraction of >30% affinity CB:CC:out\n";
   log048aff.open("s048aff.out");
   log048aff << "! Fraction of <40% affinity CB:CC:out : 40-80% : >=80%\n";

   logdiversity.open("diversity.out");
   logdiversity << "! Number of different encoded antibodies present in the GC\n"
                << "! time : total in GC : CB : CC : output in GC : output external\n";

   // additional files with emphasis on multiple antigens:
   ofstream saffin_ag("saffin_ags.out");
   saffin_ag << "! time :  mean affinity of all CB+CC to each antigen(columns)\n";
   saffin_ag.close();
   ofstream saffin_t10("saffin_ags_t10.out");
   saffin_t10 << "! time : mean affinity of CB+CC with aff>0.1 to each antigen(columns)\n";
   saffin_t10.close();
   ofstream saffin_t20("saffin_ags_t20.out");
   saffin_t20 << "! time : mean affinity of CB+CC with aff>0.2 to each antigen(columns)\n";
   saffin_t20.close();
   ofstream saffin_out_ag("saffin_out_ags.out");
   saffin_out_ag << "! time :  mean affinity of accumulated output"
                 << " to each antigen(columns)\n";
   saffin_out_ag.close();
   ofstream saffin_out_t10("saffin_out_ags_t10.out");
   saffin_out_t10 << "! time : mean affinity of accumulated output"
                  << " with aff>0.1 to each antigen(columns)\n";
   saffin_out_t10.close();
   ofstream saffin_out_t20("saffin_out_ags_t20.out");
   saffin_out_t20 << "! time : mean affinity of accumulated output"
                  << " with aff>0.2 to each antigen(columns)\n";
   saffin_out_t20.close();
   ofstream shaffin_ag("shaffin_ags.out");
   shaffin_ag << "! time : fraction of all CB+CC with affinity >0.3 to each antigen(columns)"
              << " : sum of all fractions\n";
   shaffin_ag.close();
   ofstream shaffin_ag_t80("shaffin_ags_t80.out");
   shaffin_ag_t80 << "! time : fraction of all and number of CB+CC with affinity>0.8"
		  <<" to each antigen(columns) : sum of all fractions\n";
   shaffin_ag_t80.close();
   ofstream shaffin_out_ag("shaffin_out_ags.out");
   shaffin_out_ag << "! time : fraction of accumulated output"
                  << " with affinity >0.3 to each antigen(columns) : sum of all fractions\n";
   shaffin_out_ag.close();
   ofstream cross_gcr("cross_reactivity_gcr.out");
   cross_gcr << "! time : cross-reactivity = "
             << "[mean affinity of CB+CC (all : aff>0.1 : >0.2) to all antigens]\n";
   cross_gcr.close();
   ofstream cross_out("cross_reactivity_out.out");
   cross_out << "! time : cross-reactivity = "
             << "[mean affinity of accumulated output (all : aff>0.1 : >0.2) to all antigens]\n";
   cross_out.close();
   // add a file for a histogram of GC-BCs versus Hamming distance to optimal clone
   ofstream gcbc_hamming("gcbc_hamming.out");
   gcbc_hamming << "! time[days] : time[hr] : Hamming distance : # of GC-BC nearest Ag "
                << ": # of GC-BC to mean Ag\n";
   gcbc_hamming.close();
   ofstream gcbc_affinity("gcbc_affinity.out");
   gcbc_affinity << "! time[days] : time[hr] : affinity-min : # of GC-BC nearest Ag "
                 << ": # of GC-BC to mean Ag\n";
   gcbc_affinity.close();

   ana << " done\n";
   ana << "Shape space READY!\n";
}
SS::~SS() {
   // cout<<"in ~SS()\n";
   delete[] ssp;
   delete[] Antigen;
   delete[] Seeder;
}
void SS::close_files() {
   cout << "Close files in shape space ... ";
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) {
         logdata[i].close();
      }
   }
   logmeanaff.close();
   loghighaff.close();
   logdiversity.close();
   log048aff.close();
   cout << "done.\n";
}
// =================================================================
// Routinen fuer Nummerierung, Initialisierung:
// =================================================================

int SS::get_n_Antigen() {
   return n_Antigen;
}
long int SS::Index(long int * k) {
   long int i = 0;
   for (int n = 0; n < dim; n++) {
      i += k[n] * long (pow(double (PointsPerDim), n));
   }
   return i;
}
void SS::get_koord(long int wo, long int * k) {
   long int reduce, p;
   for (int n = dim - 1; n >= 0; n--) {
      p = long (pow(double (PointsPerDim), n));
      reduce = wo / p;
      k[n] = reduce;
      wo -= reduce * p;
   }
}
double SS::euklid(long int * k, long int * l) {
   // double middle=PointsPerDim/2.;
   double distance = 0.;
   double xdist = 0.;
   for (int i = 0; i < dim; i++) {
      // xdist=double(fabs(double(k[i]-l[i]))); // fabs is double and fabs is not necessary
      xdist = k[i] - l[i];
      // if (par.Value.CyclicBorder==1)
      //   if (xdist>middle) xdist=PointsPerDim-xdist;
      distance += xdist * xdist;
   }
   return distance;
}
int SS::intN_mutation(long * k, long * l) {
   int distance = 0.;
   for (int i = 0; i < dim; i++) {
      distance += fabs(k[i] - l[i]);
   }
   return distance;
}
int SS::best_hamming_distance(long &n) {
   /* returns the distance of shape space point <n> to the nearest antigen
    * always uses the Hamming distance (irrespective of <measure>)
    */
   // First determine the nearest antigen:
   int distance = 1.0e+09; // stupid large number
   int thisdistance = 0;
   // +++++++++++++++++++ OPTION ++++++++++++++++++++++++++++++++++++++++++++
   // This line is wrong (see WARNING in write_gcbc_hamming(...))
   // and is used in order to show the affinity to the first antigen only
   // which might make sense when in an experiment only one antigen was tested in ELISA.
   //for (int nAg = 0; nAg < 1; nAg++) {
   // In the general case use this one:
   //cerr << "[n=" << n << ":ssp[n].k=(";
   //for (int jj=0; jj<dim; jj++) { cerr << ssp[n].k[jj] << ","; }
   //cerr << ");";
   for (int nAg = 0; nAg < n_Antigen; nAg++) {
      // +++++++++++++++ end OPTION ++++++++++++++++++++++++++++++++++++++++++++
     //cerr << "[nAg=" << nAg 
     //  << ",Antigen["<<nAg<<"]=" << Antigen[nAg] << ",ssp["<<Antigen[nAg]<<"].k=(";
     //for (int jj=0; jj<dim; jj++) { cerr << ssp[Antigen[nAg]].k[jj] << ","; }
     //cerr << ")";
      thisdistance = intN_mutation(ssp[n].k,ssp[Antigen[nAg]].k);
      //cerr << "this="<<thisdistance<<",";
      if (thisdistance < distance) {
         distance = thisdistance;
      }
      //cerr << "dist="<<distance<<"]";
   }
   //cerr << "]; ";
   return distance;
}
double SS::mean_hamming_distance(long &n) {
   /* returns the mean distance of shape space point <n> to all antigens/epitopes
    * always uses the Hamming distance (irrespective of <measure>)
    */
   // First determine the nearest antigen:
   int sumdistance = 0;
   for (int nAg = 0; nAg < n_Antigen; nAg++) {
      sumdistance += intN_mutation(ssp[n].k,ssp[Antigen[nAg]].k);
   }
   double meandistance = 0;
   if (n_Antigen > 0) {
      meandistance = double (sumdistance) / double (n_Antigen);
   }
   return meandistance;
}
double SS::N_mutation(long * k, long * l) {
   // double middle=PointsPerDim/2.;
   double distance = 0.;
   for (int i = 0; i < dim; i++) {
      // if (par.Value.CyclicBorder==1)
      // if (xdist>middle) xdist=par.Value.RangePerDim-xdist;
      distance += fabs(k[i] - l[i]);
   }
   return distance * distance;
}
double SS::Abstandquad(long &n, long &n2) {
   /* calculates the distance between shape space points
    * denoted by their indices.
    * Uses the Hamming or euklidean distance depending on <measure>
    */
  if (measure == 1) {
    return N_mutation(ssp[n].k, ssp[n2].k);
  } else if (measure == 2) {
    return euklid(ssp[n].k, ssp[n2].k);
  }
  return -1.0;
}
void SS::fix_neighbors(const long int &wo, long * nextn) {
   // Temporaeres
   long int zahl;
   long int tmp[dim];
   // Go through the dimensions
   int i;
   for (int n = 0; n < dim; n++) {
      // Linke Nachbarn
      zahl = ssp[wo].k[n] - 1;
      for (i = 0; i < dim; i++) {
         tmp[i] = ssp[wo].k[i];
      }
      if (zahl < 0) {
         zahl = -1;
      }
      tmp[n] = zahl;
      if (zahl >= 0) {
         nextn[n] = Index(tmp);
      } else {
         nextn[n] = -1;
      }
      // Rechte Nachbarn
      zahl = ssp[wo].k[n] + 1;
      for (i = 0; i < dim; i++) {
         tmp[i] = ssp[wo].k[i];
      }
      if (zahl >= PointsPerDim) {
         zahl = -1;
      }
      tmp[n] = zahl;
      if (zahl >= 0) {
         nextn[n + dim] = Index(tmp);
      } else {
         nextn[n + dim] = -1;
      }
   }
}
void SS::fix_neighbors(const long int &wo) {
   // Temporaeres
   long int zahl;
   long int tmp[dim];
   // Go through the dimensions
   int i;
   for (int n = 0; n < dim; n++) {
      // Linke Nachbarn
      zahl = ssp[wo].k[n] - 1;
      for (i = 0; i < dim; i++) {
         tmp[i] = ssp[wo].k[i];
      }
      if (zahl < 0) {
         zahl = -1;
      }
      tmp[n] = zahl;
      if (zahl >= 0) {
         ssp[wo].nn[n] = Index(tmp);
      } else {
         ssp[wo].nn[n] = -1;
      }
      // Rechte Nachbarn
      zahl = ssp[wo].k[n] + 1;
      for (i = 0; i < dim; i++) {
         tmp[i] = ssp[wo].k[i];
      }
      if (zahl >= PointsPerDim) {
         zahl = -1;
      }
      tmp[n] = zahl;
      if (zahl >= 0) {
         ssp[wo].nn[n + dim] = Index(tmp);
      } else {
         ssp[wo].nn[n + dim] = -1;
      }
   }
}
double SS::affinity(long int n, long int n2) {
   // Berechne das Betragsquadrat des Abstands der Punkte im Shapespace:
   // double distquad=Abstandquad(n,n2);
   // return exp(-1.*distquad/width2);
   return amplitude * exp(-1. * Abstandquad(n, n2) / width2);
}
double SS::affinity(long int n, long int n2, double &tr) {
   // Passe die Breite durch den threshold an:
   double w = width2 * (1.0 - tr) * (1.0 - tr);
   // Berechne das Betragsquadrat des Abstands der Punkte im Shapespace:
   return amplitude * exp(-1. * Abstandquad(n, n2) / w);
}
double SS::affinity_norm(long int n, long int n2) {
   // Berechne das Betragsquadrat des Abstands der Punkte im Shapespace:
   // double distquad=Abstandquad(n,n2);
   // return exp(-1.*distquad/width2);
   return exp(-1. * Abstandquad(n, n2) / width2);
}
void SS::correct_average_affinity(cells celltyp, long &pos, double &average) {
   // Never call this routine before the first cell of type celltyp was created in SS!
   if (sum_cell[celltyp] == 1) {
      average = 0;
   }
   // cout<<"average vorher="<<average<<"; add aff="<<best_affinity(pos)<<"; ";
   average = ((sum_cell[celltyp] - 1) * average + best_affinity_norm(pos)) / (sum_cell[celltyp]);
   // cout<<"sum="<<sum_cell[celltyp]<<"; threshold="<<average<<"\n";
}
// ============================================================
// Public routines:
// ============================================================

long int SS::get_random_SS_position() {
   long x[dim];
   for (int j = 0; j < dim; j++) {
      x[j] = irandom(int (PointsPerDim));
   }
   return Index(x);
}
long int SS::get_random_SS_position(double mindist, double maxdist) {
   /* Routine for the definition of a pool of seeder cells.
    * Returns a random shape space position in a defined range.
    * The range is defined by the parameters <mindist> and <maxdist>,
    * which are distances to the nearest antigen/best clone.
    * The distances are interpreted with the same metric as for
    * affinity calculations.
    * When maxdist<mindist or when at least one value is negative,
    * a random position without restriction is returned.
    * This allows to always call this routine from the constructor.
    */
   if (maxdist < mindist) {
      return get_random_SS_position();
   }
   long rseed = 0;
   long bestclone = 0;
   double dist = 0.;
   do {
      // get a random position in shape space
      rseed = get_random_SS_position();
      // determine the nearest best clone
      bestclone = Antigen[get_nearest_Antigen(rseed)];
      // get the distance between clone and best clone
      dist = sqrt(Abstandquad(rseed, bestclone));
      //    cout<<mindist<<" "<<dist<<" "<<maxdist<<"--> "
      // <<((mindist>0 && dist<mindist) || (maxdist>0 && dist>maxdist))<<"\n";
   } while ((mindist > 0 && dist < mindist) || (maxdist > 0 && dist > maxdist));
   return rseed;
}
bool SS::add_new_Antigen() {
   // ### to be programmed
   return 0;
}
bool SS::add_new_Antigen(long ag_index) {
   // ### to be programmed
   return 0;
}
long int SS::get_Seeder() {
   /* Returns a random position for a seeder cell within the pool of
    * predefined seeder cells <Seeder[n_Seeder]>.
    * <Seeder[n_Seeder]> depends on parameter values and is defined
    * in the SS-constructor.
    */
   int p = irandom(n_Seeder);
   return Seeder[p];
}
long int SS::get_Antigen() {
   // Zufaellige Position im Pool ermitteln
   int p = irandom(n_Antigen);
   return Antigen[p];
}
long int SS::get_Antigen(int i) {
   if (i < n_Antigen) {
      return Antigen[i];
   } else {
      return -1;
   }
}
int SS::get_nearest_Antigen(long n) {
   int best_antigen = 0;
   double bestdist = Abstandquad(n, Antigen[best_antigen]);
   double tmp = 0;
   for (int i = 1; i < n_Antigen; i++) {
      tmp = Abstandquad(n, Antigen[i]);
      if (tmp < bestdist) {
         bestdist = tmp;
         best_antigen = i;
      }
   }
   return best_antigen;
}
long int SS::get_Seeder(int i) {
   if (i < n_Seeder) {
      return Seeder[i];
   } else {
      return get_Seeder();
   }
}
double SS::best_affinity(long pos) {
   double tmp;
   double aff = 0.;
   for (int j = 0; j < n_Antigen; j++) {
      tmp = affinity(pos, Antigen[j]);
      if (tmp > aff) {
         aff = tmp;
      }
   }
   // if (aff>1.) { cout<<"Affinity="<<aff<<" at SSIndex="<<pos<<"\n"; exit(1); }
   return aff;
}
double SS::best_affinity_norm(long pos) {
   double tmp;
   double aff = 0.;
   for (int j = 0; j < n_Antigen; j++) {
      tmp = affinity_norm(pos, Antigen[j]);
      if (tmp > aff) {
         aff = tmp;
      }
   }
   // if (aff>1.) { cout<<"Affinity="<<aff<<" at SSIndex="<<pos<<"\n"; exit(1); }
   return aff;
}
double SS::best_distance(long pos) {
    /* returns the distance of shape space point <pos> to the nearest antigen
     */
    // First determine the nearest antigen:
    int distance = 1.0e+09; // stupid large number
    int thisdistance = 0;

    for (int nAg = 0; nAg < n_Antigen; nAg++) {
        if (measure==1) {thisdistance = intN_mutation(ssp[pos].k,ssp[Antigen[nAg]].k);}
        else if (measure==2) {thisdistance = sqrt(euklid(ssp[pos].k,ssp[Antigen[nAg]].k));}
        if (thisdistance < distance) {
            distance = thisdistance;
        }
    }
    return distance;
}
double SS::wrong_best_affinity_norm(long pos) {
   double tmp;
   double aff = 0.;
   // ++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++++++
   double wrongwidth2=2.8*2.8;
   // ++++++++++++++++ end OPTION +++++++++++++++++++++++++++++++++++
   for (int j = 0; j < n_Antigen; j++) {
      tmp = exp(-1. * Abstandquad(pos, Antigen[j]) / wrongwidth2);
      if (tmp > aff) {
         aff = tmp;
      }
   }
   // if (aff>1.) { cout<<"Affinity="<<aff<<" at SSIndex="<<pos<<"\n"; exit(1); }
   return aff;
}
double SS::mean_affinity_norm(long pos) {
   double tmp = 0.;
   for (int j = 0; j < n_Antigen; j++) {
      tmp += affinity_norm(pos, Antigen[j]);
   }
   tmp /= double (n_Antigen);
   return tmp;
}
/*double SS::get_affinity2ag(long pos, int i) {
 *  return affinity_norm(pos, Antigen[i]);
 *  }*/
// ============================================================
// Dateien schreiben
// ============================================================
void SS::mean_affinities_ag(double * mean_affinities, const int larr, int ag_index) {
  double aff, nBCs, nOUT;
  double totalBC[larr];
  // intialise affinities
  for (int i = 0; i < larr; i++) {
    mean_affinities[i] = 0.;
    totalBC[i] = 0.;
  }
  // go through all points of the AffinitySpace
  for (long i = 0; i < PointsTotal; i++) {
    // load affinities[] with corresponding values
    /* [0]: mean affinity of CB+CC to Antigen[ag_index]
     *  [1]: as [0] but only including cells with aff>0.1
     *  [2]: as [0] but only including cells with aff>0.2
     *  [3-5]: as [0-2] for accummulated output
     *  [6]: fraction of CB+CC cells with aff>0.3
     *  [7]: fraction of OUT cells with aff>0.3
     *  [8]: fraction of CB+CC cells with aff>0.8
     *  [9]: fraction of OUT cells with aff>0.8 (not used)
     *  [10]: number of CB+CC with aff>0.8
     *  [11]: number of OUT with aff>0.8 (not used)
     */
    // determine the affinity of this point to the antigen denoted by ag_index
    aff = affinity_norm(i, Antigen[ag_index]);
    // get the total number of BCs at this point
    nBCs = ssp[i].n_cell[sCB] + ssp[i].n_cell[sCC];
    // get the total number of accumulated output at this point
    nOUT = ssp[i].n_cell[sout];
    // add this to the mean with all BC and OUT
    mean_affinities[0] += aff * nBCs;
    mean_affinities[3] += aff * nOUT;
    // add this to the total cell number
    totalBC[0] += nBCs;
    totalBC[3] += nOUT;
    // do the same for the subgroup of points with a threshold affinity:
    if (aff > 0.1) {
      mean_affinities[1] += aff * nBCs;
      mean_affinities[4] += aff * nOUT;
      totalBC[1] += nBCs;
      totalBC[4] += nOUT;
      }
    if (aff > 0.2) {
      mean_affinities[2] += aff * nBCs;
      mean_affinities[5] += aff * nOUT;
      totalBC[2] += nBCs;
      totalBC[5] += nOUT;
    }
    if (aff > 0.3) {
      totalBC[6] += nBCs;
      totalBC[7] += nOUT;
    }
    if (aff > 0.8) {
      totalBC[8] += nBCs;
      totalBC[9] += nOUT;
    }
  }
  // normalise to get the mean affinities
  for (int i = 0; i < 6; i++) {
    if (totalBC[i] > 0) {
      mean_affinities[i] /= totalBC[i];
    } else {
      mean_affinities[i] = 0.;
    }
  }
  // get the fraction of aff>0.3 >0.8 for BC
  mean_affinities[10] = totalBC[8];
  if (totalBC[0] > 0) {
    mean_affinities[6] = totalBC[6] / totalBC[0];
    mean_affinities[8] = totalBC[8] / totalBC[0];
  } else {
    mean_affinities[6] = 0.;
    mean_affinities[8] = 0.;
  }
  // ... and for accumulated output
  mean_affinities[11] = totalBC[9];
  if (totalBC[3] > 0) {
    mean_affinities[7] = totalBC[7] / totalBC[3];
    mean_affinities[9] = totalBC[9] / totalBC[3];
  } else {
    mean_affinities[7] = 0.;
    mean_affinities[9] = 0.;
  }
}
void SS::write_gcbc_hamming(double time) {
  //  cerr << "Start SS:write_gcbc_hamming ...\n";
   // +++++++++++++++ OPTION ++++++++++++++++
   int max_hamming = 12;
   // +++++++++++ end OPTION ++++++++++++++++
   // define a set of counters for each Hamming distance
   int nbc[max_hamming + 1];
   // .. and for the mean Hamming distance for all antigens/epitopes
   int nbc_mean_ag[max_hamming + 1];
   int nbchere = 0;
   for (int n = 0; n <= max_hamming; n++) {
      nbc[n] = 0;
      nbc_mean_ag[n] = 0;
   }
   int hamming = 0;
   double hamming_mean = 0.;
   //cerr << "Before loop; ";
   // remove this line in the general case:
   //cout << "WARNING! Modified wrong code in SS::best_hamming_distance(long&)\n";
   // Go through all points of the shape space and add to the counters:
   for (long i = 0; i < PointsTotal; i++) {
     //cerr << "i=" << i << ": ";
     //cerr << "(" << ssp[i].n_cell[sCB] << ",";
     //cerr << ssp[i].n_cell[sCC] << ")";
      // get the total number of GC-BCs at this point
      nbchere = ssp[i].n_cell[sCB] + ssp[i].n_cell[sCC];
      //cerr << "A";
      // ##### Note that this may be speeded up by only continuing when nbchere > 0!!!
      // determine the Hamming distance of this point to the nearest antigen
      hamming = best_hamming_distance(i);
      //cerr << "B(" << hamming << ")";
      hamming_mean = mean_hamming_distance(i);
      //cerr << "C(" << hamming_mean << ")";
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      if (hamming >= max_hamming) {
         nbc[max_hamming] += nbchere;
      } else {
         nbc[hamming] += nbchere;
      }
      //cerr << "D(";
      //for (int jj=0; jj<=max_hamming; jj++) { cerr << nbc[jj] << ","; }
      //cerr << ")";
      // now do rounding of the mean hamming distance of SS-point i:
      hamming = int (hamming_mean + 0.5);
      //cerr << "E(" << hamming << ")";
      // ... and attribute to the
      if (hamming >= max_hamming) {
         nbc_mean_ag[max_hamming] += nbchere;
      } else {
         nbc_mean_ag[hamming] += nbchere;
      }
      //cerr << "F(";
      //for (int jj=0; jj<=max_hamming; jj++) { cerr << nbc_mean_ag[jj] << ","; }
      //cerr << ")";
   }
   //cerr << "loop finished.\n";
   // open the output file (append)
   ofstream gcbc_hamming;
   gcbc_hamming.open("gcbc_hamming.out", ofstream::app);
   for (int n = 0; n <= max_hamming; n++) {
     //cerr << "n=" << n << ": (" << nbc[n] << "," << nbc_mean_ag[n] << "); ";
      gcbc_hamming << time / 24. << "   " << time << "   "
                   << n << "   "
                   << nbc[n] << "   "
                   << nbc_mean_ag[n]
                   << "\n";
   }
   gcbc_hamming.close();
   //cerr << "Finished SS::write_gcbc_hamming.\n";
}
void SS::write_gcbc_affinity(double time) {
  //cerr << "Start SS::write_gcbc_affinity ...\n";
   // +++++++++++++++ OPTION ++++++++++++++++
   double logfactor = 3.0;
   int n_bins = 10;
   double affarray[n_bins];
   affarray[0] = 0;
   affarray[1] = 0.0002 / logfactor;
   for (int b = 2; b < n_bins; b++) {
      affarray[b] = logfactor * affarray[b - 1];
   }
   /*
    *  for (int b = 0; b < n_bins; b++) {
    *  cout << affarray[b] << ",";
    *  }
    *  cout << "\n";
    */
   // +++++++++++ end OPTION ++++++++++++++++
   // double check that the whole range of binding probabilities is covered
   if (affarray[n_bins - 1] * logfactor < 1.0) {
      cerr << "In SS::write_gcbc_affinity(double): largest affinity is smaller than 1.\n"
           << "Abort simulation.\n";
      exit(1);
   }
   // define a set of counters for each Hamming distance
   int nbc[n_bins];
   // .. and for the mean Hamming distance for all antigens/epitopes
   int nbc_mean_ag[n_bins];
   int nbchere = 0;
   for (int n = 0; n < n_bins; n++) {
      nbc[n] = 0;
      nbc_mean_ag[n] = 0;
   }
   double affbest = 0;
   double affmean = 0.;
   // Go through all points of the shape space and add to the counters:
   for (long i = 0; i < PointsTotal; i++) {
      // get the total number of GC-BCs at this point
      nbchere = ssp[i].n_cell[sCB] + ssp[i].n_cell[sCC];
      // determine the Hamming distance of this point to the nearest antigen
      // ++++++++++++++++++++++++++++ OPTION +++++++++++++++++++++++++++++++++++++++++++++
      /* Here, either the normal affinity as it is used for the binding probability in the
       * simulation itself is used to plot the affinity distribution (best_affinity_norm), 
       * or the same but with a locally fixed width of the binding probability is used
       * (wrong_best_affinity_norm). The latter can be useful when the width of the
       * binding probability is changed in order to modulate the binding probability,
       * but still realistic affinities should be plotted.
       */
      affbest = best_affinity_norm(i);
      //affbest = wrong_best_affinity_norm(i);
      // ++++++++++++++++++++++++ end OPTION +++++++++++++++++++++++++++++++++++++++++++++
      affmean = mean_affinity_norm(i);
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      int binbest = 0, binmean = 0;
      while (binbest<n_bins&&affbest> affarray[binbest]) {
         ++binbest;
      } // now binbest is one too far, as affbest>0, binbest>0 as well
      --binbest;
      // this is the array position for nbc
      nbc[binbest] += nbchere;
      // repeat the same for affmean
      while (binmean<n_bins&&affmean> affarray[binmean]) {
         ++binmean;
      }
      --binmean;
      nbc_mean_ag[binmean] += nbchere;
   }
   // open the output file (append)
   ofstream gcbc_affinity;
   gcbc_affinity.open("gcbc_affinity.out", ofstream::app);
   for (int n = 0; n < n_bins; n++) {
      gcbc_affinity << time / 24. << "   "
                    << time << "   "
                    << affarray[n] << "   "
                    << nbc[n] << "   "
                    << nbc_mean_ag[n]
                    << "\n";
   }
   gcbc_affinity.close();
   //cerr << "Finished SS::write_gcbc_affinity.\n";

}
void SS::to_multiag_files(double time) {
  ofstream
    saffin_ags, saffin_t10, saffin_t20,
    saffin_out_ags, saffin_out_t10, saffin_out_t20,
    shaffin_ag, shaffin_ag_t80, shaffin_out_ag,
    cross_gcr, cross_out;
  saffin_ags.open("saffin_ags.out", ofstream::app);
  saffin_ags << time << "  ";
  saffin_t10.open("saffin_ags_t10.out", ofstream::app);
  saffin_t10 << time << "  ";
  saffin_t20.open("saffin_ags_t20.out", ofstream::app);
  saffin_t20 << time << "  ";
  saffin_out_ags.open("saffin_out_ags.out", ofstream::app);
  saffin_out_ags << time << "  ";
  saffin_out_t10.open("saffin_out_ags_t10.out", ofstream::app);
  saffin_out_t10 << time << "  ";
  saffin_out_t20.open("saffin_out_ags_t20.out", ofstream::app);
  saffin_out_t20 << time << "  ";
  shaffin_ag.open("shaffin_ags.out", ofstream::app);
  shaffin_ag << time << "  ";
  shaffin_ag_t80.open("shaffin_ags_t80.out", ofstream::app);
  shaffin_ag_t80 << time << "  ";
  shaffin_out_ag.open("shaffin_out_ags.out", ofstream::app);
  shaffin_out_ag << time << "  ";
  cross_gcr.open("cross_reactivity_gcr.out", ofstream::app);
  cross_gcr << time << "  ";
  cross_out.open("cross_reactivity_out.out", ofstream::app);
  cross_out << time << "  ";
  double fracsum = 0., fracOUTsum = 0., fracsumt80 = 0.;
  // variables for cross reactivity
  double cross_reactivity[6];
  for (int i = 0; i < 6; i++) {
    cross_reactivity[i] = 0.;
  }
  const int larr = 12;
  for (int a = 0; a < n_Antigen; a++) {
    double mean_affinities[larr];
    mean_affinities_ag(mean_affinities, larr, a);
    saffin_ags << mean_affinities[0] << "  ";
    saffin_t10 << mean_affinities[1] << "  ";
    saffin_t20 << mean_affinities[2] << "  ";
    saffin_out_ags << mean_affinities[3] << "  ";
    saffin_out_t10 << mean_affinities[4] << "  ";
    saffin_out_t20 << mean_affinities[5] << "  ";
    shaffin_ag << mean_affinities[6] << "  ";
    fracsum += mean_affinities[6];
    shaffin_ag_t80 << mean_affinities[8] << "  "     // fraction with aff>0.8 for each ag
		   << mean_affinities[10] << "    "; // number of BC with aff>0.8 for each ag
    fracsumt80 += mean_affinities[8];
    shaffin_out_ag << mean_affinities[7] << "  ";
    fracOUTsum += mean_affinities[7];
    // make the sum of all antigens to get cross-reactivity (with thresholds) for BC and out
    for (int i = 0; i < 6; i++) {
      cross_reactivity[i] += mean_affinities[i];
    }
  }
  // add the sum of the fractions, which is a measure of cross-reactivity, to the fraction-file
  shaffin_ag << fracsum << "\n";
  shaffin_ag.close();
  shaffin_ag_t80 << fracsumt80 << "\n";
  shaffin_ag_t80.close();
  shaffin_out_ag << fracOUTsum << "\n";
  shaffin_out_ag.close();
  saffin_ags << "\n";
  saffin_t10 << "\n";
  saffin_t20 << "\n";
  saffin_out_ags << "\n";
  saffin_out_t10 << "\n";
  saffin_out_t20 << "\n";
  saffin_ags.close();
  saffin_t10.close();
  saffin_t20.close();
  saffin_out_ags.close();
  saffin_out_t10.close();
  saffin_out_t20.close();
  // normalise cross-reactivities with the number of antigens
  for (int i = 0; i < 3; i++) {
    if (n_Antigen > 0) {
      cross_reactivity[i] /= n_Antigen;
      cross_reactivity[i + 3] /= n_Antigen;
    } else {
      cross_reactivity[i] = 0.;
      cross_reactivity[i + 3] = 0.;
    }
    cross_gcr << cross_reactivity[i] << "  ";
    cross_out << cross_reactivity[i + 3] << "  ";
  }
  cross_gcr << "\n";
  cross_out << "\n";
  cross_gcr.close();
  cross_out.close();
}
double SS::mean_affinity(double * affinities) {
   long i;
   double aff, all;
   for (i = 0; i < 15; i++) {
      affinities[i] = 0.;
   }
   for (i = 0; i < PointsTotal; i++) {
      // Determine affinity between point i and antigen
      // Find nearest antigen-type
      aff = best_affinity_norm(i);
      if (aff > 0.3) {
         affinities[3] += double (ssp[i].n_cell[sCB]);
         affinities[4] += double (ssp[i].n_cell[sCC]);
         affinities[5] += double (ssp[i].n_cell[sout]);
      }
      if (aff < 0.4) {
         affinities[6] += double (ssp[i].n_cell[sCB]);
         affinities[7] += double (ssp[i].n_cell[sCC]);
         affinities[8] += double (ssp[i].n_cell[sout]);
      }
      if ((aff >= 0.4) && (aff < 0.8)) {
         affinities[9] += double (ssp[i].n_cell[sCB]);
         affinities[10] += double (ssp[i].n_cell[sCC]);
         affinities[11] += double (ssp[i].n_cell[sout]);
      }
      if (aff >= 0.8) {
         affinities[12] += double (ssp[i].n_cell[sCB]);
         affinities[13] += double (ssp[i].n_cell[sCC]);
         affinities[14] += double (ssp[i].n_cell[sout]);
      }
      // add corresponding affinity with weight #cb and #cc to ccaff and cbaff
      affinities[0] += aff * double (ssp[i].n_cell[sCB]);
      affinities[1] += aff * double (ssp[i].n_cell[sCC]);
      affinities[2] += aff * double (ssp[i].n_cell[sout]);
   }
   // Divide through total cell numbers:
   all = affinities[0] + affinities[1];
   if (sum_cell[sCB] > 0) {
      affinities[0] /= sum_cell[sCB];
      affinities[3] /= sum_cell[sCB];
      affinities[6] /= sum_cell[sCB];
      affinities[9] /= sum_cell[sCB];
      affinities[12] /= sum_cell[sCB];
   } else {
      affinities[0] = 0.;
   }
   if (sum_cell[sCC] > 0) {
      affinities[1] /= sum_cell[sCC];
      affinities[4] /= sum_cell[sCC];
      affinities[7] /= sum_cell[sCC];
      affinities[10] /= sum_cell[sCC];
      affinities[13] /= sum_cell[sCC];
   } else {
      affinities[1] = 0.;
   }
   if (sum_cell[sout] > 0) {
      // cout<<outs<<" "<<ha_out<<" "<<sum_cell[sout]<<"\n";
      affinities[2] /= sum_cell[sout];
      affinities[5] /= sum_cell[sout];
      affinities[8] /= sum_cell[sout];
      affinities[11] /= sum_cell[sout];
      affinities[14] /= sum_cell[sout];
      // cout<<outs<<" "<<ha_out<<"\n";
   } else {
      affinities[2] = 0.;
   }
   if (sum_cell[sCB] + sum_cell[sCC] > 0) {
      all /= (sum_cell[sCB] + sum_cell[sCC]);
   } else {
      all = 0.;
   }
   return all;
}
//MSchips
void SS::meanSD_out_affinity(double * aff_mean, double * aff_sd) {
   long i;
   double aff;
//   double aff_mean=0, aff_sd=0;
   //0 : 1 : 2 --> sout (self : total : nonSelf)
   //3 --> total bcs

   for (i = 0; i < PointsTotal; i++) {
      // Determine affinity between point i and antigen
      // Find nearest antigen-type
      aff = best_affinity_norm(i);
      aff_mean[0] += (double(ssp[i].n_cell[soutSelf])*aff);
      aff_mean[1] += (double(ssp[i].n_cell[sout])*aff);
      aff_mean[2] += (double(ssp[i].n_cell[soutNonSelf])*aff);
      double liveBCs=double(ssp[i].n_cell[sCB]+ssp[i].n_cell[sCC]-ssp[i].n_cell[sCCapoptosis]);
      aff_mean[3] += liveBCs*aff;
   }
   // Divide through total cell numbers:
   if (sum_cell[soutSelf] > 0) {
      aff_mean[0] /= sum_cell[soutSelf];
   } else {
      aff_mean[0] = 0.;
   }
   if (sum_cell[sout] > 0) {
      aff_mean[1] /= sum_cell[sout];
   } else {
      aff_mean[1] = 0.;
   }
   if (sum_cell[soutNonSelf] > 0) {
      aff_mean[2] /= sum_cell[soutNonSelf];
   } else {
      aff_mean[2] = 0.;
   }
   if (sum_cell[sCB]>0 || sum_cell[sCC]>0) {
       aff_mean[3] /= (sum_cell[sCB]+sum_cell[sCC]-sum_cell[sCCapoptosis]);
   } else {
       aff_mean[3] = 0;
   }
   for (i = 0; i < PointsTotal; i++) {
      // Determine affinity between point i and antigen
      // Find nearest antigen-type
      aff = best_affinity_norm(i);
      aff_sd[0] += ((aff-aff_mean[0])*(aff-aff_mean[0])*double(ssp[i].n_cell[soutSelf]));
      aff_sd[1] += ((aff-aff_mean[1])*(aff-aff_mean[1])*double(ssp[i].n_cell[sout]));
      aff_sd[2] += ((aff-aff_mean[2])*(aff-aff_mean[2])*double(ssp[i].n_cell[soutNonSelf]));
      double liveBCs=double(ssp[i].n_cell[sCB]+ssp[i].n_cell[sCC]-ssp[i].n_cell[sCCapoptosis]);
      aff_sd[3] += ((aff-aff_mean[3])*(aff-aff_mean[3])*liveBCs);
   }
   if (sum_cell[soutSelf] > 0) {
      aff_sd[0] /= sum_cell[soutSelf];
      aff_sd[0] = sqrt(aff_sd[0]);
   } else {
      aff_sd[0] = 0.;
   }
   if (sum_cell[sout] > 0) {
       aff_sd[1] /= sum_cell[sout];
       aff_sd[1] = sqrt(aff_sd[1]);
   } else {
      aff_sd[1] = 0.;
   }
   if (sum_cell[soutNonSelf] > 0) {
       aff_sd[2] /= sum_cell[soutNonSelf];
       aff_sd[2] = sqrt(aff_sd[2]);
   } else {
      aff_sd[2] = 0.;
   }
   if (sum_cell[sCB]>0 || sum_cell[sCC]>0) {
       aff_sd[3] /= (sum_cell[sCB]+sum_cell[sCC]-sum_cell[sCCapoptosis]);
       aff_sd[3] = sqrt(aff_sd[3]);
   } else {
      aff_sd[3] = 0.;
   }
}
void SS::get_diversity(double time) {
   logdiversity << time << "   ";
   int diversity[number];
   for (int j = 0; j < number; j++) {
      diversity[j] = 0;
   }
   for (long i = 0; i < PointsTotal; i++) {
      if (ssp[i].n_cell[sCB] > 0.5) {
         ++diversity[sCB];
      }
      if (ssp[i].n_cell[sCC] > 0.5) {
         ++diversity[sCC];
      }
      if (ssp[i].n_cell[sout] > 0.5) {
         ++diversity[sout];
      }
      if (ssp[i].n_cell[sCB] + ssp[i].n_cell[sCC] + ssp[i].n_cell[sout] > 0.5) {
         ++diversity[total];
      }
      if (ssp[i].n_cell[soutext] + ssp[i].n_cell[soutextproduce] > 0.5) {
         ++diversity[soutext];
      }
   }
   logdiversity << diversity[total] << "   " << diversity[sCB] << "   " << diversity[sCC]
                << "   " << diversity[sout]
                << "   " << diversity[soutext] << "\n";
}
void SS::to_ssfiles(double time) {
   // sum_check();
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) {
         if (i != soutdec) {
            logdata[i] << time << "     " << sum_cell[i];
         }
         if (i == sout) {
            logdata[sout] << "    " << sum_cell[sout] - oldsum_cell[sout];
         } else if (i == soutdec) {
            logdata[sout] << "   " << sum_cell[soutdec] << "   ";
            if (sum_cell[sout] > 0) {
               logdata[sout] << sum_cell[soutdec] / sum_cell[sout] << "\n";
            } else {
               logdata[sout] << "0\n";
            }
         } else if (i == sallapoptosis) {
            logdata[sallapoptosis] << "    " << sum_cell[sallapoptosis]
               - oldsum_cell[sallapoptosis] << "\n";
         } else {
            logdata[i] << "\n";
         }
      }
      oldsum_cell[i] = sum_cell[i];
   }

   double affinities[15];
   double a;
   a = mean_affinity(affinities);
   logmeanaff << time << "   " << a << "   " << affinities[0] << "   " << affinities[1] << "   "
              << affinities[2]
              << "\n";
   loghighaff << time << "   " << affinities[3] << "   " << affinities[4] << "   "
              << affinities[5] << "\n";
   log048aff << time;
   for (short int i = 6; i < 15; i++) {
      log048aff << "   " << affinities[i];
   }
   log048aff << "\n";
   // for analysis
   CB_haffinity = affinities[3];
   CC_haffinity = affinities[4];
   OUT_haffinity = affinities[5];

   get_diversity(time);

   if ((time > 144. - 1.e-08) && (time < 144. + 1.e-08)) {
      OUT_steepness = double (sum_cell[sout]);
   }
   if ((time > 288. - 1.e-08) && (time < 288. + 1.e-08)) {
      OUT_steepness = double (sum_cell[sout]) / OUT_steepness;
   }
}
// ============================================================
// Verwaltung von Zellen im shape space
// ============================================================

void SS::add_cell(cells typ, long int pos) {
   ++ssp[pos].n_cell[typ];
   ++sum_cell[typ];
   if (typ == soutext) { set_external_cell(pos); }
}
void SS::rem_cell(cells typ, long int pos) {
   --ssp[pos].n_cell[typ];
   --sum_cell[typ];
   // ##### remove this when stable for a while (introduced on 2016nov14)
   if (ssp[pos].n_cell[typ] < 0) {
     cerr<<"sspos="<<pos<<", typ="<<typ<<", ss_n_cell_value="<<ssp[pos].n_cell[typ]<<" !!!!\n";
     exit(1);
   }
}
void SS::set_external_cell(long int pos) {
   // set flags to save that an external cell is at this place:
   if (ssp[pos].external_cell == false) {
      ssp[pos].external_cell = true;
      external_cells.push_back(pos);
   }
}
void SS::PM_differentiate(cells a, cells b, double dt) {
   int n_external = external_cells.size();
   for (int n = 0; n < n_external; n++) {
      long i = external_cells[n];
      double transit = ssp[i].n_cell[a] * pm_differentiation_rate * dt;
      /*
       * if (ssp[n].n_cell[a]>0.)
       * cout<<"a="<<a<<"; b="<<b<<"; dt="<<dt<<"; r="<<pm_differentiation_rate
       * <<"; n_cell[a]="<<ssp[n].n_cell[a]<<"; transit="<<transit
       * <<"\n";
       */
      if (transit > 0) {
         ssp[i].n_cell[a] -= transit;
         sum_cell[a] -= transit;
         ssp[i].n_cell[b] += transit;
         sum_cell[b] += transit;
         /*    double transit=ssp[n].d_cell[a]*pm_differentiation_rate*dt;
          * ssp[n].d_cell[a]-=transit;
          * ssp[n].d_cell[b]+=transit;
          * ssp[n].n_cell[a]=long(ssp[n].d_cell[a]+0.5);
          * ssp[n].n_cell[b]=long(ssp[n].d_cell[b]+0.5);
          */
         // set flags to save that an antibody producer is at this place:
         if (ssp[i].ab_producer == false) {
            ssp[i].ab_producer = true;
            ab_producers.push_back(i);
         }
      }
   }
}
/* If antibody producers would also die, this would impact on the class AntibodyDyn */

double SS::get_AbProducingCells(int n) {
  //  cerr << "Ab-producing cells = " << ssp[ab_producers[n]].n_cell[soutextproduce] << endl;
  return ssp[ab_producers[n]].n_cell[soutextproduce];
}
long SS::get_AbProducingIndex(int n) {
   return ab_producers[n];
}
int SS::get_n_ab_producers() {
   return ab_producers.size();
}
void SS::put_Ab(long index,double d_ab) {
   ssp[index].antibody += d_ab;
   if (ssp[index].antibody < 0) {
      ssp[index].antibody = 0;
      cerr << "WARNING! antibody got negative at index " << index << " in shape space.\n";
   }
}
double SS::get_AbAmount(long index) {
   return ssp[index].antibody;
}
short int SS::sum_check() {
   long int i;
   double dtmp;
   short int err = 0;
   cout << "Shape space sum check ... ";
   for (int n = 0; n < number; n++) {
      // Durchlaufe Typen
      dtmp = 0.;
      // Durchlaufe shape space Punkte
      for (i = 0; i < PointsTotal; i++) {
         if (ssp[i].n_cell[n] < 0) {
            cout << "Zellzahl im SS[" << i << "].[" << n << "] = " << ssp[i].n_cell[n]
                 << " <0\n";
            exit(1);
         }
         // if (n!=soutext && n!= soutextproduce)
         dtmp += ssp[i].n_cell[n];
         // else if (n==soutext) { dtmp+=ssp[i].d_cell[n]; dtmp+=ssp[i].d_cell[n+1]; }
         // else if (n==soutext) { dtmp+=ssp[i].n_cell[n]; dtmp+=ssp[i].n_cell[n+1]; }
      }
      if ((dtmp - sum_cell[n] > 1.e-06) || (dtmp - sum_cell[n] < -1.e-06)) {
         cout << "Wrong ss sum for cell type " << n << " !\n";
         cout << "sum=" << sum_cell[n] << " but real sum=" << dtmp << "\n";
         // exit(1);
         err = 1;
      }
   }
   cout << "done\n";
   return err;
}
