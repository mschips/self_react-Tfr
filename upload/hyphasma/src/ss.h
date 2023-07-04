#ifndef i_SS
#define i_SS
#include "setparam.h"
#include "affinityspace.h"

/** @brief Class to store one point in the shape space, with its index and neighbors, */
// Ein Punkt im shapespace besteht aus seinem Index und den Indizes aller Nachbarn.
class SSpoint {
  public:
   /** @brief Konstructor  */
   SSpoint();           // SSpoint(Parameter& par);
   /** @brief Resizes fields according to par.Value.DimShapeSpace */
   void init(Parameter &par);
   /** @brief Destructor    */
   ~SSpoint();

   /** @brief Index of this point in the shape space    */
   long int index;

   /** @brief vector of coordinates
    * it will be of size 'dim' from the SS class it belongs   */
   long * k;
   /** @brief vector with the index of the neighboring points of the shape space*/
   long * nn;

   /** @brief Stores the number of cells that are in this point of the shape space
    * @param type of cell (see 'enum cells', sCB = 0 ...) */
   double n_cell[number];

   /** @brief switch for the existence of external cells and antibody producers at this position */
   bool external_cell,ab_producer;

   /** @brief to store the amount of antigen at this position */
   double antigen;
   /** @brief Amount of antibody of this type */
   double antibody;
};

/** @brief Operators to say if two points are equal or unequal in the shape space   */
char operator ==(const SSpoint &a, const SSpoint &b);
char operator !=(const SSpoint &a, const SSpoint &b);

/** @brief A shape space class to link a position in shape space in to affinity and
 *  to follow over time (output) where the cells are located in the shape space
 *
 *  Layout :
 *  - SS it is a multidimensional grid of points of type SSpoint,
 *  - the size of the grid can change dynamically
 *  - individual points can be accessed using SSpoint[index]
 *    or by using their coordinates (how ?)
 *  - the antigens are the positions in the shape space where the affinity is best.
 *    in most simulations there is only one antigen (the 'optimal' position in the shape space).
 *  - the affinity of a cell will be given by its distance to the closest antigen using
 *    a gaussian function. Two distances are possible : hamming (SS::N_mutation(pt1, pt2)) or
 * euclidian (SS::euklid(pt1,pt2)).
 *  - there are functions to pre-define a pool of 'seeder cells' (= initial positions of cells in
 * the shape space)
 *    and pools of antigens (places with best affinity).
 *    they are caracterized by their index.
 *
 * Implicit (fields) output ofstreams filled by this class :
 *    logmeanaff
 *    loghighaff
 *    logdiversity
 *    log048aff
 *    logdata[each cell type]
 * output files generated :
 *      logdata[sCB].open("ssumcb.out");
 *      logdata[sCC].open("ssumcc.out");
 *      logdata[sCCunselected].open("ssumccun.out");
 *      logdata[sCCcontact].open("ssumccco.out");
 *      logdata[sCCselected].open("ssumccse.out");
 *      logdata[sCCapoptosis].open("ssumapo1.out");
 *      logdata[sallapoptosis].open("ssumapo2.out");
 *      ... FDC and Tcell are ignored !
 *      logdata[sout].open("ssumout.out");
 *      logdata[soutext].open("ssumoute.out");
 *      logdata[soutextproduce].open("ssumoutep.out");
 *      logdata[total].open("ssum.out");
 *
 *      logmeanaff.open("saffin.out"):  Mean affinity of CB+CC : CB : CC : out\n";
 *      loghighaff.open("shaffin.out"): Fraction of >30% affinity CB:CC:out\n";
 *      log048aff.open("s048aff.out"):  Fraction of <40% affinity CB:CC:out : 40-80% : >=80%\n";
 *      logdiversity.open("diversity.out"): Number of different encoded antibodies present in the
 * GC\n"  */
/* ist ein Mehrdimensionales Gitter von Punkten des Typs SSpoint.
 * Die Groesse des Gitters ist dynamisch.
 * Die einzelnen Punkte sind point[..]
 * Die Indizes der Punkte werden durchnummeriert koennen aber auch
 * in Koordinatenform ausgedrueckt werden.
 * Es stehen Routinen zur Verfuegung, um Pools von Antigenen und Seeder cells
 * zu definieren und von aussen abzurufen. Dabei sind diese durch ihren Index
 * im shapespace charakterisiert.
 * Routinen fuer den Abstand zwischen zwei Punkten und fuer die Affinitaet
 * von Antigen und antibody sind vorhanden.*/

class SS: public AffinitySpace {
  public:
   /** @brief Useless constructor */
   SS();
   /** @brief Useful constructor
    *  @param par   (Parameter) global list of parameters
    *  @param ana   (ofstream)  to output what's happening to the shape space
    *
    *  important parameters, inside par are (kenntext, access, meaning :)
    *
    *  From the shape space paragraph
    *       - Dimension of Shapespace (<11):          par.value.DimShapeSpace
    *           --> stored in dim
    *       - Use metric (1:N_mutation; 2:euclidian): par.Value.metrik
    *           --> stored in measure
    *       - Number of antibody types:               par.Value.SSStates      
    *           --> total number of points in the shape space
    *           --> stored in PointsTotal
    *       - Number per Dimension (Fixed,int-type):  par.Value.SSRangePerDim 
    *           --> size of the grid.
    *           --> stored in PointsPerDim
    *               note that these parameters should fit : PointsTotal = PointsPerDim ^
    * DimShapeSpace
    *
    *       - Number of initial Anitgen Peaks in its Shapespace (int-type): par.Value.APeakNumber
    *            --> stored in n_Antigen;
    *       - Fix Epitop presentation (max 10 values):  par.Values.takeA[1...] (dynArray)
    *            --> the first n_Antigen ones are taken, and random ones 
    *                (random index/position in the shape space) 
    *                are created if nAntigen > takeA.size
    *
    *       - Width of gaussian affinity weight function:     par.Value.GammaGauss
    *           --> stored in width2
    *       - Amplitude of Gauss affinity weight function (0<a<=1): par.Value.amplitudeGauss
    *           --> stored in amplitude
    *
    *  From the centroblast paragraph
    *       - Total Number of initial B-Cells types: par.Value.totalBss --> initial number of
    *seeders
    *           --> saved in n_Seeder
    *       - Fix initial Centroblast distribution (max 10 values): par.Value.takeB[]
    *           --> the first n_Seeder ones are taken, and random ones 
    *               (random index/position in the shape space) 
    *               are created if n_Seeders > takeB.size
    *
    *  NOT USED PARAMETERS :
    *       - Total initial Number of presented Antigen Epitops (int-type):totalA
    *       - Use relative epitop weights in selection (=1): EpitopWeight
    *
    *  Other parameters :
    *       pm_differentiation_rate = log(2.)/par.Value.pm_differentiation_time; --> rate of
    * differentiation into plasmablasts       */

   SS(Parameter &par, ofstream &ana);
   /** @brief Destruktor */
   ~SS();
   
   /** @brief Givesthe number of antigene, */
   int get_n_Antigen();

   bool add_new_Antigen(); // ### to be programmed
   bool add_new_Antigen(long ag_index); // ### to be programmed

   /** @brief Gives the index of one of the antiges,
    *  randomly chosen from the predefined pool*/
   // Liefere einen zufaelligen Index aus dem Antigen Pool
   long int get_Antigen();

   /** @brief Gives the index of antigen nr i */
   // Liefere einen bestimmten Index aus dem Antigen Pool
   long int get_Antigen(int i);

   /** @brief Return the index of the antigen nearest to the point <n>: */
   int get_nearest_Antigen(long n);

   /** @brief Gives the index of one of the seeding positions for cells,
    *  randomly chosen from the pre-defined pool*/
   // Liefere einen zufaelligen Index aus dem seeder-cell Pool
   long int get_Seeder();

   /** @brief Gives the index of seeding position nr i */
   // Liefere einen bestimmten Index aus dem seeder-cell Pool
   long int get_Seeder(int i);

   /** Interface functions with the simulation */

   /** @brief Mutates a sequence and gets the new position in the shape space */
   long getMutation(long from_position);

   /** @brief  Get binding affinity between two points */
   /** variant1 : using gaussian function (amplitude, width) */
   double affinity(long int n, long int m);          // ### only used in cellTC::make_tx_cc_link.
                                                     // should be affinity_norm?
   /** variant2 : using gaussian (amplitude, width), with a correction from the threshold */
   double affinity(long int n, long int n2, double &tr);
   /** variant3 : using gaussian (width only, with amplitude = 1)*/
   double affinity_norm(long int n, long int m);

   /** @brief  From a simulation, ask to change the numbers of cells at different SS points */
   // Variation of numbers at different SS points
   void add_cell(cells typ, long int pos);
   void rem_cell(cells typ, long int pos);

   /** side functions to observe what is happening inside the shape space */

   // average affinity (made for output cells)
   void correct_average_affinity(cells celltyp, long &pos, double &average);
   
   /** @brief Gives the best affinity from position <pos> to the Ags, uses best_affinity_norm*/
   double best_affinity(long pos);  // deprecated version
   double best_affinity_norm(long pos);

   /** @brief Gives the affinity of position <pos> to antigen with list-index <i> */
   // double get_affinity2ag(long pos, int i); --> replaced by affinity_norm directly in the code

   /** @brief square distance between two points in shape space */
   double Abstandquad(long int &n, long int &n2);
   
   double get_sum_cell(cells celltype);
   double get_oldsum_cell(cells celltype);
   /** @brief Checks if the sum of all the cells in sum_cell were correctly computed */
   short int sum_check();

   /** @brief  writes a sum up in the ofstreams */
   // general for single antigen
   void to_ssfiles(double time);
   // and for multiple antigen
   void to_multiag_files(double time);

   void write_gcbc_hamming(double time);
   void write_gcbc_affinity(double time);

   /** @brief To close the ofstreams */
   void close_files();

   /** @brief  Indicators that are updated when to_ss files is called */
   // Charakteristische Daten zum Output
   double OUT_haffinity, OUT_steepness, CB_haffinity, CC_haffinity;

   /** @brief plasma cell differentiation : 
    * converts soutext cells to soutextproduce cells with the
    * appropriate rate */
   void PM_differentiate(cells a, cells b, double dt);

   // returns the number of producing cells at the point of element n of vector ab_producers
   double get_AbProducingCells(int n);
   // returns the index saved in element n of vector ab_producers
   long get_AbProducingIndex(int n);
   // get total number of antibody producers
   int get_n_ab_producers();

   /** @brief handling of antibodies associated with points in AffinitySpace */
   void put_Ab(long index,double d_ab);
   double get_AbAmount(long index);

  private:
   /** @brief Containers for Shape Space points,
    *  size : PointsTotal
    *  access : ssp[i] = point which index i */
   SSpoint * ssp;

   /** @brief Number of Dimensions */
   short int dim;

   /** @brief Size of each dimension */
   long int PointsPerDim;
   /** @brief Total Number of Points */
   long int PointsTotal;

   /** @brief Antigen Pool */
   long * Antigen;
   int n_Antigen;
   /** @brief  Pool of seeder cells */
   long * Seeder;
   int n_Seeder;

   // Mittlere Affinitaet der Zellen zum Antigen
   /** @brief returns average affinity of each cells, and fills the table affinities with :
    *   affinity[0,1,2] = average affinities of cells from [sCB, sCC, sout] to have aff > 0.3
    *   affinity[3,4,5] = % of cells from [sCB, sCC, sout] to have aff > 0.3
    *   affinity[6,7,8] = % of cells from [sCB, sCC, sout] to have aff < 0.4
    *   affinity[9,10,11] = % of cells from [sCB, sCC, sout] to have 0.4 =< aff < 0.8
    *   affinity[12,13,14] = % of cells from [sCB, sCC, sout] to have aff >= 0.8  */
   double mean_affinity(double * affinities);
   //Mschips
   void meanSD_out_affinity(double * aff_mean, double * aff_sd);
   double best_distance(long pos);

   /** @brief Gives back the affinities (mean and sd) of the required type */
   /* for multiple antigens this returns corresponding values differentiated for each antigen
    * [0]: mean affinity of CB+CC to Antigen[ag_index]
    * [1]: as [0] but only including cells with aff>0.1
    * [2]: as [0] but only including cells with aff>0.2
    * [3-5]: as [0-2] for accummulated output
    * [6]: fraction of CB+CC cells with aff>0.3
    * [7]: fraction of OUT cells with aff>0.3 */
   void mean_affinities_ag(double * mean_affinities, const int larr, int ag_index);
   /** @brief Gives back the Hamming distance (not squared) between to shape space vectors */
   int intN_mutation(long * k, long * l);
   /** @brief Gives back the best Hamming distance of SS-position n to any of the antigens */
   int best_hamming_distance(long &n);
   /** @brief Gives back the mean Hamming distance of SS-position n to all antigens */
   double mean_hamming_distance(long &n);
   /** same as best_affinit_norm(long) but returning the mean of affinities to all antigens */
   double mean_affinity_norm(long pos);
   /** @brief Returns the same as best_affinity_norm(long pos) but using a locally fixed width */
   double wrong_best_affinity_norm(long pos);
  
   /** @brief Sum of each kind of cells over the shape space */
   double sum_cell[number];
   /** @brief  saves the last value of last output for logdata output to files */
   double oldsum_cell[number];

   // vector of inidices at which external output cells do exist
   vector<long int> external_cells;
   // Asks whether this point in AffinitySpace is getting soutext cells for the first time. If yes,
   // add to external_cells
   void set_external_cell(long int pos);

   // vector of inidices at which antibody producers do exist
   vector<long int> ab_producers;

   // Return a random position in shape space
   long int get_random_SS_position();
   // or restricted to specified distances to the best clone
   long int get_random_SS_position(double mindist, double maxdist);

   // Metric
   short int measure;
   // Width of Gauss-Affinity squared
   double width2, amplitude;

   // Berechne den Index im Feld pop zu den Koordinaten k
   // Funktioniert fuer Gitter mit gleichen Ausmassen pro Dimension
   long int Index(long int * k);

   // Berechne die Koordinaten zu index und speichere sie in k
   // Funktioniert fuer Gitter mit gleichen Ausmassen pro Dimension
   void get_koord(long int index, long int * k);

   /* Berechnet den euklidischen Abstand der beiden Punkte k,l (n,n2) im
    * Shapespace (eventuell mit Beruecksichtigung periodischer
    * Randbedingungen) */
   double euklid(long int * k, long int * l);
   double N_mutation(long int * k, long int * l);

   // Berechne die Indizes der Nachbarn zu k und speichere sie in wert.neighbor
   /* Dabei werden die "linken" Nachbarn in 0,..,d-1 und die "rechten"
    * in d,..,2d-1 gespeichert. */
   /* Funktioniert fuer Gitter mit gleichen Ausmassen pro Dimension */
   void fix_neighbors(const long int &wo);
   void fix_neighbors(const long int &wo, long * nextn);

   // Ausgabe Dateien
   ofstream logdata[SSlogs];
   ofstream logmeanaff, loghighaff, log048aff, logdiversity;
   void get_diversity(double time);

   double pm_differentiation_rate;

   void getStatistics(double &_OUT_haffinity,
                      double &_OUT_steepness,
                      double &_CB_haffinity,
                      double &_CC_haffinity) {
      _OUT_haffinity = OUT_haffinity;
      _OUT_steepness = OUT_steepness;
      _CB_haffinity = CB_haffinity;
      _CC_haffinity = CC_haffinity;
   }
};

#endif
