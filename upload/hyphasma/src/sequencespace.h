#ifndef SEQUENCESPACE_H
#define SEQUENCESPACE_H

#include "affinityspace.h"
#include <sstream>

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         1 - A general class extending AffinitySpace, to handle <sequences> of BCRs, antigens and
// TCRs
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief A subclass of AffinitySpace to handle lists of sequences (boolean or Arup),
 *  and is therefore shared between ArupSpace and SequenceSpace (as mother class).
 *  The functions handling communication with affinity are however specific to  Sequence Space vs
 *ArupSpace.*/
/** Constraints that a sequence (of type T = whatever) should have to be able to use this class :
 *  - a field n_cells[]
 *  - a field antibody */
template <typename T>
class generalSequenceContainer: public AffinitySpace {
  public:
   generalSequenceContainer();
   T* getSequence(long id);
   int get_n_Seeder();
   long int get_Seeder();
   long int get_Seeder(int i);
   int get_n_Antigen();
   long int get_Antigen();                 // cannot change it, reimplemented for mother class
   long int get_Antigen(int i);            // cannot change it, reimplemented for mother class
   int get_n_TCR();
   long int get_TCR();
   long int get_TCR(int i);

   /** @brief To put or read antibody amounts for each position in AffinitySpace */
   void put_Ab(long index,double d_ab);
   double get_AbAmount(long index);

   /** @brief Updating the field .n_cell[] of a sequence stored. Calls set_external_cell
    * automatically */
   void add_cell(cells typ, long int pos);
   void rem_cell(cells typ, long int pos);
   
   int get_n_ab_producers();
   double get_AbProducingCells(int n);
   long get_AbProducingIndex(int n);
   void PM_differentiate(cells a, cells b, double dt);
   
   double get_sum_cell(cells celltype);
   double get_oldsum_cell(cells celltype);
   short int sum_check();
   
  protected:
   vector<T*> poolSequences;
   vector<long> indexParentSeq;
   long n_Sequences;
   vector<long> Seeders;
   int n_Seeder;
   vector<long> Antigens;
   int n_Antigen;
   vector<long> TCRs;
   int n_TCRs;

   vector<long int> external_cells;
   vector<long int> ab_producers;
   void set_external_cell(long int pos);

   static double pm_differentiation_rate;

   double sum_cell[number];
   double oldsum_cell[number];

   /// Remaining virtual functions */
   /*virtual int get_nearest_Antigen(long n) = 0;
    *  virtual long getMutation(long from_position) = 0;
    *  virtual double affinity(long int n, long int m) = 0;
    *  virtual double best_affinity(long pos) = 0;
    *  virtual double best_affinity_norm(long pos) = 0;
    *  virtual double get_affinity2ag(long pos, int i) = 0; // No !
    *  virtual double Abstandquad(long int &n, long int &n2) = 0;
    *  virtual void correct_average_affinity(cells celltyp, long &pos, double &average) = 0;
    *  //virtual double mean_affinity(double * affinities) = 0;
    *  virtual void to_ssfiles(double time) = 0;
    *  virtual void to_multiag_files(double time) = 0;
    *  virtual void write_gcbc_hamming(double time) = 0;
    *  virtual void write_gcbc_affinity(double time) = 0;
    *  virtual void close_files() = 0;
    *  virtual void getStatistics(double &_OUT_haffinity,
    *                          double &_OUT_steepness,
    *                          double &_CB_haffinity,
    *                          double &_CC_haffinity) = 0;*/
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         2 - SEQUENCES IN SAHAMS AFFINITY MODEL
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum typeSeq {
   typeBCR, typeTCR, typeAG, numberTypesSeq
};
enum typeAffFunction {
   seqAff, seqAffNorm, seqAffWindow, numberTypesAffFunctions
};

/** @brief A common mother class to all boolean sequences.
 *  The subclasses will be :
 *     - TCRs,
 *     - BCRs,          with the additional property: nb_AB_producers()
 *     - Antigens,      */
struct sequence {
   /** @brief Constructor : sequence of '0's of this size */
   sequence(int _size);
   /** @brief Constructor : from a string of '0' and '1'. Other characters become '0' */
   sequence(string boolText);
   /** @brief Constructor : copies; the other sequence can be deleted without problem) */
   sequence(sequence * toCopy);
   /** @brief Constructor : the long number that will be converted to binary. Max 2^63 */
   sequence(int _size, long ID);
   virtual ~sequence() { }
   void operator =(sequence B) { content = B.content; size = B.size; }

   /** @brief boolean storage. Could not use bitset : size unknown before compiling */
   std::vector<bool> content;           // alternately, vector<short>
   int size;                            // alternately, could be static

   /** @brief Reading sequence at position i via [i]. seq[i] = xxx is not possible. */
   int operator [](int i);

   /** @brief Number of cells of each type carrying this sequence. 0 for antigens*/
   double n_cell[number];
   /** @brief the best affinity of BCRs and TCRs to antigens. Will be calculated once when a
    *  sequence is created and updated each time the antigen pool is changed. 0 For antigens.*/
   double max_affinity_to_antigens;
   /** @brief Amount of antibodies produced for this sequence. 0 for TCR or antigens. */
   double antibody;

   /** @brief flips every position with a probability 0.5 */
   void randomize();

   /** @brief Pick randomly a position and flips it */
   void mutateOnePosition();

   /** @brief compare content[..] of two sequences ONLY. (use : sequence::compSeq(a,b)) */
   static bool compSeq(sequence * a, sequence * b);

   /** @brief affinity between two sequences (from their pointers), with exponent r
    *  formula = [sum((matches_sizes)^r)] / total_size^r.  (use : sequence::seq_affinity(a,b))  */
   static double seq_affinity(sequence * x, sequence * y, double r,
                              int maxSizeClusters, int type_affinity_function);

   /** @brief hamming distance (nb of mutations) between two sequences from their pointers */
   static double hamming(sequence * x, sequence * y);

   // ===== Functions reimplemented by subclasses(BCR, TCR or Antigen) =====

   /** @brief type of the subclass */
   virtual typeSeq getType() { return numberTypesSeq; }

   virtual double nb_AB_producers() {
      cerr << "ERR : 'add_producer' called from non antibody : " << typeInString() << endl;
      return 0;
   }
   virtual bool hasExternalCell() {
      cerr << "ERR : 'hasExternalCell' called from non antibody : " << typeInString() << endl;
      return 0;
   }
   // ===== Tool functions =====

   static string testeAffinityFunctions(double L,
                                        double R,
                                        int maxClusters,
                                        int typeAffinityFunction);
   static string typeCellInString(cells index);
   string printNbCells();
   string typeInString();
   virtual std::string print();
};

struct BCR: public sequence {
   BCR(int size) : sequence(size) { }
   BCR(string txt) : sequence(txt) { }
   BCR(BCR * anotherBCR) : sequence((sequence*) anotherBCR) { }
   BCR(sequence * anotherSeq) : sequence(anotherSeq) { }
   double nb_AB_producers() { return n_cell[soutextproduce]; }
   bool hasExternalCell() { return (n_cell[soutext] + n_cell[soutextproduce]) > 0; }
   string print() {
      stringstream res;
      res << "BCR:" << sequence::print() << "\tNbProducers:\t" << nb_AB_producers();
      return res.str();
   }
   typeSeq getType() { return typeBCR; }
   ~BCR() { }
};

struct TCR: public sequence {
   TCR(int size) : sequence(size) { }
   TCR(string txt) : sequence(txt) { }
   TCR(TCR * anotherTCR) : sequence((sequence*) anotherTCR) { }
   TCR(sequence * anotherSeq) : sequence(anotherSeq) { }
   string print() {
      stringstream res;
      res << "TCR:" << sequence::print();
      return res.str();
   }
   typeSeq getType() { return typeTCR; }
   ~TCR() { }
};

struct Antigen: public sequence {
   // static int number_of_bins;
   Antigen(int size) : sequence(size) { }
   Antigen(string txt) : sequence(txt) { }
   Antigen(Antigen * anotherAntigen) : sequence((sequence*) anotherAntigen) { }
   Antigen(sequence * anotherSeq) : sequence(anotherSeq) { }
   string print() {
      stringstream res;
      res << "Antigen:" << sequence::print();
      return res.str();
   }
   typeSeq getType() { return typeAG; }
   ~Antigen() { }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         3 - SEQUENCES IN ARUPS AFFINITY MODEL
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @Brief Class to generate random numbers from an arbitrary distribution,
 *  defined by a liste of values (classes) and their densities
 *  THE VALUES OF THE CLASSES SHOULD BE ORDERED INCREASING */
struct probaLawFromTable {
   int size;
   vector<double> Xs;
   vector<double> NormalizedDensities;      // to probabilities
   vector<double> cumulatedValuesLow;       // cumulated probabilities
   double normalizingCoeff;                 // just for remembering the normalisation

   bool uniformizeInsideClasses;
   vector<double> lowBoundsXs;
   vector<double> highBoundsXs;

   virtual double getRandValue();
   int getRandClass();
   probaLawFromTable(vector<double> &_Xs, vector<double> &densities, bool _uniformizeInsideClasses);
   static double IndependentRandRealGen();  // If later someone wants to separate generators for
                                            // seeding.
   virtual string print();
};

/** @Brief Class to test/control the behavior of probaLawFromTable.
 *  by being daughter cell, reimplements the getRandomValue function and
 *  does statistics on how the mother class behaves */
struct probaLawFromTableStoringResults: public probaLawFromTable {
   vector<int> nbOfTimesPerClass;
   int totalNbEvents;

   probaLawFromTableStoringResults(vector<double> &_Xs,
                                   vector<double> &densities,
                                   bool _uniformizeInsideClasses);
   double getRandValue();
   double getFrequency(int classIndex);
   int fitDoubleInAClass(double val);
   void clearRecord();
   string print();
   static string TestprobaLawFromTable();
};

/* List of parameters for using Arup sequences
 *
 *  Use arup space (1/0) [use_arup_space]:
 *  0
 *  Length of sequences [arup_length_sequences]:
 *  46
 *  Nb conserved residues [arup_N_conserved]:
 *  18
 *  Nb mutated residues [arup_N_mutates]:
 *  22
 *  Nb shielded residues (the rest) [arup_N_shielded]:
 *  6
 *  Number of Initial Antigen sequences (int-type) [arup_nb_ini_antigens]:
 *  3
 *  Fix Antigen Sequence presentation. Order is Conserved(should be 1), Variable, shielded (should
 * be 1). max 1000 values [arup_ini_antigens[...]]:
 *  1111111111111111111111111111111111111111111111
 *  1111111111111111110000000000011111111111111111
 *  1111111111111111111111111111100000000000111111
 *  -1
 *  Fraction of Arup Ags (non-fixed Ag enter with same fraction) [arup_ag_fraction[...]]:
 *  -1
 *  Number of mutations for mutated strains (if more initial Ag have to be generated)
 * [arup_nb_mutations_gen_strains]:
 *  11
 *  Activation threshold (kcal/mol) [arup_threshold_activation]:
 *  10.8
 *  Initial interval for BCR sequences [arup_h_min, arup_h_max]:
 *  -0.18
 *  0.9
 *  Fix initial Repertoire distribution (max 1000 values) [arup_ini_bcrs[...]]:
 *  -1
 #would be 0.1 0.8 0.5 -0.2 0.1 ... with values
 *  Mutation rate per sequence per division per residue [arup_mutation]:
 *  0.003
 *  probability of a mutation being lethal [arup_proba_lethal_mut]:
 *  0.3
 *  probability of a mutation being affecting [arup_proba_affecting_mut]:
 *  0.2
 *  probability silent [arup_proba_silent_mut]:
 *  0.5
 *  distribution of affinity changes by mutation [nbLinesToRead, arup_law_mut_Xs[...],
 * arup_law_mut_Densities[...]]:
 *  40
 *  -6,4	0,019047619
 *  -6,2	0,003809524
 *  -6	0,003809524
 *  -5,8	0,003809524
 *  -5,6	0,003809524
 *  -5,4	0,024380952
 *  -5,2	0,042666667
 *  -5	0,028952381
 *  -4,8	0,03352381
 *  -4,6	0,03352381
 *  -4,4	0,058666667
 *  -4,2	0,084571429
 *  -4	0,089142857
 *  -3,8	0,054857143
 *  -3,6	0,09447619
 *  -3,4	0,095238095
 *  -3,2	0,079238095
 *  -3	0,14552381
 *  -2,8	0,11352381
 *  -2,6	0,134095238
 *  -2,4	0,095238095
 *  -2,2	0,179809524
 *  -2	0,234666667
 *  -1,8	0,089142857
 *  -1,6	0,134095238
 *  -1,4	0,204190476
 *  -1,2	0,290285714
 *  -1	0,249904762
 *  -0,8	0,215619048
 *  -0,6	0,170666667
 *  -0,4	0,340571429
 *  -0,2	0,361142857
 *  0	0,365714286
 *  0,2	0,316190476
 *  0,4	0,324571429
 *  0,6	0,043428571
 *  0,8	0,084571429
 *  1	0,048761905
 *  1,2	0,04952381
 *  1,4	0,03352381
 *  Coefficient of unshielding [arup_alpha]:
 *  2
 *  Bounding of the shielding h' [arup_hprime_min, arup_hprime_max]:
 *  -1.5
 *  1.5
 *  Bounding of the residues facing mutated h' [arup_hmut_min, arup_hmut_max]:
 *  -1.5
 *  1e6 */

// an antigen is a sequence of +1/-1. A nn antiboy is a sequence of double (h values)
enum typeArupSeq {
   typeArupBCR, typeArupAG, numberTypesArup
};

/** @brief A class for Arup arupProteins (double vectors).*/
struct arupProtein {
   static int arup_length_sequences, arup_N_conserved, arup_N_mutates, arup_N_shielded;
   static double arup_alpha, arup_h_min, arup_h_max, arup_hprime_min, arup_hprime_max,
                 arup_hmut_min, arup_hmut_max;
   /** @brief to store the distribution of affinities following mutation */
   static probaLawFromTable * lawMutations;
   static void set_statics(Parameter &par, ofstream &ana);

   /** @brief Constructor : arupProtein of '1's of this size */
   arupProtein();
   /** @brief Constructor : arupProtein with the values from this vector */
   arupProtein(vector<double> numbers);
   /** @brief Constructor : arupProtein from a string of double values */
   arupProtein(string numbersInText);
   /** @brief Constructor : copies the content. the other arupProtein can be deleted */
   arupProtein(arupProtein * toCopy);
   /** @brief to copy from another arupProtein */
   void copy(arupProtein * toCopy);
   /** @brief Assigning arupProtein content using '=' */
   void operator =(arupProtein B) { content = B.content; size = B.size; }
   /** @brief Loses memory. Similar to Alzheimer. */
   void clear();

   virtual ~arupProtein();

   /** @brief Storage of the arupProtein */
   std::vector<double> content;
   int size;                            // alternately, put it static,
   /** @brief Accessing a position of the arupProtein via [i]. Note that seq[i] = xxx is not
    * possible. need to say seq.content[i] = xxx */
   int operator [](int i);

   /** @brief compare content[..] of two arupPoteins ONLY. (use : sequence::compSeq(a,b)) */
   static bool compSeq(arupProtein * a, arupProtein * b);

   /** @brief function to give the affinity between two arupProteins (from their pointers) */
   static double arup_affinity(arupProtein * BCR, arupProtein * Antigen);

   /** @brief sum of distances between two arupProteins from their pointers */
   static double continuoushamming(arupProtein * x, arupProtein * y);

   /** @brief Number of cells of each type carrying this sequence. 0 for antigens*/
   double n_cell[number];
   /** @brief the best affinity of BCRs and TCRs to antigens. Will be calculated once when a
    *  sequence is created and updated each time the antigen pool is changed. 0 For antigens.*/
   double max_affinity_to_antigens;
   /** @brief Amount of antibodies produced for this sequence. 0 for TCR or antigens. */
   double antibody;

   // ===== Functions reimplemented by subclasses(BCR, TCR or Antigen) =====

   /** @brief type of the subclass */
   virtual typeArupSeq getType();

   virtual double nb_AB_producers() {
      cerr
      << "ERR : you are calling 'add_producer' from a non antibody arupProtein : "
      << typeInString() << endl;
      return 0;
   }

   virtual bool hasExternalCell() {
      cerr
      << "ERR : calling hasExternalCell from a non-antibody : " << typeInString() << endl;
      return false;
   }

   // ===== Tool functions =====

   // functions reserved for antigens. (implemented but not used)
   // virtual void add_producer(double affinity, int nbToAdd = 1){cerr << "ERR : you are calling
   // 'add_producer' from a non antigen arupProtein" << typeInString() << endl;} // because is
   // reimplemeneted for Antigens only
   // virtual void rem_producer(double affinity, int nbToRemove = 1){cerr << "ERR : you are calling
   // 'rem_producer' from a non antigen arupProtein" << typeInString() << endl;} // because is
   // reimplemeneted for Antigens only

   /** @brief Test Function to show the properties of the different affinity functions */
   static string testeAffinityFunctions();

   // function to print the types of cells as string (used for printing in the next function)
   string typeCellInString(cells index);
   string printNbCells();

   /** @brief For printing / debugging */
   string typeInString();
   /** @brief Returns arupProtein as a string of '0' and '1'*/
   virtual std::string print();
};

struct ArupBCR: public arupProtein {
   ArupBCR() : arupProtein() { }
   ArupBCR(string txt) : arupProtein(txt) { }
   ArupBCR(ArupBCR * anotherBCR) : arupProtein((arupProtein*) anotherBCR) { }
   ArupBCR(arupProtein * anotherSeq) : arupProtein(anotherSeq) { }
   double nb_AB_producers() { return n_cell[soutextproduce]; }
   bool hasExternalCell() { return (n_cell[soutext] + n_cell[soutextproduce]) > 0; }
   string print() {
      stringstream res;
      res << "ArupBCR:" << arupProtein::print() << "\tNbProducers:\t" << nb_AB_producers();
      return res.str();
   }
   typeArupSeq getType() { return typeArupBCR; }
   void mutate(double targetChangeEnergy);
   void randomize(double targetEnergy, bool beExact = false);
};

struct ArupAntigen: public arupProtein {
   ArupAntigen() : arupProtein() { }
   ArupAntigen(string txt) : arupProtein(txt) { }
   ArupAntigen(ArupAntigen * anotherAntigen) : arupProtein((arupProtein*) anotherAntigen) { }
   ArupAntigen(arupProtein * anotherSeq) : arupProtein(anotherSeq) { }
   string print() {
      stringstream res;
      res << "ArupAntigen:" << arupProtein::print();
      return res.str();
   }
   typeArupSeq getType() { return typeArupAG; }
   void randomize(int nbResiduesToMutate = -1);
   void mutate();  // just one position
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         4 - (BOOLEAN/SAHAMs) SEQUENCE SPACE
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/** @brief Parameters necessary for a sequence space :
 *  @param par   (Parameter) global list of parameters
 *  @param ana   (ofstream)  to output what's happening to the shape space
 *
 *  important parameters, inside 'par' are (kenntext, access, meaning :)
 *
 *  From the sequence space paragraph
 *       - Length of sequences:                                   par.value.size_sequences
 *       - Number of Initial Antigen sequences (int-type):
 *         par.Value.init_antigen_sequences
 *       - Fix Antigen Sequence presentation (max 1000 values):   string par.Value.initAntigenSeqs[]
 *       - Maximum Hamming distance between antigens:             par.Value.max_hamming_antigens
 *       - Minimum Hamming distance between antigens:             par.Value.min_hamming_antigens
 *
 *  From the CB antibody type paragraph
 *       - Total Number of initial B-Cells types:                 par.Value.totalBss
 *
 *  From the sequence space paragraph (again)
 *       - Fix initial Repertoire distribution (max 1000 values): string par.Value.initBCRSeqs[]
 *       - Maximum Initial Hamming distance between BCRs:         par.Value.max_hamming_BCRs
 *       - Minimum affinity of initial BCRs to Antigens:
 *            par.Value.min_initial_affinity_BCRs
 *       - Maximum affinity of initial BCRs to Antigens:
 *            par.Value.max_initial_affinity_BCRs
 *       - Fix initial Repertoire distribution for T cells:       string par.Value.initTCRSeqs[]
 *       - Maximum Initial Hamming distance between TCRs:         par.Value.max_hamming_TCRs
 *       - Minimum affinity of initial TCRs to Antigens:
 *            par.Value.min_initial_affinity_TCRs
 *       - Maximum affinity of initial TCRs to Antigens:
 *            par.Value.max_initial_affinity_TCRs
 *       - Specifity of sequences affinity (double R):            par.Value.R_affinity
 *
 *  Other parameters :
 *       pm_differentiation_rate                     --> the rate of cells converting into AB
 *producers
 */

class sequenceSpace: public generalSequenceContainer<sequence> {
  public:
   sequenceSpace(Parameter &par, ofstream &ana);
   ~sequenceSpace();

   /** @brief Add a new sequence (its pointer) to the pool and gets its new ID (position in
    * poolSequences). this sequence should already be of type BCR, TCR or Antigen.
    * NOTE :  - when a BCR or TCR is added, its affinity to antigens is automatically computed
    *         - when an antigen is added, the affinities of all BCRs for antigens are updated */
   long index_adding_sequence(sequence * toAdd, long int ID_of_mother_sequence = -1);

   /** @brief dynamically add new seeding sequences from a string or a BCR* (will update n_Seeder).
    * */
   long int add_Seeder(string newSequence);
   long int add_Seeder(BCR * newSequence);

   /** @brief dynamically add new Antigen sequences from a string or an Antigen* (will update
    * n_Antigen).
    * note : will lead to the re-computation of affinities for all BCRs and TCRs to antigens */
   long int add_Antigen(string newSequence);
   long int add_Antigen(Antigen * newSequence);
   bool add_new_Antigen(); // ### to be programmed   // Michael: This might be redundant?
   bool add_new_Antigen(long ag_index); // ### to be programmed

   /** @brief dynamically add new TCR sequences from a string or a TCR* (will update n_TCRs). */
   long int add_TCR(string newSequence);
   long int add_TCR(TCR * newSequence); // {return generalSequenceContainer::add_TCR((sequence*)
                                        // newSequence);}

   /** @brief Generates a new BCR sequence in poolSequences, and returns its ID.
    * (note : the mother of this new sequence is automatically updated as to from_id_sequence) */
   long getMutation(long from_ID_sequence);

   // PLEASE ALWAYS CALL AFFINITY(ANTIBODY, ANTIGEN) in this order

   /** @brief basic affinity. For boolean sequences, it is symmetric. */
   double affinity(long int seqId1, long int SeqId2);
   double affinity_norm(long int seqIAntibody, long int SeqIdAntigen);
   /** @brief affinity normalized by a threshold */
   double affinity(long int seqId1, long int seqId2, double &tr);

   /** @brief gives the highest affinity (in both cases)*/
   double best_affinity(long IDSeq);
   double best_affinity_norm(long pos);
   /** @brief computes the average, using affinity_norm */
   void correct_average_affinity(cells celltyp, long &pos, double &average);

   /** @brief affinity of a sequence to all antigens and returns the maximum */
   double maxAffinitySequenceToAntigens(sequence * seq);

   /** @brief returns square distance between two points/sequences (hamming in case of sequences */
   double Abstandquad(long int &n, long int &n2);

   /** @brief Gives back the Hamming distance (not squared) between two shape space vectors */
   int intN_mutation(sequence *x, sequence *y);
   /** @brief Gives back the best Hamming distance of SS-position n to any of the antigens */
   int best_hamming_distance(long n);
   /** @brief Gives back the mean Hamming distance of SS-position n to all antigens */
   double mean_hamming_distance(long n);
   /** same as best_affinit_norm(long) but returning the mean of affinities to all antigens */
   double mean_affinity_norm(long pos);

   double N_mutation(sequence * k, sequence * l);

   /** @brief (For Analysis) returns the number of the closest antigen (0..n_Antigen-1).
    * To get its sequence space index, use get_Antigen(get_nearest_antigen( )). */
   int get_nearest_Antigen(long n);

   string printSequences(bool showSequencesWithNoAliveCells = false, bool showTree = false);
   static void testSequenceSpace();
   void to_ssfiles(double time);
   void to_multiag_files(double time);
   void write_gcbc_hamming(double time);
   void write_gcbc_affinity(double time);
   void close_files();
   void meanSD_out_affinity(double * aff_mean, double * aff_sd);
   double best_distance(long pos);
   
  private:
   int size_sequences;
   double R_affinity;
   int maxSizeClusters;
   int typeAffinityFunction;
   int use_logarithmic_seq_affinity;
   double amplitude;
   /** @brief  Indicators that are updated when to_ss files is called */
   double OUT_haffinity,OUT_steepness,CB_haffinity,CC_haffinity;
   /**  */

   // double pm_differentiation_rate; // from sequenceContainer now
   /** @brief fields for writing data into files */
   ofstream logdata[number];
   ofstream logmeanaff,loghighaff,log048aff,logdiversity;

   double mean_affinity(double * affinities); // ### Michael made this private
   void mean_affinities_ag(double * mean_affinities, int ag_index);
   void get_diversity(double time);
   void getStatistics(double &_OUT_haffinity, double &_OUT_steepness,
                      double &_CB_haffinity, double &_CC_haffinity) {
      _OUT_haffinity = OUT_haffinity;
      _OUT_steepness = OUT_steepness;
      _CB_haffinity = CB_haffinity;
      _CC_haffinity = CC_haffinity;
   }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         5 - ARUP SPACE
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*Potentially additional parameters :
 *   Probability of survival is threshold is reached     arup_Pa
 *   Probability get Th help                             arup_Pth
 *   Recycling probability                               0.9     */

class arupSpace: public generalSequenceContainer<arupProtein> {
  public:
   arupSpace(Parameter &par, ofstream &ana);
   ~arupSpace();

   long index_adding_arupProtein(arupProtein * toAdd, long int ID_of_mother_arupProtein = -1);

   long int add_Seeder(string newarupProtein);
   long int add_Seeder(ArupBCR * newarupProtein);

   /** @brief ONLY APPLIES TO BCRS */
   long getMutation(long from_ID_arupProtein);

   /** @brief returns the number of the closest antigen (0..n_Antigen-1). To get its arupProtein
    * space index, use get_Antigen(get_nearest_antigen( )). */
   int get_nearest_Antigen(long n);

   /** @brief possible to dynamically add new Antigen arupProteins from a string or a Antigen* (will
    * update n_Antigen).
    * note : adding an antigen will lead to the re-computation of affinities for all BCRs and TCRs
    *to antigens */
   long int add_Antigen(string newarupProtein);
   long int add_Antigen(ArupAntigen * newarupProtein);
   bool add_new_Antigen();                // antigen chosen by AffinitySpace
   bool add_new_Antigen(long ag_index);   // antigen index provided by the calling routine

   // PLEASE ALWAYS CALL AFFINITIES(ANTIBODY, ANTIGEN) in this order
   /** @brief basic affinity. NOT symmetric. */
   double affinity(long int seqIAntibody, long int SeqIdAntigen);
   double affinity_norm(long int seqIAntibody, long int SeqIdAntigen);
   /** @brief affinity normalized by a threshold */
   double affinity(long int seqId1, long int seqId2, double &tr);

   /** @brief gives the highest affinity (in both cases)*/
   double best_affinity(long IDSeq);
   double best_affinity_norm(long pos);

   /** @brief Gives back the Hamming distance (not squared) between two shape space vectors */
   int intN_mutation(arupProtein * k, arupProtein * l);
   /** @brief Gives back the best Hamming distance of SS-position n to any of the antigens */
   int best_hamming_distance(long n);
   /** @brief Gives back the mean Hamming distance of SS-position n to all antigens */
   double mean_hamming_distance(long n);
   /** same as best_affinit_norm(long) but returning the mean of affinities to all antigens */
   double mean_affinity_norm(long pos);

   double N_mutation(arupProtein * k, arupProtein * l);

   /** @brief computes the average, using affinity_norm */
   void correct_average_affinity(cells celltyp, long &pos, double &average);

   /** @brief affinity of a sequence to all antigens and returns the maximum */
   double maxAffinityarupProteinToAntigens(arupProtein * seq);

   /** @brief returns square distance between two points/sequences (hamming in case of sequences */
   double Abstandquad(long int &n, long int &n2);

   // ** @brief  Indicators that are updated when to_ss files is called */
   double OUT_haffinity,OUT_steepness,CB_haffinity,CC_haffinity;
   
   string printarupProteins(bool showarupProteinsWithNoAliveCells = false, bool showTree = false);
   void testarupSpace();
   void to_ssfiles(double time);
   void to_multiag_files(double time);
   void write_gcbc_hamming(double time);
   void write_gcbc_affinity(double time);
   void meanSD_out_affinity(double * aff_mean, double * aff_sd);
   double best_distance(long pos);
   void close_files();
   
  private:
   double mean_affinity(double * affinities);  // ### Michael made this private
   void mean_affinities_ag(double * mean_affinities, int ag_index);
   void get_diversity(double time);
   void getStatistics(double &_OUT_haffinity,
                      double &_OUT_steepness,
                      double &_CB_haffinity,
                      double &_CC_haffinity) {
      _OUT_haffinity = OUT_haffinity;
      _OUT_steepness = OUT_steepness;
      _CB_haffinity = CB_haffinity;
      _CC_haffinity = CC_haffinity;
   }
   
   static int arup_nb_ini_antigens;
   static double arup_threshold_activation;
   static double arup_mutation, arup_proba_lethal_mut, arup_proba_affecting_mut,
                 arup_proba_silent_mut;
   double amplitude;

   // ** @brief fields for writing data into files */
   ofstream logdata[number];  /// Philippe, use number instead of SSlogs
   ofstream logmeanaff,loghighaff,log048aff,logdiversity;

};

#endif // arupProteinSPACE_H
