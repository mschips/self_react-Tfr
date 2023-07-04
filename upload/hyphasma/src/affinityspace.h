#ifndef AFFINITYSPACE_H
#define AFFINITYSPACE_H

#include <string>
#include <vector>
#include <fstream>
#include "setparam.h"

/** @brief type of informations to store for each point in the affinity space :
 *  Note : the states sCCxxx correspond to the states defined in cellthis.h:
 *  enum centrocytes {unselected,contact,FDCselected,TCcontact,selected,apoptosis}; 
 *  The sequence of the states is essential because there are +-operations done
 *  in the code. If you add, add at the end before "number".
 */
enum cells {
   sCB,              // Number of CB
   sCC,              // Sum of all the CC
   sCCunselected,       // in particular, number of unselected CC
   sCCcontact,          // in particular, number of CC that contact a FDC but didn't take antigen
                        // yet
   sCCFDCselected,      // in particular, CC selected from the FDC and waiting for T cell help
   sCCTCcontact,        // in particular, CC in contact with a TC but not recycled yet ?
   sCCselected,         // in particular, CC selected for recycling ?
   sCCapoptosis,     // the total number of cells on the grid (not all dead)
   sallapoptosis,    // all cells that died
   sFDC,
   sTcell,
   sout,             // Zahl der output-Zellen die erzeugt wurden (alle)
   soutext,          // Zahl der ausgeworfenen output-Zellen (ohne Raumgitter)
   soutextproduce,   // Zahl der ausgeworfenen output-Zellen, die Antikoerper produzieren
   total,            // Zahl lebender CB,CC,out im AffinitySpace auf und ausserhalb des Gitters --> I DONT GET THIS ONE: CB-->CC-->OUT no total is removed(multiple counting?)
   soutdec,          // Zahl der DEC205+ output-Zellen, die erzeugt wurden (wie sout nur dec)
    //MSchips
    sCCTFRcontact,
    soutSelf,
    soutNonSelf,
//    sCCself, // number of self reactive CC --non apoptotic --> sthg wrong check before using!!
   number            // Number of cell types
};

const short int SSlogs = number;
// (Seb) Having public fields in abstract classes kind of defeats the purpose !

class AffinitySpace {
  public:
   AffinitySpace() { }
   virtual ~AffinitySpace() { }

   // ================== Getting or adding sequences =================

   /** @brief Adds a new antigen to a running system; returns success ### to be done:
    *  Note that adding a new antigen on AffinitySpace is not sufficient.
    *  It has also to be distributed on the FDCs  */
   virtual bool add_new_Antigen() = 0;              // antigen chosen by AffinitySpace
   virtual bool add_new_Antigen(long ag_index) = 0; // antigen index provided by the calling routine

   /** @brief Return the number of currently active antigens: */
   virtual int get_n_Antigen() = 0;

   /** @brief Gives the index of one of the antigens, randomly chosen from the predefined pool*/
   virtual long int get_Antigen() = 0;

   /** @brief Gives the index of antigen nr i */
   virtual long int get_Antigen(int i) = 0;

   /** @brief Return the index of the antigen nearest to the point/sequence <n>:*/
   virtual int get_nearest_Antigen(long n) = 0;

   /** @brief Gives the index of one of the seeding positions for cells,
    *  randomly chosen from the pre-defined pool*/
   virtual long int get_Seeder() = 0;

   /** @brief Gives the index of seeding position nr i among n_Seeders */
   virtual long int get_Seeder(int i) = 0;

   // ================== Mutation ! =================

   /** @brief returns the ID of a new, mutated position/sequence */
   virtual long getMutation(long from_position) = 0;

   // ================== Affinities =================
   // FOR LATER EXTENSIONS, PLEASE ALWAYS CALL AFFINITIES(ANTIBODY, ANTIGEN) in this order

   /** @brief  Get binding affinity between two points. variant1 : using gaussian function
    * (amplitude, width) */
   virtual double affinity(long int n, long int m) = 0;
   /** @brief Affinity, with a threshold */
   virtual double affinity(long int n, long int n2, double &tr) = 0;
   /** @brief Normalized affinity to maximum 1 (in case amplitude is not 1), for analysis */
   virtual double affinity_norm(long int n, long int m) = 0;

   /** @brief Gives the best affinity from position <pos> to the Ags, uses best_affinity_norm*/
   virtual double best_affinity(long pos) = 0;
   virtual double best_affinity_norm(long pos) = 0;

   /** @brief returns square distance between two points/sequences (hamming in case of sequences */
   virtual double Abstandquad(long int &n, long int &n2) = 0;

   //Mschips
   virtual void meanSD_out_affinity(double * aff_mean, double * aff_sd) =0;
   /** @brief returns the best distance fromposition <pos> to the Ags */
   virtual double best_distance(long pos) = 0;

   /** @brief Gives the affinity of position <pos> to antigen with list-index <i> */
   // virtual double get_affinity2ag(long pos, int i) = 0; --> replaced by affinity_norm directly in
   // the code

   // ================== Nb of cells and antibodies per point =================

   /** @brief  From a simulation, ask to change the numbers of cells at different AffinitySpace
    * points */
   virtual void add_cell(cells typ, long int pos) = 0;
   virtual void rem_cell(cells typ, long int pos) = 0;

   /** @brief plasma cell differentiation and antibody production */
   virtual void PM_differentiate(cells a, cells b, double dt) = 0;

   /** @brief To put or read antibody amounts for each position in AffinitySpace
    * This uses the affinity space as a container to store the amount of produced antibody
    * for each sequences. Note: the same sequence might appear multiple times in the list.*/
   virtual void put_Ab(long index,double d_ab) = 0;
   virtual double get_AbAmount(long index) = 0;

   /** @brief Sum of each kind of cells over the AffinitySpace */
   virtual double get_sum_cell(cells celltype) = 0;

   /** @brief  saves the last value of last output for logdata output to files */
   virtual double get_oldsum_cell(cells celltype) = 0;

   // ================== List of points matching to AB producing cells =================

   /** @brief asks to add this sequence ID in the list of external cells.
    * If it is not already in, add it */
   virtual void set_external_cell(long int pos) = 0;
   /** @brief Number of positions holding producers */
   virtual int get_n_ab_producers() = 0;
   /** @brief Index of the position nr n in the list */
   virtual double get_AbProducingCells(int n) = 0;
   /** @brief Amount of antibody producing cells at this position */
   virtual long get_AbProducingIndex(int n) = 0;

   // ================== Functions for analyses  =================

   /** @brief side functions to observe what is happening inside the shape space */
   virtual void correct_average_affinity(cells celltyp, long &pos, double &average) = 0;

   /** @brief returns average affinity of each cells, and fills the table affinities with :
    *   affinity[0,1,2] = average affinities of cells from [sCB, sCC, sout] to have aff > 0.3
    *   affinity[3,4,5] = % of cells from [sCB, sCC, sout] to have aff > 0.3
    *   affinity[6,7,8] = % of cells from [sCB, sCC, sout] to have aff < 0.4
    *   affinity[9,10,11] = % of cells from [sCB, sCC, sout] to have 0.4 =< aff < 0.8
    *   affinity[12,13,14] = % of cells from [sCB, sCC, sout] to have aff >= 0.8  */
   // virtual double mean_affinity(double * affinities) = 0; --> now private

   /** @brief  to check the sum of cells in each space position is consistent with the total */
   virtual short int sum_check() = 0;

   /** @brief  writes a sum up in the ofstreams */
   virtual void to_ssfiles(double time) = 0;
   virtual void to_multiag_files(double time) = 0;

   /** @brief Functions to export affinity / hamming profiles */
   virtual void write_gcbc_hamming(double time) = 0;
   virtual void write_gcbc_affinity(double time) = 0;

   /** @brief To close the ofstreams */
   virtual void close_files() = 0;

   /** @brief to access the private fields OUT_ ... by asking a copy */
   virtual void getStatistics(double &_OUT_haffinity,
                              double &_OUT_steepness,
                              double &_CB_haffinity,
                              double &_CC_haffinity) = 0;
};
#endif
