#include "sequencespace.h"
#include <cmath>                // for std::pow, seems that math.h doesn't have it
#include "random.h"
#include <sstream>
#include <fstream>
using namespace std;

// necessary for tests:
#include <bitset>
#include <algorithm>
#include <iomanip>

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         1 - A general class extending AffinitySpace, to handle sequences of BCRs, antigens and
// TCRs
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Pre-declaration of the available sequence classes :
struct sequence;
struct arupProtein;

/** @brief Generally, templates can not be implemented in the .cpp because everything has to be
 * known from the .h.
 *  Found two ways to go around :
 *  1/ explicit say here which T types are allowed, then the implementation can be kept in the .cpp
 **/
template class generalSequenceContainer<sequence>;
template class generalSequenceContainer<arupProtein>;
/** 2/ (that could be done) : put the implementation of the template class in a xxx.tpp file and do
 * #include "xxx.tpp" at the end of the .h file */

template <typename T>
double generalSequenceContainer<T>::pm_differentiation_rate = NAN; // to raise errors if not
                                                                   // instanciated ...

template <typename T>
generalSequenceContainer<T>::generalSequenceContainer() {
   n_Sequences = 0;
   n_Seeder = 0;
   n_Antigen = 0;
   n_TCRs = 0;
   // pm_differentiation_rate = 0; // Has to be initialized by the subclass constructor
}
template <typename T>
T* generalSequenceContainer<T>::getSequence(long id) {
   if ((id < 0) || (id >= n_Sequences)) {
      cerr << "ERR: sequenceSpace::getSequence(id=" << id << "), incorrect id, only "
           << n_Sequences << " definde\n";
      return NULL;
   }
   return poolSequences[id];
}
template <typename T>
int generalSequenceContainer<T>::get_n_Antigen() { return n_Antigen; }
template <typename T>
int generalSequenceContainer<T>::get_n_Seeder() { return n_Seeder; }
template <typename T>
int generalSequenceContainer<T>::get_n_TCR() { return n_TCRs; }
template <typename T>
long int generalSequenceContainer<T>::get_Antigen() {
   return Antigens[irandom(n_Antigen)];
}
template <typename T>
long int generalSequenceContainer<T>::get_Seeder() {
   return Seeders[irandom(n_Seeder)];
}
template <typename T>
long int generalSequenceContainer<T>::get_TCR() {
   return TCRs[irandom(n_Seeder)];
}
template <typename T>
long int generalSequenceContainer<T>::get_Antigen(int i) {
   if ((i >= 0) && (i < n_Antigen)) { return Antigens[i]; } else { return -1; }
}
template <typename T>
long int generalSequenceContainer<T>::get_Seeder(int i) {
   if ((i >= 0) && (i < n_Seeder)) { return Seeders[i]; } else { return get_Seeder(); }
}
template <typename T>
long int generalSequenceContainer<T>::get_TCR(int i) {
   if ((i >= 0) && (i < n_TCRs)) { return TCRs[i]; } else { return get_TCR(); }
}
template <typename T>
void generalSequenceContainer<T>::put_Ab(long index,double d_ab) {
   getSequence(index)->antibody += d_ab;
   if (getSequence(index)->antibody < 0) {
      getSequence(index)->antibody = 0;
      cerr << "WARNING! antibody got negative at index " << index << " in Sequence/Arup space.\n";
   }
}
template <typename T>
double generalSequenceContainer<T>::get_AbAmount(long index) {
   return getSequence(index)->antibody;
}
template <typename T>
void generalSequenceContainer<T>::add_cell(cells typ, long int pos) {
   if ((pos < 0) || (pos >= n_Sequences)) {
      cerr << "ERR: sequenceSpace::add_cell(typ, IDseq), the IDseq is invalid. Only "
           << n_Sequences
           << " defined. You should make sure a sequence exists before adding a cell...\n";
      return;
   }
   if (typ == soutext) {
      set_external_cell(pos);    // has to be done before n_cell[soutext] is increased,
   }                           // to detect when it goes from 0 to another number

   if (typ == soutextproduce) {
      cerr
         <<
         "WRN: sequenceSpace::add_cell(soutextproduce), you are not supposed to add soutextproduce from the simulation. It should happen by differentiation by calling AffinitySpace::pm_differentiate()"
         << endl;
   }
   poolSequences[pos]->n_cell[typ] += 1;
   ++sum_cell[typ];
}
template <typename T>
void generalSequenceContainer<T>::rem_cell(cells typ, long int pos) {
   if ((pos < 0) || (pos >= n_Sequences)) {
      cerr << "ERR: sequenceSpace::rem_cell(typ, IDseq =" << pos << "), invalid ID. [0.."
           << n_Sequences - 1 << "]" << endl;
      return;
   }
   poolSequences[pos]->n_cell[typ] -= 1;
   --sum_cell[typ];
   // note :  external cells is not updated if soutext/soutextproduce are removed
}
template <typename T>
double generalSequenceContainer<T>::get_AbProducingCells(int n) {
   if ((n < 0) || (n >= (int) ab_producers.size())) {
      cerr << "ERR:sequenceSpace::getAbProducingCells(" << n
           << "), index out of bounds [0.." << ab_producers.size() - 1 << "]" << endl;
      return 0;
   }
   T * bseq = getSequence(ab_producers[n]);
   return bseq->nb_AB_producers();      // (returns ->n_cell[soutextproduce];
}
template <typename T>
long int generalSequenceContainer<T>::get_AbProducingIndex(int n) {
   if ((n < 0) || (n >= (int) ab_producers.size())) {
      cerr << "ERR:sequenceSpace::getAbProducingIndex(" << n
           << "), index out of bounds [0.." << ab_producers.size() - 1 << "]" << endl;
      return 0;
   }
   return ab_producers[n];
}
template <typename T>
int generalSequenceContainer<T>::get_n_ab_producers() {
   return ab_producers.size();
}
template <typename T>
void generalSequenceContainer<T>::set_external_cell(long int pos) {
   T * seq = getSequence(pos);
   if (!seq->hasExternalCell()) {
      // it is made such that, for sequences or ArupProteins, if this is called from a non-BCR, an
      // error will be launched
      external_cells.push_back(pos);    // the first time that a sequence gets producers (meaning
                                        // n_cells[soutext] gets > 0),
   }
}
// thie function is a bit of cheating because it only works for a=soutext and b=soutextproduce
template <typename T>
void generalSequenceContainer<T>::PM_differentiate(cells a, cells b, double dt) {
   int n_external = external_cells.size();
   for (int i = 0; i < n_external; i++) {
      long n = external_cells[i];
      // if (poolSequences[n]->getType() == typeBCR) { // is sequence type dependent. external cells
      // should only contain BCR ids so it's ok
      double transit = poolSequences[n]->n_cell[a] * pm_differentiation_rate * dt;
      if ((transit > 0) && (poolSequences[n]->n_cell[b] == 0)) {
         // the first time there is a producer at this position
         ab_producers.push_back(n);
      }
      poolSequences[n]->n_cell[a] -= transit;
      sum_cell[a] -= transit;
      poolSequences[n]->n_cell[b] += transit;
      sum_cell[b] += transit;
      // }
   }
}
template <typename T>
double generalSequenceContainer<T>::get_sum_cell(cells celltype) { return sum_cell[celltype]; }
template <typename T>
double generalSequenceContainer<T>::get_oldsum_cell(cells celltype) {
   return oldsum_cell[celltype];
}
template <typename T>
short int generalSequenceContainer<T>::sum_check() {
   long int i;
   double dtmp;
   short int err = 0;
   cout << "Sequence/Arup space sum check ... ";
   for (int n = 0; n < number; n++) {
      // Durchlaufe Typen
      dtmp = 0.;
      for (i = 0; i < n_Sequences; i++) {
         // Durchlaufe shape space Punkte
         if (getSequence(i)->n_cell[n] < 0) {
            cout << "ERR: NbCells in sequence[ID=" << i << "].[typeCell=" << n << "] = "
                 << getSequence(i)->n_cell[n] << " <0\n";
            exit(1);
         }
         // if (n!=soutext && n!= soutextproduce)
         dtmp += getSequence(i)->n_cell[n];
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
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         2 - SEQUENCES IN SAHAMS AFFINITY MODEL
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/* ---- Constructors ---- */

int sequence::operator [](int i) {
   if ((i < 0) || (i > size)) {
      cerr << "sequence[" << i << "] out of bounds (size=" << size << ")\n";
      return 0;
   }
   return (content[i]) ? 1 : 0;
}
sequence::sequence(string boolText) {
   size = boolText.size();
   content.resize(size, false);
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   for (int i = 0; i < size; ++i) {
      // if((boolText[i] != '1') && (boolText[i] != '0')) {cerr << "ERR:sequence(string), parsing
      // sequence with un-allowed character : " << boolText[i] << endl;}
      content[i] = (boolText[i] == '1');
   }
   max_affinity_to_antigens = -1;
}
sequence::sequence(sequence * toCopy) {
   size = toCopy->size;
   content.resize(size);
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   for (int i = 0; i < size; ++i) {
      content[i] = toCopy->content[i];
   }
   max_affinity_to_antigens = toCopy->max_affinity_to_antigens;     // ??
}
sequence::sequence(int _size) {
   if (_size < 0) { _size = 0; }
   size = _size;
   content.resize(_size, false);
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   max_affinity_to_antigens = -1;
}
// A function to transform a long number into binary sequence, inside a given vector.
void recursive_binary(long ID, vector<bool> &towrite, int currentposition) {
   // to be called with : recursive_binary(2354, aboolvectortofill, sizeofthisvector - 1);
   if (ID < 0) { cerr << "ERR::sequence::binary, negative number " << ID << endl; return; }
   if (currentposition
       < 0) {
      cerr
         <<
         "Err: sequence::binary, the number was too big, exceeds the size of the output vector."
         << endl;
      return;
   }
   if (currentposition
       >= (int) towrite.size()) {
      cerr << "Err: sequence::binary, currentposition is out of bounds ("
           << currentposition << "), while vector size is " << towrite.size() << endl;
      return;
   }

   if (ID <= 1) { towrite[currentposition] = (ID != 0); return; }  // to stop recursion
   long remainder = ID % 2;
   towrite[currentposition] = (remainder != 0);
   recursive_binary(ID >> 1, towrite, currentposition - 1);
}
sequence::sequence(int _size, long ID) {
   /*: sequence(_size)*/
   // since doing : sequence(_size) leads to a warning because only works in C++11, I copy the
   // sequence constructor here :
   if (_size < 0) { _size = 0; }
   size = _size;
   content.resize(_size, false);
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   max_affinity_to_antigens = -1;

   // now, converts the ID into binary
   recursive_binary(ID, this->content, this->content.size() - 1);

   // Alternative way using bitset, but don't use that anymore because the testSize has to be known
   // at compilation
   /*string seq = std::bitset<testSize>(v).to_string();
    * for(int i = 0; i < testSize; ++i){
    *  a->content[i] = (seq[i] == '1');
    * }
    * return a;*/
}
/* ---- Affinity, hamming and comparison (static -> can be called outside any instance, by
 * 'sequence::' ---- */

double sequence::seq_affinity(sequence * x,
                              sequence * y,
                              double r,
                              int maxSizeClusters,
                              int type_affinity_function) {
   switch (type_affinity_function) {
      case seqAff: {
         // The original seq_affinity function

         // double sequence::seq_affinity(sequence* x, sequence* y, double r){
         //                        // This is slightly different than sahams model here by saying
         // there
         // is affinity if sequence are EQUAL instead of COMPLEMENTARY, in order to have BCR-AG-TCRs
         // are close if they have the same sequence (easier to interprete)
         if (x->size
             != y->size) {
            std::cout << "Err affinity: length of antigen and BCR are not identical!";
            return -1;
         }
         int l = 0;
         double sum = 0;
         for (int i = 0; i < x->size; ++i) {
            if ((*x)[i] == (*y)[i]) {
               // as long as the sequences are equal, increases the matching size l
               ++l;
            } else {
               sum += std::pow((double) l,r);      // when they are not equal anymore, close the
                                                   // match, and calculates its size power r
               l = 0;
            }
         }
         if (l > 1e-6) {
            sum += std::pow(l,r);          // if the match was not closed, close it ; l = 0 if
                                           // already closed
         }
         return sum / std::pow((double) x->size,r);
         // }           // Other Idea : implement our own my_pow that precomputes all the values of
         // pow(int i[0..length], r) to save time
      }

      case seqAffNorm: {
         // Same but instead of normalizing by L^r, normalizes to maxSizeClusters^r, and brings back
         // to
         // 1.0 if exceeds

         if (x->size
             != y->size) {
            std::cout << "Err affinity: length of antigen and BCR are not identical!";
            return -1;
         }
         int l = 0;
         double sum = 0;
         for (int i = 0; i < x->size; ++i) {
            if ((*x)[i] == (*y)[i]) {
               // as long as the sequences are equal, increases the matching size l
               ++l;
            } else {
               sum += std::pow((double) l,r);      // when they are not equal anymore, close the
                                                   // match, and calculates its size power r
               l = 0;
            }
         }
         if (l > 1e-6) {
            sum += std::pow(l,r);          // if the match was not closed, close it ; l = 0 if
                                           // already closed
         }
         return min(1.0, (sum / std::pow((double) maxSizeClusters,r)));
      }

      case seqAffWindow: {
         int sizewindow = maxSizeClusters;

         // the new affinity function with sizewindow. Call with sizewindow = -1 to get the original
         // affinity function.
         // to keep the algorithm linear and not compute an affinity for each window separately,
         // computes the
         // affinity of the sliding window in real time by adding the contribution of the right new
         // position taken
         // and removing the contribution of the left position which is left by the sliding windows.
         // double sequence::seq_affinity(sequence* x, sequence* y, double r, int sizewindow){
         if (x->size
             != y->size) {
            std::cout << "Err affinity: length of antigen and BCR are not identical!";
            return -1;
         }
         if (sizewindow < 1) { sizewindow = x->size; }

         // Step 1 : Counting inside the matches, from the left and from the right. Example, if the
         // match is 011111000111,
         vector<int> left(x->size, -1);     // -1 1 2 3 4 5 6 -1 -1 -1 1 2 3       -> remaining of
                                            // the cluster at this position on the right
         vector<int> right(x->size, -1);    // -1 6 5 4 3 2 1 -1 -1 -1 3 2 1       -> size of the
                                            // cluster if taking one more position on the right
         int l = 0;
         for (int i = 0; i < x->size + 1; ++i) {
            // the +1 allows to close the last match if it touches the end of the sequence
            if ((i < x->size) && ((*x)[i] == (*y)[i])) {
               // as long as the sequences are equal, increases the matching size l. For the last
               // check
               // (i = x->size, directly go to the 'else' to close the last match)
               ++l;
               left[i] = l;                         // 1 2 3 4 5 ... inside the match/cluster
            } else {
               for (int j = 1; j <= l; ++j) {
                  // goes backwards from the end of the cluster (i-1) to its beginning (i-l)
                  right[i - j] = j;                 // ... 5 4 3 2 1  inside the match/cluster
               }
               l = 0;
            }
         }

         // Step2 : moving a sliding window and getting its affinity (windowsum) each time.
         double windowsum = 0;          // affinity of the sliding window
         double maxWindowSum = 0;       // maximum encountered
         for (int i = 0; i < x->size; ++i) {
            // the sliding window now 'eat' the position i and 'leaves' the position i - sizewindow
            int posToRemove = i - sizewindow;

            // Eating : when the sliding windows take a new position of a match/cluster
            if (left[i] > 0) {
               // if new position taken by the window is in a match
               windowsum += std::pow(left[i], r);                                       // then add
                                                                                        // the macth
               if ((i > 0) && (left[i - 1] > 0)) {
                  windowsum -= std::pow(left[i - 1], r);                                // but if
                                                                                        // the last
                                                                                        // position
                                                                                        // was
                                                                                        // already
                                                                                        // in the
                                                                                        // macth,
                                                                                        // then
                                                                                        // remove
                                                                                        // match-1
               }
            }

            // Pooing : when the sliding windows leaves a position that was in a match
            if ((posToRemove >= 0)          // if the sliding window started leaving positions
                && (right[posToRemove] > 0)) {
               // if the position that the window is leaving was in a match,
               windowsum -= std::pow(right[posToRemove], r);                                    // then
                                                                                                // remove
                                                                                                // that
                                                                                                // match
               if (right[posToRemove + 1] > 0) {
                  windowsum += std::pow(right[posToRemove + 1], r);                             // but
                                                                                                // in
                                                                                                // case
                                                                                                // it
                                                                                                // stays
                                                                                                // in
                                                                                                // the
                                                                                                // match,
                                                                                                // restore
                                                                                                // the
                                                                                                // remaining:
                                                                                                // match-1
               }
            }

            maxWindowSum = max(maxWindowSum, windowsum);        // Note : do not put the condition
                                                                // if(i >= sizewindow-1)  in case
                                                                // the sequence is smaller than the
                                                                // window. Anyway, as long as i <
                                                                // sizewindow, the sum only
                                                                // increases so the maximum will
                                                                // always increase and it doesn't
                                                                // make a difference to put this
                                                                // condition or not.

            // cout << ((i < sizewindow-1)? 0 : windowsum) << endl;
            // cout << "left [" << i << "] = " << left[i] << "\t";
            // cout << "right[" << i << "] = " << right[i] << "\t";
         }
         return maxWindowSum / std::pow((double) sizewindow, r);
      }
   }  // end switch
   cerr << "sequence::seq_affinity(), unknown type of affinity function "
        << type_affinity_function << endl;
   return 0.0;
}
double sequence::hamming(sequence * x, sequence * y) {
   if (x->size
       != y->size) {
      std::cout << "Err affinity: length of antigen and BCR are not identical!";
      return -1;
   }
   int length = x->size;
   int sum = 0;
   for (int i = 0; i < length; ++i) {
      sum += ((*x)[i] != (*y)[i]) ? 1 : 0;
   }
   return sum;
}
bool sequence::compSeq(sequence * a, sequence * b) {
   bool res = (a->size == b->size);
   if (res) {
      for (int i = 0; i < a->size; ++i) {
         res = res && ((*a)[i] == (*b)[i]);
      }
   }
   return res;
}
string sequence::print() {
   stringstream out;
   for (int i = 0; i < size; ++i) {
      out << content[i];
   }
   return out.str();
}
/* ---- Asking a sequence to mutate ---- */

void sequence::mutateOnePosition() {
   int position = irandom(size);
   if ((position >= size)
       || (position < 0)) {
      cerr << "Random does shit, sequencespace.cpp::mutateOnePosition";
      return;
   }
   content[position] = !content[position];
}
void sequence::randomize() {
   for (int i = 0; i < size; ++i) {
      if (drandom() < 0.5) {
         content[i] = !content[i];
      }
   }
}
string sequence::typeInString() {
   switch (getType()) {
      case typeAG:
         return string("Antigene");

      case typeTCR:
         return string("TCR");

      case typeBCR:
         return string("BCR/Antibody");

      case numberTypesSeq:
         return string("Unknown");
   }
   return string("Incorrect type");
}
string sequence::printNbCells() {
   stringstream res1, res2;
   res1 << "type of cells :";
   res2 << "Number";
   for (int i = 0; i < number; ++i) {
      res1 << "\t" << typeCellInString((cells) i);
      res2 << "\t\t" << n_cell[i];
   }
   return res1.str() + string("\n") + res2.str();
}
/* @brief For each position, mutate it with the probability probPerBaseCorrected
 *  and returns true if the mutation occured.
 *  Note that the probability probPerBaseCorrected is the 'observed mutation rate' per position
 *  and has to be corrected from the mechanistic mutation rate in order to deduce the probability
 *  of a position mutating twice (i.e. not seen), three times (seen), four times (not seen) etc...
 *  therefore probPerBaseCorrected = sequence::correctMut(mechanisticProbaPerBase).
 *  bool mutate(double probPerBaseCorrected);
 *
 *  @brief see mutate(proba)
 *  static double correctMut(double p);
 *  // in case needed, transforms a proba of mutation per base to a proba of being really mutated,
 *per
 *  // base (removing the cases or mutation and reverse mutation happening at the same place)
 *  double sequence::correctMut(double p) {
 *  return p - pow(p,2) + pow(p,3) - pow(p,4) + pow(p,5);
 *  }
 *  bool sequence::mutate(double probPerBaseCorrected) {
 *  // int events = binomial(size, probPerBaseCorrected);
 *  bool isMutated = false;
 *  for (int i = 0; i < size; ++i) {
 *   if (drandom() < probPerBaseCorrected) {
 *     content[i] = !content[i];
 *     isMutated = true;
 *   }
 *  }
 *  return isMutated;
 *  }
 */

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         3 - SEQUENCES IN ARUPS AFFINITY MODEL
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double probaLawFromTable::getRandValue() {
   int classTaken = getRandClass();
   if (uniformizeInsideClasses) {
      return lowBoundsXs[classTaken] + IndependentRandRealGen()
             * (highBoundsXs[classTaken] - lowBoundsXs[classTaken]);
   } else {
      return Xs[classTaken];
   }
}
int probaLawFromTable::getRandClass() {
   if (size == 0) { cerr << "Err : probaLawFromTable::getRandClass(), NO DATA" << endl; }
   double rdreal = IndependentRandRealGen();
   for (int i = 0; i < size; ++i) {
      if (rdreal > cumulatedValuesLow[i]) { return i; }
   }
   return size - 1;
}
probaLawFromTable::probaLawFromTable(vector<double> &_Xs,
                                     vector<double> &densities,
                                     bool _uniformizeInsideClasses) {
   if (_Xs.size()
       != densities.size()) {
      cerr
         << "ERR: probaLawFromTable, Xs and Densities vectors don't have the same size." << endl;
      return;
   }
   size = _Xs.size();
   Xs.resize(size, 0.0);
   NormalizedDensities.resize(size, 0.0);
   cumulatedValuesLow.resize(size, 0.0);
   if (size > 0) {
      double sumDensities = 0.0;
      for (int i = 0; i < size; ++i) {
         Xs[i] = _Xs[i];
         if ((i > 0)
             && (Xs[i]
                 < Xs[i - 1])) {
            cerr << "ERR : probaLawFromTable, the list of X values should be increasing only."
                 << endl;
         }
         sumDensities += densities[i];
      }
      normalizingCoeff = sumDensities;
      for (int i = 0; i < size; ++i) {
         NormalizedDensities[i] = densities[i] / sumDensities;
         if (i
             > 0) { cumulatedValuesLow[i] = cumulatedValuesLow[i - 1] + NormalizedDensities[i - 1];
         }
      }
      if (fabs((NormalizedDensities[size - 1] + cumulatedValuesLow[size - 1]) - 1.0)
          > 1e-6) {
         cerr
            << "ERR : probaLawFromTable, sum of probabilities is not 1. Should not happen. Why ??"
            << endl;
      }

      // In case one wants all the possible values inside each class,
      if (uniformizeInsideClasses) {
         lowBoundsXs.resize(size);
         highBoundsXs.resize(size);
         for (int i = 0; i < size; ++i) {
            if (i > 0) { lowBoundsXs[i] = (Xs[i] + Xs[i - 1]) / 2.0; }
            if (i < size - 1) { highBoundsXs[i] = (Xs[i + 1] + Xs[i]) / 2.0; }
         }
         lowBoundsXs[0] = Xs[0] - (highBoundsXs[0] - Xs[0]);
         highBoundsXs[size - 1] = Xs[size - 1] + (Xs[size - 1] - lowBoundsXs[size - 1]);
      }
   }
}
double probaLawFromTable::IndependentRandRealGen() {
   // If later someone wants to separate generators for seeding.
   /* needs C++11 and #include <random>
    *  static std::random_device *rd                       = new std::random_device();;
    *  static std::mt19937 *gen                            = new std::mt19937 ((*rd)());
    *  static std::uniform_real_distribution<> *RealDistrib= new std::uniform_real_distribution<>
    * (0,1);
    *  return (*RealDistrib)(*gen);  */
   return drandom();
}
string probaLawFromTable::print() {
   if (size == 0) { return string(); }
   stringstream res;
   res << "probaLawFromTable with following classes : " << endl;
   if (uniformizeInsideClasses) {
      res
         << " --> Random values are taken uniformly inside each class\n";
   } else { res << " --> Random values are always the exact value of each class\n"; }
   res << "i\tX[i]\tproba\tcumulLow\toriginalDensity"
       << ((uniformizeInsideClasses) ? "\tLowBound\tHighBound" : "") << endl;
   for (int i = 0; i < size; ++i) {
      res << i << "\t" << Xs[i] << "\t"
          << ((i
           < size - 1) ? cumulatedValuesLow[i + 1] - cumulatedValuesLow[i] : 1
          - cumulatedValuesLow[size - 1]) << "\t" << cumulatedValuesLow[i] << "\t"
          << (normalizingCoeff
          * ((i < size - 1) ? cumulatedValuesLow[i + 1] : 1 - cumulatedValuesLow[size - 1]));
      if (uniformizeInsideClasses) { res << "\t" << lowBoundsXs[i] << "\t" << highBoundsXs[i]; }
      res << endl;
   }
   res << "... The initial density was of total area " << normalizingCoeff << endl;
   return res.str();
}
probaLawFromTableStoringResults::probaLawFromTableStoringResults(vector<double> &_Xs,
                                                                 vector<double> &densities,
                                                                 bool _uniformizeInsideClasses) :
   probaLawFromTable(_Xs, densities,  _uniformizeInsideClasses) {
   nbOfTimesPerClass.resize(size, 0);
   totalNbEvents = 0;
}
double probaLawFromTableStoringResults::getRandValue() {
   double valueTaken = probaLawFromTable::getRandValue();  // the function from the mother class
   int correspondingClass = fitDoubleInAClass(valueTaken);
   if ((correspondingClass < 0)
       || (correspondingClass
           >= size)) {
      cerr
           <<
      "Err :  probaLawFromTableStoringResults, the getRandValue from probaLawFromTable was giving a value which is out of bounds : "
           << valueTaken << endl;
   } else {
      nbOfTimesPerClass[correspondingClass]++;
      totalNbEvents++;
   }
   return valueTaken;
}
double probaLawFromTableStoringResults::getFrequency(int classIndex) {
   if ((classIndex < 0) || (classIndex >= size)) { return -1; }
   return ((double) nbOfTimesPerClass[classIndex]) / (totalNbEvents);
}
int probaLawFromTableStoringResults::fitDoubleInAClass(double val) {
   if (uniformizeInsideClasses) {
      for (int i = 0; i < size; ++i) {
         if ((val > lowBoundsXs[i]) && (val < highBoundsXs[i])) {
            return i;
         }
      }
   } else {
      for (int i = 0; i < size; ++i) {
         if (fabs(val - Xs[i]) < 1e-6) {
            return i;
         }
      }
   }
   return -1;
}
void probaLawFromTableStoringResults::clearRecord() {
   nbOfTimesPerClass.clear();
   nbOfTimesPerClass.resize(size, 0);
   totalNbEvents = 0;
}
string probaLawFromTableStoringResults::print() {
   if (size == 0) { return string(); }
   stringstream res;
   res << "probaLawFromTableStoringResults, extending ";
   res << probaLawFromTable::print() << endl;
   res << totalNbEvents << " events were recorded\n";
   res << "i\tXs[i]\tFrequency\tRequestedDensity\n";
   for (int i = 0; i < size; ++i) {
      res << i << "\t" << Xs[i] << "\t" << getFrequency(i) << "\t"
          << ((i
           < size - 1) ? cumulatedValuesLow[i + 1] : 1 - cumulatedValuesLow[size - 1]) << endl;
   }
   return res.str();
}
string probaLawFromTableStoringResults::TestprobaLawFromTable() {
   stringstream res;
   vector<double> Xs;   // in case people doesn't use C++11, = {bla, bla, bla ...} doesn't work
   Xs.push_back(0.5);
   Xs.push_back(1);
   Xs.push_back(1.5);
   Xs.push_back(2);
   vector<double> densities;
   densities.push_back(0.4);
   densities.push_back(0.8);
   densities.push_back(1.2);
   densities.push_back(0.4);
   probaLawFromTable Law = probaLawFromTable(Xs, densities, false);
   cout << Law.print() << endl;
   for (int i = 0; i < 50; ++i) {
      res << Law.getRandValue() << endl;
   }

   probaLawFromTableStoringResults Law2 = probaLawFromTableStoringResults(Xs, densities, true);
   cout << Law2.print() << endl;
   for (int i = 0; i < 10000; ++i) {
      Law2.getRandValue();
   }
   cout << Law2.print() << endl;
   return res.str();
}
/* ================================= Manipulation of arupProteins
 * ====================================== */

/* ---- Constructors ---- */

int arupProtein::arup_length_sequences = 46;
int arupProtein::arup_N_conserved = 18;
int arupProtein::arup_N_mutates = 22;
int arupProtein::arup_N_shielded = 6;
double arupProtein::arup_alpha = 2.0;
double arupProtein::arup_h_min = -0.18;
double arupProtein::arup_h_max = 0.9;
double arupProtein::arup_hprime_min = -1.5;
double arupProtein::arup_hprime_max = 1.5;
double arupProtein::arup_hmut_min = -1.5;
double arupProtein::arup_hmut_max = 1e6;
probaLawFromTable* arupProtein::lawMutations = NULL;

void arupProtein::set_statics(Parameter &par, ofstream &ana) {
   arup_length_sequences = par.Value.arup_length_sequences;
   arup_N_conserved = par.Value.arup_N_conserved;
   arup_N_mutates = par.Value.arup_N_mutates;
   arup_N_shielded = par.Value.arup_N_shielded;
   arup_alpha = par.Value.arup_alpha;
   arup_hprime_min = par.Value.arup_hprime_min;
   arup_hprime_max = par.Value.arup_hprime_max;
   arup_hmut_min = par.Value.arup_hmut_min;
   arup_hmut_max = par.Value.arup_hmut_max;
   arup_h_min = par.Value.arup_h_min;
   arup_h_max = par.Value.arup_h_max;
   lawMutations = new probaLawFromTable(par.Value.arup_law_mut_Xs,
                                        par.Value.arup_law_mut_Densities,
                                        true);
   ana << "Initializing parameters for Arup Sequences : \n";
   ana << "   Number of residues : " << arup_length_sequences << " among which \n";
   ana << "    - " << arup_N_conserved << " conserved\n";
   ana << "    - " << arup_N_mutates << " mutables\n";
   ana << "    - " << arup_N_shielded << " shielded\n";
   ana << "   Unshielding coefficient (alpha) = " << arup_alpha << "\n";
   ana << "   Initial boundaries for each position in a sequence : [" << arup_h_min << " ; "
       << arup_h_max << "]\n";
   ana << "   Boundaries for mutable positions during mutation   : [" << arup_hmut_min << " ; "
       << arup_hmut_max << "]\n";
   ana << "   Boundaries for shielded positions during mutation  : [" << arup_hprime_min << " ; "
       << arup_hprime_max << "]\n";
   ana << "   Law of probabilities to generate mutations : " << lawMutations->print() << "\n\n";
}
// I DO NOT UNDERSTAND THE FUCKING TARGET HERE
// check if irandom reaches the upper boundary
/* ---- Asking a arupProtein to mutate ---- */

void ArupBCR::mutate(double targetChangeEnergy) {
   int posToMutate = irandom(arup_length_sequences - 1);
   double deltaE = lawMutations->getRandValue();
   double deltaHK = deltaE;   // because target is wild-type ???
   content[posToMutate] += deltaHK;
   if (posToMutate >= arup_N_conserved + arup_N_mutates) {
      // if in the shielded area,
      int posCompensation = irandom(arup_N_conserved - 1);
      double hprime = content[posCompensation] - arup_alpha * (deltaHK);
      hprime = max(arup_hprime_min, hprime);
      hprime = min(arup_hprime_max, hprime);
      content[posCompensation] = hprime;
   }
   if ((posToMutate >= arup_N_conserved) && (posToMutate < (arup_N_conserved + arup_N_mutates))) {
      // if mutates in the variable area
      content[posToMutate] = max(arup_hmut_min, content[posToMutate]);
      content[posToMutate] = min(arup_hmut_max, content[posToMutate]);
   }
}
void ArupBCR::randomize(double targetEnergy, bool beExact) {
   arupProtein WildType = arupProtein();  // by default, only 1s (see empty constructor)
   // strategy 1 : generate randomly up to a good enough energy. In case it doesn't work, keeps the
   // best;
   arupProtein bestSoFar = arupProtein();
   double bestAffinity = -1e6;
   int count = 0;
   while ((count < 10000) && (bestAffinity < targetEnergy)) {
      for (int i = 0; i < size; ++i) {
         content[i] = arup_h_min + drandom() * (arup_h_max - arup_h_min);
      }
      double newAffinity = (arup_affinity(this, &WildType));
      if (newAffinity > bestAffinity) {
         bestAffinity = newAffinity;
         bestSoFar.copy(this);
      }
      ++count;
   }
   if ((!beExact) && (bestAffinity > targetEnergy)) {
      arupProtein::copy(&bestSoFar);   // this is a bit weird because you are a BCR , but you just
                                       // copy the ArupProtein part ...
      return;
   }
   // strategy 2 : Now mutates to reach the target (being exact and/or being high enough);
   double targetChange = targetEnergy - bestAffinity;
   int nbTries = 0;
   while ((nbTries < 10)
          && ((beExact && (fabs(targetChange) > 1e-4)) || ((!beExact) && (targetChange > 0)))) {
      mutate(targetChange);
      targetChange = targetEnergy - arup_affinity(this, &WildType);
      nbTries++;
   }
}
void ArupAntigen::randomize(int nbResiduesToMutate) {
   if (nbResiduesToMutate < 0) {
      for (int i = arup_N_conserved; i < arup_N_conserved + arup_N_mutates; ++i) {
         if (drandom() < 0.5) { content[i] = -content[i]; }
      }
   } else {
      vector<bool> mutatedPos = vector<bool> (arup_length_sequences, false);
      int nbTries = 0;
      int nbMuts = 0;
      while ((nbTries < 100) && (nbMuts < nbResiduesToMutate)) {
         int pos = arup_N_conserved + irandom(arup_N_mutates - 1);
         if (!mutatedPos[pos]) {
            content[pos] = -content[pos];
            mutatedPos[pos] = true;
            nbMuts++;
         }
      }
   }
   // a stupid check, in case mutations occurs outside the area ...
   for (int i = 0; i < arup_N_conserved; ++i) {
      if (fabs(content[i] - 1.0)
          > 1e-9) {
         cerr
            <<
            "ERR : ArupAntigen::randomize, this antigen has problems with conserved residues that are not 1"
            << endl;
      }
   }
   for (int i = arup_N_conserved + arup_N_mutates; i < arup_length_sequences; ++i) {
      if (fabs(content[i] - 1.0)
          > 1e-9) {
         cerr
            <<
            "ERR : ArupAntigen::randomize, this antigen has problems with shielded residues that are not 1"
            << endl;
      }
   }
}
void ArupAntigen::mutate() {
   int pos = arup_N_conserved + irandom(arup_N_mutates - 1);
   content[pos] = -content[pos];
}
void arupProtein::clear() {
   size = arup_length_sequences;
   content.resize(size, 1.0);   // this 1 is very important because the wild-type antigen is 111111
                                // by default.
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   max_affinity_to_antigens = NAN;
}
int arupProtein::operator [](int i) {
   if ((i < 0)
       || (i > size)) {
      cerr << "arupProtein[" << i << "] out of bounds (size=" << size << ")\n";
      return 0;
   }
   return content[i];
}
arupProtein::arupProtein(vector<double> numbers) {
   clear();
   if ((int) numbers.size()
       != arup_length_sequences) {
      cerr
         << "ERR: arupProtein::arupProtein(vector), sequences should be of size "
         << arup_length_sequences << " while the input vector has " << numbers.size() << endl;
      return;
   }
   for (int i = 0; i < size; ++i) {
      content[i] = numbers[i];
   }
}
arupProtein::arupProtein(string numbersInText) {
   clear();
   stringstream transform(numbersInText);
   double buffer;
   content.clear();
   while ((transform >> buffer) && (transform.str().size() > 0)) {
      content.push_back(buffer);
   }
   if ((int) content.size() != arup_length_sequences) {
      cerr << "ERR: arupProtein::arupProtein(vector), sequences should be of size "
           << arup_length_sequences << " while the input string of numbers has " << content.size()
           << endl;
      content.resize(arup_length_sequences);
      return;
   }
}
arupProtein::arupProtein(arupProtein * toCopy) {
   size = toCopy->size;
   content.resize(size);
   for (int i = 0; i < number; ++i) {
      n_cell[i] = 0;
   }
   for (int i = 0; i < size; ++i) {
      content[i] = toCopy->content[i];
   }
   max_affinity_to_antigens = toCopy->max_affinity_to_antigens;     // ??
}
void arupProtein::copy(arupProtein * toCopy) {
   if (toCopy->size
       != size) { cerr << "ERR:  arupProtein::copy, copying sequences of different size" << endl; }
   for (int i = 0; i < size; ++i) {
      content[i] = toCopy->content[i];
   }
   max_affinity_to_antigens = toCopy->max_affinity_to_antigens;
}
arupProtein::arupProtein() {
   clear();
}
arupProtein::~arupProtein() {
   clear();
   cerr << "~arupProtein() Should be reimplemented by daughter classes" << endl;
}
double arupProtein::arup_affinity(arupProtein * BCR, arupProtein * Antigen) {
   if (BCR->size
       != Antigen->size) {
      cerr
         << "ERR: arupProtein::arup_affinity, comparing sequences of different sizes" << endl;
   }
   double res = 0;
   for (int i = 0; i < BCR->size; ++i) {
      res += (*BCR)[i] * (*Antigen)[i];
      if ((i < arup_N_conserved)
          && (fabs(((*Antigen)[i]) - 1.0)
              > 1e-12)) {
         cerr
              <<
         "ERR: arupProtein::arup_affinity, conserved positions in the antigen should be 1. Antigen = "
              << Antigen->print() << endl;
      }
      if ((i >= (arup_N_conserved + arup_N_mutates))
          && (fabs(((*Antigen)[i]) - 1.0)
              > 1e-12)) {
         cerr
            <<
            "ERR: arupProtein::arup_affinity, shielded positions in the antigen should be 1. Antigen = "
            << Antigen->print() << endl;
      }
   }
   return res;
}
double arupProtein::continuoushamming(arupProtein * x, arupProtein * y) {
   if (x->size
       != y->size) {
      std::cout << "Err affinity: length of antigen and BCR are not identical!";
      return -1;
   }
   int length = x->size;
   int sum = 0;
   for (int i = 0; i < length; ++i) {
      sum += fabs((*x)[i] - (*y)[i]);
   }
   return sum;
}
bool arupProtein::compSeq(arupProtein * a, arupProtein * b) {
   bool res = (a->size == b->size);
   if (res) {
      for (int i = 0; i < a->size; ++i) {
         res = res && ((*a)[i] == (*b)[i]);
      }
   }
   return res;
}
string arupProtein::print() {
   stringstream out;
   for (int i = 0; i < size; ++i) {
      out << content[i] << " ";
   }
   return out.str();
}
string arupProtein::printNbCells() {
   stringstream res1, res2;
   res1 << "type of cells :";
   res2 << "Number";
   for (int i = 0; i < number; ++i) {
      res1 << "\t" << typeCellInString((cells) i);
      res2 << "\t\t" << n_cell[i];
   }
   return res1.str() + string("\n") + res2.str();
}
typeArupSeq arupProtein::getType() {
   // cerr << "ERR: arupProtein::getType(), you are calling getType from a arupProtein without type"
   // << endl;
   return numberTypesArup;
}
string arupProtein::typeInString() {
   switch (getType()) {
      case typeArupAG:
         return string("ArupAntigene");

      case typeArupBCR:
         return string("ArupBCR/Antibody");

      case numberTypesArup:
         return string("ArupUnknown");
   }
   return string("Incorrect type");
}
// function to print the types of cells as string (used for printing in the next function)
string arupProtein::typeCellInString(cells index) { return sequence::typeCellInString(index); }

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         4 - (BOOLEAN) SEQUENCE SPACE
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

sequenceSpace::sequenceSpace(Parameter &par, ofstream &ana) :
   generalSequenceContainer<sequence> (),
   OUT_haffinity(0),
   OUT_steepness(0),
   CB_haffinity(0),
   CC_haffinity(0) {
   long int n;
   ana << "Initialize Sequence Space fields ... \n";
   typeAffinityFunction = par.Value.type_affinity_function;
   use_logarithmic_seq_affinity = par.Value.use_logarithmic_seq_affinity;
   R_affinity = par.Value.R_affinity;
   maxSizeClusters = par.Value.max_affinity_cluster;
   if (maxSizeClusters < 0) { maxSizeClusters = par.Value.size_sequences; }
   size_sequences = par.Value.size_sequences;
   pm_differentiation_rate = log(2.) / par.Value.pm_differentiation_time;
   amplitude = par.Value.amplitudeGauss;
   ana << "Size of sequences : " << par.Value.size_sequences << endl;
   ana << "Specificity R parameter for sequence affinity : R= " << R_affinity << endl;
   ana << "Maximum size of clusters for affinity computing :  " << maxSizeClusters << endl;
   ana << "Type of affinity function "
       << ((typeAffinityFunction == seqAff) ? "standard affinity/seqAff(=0)" : "")
       << ((typeAffinityFunction
        == seqAffNorm) ? "affinity normalized by cluster size/seqAffNorm(=1)" : "")
       << ((typeAffinityFunction == seqAffWindow) ? "sliding window/seqAffWindow(=2)" : "") 
       << endl;
   ana
      << ((use_logarithmic_seq_affinity) ? "Using Logarithmic of the affinity" :
          "(No special rescaling)")
      << endl;
   ana << "  pm differentiation rate = " << pm_differentiation_rate << "\n";
   ana << "  amplitude (coefficient to all affinities) : " << amplitude << "\n";

   ana << "Initialize Antigen pool ...\n ";
   ana << "Starting antigens required : " << par.Value.init_antigen_sequences << endl;
   ana << "  with hamming distance between each-other within " << par.Value.min_hamming_antigens
       << " and " << par.Value.max_hamming_antigens << endl;
   for (n = 0; n < par.Value.init_antigen_sequences; n++) {
      if ((n < (int) par.Value.initAntigenSeqs.size())
          && (par.Value.initAntigenSeqs[n].size() > 2)) {
         // because this vector also gets the '-1' as ending
         int index = add_Antigen(par.Value.initAntigenSeqs[n]);

         ana << "   ... predefined seq nr" << n << " index= " << index << "(" << getSequence(
            index)->size << "):" << getSequence(index)->print() << endl;
         if (getSequence(index)->size != par.Value.size_sequences) {
            cerr
                 <<
               "ERR: SequenceSpace constructor, the sequence given as initial antigens doesn't fit to the lenth of sequences which is "
                 << par.Value.size_sequences << endl;
         }
      } else {
         Antigen * newAG = new Antigen(par.Value.size_sequences);
         int tries = 0;
         bool correct = false;
         while ((tries < 10000) && (!correct)) {
            newAG->randomize();
            correct = true;
            // correct = ((2 * hamming( Antigen[n], seedSeq)) <= par.Value.max_hamming_antigens);
            for (int i = 1; i < n_Antigen; ++i) {
               correct = correct
                         && ((sequence::hamming(getSequence(get_Antigen(i)), newAG)
                              >= par.Value.min_hamming_antigens)
                             && (sequence::hamming(getSequence(get_Antigen(i)), newAG)
                                 <= par.Value.max_hamming_antigens));
            }
            tries++;
         }
         // Now, if random sequences didn't manage to fill the constraints, build sequence manually
         // to fill them.
         if (!correct) {
            cerr
                 <<
               "Warning : SequenceSpace, randomly generated antigens didn't fit the expected hamming distance ["
                 << par.Value.min_hamming_antigens << ";" << par.Value.max_hamming_antigens
                 << "]" << endl;
            cerr << " ... now trying to build such antigens manually" << endl;
            while ((tries < 10000) && (!correct)) {
               delete newAG;
               newAG = new Antigen(getSequence(get_Antigen())); // takes a random antigen from the
                                                                // pool
               int nbAdditionalMut = 0;
               for (int j = 0; j < par.Value.min_hamming_antigens - 1; ++j) {
                  newAG->mutateOnePosition();
                  nbAdditionalMut++;
               }
               while ((!correct) && (nbAdditionalMut <= par.Value.max_hamming_antigens)) {
                  newAG->mutateOnePosition();
                  nbAdditionalMut++;
                  for (int i = 1; i < n_Antigen; ++i) {
                     correct = correct
                               && ((sequence::hamming(getSequence(get_Antigen(i)), newAG)
                                    >= par.Value.min_hamming_antigens)
                                   && (sequence::hamming(getSequence(get_Antigen(i)), newAG)
                                       <= par.Value.max_hamming_antigens));
                  }
               }
            }
            if (!correct) {
               cerr
                    <<
                  "ERROR !! : SequenceSpace, failed to build antigens within the expected hamming distance ["
                    << par.Value.min_hamming_antigens << ";" << par.Value.max_hamming_antigens
                    << "]" << endl;
               cerr << " ... takes a random antigen sequence " << endl;
               newAG->randomize();
            }
         }
         int index = add_Antigen(newAG);
         ana << "New antigen created : " << n << " index= " << index << "(" << newAG->size
             << "):" << newAG->print() << ", with index " << index << endl;
      }
   }
   ana << "... done.\n";

   ana << "Initialize Seeder cell pool ... \n";
   ana << par.Value.totalBss << " seeders required" << endl;
   ana << " with affinity between : " << par.Value.min_initial_affinity_BCRs << " and "
       << par.Value.max_initial_affinity_BCRs << endl;
   ana << " and a similarity of less than a hamming distance of " << par.Value.max_hamming_BCRs
       << "between them" << endl;

   // creates a seed sequence with a required affinity to the antigens.
   BCR * seedSeq = new BCR(par.Value.size_sequences);
   seedSeq->randomize();
   int tries2 = 0;
   while ((tries2 < 100000)
          && ((maxAffinitySequenceToAntigens(seedSeq) < par.Value.min_initial_affinity_BCRs)
              || (maxAffinitySequenceToAntigens(seedSeq) > par.Value.max_initial_affinity_BCRs))) {
      seedSeq->randomize();
      tries2++;
   }
   if (tries2 >= 100000) {
      cerr
           <<
         "Warning: SequenceSpace, could not generate a seeding BCR sequence with required affinity to antigens"
           << endl;
   }
   ana << "  Random Seeding sequence   : " << n << " (" << seedSeq->size << "):" << ", affinity="
       << maxAffinitySequenceToAntigens(seedSeq) << endl;
   for (n = 0; n < par.Value.totalBss; n++) {
      if ((n < (int) par.Value.initBCRSeqs.size()) && (par.Value.initBCRSeqs[n].size() > 2)) {
         int index = add_Seeder(new BCR(par.Value.initBCRSeqs[n]));
         ana << "  Predefined seeder BCR seq : " << n << " index= " << index << "("
             << getSequence(index)->size << "):" << getSequence(index)->print() << ", affinity="
             << maxAffinitySequenceToAntigens(getSequence(index)) << endl;
         if (getSequence(index)->size != par.Value.size_sequences) {
            cerr << "ERR: SequenceSpace constructor, the sequence given as initial BCRs has length "
                 << getSequence(index)->size
                 << " instead of the defined size of sequences which is "
                 << par.Value.size_sequences << endl;
         }
      } else {
         // creates BCRs randomly within a half-hamming distance from a seed sequence
         // and with a specific affinity to the antigens
         BCR * newPadawanSequence = new BCR(par.Value.size_sequences);
         int tries = 0;
         bool correct = false;
         while ((tries < 10000000) && (!correct)) {
            newPadawanSequence->randomize();
            correct
               = ((2
                   * sequence::hamming(newPadawanSequence, seedSeq)) <= par.Value.max_hamming_BCRs);
            correct = correct
                      && ((maxAffinitySequenceToAntigens(newPadawanSequence)
                           > par.Value.min_initial_affinity_BCRs)
                          && (maxAffinitySequenceToAntigens(newPadawanSequence)
                              < par.Value.max_initial_affinity_BCRs));
            tries++;
         }
         int index = add_Seeder(newPadawanSequence);
         ana << "  Generated new BCR sequence: " << n << " index= " << index << "("
             << newPadawanSequence->size << "):" << newPadawanSequence->print() << ", affinity="
             << maxAffinitySequenceToAntigens(newPadawanSequence) << endl;
         if (tries >= 10000000) {
            cerr
                 <<
               "Warning : SequenceSpace, it was impossible to generate a random receptor with the expected hamming distance [<"
                 << par.Value.max_hamming_BCRs << "] and affinities : ["
                 << par.Value.min_initial_affinity_BCRs << ";"
                 << par.Value.max_initial_affinity_BCRs << "]" << endl;
         }
      }
   }
   ana << "... done.\n";

   ana << "Initialize TCR initial sequences... \n";
   ana << par.Value.totalTC << " sequences required" << endl;
   ana << " with affinity between : " << par.Value.min_initial_affinity_TCRs << " and "
       << par.Value.max_initial_affinity_TCRs << endl;
   ana << " and a similarity of less than a hamming distance of " << par.Value.max_hamming_TCRs
       << "between them" << endl;

   // creates a seed sequence with a required affinity to the antigens.
   TCR * seedSeq2 = new TCR(par.Value.size_sequences);
   seedSeq2->randomize();
   tries2 = 0;
   while ((tries2 < 100000)
          && ((maxAffinitySequenceToAntigens(seedSeq2) < par.Value.min_initial_affinity_TCRs)
              || (maxAffinitySequenceToAntigens(seedSeq2)
                  > par.Value.max_initial_affinity_TCRs))) {
      seedSeq2->randomize();
      tries2++;
   }
   if (tries2 >= 100000) {
      cerr
           <<
         "Warning: SequenceSpace, could not generate a seedinf TCR sequence with required affinity to antigens"
           << endl;
   }
   ana << "  Random Seeding sequence   : " << n << " (" << seedSeq2->size << "):"
       << seedSeq2->print() << ", affinity=" << maxAffinitySequenceToAntigens(seedSeq2) << endl;
   for (n = 0; n < par.Value.totalTC; n++) {
      if ((n < (int) par.Value.initTCRSeqs.size()) && (par.Value.initTCRSeqs[n].size() > 2)) {
         int index = add_TCR(new TCR(par.Value.initTCRSeqs[n]));
         ana << "  Predefined TCR seq : " << n << " index= " << index << "("
             << getSequence(index)->size << "):" << getSequence(index)->print() << ", affinity="
             << maxAffinitySequenceToAntigens(getSequence(index)) << endl;
         if (getSequence(index)->size
             != par.Value.size_sequences) {
            cerr << "ERR: SequenceSpace constructor, the sequence given as initial TCR has length "
                 << getSequence(index)->size
                 << " instead of the defined size of sequences which is "
                 << par.Value.size_sequences << endl;
         }
      } else {
         // creates TCRs randomly within a half-hamming distance from a seed sequence
         // and with a specific affinity to the antigens
         TCR * newPadawanSequence = new TCR(par.Value.size_sequences);
         int tries = 0;
         bool correct = false;
         while ((tries < 10000000) && (!correct)) {
            newPadawanSequence->randomize();
            correct = ((2 * sequence::hamming(newPadawanSequence,seedSeq2))
                       <= par.Value.max_hamming_TCRs);
            correct = correct
                      && ((maxAffinitySequenceToAntigens(newPadawanSequence)
                           > par.Value.min_initial_affinity_TCRs)
                          && (maxAffinitySequenceToAntigens(newPadawanSequence)
                              < par.Value.max_initial_affinity_TCRs));
            tries++;
         }
         int index = add_TCR(newPadawanSequence);

         ana << "  Generated new TCR sequence: " << n << " index= " << index << "("
             << newPadawanSequence->size << "):" << newPadawanSequence->print() << ", affinity="
             << maxAffinitySequenceToAntigens(newPadawanSequence) << endl;
         if (tries >= 10000000) {
            cerr
                 <<
               "Warning : SequenceSpace, it was impossible to generate a random TCR with the expected hamming distance [<"
                 << par.Value.max_hamming_TCRs << "] and affinities : ["
                 << par.Value.min_initial_affinity_TCRs << ";"
                 << par.Value.max_initial_affinity_TCRs << "]" << endl;
         }
      }
   }
   ana << "... done.\n";

   // Opening the file streams for output during simulation - think of calling close_files(); later
   // to finish the job at the end of simulation
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

// logdata[4].open("tbeta.out");                    logdata[4] << "! Beta fuer die Summe der
// praesentierten Epitope.\n";
// logdata[5].open("tomega.out");                   logdata[5] << "! Omega fuer die Summe der
// praesentierten Epitope.\n";
// logdata[6].open("tbpersum.out");                 logdata[6] << "! % Centroblasten an den
// praesentierten Epitopen.\n";
// logdata[7].open("topersum.out");                 logdata[7] << "! % Output-Zellen an den
// praesentierten Epitopen.\n";
// logdata[8].open("tfb.out");                      logdata[8] << "! % der selektierten
// Centroblasten.\n";
// logdata[9].open("tfo.out");                      logdata[9] << "! % der selektierten
// Output-Zellen.\n";

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
   ana << "Sequence space READY! (YEAAAAH !!)\n";
}
sequenceSpace::~sequenceSpace() {
   for (int i = 0; i < (int) poolSequences.size(); ++i) {
      if (poolSequences[i]) { delete poolSequences[i]; }
   }
}
// common with ss.cpp
void sequenceSpace::close_files() {
   cout << "Close files in shape space ... ";
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) { logdata[i].close(); }
   }
   logmeanaff.close();
   loghighaff.close();
   logdiversity.close();
   log048aff.close();

   ofstream outSeqSpace;
   outSeqSpace.open("seqspace.out");
   outSeqSpace << printSequences();
   outSeqSpace.close();
   cout << "done.\n";
}
/// the sequence should already have a type (BCR, TCR or AG), but sequence only will raise an error
long sequenceSpace::index_adding_sequence(sequence * toAdd, long int ID_of_mother_sequence) {
   if ((int) poolSequences.size() != n_Sequences) {
      cerr << "ERR: sequenceSpace::index_adding_sequence, n_Sequences management error !!" << endl;
   }
   if (toAdd == NULL) {
      cerr << "ERR: sequenceSpace::index_adding_sequence(..., sequence* = NULL !!!)" << endl;
   }
   if (toAdd->getType() == numberTypesSeq) {
      cerr
           <<
         "ERR: sequenceSpace::index_adding_sequence(), forbidden to add a sequence without type (should be TCR, Antigen or BCR)"
           << endl;
      return -1;
   }

   poolSequences.push_back(toAdd);
   if (toAdd->getType() != typeAG) {
      toAdd->max_affinity_to_antigens = maxAffinitySequenceToAntigens(toAdd);
   } else { toAdd->max_affinity_to_antigens = 0.0; }
   indexParentSeq.push_back(ID_of_mother_sequence);
   n_Sequences++;

   // NOTE : when an antigen is added, the affinities of all BCRs have to be recomputed...
   //       (not used) and the bins of this new AG has to be filled with all the producers depending
   // on their affinities
   if (toAdd->getType() == typeAG) {
      for (int pos = 0; pos < n_Sequences; ++pos) {
         if (poolSequences[pos]->getType() == typeBCR) {
            poolSequences[pos]->max_affinity_to_antigens = maxAffinitySequenceToAntigens(
               poolSequences[pos]);
            // (not used) toAdd->add_producer(sequence::seq_affinity(poolSequences[pos], toAdd,
            // R_affinity, maxSizeClusters, typeAffinityFunction),
            // poolSequences[pos]->nb_AB_producers());
         }
      }
   }
   // if you create a new BCR, nothing change because there is no yet cells with this BCR
   // if you create a new TCR, well, nobody cares !
   return n_Sequences - 1;
}
long sequenceSpace::getMutation(long from_position) {
   if ((from_position < 0) || (from_position >= n_Sequences)) {
      cerr << "sequenceSpace::getMutation(id=" << from_position
           << "), wrong id (max " << n_Sequences << " sequences" << endl;
      return -1;
   }
   sequence * theSeq = poolSequences[from_position];
   if (theSeq == NULL) {
      cerr << "sequenceSpace::getMutation, poolSequences[" << from_position
           << "] is NULL ! \n";
      return -1;
   }
   if (theSeq->getType() != typeBCR) {
      cerr << "sequenceSpace::getMutation, you are trying to mutate something else than a BCR"
           << endl;
      return -1;
   }
   BCR * newBCR = new BCR((BCR*) theSeq);
   newBCR->mutateOnePosition(); // it was an error before, it was mutating the previous sequence
   return index_adding_sequence(newBCR, from_position);
}
double sequenceSpace::affinity(long int seqId1, long int seqId2) {
   if ((seqId1 < 0) || (seqId2 < 0)) {
      cerr << "ERR : sequenceSpace::Affinity, negative indices" << endl;
      return 0;
   }
   if ((seqId1 >= n_Sequences) || (seqId2 >= n_Sequences)) {
      cerr << "ERR : sequenceSpace::Affinity(seqId1, seqId2), an index is out of scope\n";
      return 0;
   }
   double res = amplitude * sequence::seq_affinity(getSequence(seqId1), getSequence(seqId2),
                                                   R_affinity, maxSizeClusters,
                                                   typeAffinityFunction);
   if (use_logarithmic_seq_affinity) {
      return log(max(res, 1e-10));
   } else {
      return res;
   }
}
double sequenceSpace::affinity_norm(long int seqId1, long int seqId2) {
   return affinity(seqId1, seqId2) / (amplitude + 1e-9);
}
double sequenceSpace::affinity(long int seqId1, long int seqId2, double &tr) {
   double v = affinity(seqId1, seqId2);
   if ((tr < 0) || (tr >= 1)) {
      cerr << "sequenceSpace::affinity(..., ..., Threshold=" << tr
           << "), the threshold should be in [0;1]\n";
      return 0.0;
   }
   return (v - tr) / (1 - tr);
}
void sequenceSpace::correct_average_affinity(cells celltyp, long &pos, double &average) {
   // Never call this routine before the first cell of type celltyp was created in SS!
   if (sum_cell[celltyp] == 1) { average = 0; }
   // cout<<"average vorher="<<average<<"; add aff="<<best_affinity(pos)<<"; ";
   average = ((sum_cell[celltyp] - 1) * average + best_affinity_norm(pos)) / (sum_cell[celltyp]);
   // cout<<"sum="<<sum_cell[celltyp]<<"; threshold="<<average<<"\n";
}
double sequenceSpace::maxAffinitySequenceToAntigens(sequence * seq) {
   if (!seq) {
      cerr << "sequenceSpace::maxAffinitySequenceToAntigens , sequence is NULLLL\n";
      return 0;
   }
   double aff = 0;
   for (int i = 0; i < n_Antigen; ++i) {
      aff = max(aff, amplitude * sequence::seq_affinity(poolSequences[Antigens[i]], seq,
                                                        R_affinity, maxSizeClusters,
                                                        typeAffinityFunction));
   }
   if (use_logarithmic_seq_affinity) {
      return log(max(aff, 1e-10));
   } else {
      return aff;
   }
}
double sequenceSpace::best_affinity(long pos) {
   if ((pos < 0) || (pos >= n_Sequences)) {
      cerr << "ERR: sequenceSpace::bestAffinity, invalid index: " << pos
           << ". Total nb of sequences :" << n_Sequences << endl;
      return 0;
   }
   if (poolSequences[pos]->max_affinity_to_antigens <= -1) {
      poolSequences[pos]->max_affinity_to_antigens
         = maxAffinitySequenceToAntigens(poolSequences[pos]);
   }
   return poolSequences[pos]->max_affinity_to_antigens;
}
double sequenceSpace::best_affinity_norm(long pos) {
   if ((pos < 0) || (pos >= n_Sequences)) {
      cerr << "ERR: sequenceSpace::bestAffinityNorm, sequence index unknown : "
           << pos << ". Total nb of sequences :" << n_Sequences << endl;
      return 0;
   }
   return best_affinity(pos) / (amplitude + 1e-9);
}
int sequenceSpace::get_nearest_Antigen(long n) {
   double aff = -1;
   int IDclosest = -1;
   for (int i = 0; i < n_Antigen; ++i) {
      double newAff = sequence::seq_affinity(poolSequences[Antigens[i]], getSequence(
                                                n), R_affinity, maxSizeClusters,
                                             typeAffinityFunction);
      if (newAff > aff) { aff = newAff; IDclosest = i; } // do not care if it is logarithmic or not
                                                         // here
   }
   if (IDclosest < 0) {
      cerr
           <<
         "sequenceSpace::get_nearest_antigen(), could not find a closest antigen. Note : n_Antigen = "
           << n_Antigen << endl;
   }
   return IDclosest;                     // Might be more consistent to change it to return
                                         // get_antigen(IDclosest) !!!
}
/** @brief returns square distance between two points/sequences (hamming in case of sequences */
double sequenceSpace::Abstandquad(long int &seqId1, long int &seqId2) {
   if ((seqId1 < 0) || (seqId2 < 0)) {
      cerr << "ERR : sequenceSpace::Abstandquad, negative indices" << endl;
      return 0;
   }
   if ((seqId1 >= n_Sequences) || (seqId2 >= n_Sequences)) {
      cerr << "ERR : sequenceSpace::Abstandquad(seqId1, seqId2), an index is out of scope\n";
      return 0;
   }
   return std::pow(sequence::hamming(getSequence(seqId1), getSequence(seqId2)), 2.0);
}
// ============================================================
// Public routines:
// ============================================================

long int sequenceSpace::add_Antigen(string newSequence) {
   Antigen * newAG = new Antigen(newSequence);
   return add_Antigen(newAG);
}
long int sequenceSpace::add_Antigen(Antigen * newAG) {
   // problem : doesn't check if this is an antigen.
   Antigens.push_back(index_adding_sequence((sequence*) newAG, -1));
   n_Antigen++;
   return Antigens[n_Antigen - 1];
}
bool sequenceSpace::add_new_Antigen() {
   // ### to be programmed
   return 0;
}
bool sequenceSpace::add_new_Antigen(long ag_index) {
   // ### to be programmed
   return 0;
}
long int sequenceSpace::add_Seeder(string newSequence) {
   BCR * newBCR = new BCR(newSequence);
   return add_Seeder(newBCR);
}
long int sequenceSpace::add_Seeder(BCR * newBCR) {
   // problem : doesn't check this is a TCR. Should be private
   Seeders.push_back(index_adding_sequence((sequence*) newBCR, -1));
   n_Seeder++;
   return Seeders[n_Seeder - 1];
}
long int sequenceSpace::add_TCR(string newSequence) {
   TCR * newTCR = new TCR(newSequence);
   return add_TCR(newTCR);
}
long int sequenceSpace::add_TCR(TCR * newTCR) {
   TCRs.push_back(index_adding_sequence((sequence*) newTCR, -1));
   n_TCRs++;
   return TCRs[n_TCRs - 1];
}
// ============================================================
// Dateien schreiben
// ============================================================

double sequenceSpace::mean_affinity(double * affinities) {
   long i;
   double aff,all;
   for (i = 0; i < 15; i++) {
      affinities[i] = 0.;
   }
   for (i = n_Antigen; i < n_Sequences; i++) {
      // Find nearest antigen-type
      aff = best_affinity_norm(i);
      if (aff > 0.3) {
         affinities[3] += double (getSequence(i)->n_cell[sCB]);
         affinities[4] += double (getSequence(i)->n_cell[sCC]);
         affinities[5] += double (getSequence(i)->n_cell[sout]);
      }
      if (aff < 0.4) {
         affinities[6] += double (getSequence(i)->n_cell[sCB]);
         affinities[7] += double (getSequence(i)->n_cell[sCC]);
         affinities[8] += double (getSequence(i)->n_cell[sout]);
      }
      if ((aff >= 0.4) && (aff < 0.8)) {
         affinities[9] += double (getSequence(i)->n_cell[sCB]);
         affinities[10] += double (getSequence(i)->n_cell[sCC]);
         affinities[11] += double (getSequence(i)->n_cell[sout]);
      }
      if (aff >= 0.8) {
         affinities[12] += double (getSequence(i)->n_cell[sCB]);
         affinities[13] += double (getSequence(i)->n_cell[sCC]);
         affinities[14] += double (getSequence(i)->n_cell[sout]);
      }
      // add corresponding affinity with weight #cb and #cc to ccaff and cbaff
      affinities[0] += aff * double (getSequence(i)->n_cell[sCB]);
      affinities[1] += aff * double (getSequence(i)->n_cell[sCC]);
      affinities[2] += aff * double (getSequence(i)->n_cell[sout]);
   }
   // Divide through total cell numbers:
   all = affinities[0] + affinities[1];
   if (sum_cell[sCB] > 0) {
      affinities[0] /= sum_cell[sCB];
      affinities[3] /= sum_cell[sCB];
      affinities[6] /= sum_cell[sCB];
      affinities[9] /= sum_cell[sCB];
      affinities[12] /= sum_cell[sCB];
   } else { affinities[0] = 0.; }
   if (sum_cell[sCC] > 0) {
      affinities[1] /= sum_cell[sCC];
      affinities[4] /= sum_cell[sCC];
      affinities[7] /= sum_cell[sCC];
      affinities[10] /= sum_cell[sCC];
      affinities[13] /= sum_cell[sCC];
   } else { affinities[1] = 0.; }
   if (sum_cell[sout] > 0) {
      // cout<<outs<<" "<<ha_out<<" "<<sum_cell[sout]<<"\n";
      affinities[2] /= sum_cell[sout];
      affinities[5] /= sum_cell[sout];
      affinities[8] /= sum_cell[sout];
      affinities[11] /= sum_cell[sout];
      affinities[14] /= sum_cell[sout];
      // cout<<outs<<" "<<ha_out<<"\n";
   } else { affinities[2] = 0.; }
   if (sum_cell[sCB] + sum_cell[sCC] > 0) {
      all /= (sum_cell[sCB] + sum_cell[sCC]);
   } else { all = 0.; }
   return all;
}
void sequenceSpace::mean_affinities_ag(double * mean_affinities, int ag_index) {
   double aff, nBCs, nOUT;
   double totalBC[8];
   // intialise affinities
   for (int i = 0; i < 8; i++) {
      mean_affinities[i] = 0.;
      totalBC[i] = 0.;
   }
   // go through all points of the AffinitySpace
   for (long i = 0; i < n_Sequences; i++) {
      // load affinities[] with corresponding values
      /* [0]: mean affinity of CB+CC to Antigen[ag_index]
       *  [1]: as [0] but only including cells with aff>0.1
       *  [2]: as [0] but only including cells with aff>0.2
       *  [3-5]: as [0-2] for accummulated output
       *  [6]: fraction of CB+CC cells with aff>0.3
       *  [7]: fraction of OUT cells with aff>0.3
       */
      // determine the affinity of this point to the antigen denoted by ag_index
      aff = affinity_norm(i, get_Antigen(ag_index));
      // get the total number of BCs at this point
      nBCs = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // get the total number of accumulated output at this point
      nOUT = getSequence(i)->n_cell[sout];
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
   }
   // normalise to get the mean affinities
   for (int i = 0; i < 6; i++) {
      if (totalBC[i] > 0) {
         mean_affinities[i] /= totalBC[i];
      } else {
         mean_affinities[i] = 0.;
      }
   }
   // get the fraction of aff>0.3 for BC
   if (totalBC[0] > 0) {
      mean_affinities[6] = totalBC[6] / totalBC[0];
   } else {
      mean_affinities[6] = 0.;
   }
   // ... and for accumulated output
   if (totalBC[3] > 0) {
      mean_affinities[7] = totalBC[7] / totalBC[3];
   } else {
      mean_affinities[7] = 0.;
   }
}
void sequenceSpace::to_multiag_files(double time) {
   ofstream
      saffin_ags, saffin_t10, saffin_t20,
      saffin_out_ags, saffin_out_t10, saffin_out_t20,
      shaffin_ag, shaffin_out_ag,
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
   shaffin_out_ag.open("shaffin_out_ags.out", ofstream::app);
   shaffin_out_ag << time << "  ";
   cross_gcr.open("cross_reactivity_gcr.out", ofstream::app);
   cross_gcr << time << "  ";
   cross_out.open("cross_reactivity_out.out", ofstream::app);
   cross_out << time << "  ";
   double fracsum = 0., fracOUTsum = 0.;
   // variables for cross reactivity
   double cross_reactivity[6];
   for (int i = 0; i < 6; i++) {
      cross_reactivity[i] = 0.;
   }
   for (int a = 0; a < n_Antigen; a++) {
      double mean_affinities[8];
      mean_affinities_ag(mean_affinities, a);
      saffin_ags << mean_affinities[0] << "  ";
      saffin_t10 << mean_affinities[1] << "  ";
      saffin_t20 << mean_affinities[2] << "  ";
      saffin_out_ags << mean_affinities[3] << "  ";
      saffin_out_t10 << mean_affinities[4] << "  ";
      saffin_out_t20 << mean_affinities[5] << "  ";
      shaffin_ag << mean_affinities[6] << "  ";
      shaffin_out_ag << mean_affinities[7] << "  ";
      fracsum += mean_affinities[6];
      fracOUTsum += mean_affinities[7];
      // make the sum of all antigens to get cross-reactivity (with thresholds) for BC and out
      for (int i = 0; i < 6; i++) {
         cross_reactivity[i] += mean_affinities[i];
      }
   }
   // add the sum of the fractions, which is a measure of cross-reactivity, to the fraction-file
   shaffin_ag << fracsum << "\n";
   shaffin_ag.close();
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
//MS
void sequenceSpace::meanSD_out_affinity(double * aff_mean, double * aff_sd) {cerr<<"notcoded";exit(1);}
double sequenceSpace::best_distance(long pos) {cerr<<"notcoded";exit(1);}

void sequenceSpace::get_diversity(double time) {
   logdiversity << time << "   ";
   int diversity[number];
   for (int j = 0; j < number; j++) {
      diversity[j] = 0;
   }
   for (long i = 0; i < n_Sequences; i++) {
      if (getSequence(i)->n_cell[sCB] > 0.5) { ++diversity[sCB]; }
      if (getSequence(i)->n_cell[sCC] > 0.5) { ++diversity[sCC]; }
      if (getSequence(i)->n_cell[sout] > 0.5) { ++diversity[sout]; }
      if (getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC]
          + getSequence(i)->n_cell[sout] > 0.5) { ++diversity[total]; }
      if (getSequence(i)->n_cell[soutext] + getSequence(i)->n_cell[soutextproduce]
          > 0.5) { ++diversity[soutext]; }
   }
   logdiversity << diversity[total] << "   "
                << diversity[sCB] << "   "
                << diversity[sCC] << "   "
                << diversity[sout] << "   "
                << diversity[soutext] << "\n";
}
void sequenceSpace::to_ssfiles(double time) {
   // sum_check();
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) {
         if (i != soutdec) { logdata[i] << time << "     " << sum_cell[i]; }
         if (i == sout) {
            logdata[sout] << "    " << sum_cell[sout] - oldsum_cell[sout];
         } else if (i == soutdec) {
            logdata[sout] << "   " << sum_cell[soutdec] << "   ";
            if (sum_cell[sout]
                > 0) {
               logdata[sout] << sum_cell[soutdec] / sum_cell[sout] << "\n";
            } else { logdata[sout] << "0\n"; }
         } else if (i == sallapoptosis) {
            logdata[sallapoptosis] << "    " << sum_cell[sallapoptosis]
               - oldsum_cell[sallapoptosis] << "\n";
         } else { logdata[i] << "\n"; }
      }
      oldsum_cell[i] = sum_cell[i];
   }

   double affinities[15];
   double a;
   a = mean_affinity(affinities);
   logmeanaff << time << "   " << a << "   " << affinities[0] << "   "
              << affinities[1] << "   " << affinities[2] << "\n";
   loghighaff << time << "   " << affinities[3] << "   "
              << affinities[4] << "   " << affinities[5] << "\n";
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
string sequenceSpace::printSequences(bool showSequencesWithNoAliveCells, bool showTree) {
   stringstream ss;
   ss << n_Antigen << " sequences of antigens\n";
   for (int i = 0; i < n_Antigen; ++i) {
      ss << "Ag nr " << i << " ,ID=\t" << get_Antigen(i) << "\t("
         << getSequence(get_Antigen(i))->size << ")\t" << getSequence(get_Antigen(i))->print()
         << "\n";
   }
   ss << n_Seeder << " sequences of seeders\n";
   for (int i = 0; i < n_Seeder; ++i) {
      ss << "Seeder BCR nr " << i << " ,ID=\t" << get_Seeder(i) << "\t("
         << getSequence(get_Seeder(i))->size << ")\tAff=\t"
         << getSequence(get_Seeder(i))->max_affinity_to_antigens
         << "\t" << getSequence(get_Seeder(i))->print() << "\n";
   }
   ss << n_TCRs << " sequences of TCRs\n";
   for (int i = 0; i < n_TCRs; ++i) {
      ss << "TCR nr " << i << " ,ID=\t" << get_TCR(i) << "\t(" << getSequence(get_TCR(i))->size
         << ")\t" << getSequence(get_TCR(i))->print() << "\n";
   }

   vector<vector<long> > tree;
   tree.resize(n_Sequences);
   for (int i = 0; i < n_Sequences; ++i) {
      if (indexParentSeq[i] >= 0) { tree[indexParentSeq[i]].push_back(i); }
   }

   ss << n_Sequences - n_Antigen - n_TCRs << " sequences of BCRs (seeders or not)\n";
   if (!showSequencesWithNoAliveCells) {
      ss << "List of sequences with alive cells carrying them:\n";
   }
   for (int i = 0; i < n_Sequences; ++i) {
      if (getSequence(i)->getType() == typeBCR) {
         if (showSequencesWithNoAliveCells
             || ((getSequence(i)->n_cell[sCB]) + (getSequence(i)->n_cell[sCC])
                 + (getSequence(i)->n_cell[sout]) > 0)) {
            sequence * seq = getSequence(i);     // doesn't need the cast actually
            ss << "ID=\t" << i << "\tNbCellsTot(BCO)=\t" << getSequence(i)->n_cell[sCB]
               + getSequence(i)->n_cell[sCC] + getSequence(i)->n_cell[sout] << "\t("
               << getSequence(i)->n_cell[sCB] << "-" << getSequence(i)->n_cell[sCC] << "-"
               << getSequence(i)->n_cell[sout] << ")\t(" << seq->size << "),affinity=\t"
               << best_affinity(i) << "\t" << seq->print() << "\n";
            if (showTree) {
               if (tree[i].size() > 0) {
                  ss << "\tDaughters are :\n";
                  for (int j = 0; j < (int) tree[i].size(); ++j) {
                     ss << "\t-->ID=\t" << tree[i][j] << "\t,nCC=\t"
                        << getSequence(tree[i][j])->n_cell[sCC]
                        << "\taffinity=\t" << best_affinity(tree[i][j]) << "\tImproving=\t"
                        << best_affinity(tree[i][j]) - best_affinity(i) << "\n";
                  }
               }
            }
         }
      }
   }

   ss << external_cells.size()
      << " Sequences IDs where there are external cells (producing or not) \n";
   for (int i = 0; i < (int) external_cells.size(); ++i) {
      ss << "ExtCell\t" << i << " ,ID=\t" << external_cells[i] << "\t";
      if (getSequence(external_cells[i])->getType() != typeBCR) {
         ss << "ERR: this is of type "
            << getSequence(external_cells[i])->typeInString() << endl;
      } else {
         BCR * bseq = (BCR*) getSequence(external_cells[i]);
         ss << "n_soutext=\t" << bseq->n_cell[soutext] << "\tn_soutextproduce:\t"
            << bseq->n_cell[soutextproduce] << endl;
      }
   }

   ss << ab_producers.size() << " Sequences IDs where there are producing cells \n";
   for (int i = 0; i < (int) ab_producers.size(); ++i) {
      ss << "ProducerCel\t" << i << " ,ID=\t" << ab_producers[i] << "\t";
      if (getSequence(ab_producers[i])->getType() != typeBCR) {
         ss << "ERR: this is of type "
            << getSequence(ab_producers[i])->typeInString() << endl;
      } else {
         BCR * bseq = (BCR*) getSequence(ab_producers[i]);
         ss << "n_soutext=\t" << bseq->n_cell[soutext] << "\tn_soutextproduce:\t"
            << bseq->n_cell[soutextproduce] << endl;
      }
   }

   return ss.str();
}
/// Philippe 17.05.16 adapting the output functions for affinity histograms

int sequenceSpace::intN_mutation(sequence * x, sequence * y) {
   return sequence::hamming(x,y);
}
int sequenceSpace::best_hamming_distance(long n) {
   int distance = 1.0e+09;  // stupid large number
   for (int nAg = 0; nAg < get_n_Antigen(); nAg++) {
      distance = min(distance, intN_mutation(getSequence(n),getSequence(get_Antigen(nAg))));
   }
   return distance;
}
double sequenceSpace::mean_hamming_distance(long n) {
   int sumdistance = 0;
   for (int nAg = 0; nAg < get_n_Antigen(); nAg++) {
      sumdistance += intN_mutation(getSequence(n),getSequence(get_Antigen(nAg)));
   }
   if (get_n_Antigen() > 0) {
      return double (sumdistance) / double (get_n_Antigen());
   }
   return 0;
}
double sequenceSpace::N_mutation(sequence * k, sequence * l) {
   double hamm = sequence::hamming(k,l);
   return hamm * hamm;
}
double sequenceSpace::mean_affinity_norm(long pos) {
   double tmp = 0.;
   for (int j = 0; j < get_n_Antigen(); j++) {
      tmp += affinity_norm(pos, get_Antigen(j));
   }
   tmp /= double (get_n_Antigen());
   return tmp;
}
void sequenceSpace::write_gcbc_hamming(double time) {
   int max_hamming = size_sequences;
   vector<int> nbc(max_hamming + 1, 0);          // histogram counters for each Hamming distance,
   vector<int> nbc_mean_ag(max_hamming + 1, 0);  // .. and for the mean Hamming distance for all
                                                 // antigens/epitopes

   // remove this line in the general case:
   /// Philippe : don't understand this warning
   cout << "WARNING! Modified wrong code in SS::best_hamming_distance(long&)\n";

   int nbchere = 0;
   int hamming = 0;
   double hamming_mean = 0.;
   for (long i = 0; i < n_Sequences; i++) {
      // get the total number of GC-BCs at this point
      nbchere = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // determine the Hamming distance of this point to the nearest antigen
      hamming = best_hamming_distance(i);
      hamming_mean = mean_hamming_distance(i);
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      if (hamming >= max_hamming) {
         nbc[max_hamming] += nbchere;
      } else {
         nbc[hamming] += nbchere;
      }
      // now do rounding of the mean hamming distance of sequence i:
      hamming = int (hamming_mean + 0.5);
      // ... and attribute to the
      if (hamming >= max_hamming) {
         nbc_mean_ag[max_hamming] += nbchere;
      } else {
         nbc_mean_ag[hamming] += nbchere;
      }
   }
   // open the output file (append)
   ofstream gcbc_hamming;
   gcbc_hamming.open("gcbc_hamming.out", ofstream::app);
   for (int n = 0; n <= max_hamming; n++) {
      gcbc_hamming << time / 24. << "   " << time << "   " << n
                   << "   " << nbc[n] << "   " << nbc_mean_ag[n] << "\n";
   }
   gcbc_hamming.close();
}
void sequenceSpace::write_gcbc_affinity(double time) {
   // +++++++++++++++ OPTION ++++++++++++++++
   double logfactor = 3.0;
   int n_bins = 10;
   double affarray[n_bins];
   affarray[0] = 0;
   affarray[1] = 0.0002 / logfactor;
   for (int b = 2; b < n_bins; b++) {
      affarray[b] = logfactor * affarray[b - 1];
   }
   /* for (int b = 0; b < n_bins; b++) {
    *  cout << affarray[b] << ",";
    *  } cout << "\n";  */
   // +++++++++++ end OPTION ++++++++++++++++

   // double check that the whole range of binding probabilities is covered
   if (affarray[n_bins - 1] * logfactor < 1.0) {
      cerr << "In sequenceSpace::write_gcbc_affinity(double): largest affinity is smaller than 1.\n"
           << "Abort simulation.\n";
      exit(1);
   }
   // define a set of counters for each Hamming distance
   vector<int> nbc(n_bins, 0);
   // .. and for the mean Hamming distance for all antigens/epitopes
   vector<int> nbc_mean_ag(n_bins, 0);
   int nbchere = 0;
   double affbest = 0;
   double affmean = 0.;
   // Go through all points of the shape space and add to the counters:
   for (long i = 0; i < n_Sequences; i++) {
      // get the total number of GC-BCs at this point
      nbchere = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // determine the Hamming distance of this point to the nearest antigen
      affbest = best_affinity_norm(i);
      affmean = mean_affinity_norm(i);
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      int binbest = 0, binmean = 0;
      while (binbest<n_bins&&affbest> affarray[binbest]) {
         ++binbest;
      }   // now binbest is one too far, as affbest>0, binbest>0 as well
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
      gcbc_affinity << time / 24. << "   " << time << "   " << affarray[n] << "   "
                    << nbc[n] << "   " << nbc_mean_ag[n] << "\n";
   }
   gcbc_affinity.close();
}
/// End 17.05.16

bool compSequences(pair<double, sequence*> a, pair<double, sequence*> b) {
   return a.first > b.first;
}
string sequence::testeAffinityFunctions(double L, double R, int maxClusters,
                                        int typeAffinityFunction) {
   stringstream out;
   out << "Testing the properties ot the affinity function for the following parameters : \n";
   out << "   ->     L= " << L << "\t(Size of sequences)" << endl;
   out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
   out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
   switch (typeAffinityFunction) {
      case seqAff: {
         out << "   -> Using standard affinity (Saham's)\n";
         break;
      }

      case seqAffNorm: {
         out << "   -> Using standard affinity normalized by maxCl^r\n";
         break;
      }

      case seqAffWindow: {
         out << "   -> Using the maximum affinity of a sliding window\n";
         break;
      }
   }

   out
       <<
      "==== Part 1 : enumerates all (if possible), or a lot of sequences and sort them by affinity to get the best ones : ===="
       << endl;

    #define resolutiondistrib 100
    #define maxSequencesToEnumeate 1100000

   vector<pair<double, sequence*> > store;
   sequence * ref = new sequence(L);
   vector<double> distribution;
   distribution.resize(resolutiondistrib + 1);

   int total = 0;
   int maxim = pow(2, L);
   if (L > 26) {
      maxim = maxSequencesToEnumeate + 1;           // to avoid it to become negative ...
   }
   bool enumerateAll = (maxim < maxSequencesToEnumeate);
   if (enumerateAll) {
      for (int i = 0; i < maxim; ++i) {
         sequence * a = new sequence(L, (long) i);
         double affi = sequence::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
         distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0; // put into the histogram
         store.push_back(pair<double, sequence*> (affi,a));
         total++;
      }
   } else {
      for (int i = 0; i < maxSequencesToEnumeate; ++i) {
         sequence * a = new sequence(L);
         a->randomize();
         double affi = sequence::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
         distribution[(int) (((double) resolutiondistrib) * affi)] += 1.0;
         store.push_back(pair<double, sequence*> (affi,a));
         total++;
      }
   }

   out << "Distribution of affinities\n";
   for (int i = 0; i < (int) distribution.size(); ++i) {
      distribution[i] /= (double) total;
      out << i << "\t" << distribution[i] << "\t" << double (i)
         * (1.0 / (double) resolutiondistrib) << "\t" << double (i + 1)
         * (1.0 / (double) resolutiondistrib) << endl;
   }

   out << "\nSequences and affinity, "
       << ((enumerateAll) ? " in the order of ID\n" : " randomly generated\n");
   for (int i = 0; i < 200; ++i) {
      out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
   }

   out << "\nSequences, sorted from the best, out of the "
       << min(maxim,(int) maxSequencesToEnumeate) << " evaluated sequences\n";
   std::sort(store.begin(), store.end(), compSequences);
   for (int i = 0; i < 200; ++i) {
      out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
   }
   if (enumerateAll) {
      out << "\nAffinity of sequences taken randomly\n";
      for (int i = 0; i < 100; ++i) {
         sequence * seqtmp = new sequence(L);
         seqtmp->randomize();
         out << i << "\t" << sequence::seq_affinity(seqtmp,
                                                    ref,
                                                    R,
                                                    maxClusters,
                                                    typeAffinityFunction) << "\t"
             << seqtmp->print() << "\n";
      }
   }
   for (int i = 0; i < (int) store.size(); ++i) {
      delete store[i].second;
   }

   out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;

   int nbAntigens = 10;
   out << "Generating randomly " << nbAntigens << " antigens " << endl;
   vector<sequence*> ags;
   for (int i = 0; i < nbAntigens; ++i) {
      sequence * seq = new sequence(L);
      seq->randomize();
      ags.push_back(seq);
      out << "\tAg nr " << i << "\t" << seq->print() << endl;
   }
   out << "\nNumber of antigens recognized by randomly generated sequences, based on threshold\n";

   out << "  -> (for the first 100 sequences : ) In the case of random sequences" << endl;
   total = 0;
    #define thresholdRecoAg 0.1
   int nbDiscardedSeq = 0;  // sequences that don't recognize anything
   int countprint = 0;
   for (int k = 0; k < min(maxim, (int) maxSequencesToEnumeate); ++k) {
      if (k == 100) {
         out
             <<
            "  -> (for the remaining sequences) for sequences recognizing at least an antigen with affinity 0.1"
             << endl;
      }
      total++;

      // for each sequence,
      bool recoAtLeastOne = false;
      vector<double> nbRecoDepThresh(10, 0.0);
      vector<double> affinityEach(nbAntigens, 0.0);
      sequence * seqtmp = new sequence(L);
      seqtmp->randomize();
      for (int j = 0; j < nbAntigens; ++j) {
         double thisAff = sequence::seq_affinity(seqtmp,
                                                 ags[j],
                                                 R,
                                                 maxClusters,
                                                 typeAffinityFunction);
         if ((thisAff > thresholdRecoAg) || (k < 100)) {
            recoAtLeastOne = true;
         } else { nbDiscardedSeq++; }
         affinityEach[j] = thisAff;
         for (int i = 0; i <= (int) (9.99 * thisAff); ++i) {
            if (i < 10) { nbRecoDepThresh[i]++; }
         }
      }
      if (recoAtLeastOne && (countprint < 5000)) {
         countprint++;
         out << "RandSeq " << k << ", " << seqtmp->print() << " ";
         out << "nbAgPerThreshold:";
         for (int i = 0; i < 10; ++i) {
            out << "\t" << nbRecoDepThresh[i];
         }
         out << "\taffPerAg:";
         for (int i = 0; i < nbAntigens; ++i) {
            out << "\t" << affinityEach[i];
         }
         out << endl;
      }
      delete seqtmp;
   }
   out << "   ... Nb of sequences analyzed: " << total << endl;
   out << "   ... Nb of sequences discarded: " << nbDiscardedSeq
       << "(except the 100 first ones, i.e. among the :" << total - 100 << " remaining)" << endl;

   out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;

   sequence * start = new sequence(L);    // starting by '0000' : the best sequence
   out << "NbMut\tsequence\taffinity\n";
   for (int i = 0; i < 2 * L; ++i) {
      out << i << "\t" << start->print() << "\t" << sequence::seq_affinity(start,
                                                                           ref,
                                                                           R,
                                                                           maxClusters,
                                                                           typeAffinityFunction)
          << endl;
      start->mutateOnePosition();
   }

   out << "\tReaching a good affinity\t" << endl;
   start->randomize();
   double prevaff = sequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);

   bool stop = false;
   for (int i = 0; (i < L) && (!stop); ++i) {
      out << "sequence : " << start->print() << "\tAff:\t" << prevaff << "\t";
      out << "PossibleMut:";
      vector<int> posGoodMutations;
      for (int i = 0; i < L; ++i) {
         sequence stmp = sequence(start);
         stmp.content[i] = !stmp.content[i];
         double newaff
            = sequence::seq_affinity(&stmp, ref, R, maxClusters, typeAffinityFunction);
         out << "\t" << newaff;
         if (newaff > prevaff) { posGoodMutations.push_back(i); }
      }
      out << endl;
      if (posGoodMutations.size() > 0) {
         int nextmut = irandom(posGoodMutations.size() - 1);
         start->content[posGoodMutations[nextmut]] = !start->content[posGoodMutations[nextmut]];
         prevaff = sequence::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
      } else {
         stop = true;
      }
   }

   return out.str();
}
void sequenceSpace::testSequenceSpace() {
   // use cerr to catch errors and avoid delays between cout/cerr
    #define OutputTest cerr
   OutputTest << "============ Testing sequences ... ============\n" << endl;
   sequence * a = new sequence(10);
   OutputTest << "Empty sequence  a :" << a->print() << endl;
   a->mutateOnePosition();
   OutputTest << "One Mutation    a :" << a->print() << endl;
   a->mutateOnePosition();
   OutputTest << "Another one     a :" << a->print() << endl;
   a->randomize();
   OutputTest << "Randomized :    a :" << a->print() << endl;
   sequence * b = new sequence(string("01111111111110"));
   OutputTest << "New sequence    b :" << b->print() << endl;
   // b->mutate(0.5);
   OutputTest << "Mutate 50%/base b:" << b->print() << endl;
   sequence * c = new sequence(b);
   OutputTest << "New sequence  c=b :" << c->print() << endl;
   OutputTest << ((sequence::compSeq(b,c)) ? "c equals b" : "problem : c != b") << endl;
   OutputTest << "affinity b-c (r=2): " << sequence::seq_affinity(b, c, 2.0, -1, seqAff) << endl;
   sequence * d = new sequence(string("1111100000"));
   sequence * e = new sequence(string("1111011111"));
   OutputTest << "affinity d-e (r=3): " << sequence::seq_affinity(d, e, 3.0, -1, seqAff)
              << " between " << d->print() << "\t" << e->print() << endl;
   OutputTest << "hamming(d,e) = " << sequence::hamming(d,e) << endl;

   // Antigen::number_of_bins = 10;
   OutputTest << "Getting type of a sequence : in enum, " << e->getType() << " and as string: "
              << e->typeInString() << "\t" << e->print() << endl;
   BCR * s1 = new BCR(10);
   OutputTest << "Getting type of a BCR      : in enum, " << s1->getType() << " and as string: "
              << s1->typeInString() << "\t" << s1->print() << endl;
   TCR * s2 = new TCR(10);
   OutputTest << "Getting type of a TCR      : in enum, " << s2->getType() << " and as string: "
              << s2->typeInString() << "\t" << s2->print() << endl;
   Antigen * s3 = new Antigen(10);
   OutputTest << "Getting type of an Antigen : in enum, " << s3->getType() << " and as string: "
              << s3->typeInString() << "\t" << s3->print() << endl;
   /*s1->add_producer(0.15, 4);
    * s1->add_producer(0.35, 5);
    * s1->rem_producer(0.31, 2);
    * s1->rem_producer(0.0, 1);
    * s1->add_producer(0.0, 2);
    * s1->printNbCells();
    * s1->print();*/

   OutputTest << "============ Testing sequenceSpace ... ============\n" << endl;

   ofstream ana("testFile.out");
   Parameter par;
   par.Value.size_sequences = 10;
   par.Value.init_antigen_sequences = 8;
   par.Value.initAntigenSeqs.resize(100, string(""));
   par.Value.initAntigenSeqs[0] = string("0000000001");
   par.Value.initAntigenSeqs[0] = string("-1");
   par.Value.max_hamming_antigens = 2;
   par.Value.min_hamming_antigens = 1;
   par.Value.totalBss = 5;  // nb seeder cells
   par.Value.initBCRSeqs.resize(100, string(""));
   par.Value.initBCRSeqs[0] = string("1111111110");
   par.Value.initBCRSeqs[0] = string("-1");
   par.Value.max_hamming_BCRs = 2;
   par.Value.min_initial_affinity_BCRs = 0.1;
   par.Value.max_initial_affinity_BCRs = 1.0;
   par.Value.initTCRSeqs.resize(100, string(""));
   par.Value.initTCRSeqs[0] = string("1111111110");
   par.Value.initTCRSeqs[0] = string("-1");
   par.Value.max_hamming_TCRs = 2;
   par.Value.min_initial_affinity_TCRs = 0.1;
   par.Value.max_initial_affinity_TCRs = 1.0;
   par.Value.R_affinity = 2.0;
   par.Value.pm_differentiation_time = 0.5;
   ana << "hi" << endl;
   sequenceSpace sp(par, ana);

   OutputTest << "Sequence space at initialisation : " << endl;
   OutputTest << sp.printSequences() << endl;

   sequence * newseq = new sequence("0001111110");
   BCR * aBCR = new BCR(newseq);
   BCR * bBCR = new BCR("1100111000");

   OutputTest
              <<
      "when trying to add a sequence (without type) to the sequencespace, should raise an error. "
              << endl;
   sp.index_adding_sequence(newseq);   // should raise an error

   OutputTest << endl;
   long newIda = sp.index_adding_sequence(aBCR);   // should raise an error
   long newIdb = sp.index_adding_sequence(bBCR);   // should raise an error
   OutputTest
      <<
      "inserting two BCRs to the space, and got the IDs, and best affinity to antigens (automatically updated :"
      << endl;
   OutputTest << sp.getSequence(newIda)->print() << " with ID " << newIda << "\t"
              << sp.best_affinity(newIda) << endl;
   OutputTest << sp.getSequence(newIdb)->print() << " with ID " << newIdb << "\t"
              << sp.best_affinity(newIdb) << endl;

   sp.add_cell(soutext, newIda);
   sp.add_cell(soutext, newIda);
   sp.add_cell(soutext, newIda);
   OutputTest
              <<
      "After adding three cell to the first BCR sequence, now list of cells (n_cell) for this sequence:"
              << endl;
   OutputTest << sp.getSequence(newIda)->printNbCells() << endl;
   cout << "differentiating the three cells for dt = 0.5 with pm_differentiation_time = 0.5\n";
   sp.PM_differentiate(soutext, soutextproduce, 0.5);
   OutputTest << "new state of the sequence space : " << endl;
   OutputTest << sp.printSequences(true) << endl;

   Antigen * newA = new Antigen(string("01000001010"));
   Antigen * newB = new Antigen(string("01111001010"));
   long idNewA = sp.add_Antigen(newA);
   long idNewB = sp.add_Antigen(newB);
   OutputTest << "Now, adding two antigens : " << newA->print() << " ID=" << idNewA << "\t"
              << newB->print() << " ID=" << idNewB << endl;
   OutputTest << sp.printSequences(true) << endl;
}
string sequence::typeCellInString(cells index) {
   switch (index) {
      case sCB:
         return string("Centroblasts  ");   // Centroblasts

      case sCC:
         return string("TotalCentroc. ");   // Centrocytes (Sum of all following sCC)

      case sCCunselected:
         return string("Centr.Unselec.");

      case sCCcontact:
         return string("Centr.ContaAG.");

      case sCCFDCselected:
         return string("Centr.FDCselec");

      case sCCTCcontact:
         return string("Centr.TCcontac");

      case sCCselected:
         return string("Centr.Selected");

      case sCCapoptosis:
         return string("Centr.Apoptosi");   // die toten Zellen auf dem Gitter (nicht alle toten)

      case sallapoptosis:
         return string("AllDeadCells. ");   // alle jemals gestorbenen Zellen

      case sFDC:
         return string("FDCs.         ");

      case sTcell:
         return string("T cells       ");

      case sout:
         return string("Cells to goout");   // Zahl der output-Zellen die erzeugt wurden (alle)

      case soutext:
         return string("Cells thatLeft");   // Zahl der ausgeworfenen output-Zellen (ohne
                                            // Raumgitter)

      case soutextproduce:
         return string("CellsExtProdAB");   // Zahl der ausgeworfenen output-Zellen, die Antikoerper
                                            // produzieren

      case total:
         return string("TotalNbCells  ");   // Zahl lebender CB,CC,out im AffinitySpace auf und
                                            // ausserhalb des Gitters

      case soutdec:
         return string("DEC205+ outCel");   // Zahl der DEC205+ output-Zellen, die erzeugt wurden
                                            // (wie sout nur dec)

      case number:
         return string("not a type !  ");
   }
   return string("Error");
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///
///
///
///         5 - ARUP SPACE
///
///
///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int arupSpace::arup_nb_ini_antigens;
double arupSpace::arup_threshold_activation;
double arupSpace::arup_mutation;
double arupSpace::arup_proba_lethal_mut;
double arupSpace::arup_proba_affecting_mut;
double arupSpace::arup_proba_silent_mut;

arupSpace::arupSpace(Parameter &par, ofstream &ana) : generalSequenceContainer<arupProtein> (),
   OUT_haffinity(0),
   OUT_steepness(0),
   CB_haffinity(0),
   CC_haffinity(0) {
   // Antigen::number_of_bins = par.Value.antibodies_resolution + 2;
   long int n;
   ana << "Initialize arupProtein Space fields ... \n";
   pm_differentiation_rate = log(2.) / par.Value.pm_differentiation_time;
   amplitude = par.Value.amplitudeGauss;
   ana << "   pm differentiation rate = " << pm_differentiation_rate << "\n";
   ana << "   amplitude (coefficient to affinities) = " << amplitude << "\n";

   arup_threshold_activation = par.Value.arup_threshold_activation;
   arup_nb_ini_antigens = par.Value.arup_nb_ini_antigens;
   arup_mutation = par.Value.arup_mutation;
   arup_proba_lethal_mut = par.Value.arup_proba_lethal_mut;
   arup_proba_affecting_mut = par.Value.arup_proba_affecting_mut;
   arup_proba_silent_mut = par.Value.arup_proba_silent_mut;
   ana << "   Threshold of activation = " << arup_threshold_activation << "\n";
   ana << "   Probability of mutation = " << arup_mutation << endl;
   ana << "     - probability of it being lethal    = " << arup_proba_lethal_mut << endl;
   ana << "     - probability of it to be affecting = " << arup_proba_affecting_mut << endl;
   ana << "     - probability of it to be silent    = " << arup_proba_silent_mut << endl;
   ana << "   Initial number of antigens defined = " << arup_nb_ini_antigens << "\n";
   ana << "   In case antigens have to be generated automatically, they will have "
       << par.Value.arup_nb_mutations_gen_strains << " mutations " << endl;

   // vector<double> arup_ag_fraction;

   ana << "Initialize Antigen pool (Arup)...\n ";
   int nbAg = par.Value.arup_ini_antigens.size();
   for (n = 0; n < min(nbAg, arup_nb_ini_antigens); n++) {
      int index = add_Antigen(par.Value.arup_ini_antigens[n]);
      ana << "   ... predefined seq nr" << n << " index= " << index << "("
          << getSequence(index)->size << "):" << getSequence(index)->print() << endl;
   }
   for (n = nbAg; n < arup_nb_ini_antigens; n++) {
      ArupAntigen * AAG = new ArupAntigen();
      AAG->randomize(par.Value.arup_nb_mutations_gen_strains);
      int index = add_Antigen(AAG);
      ana << "   ... Generated antigen " << n << " index= " << index << "("
          << getSequence(index)->size << "):" << getSequence(index)->print() << endl;
   }
   ana << "... done.\n";

   ana << "Initialize Seeder cell pool ... \n";
   ana << par.Value.totalBss << " seeders required" << endl;
   int nbBCRs = par.Value.arup_ini_bcrs.size();
   for (n = 0; n < min(nbBCRs, (int) par.Value.totalBss); n++) {
      int index = add_Seeder(par.Value.arup_ini_bcrs[n]);
      ana << "   ... predefined seq nr" << n << " index= " << index << "("
          << getSequence(index)->size << "):" << getSequence(index)->print() << endl;
   }
   for (n = nbBCRs; n < (int) par.Value.totalBss; n++) {
      ArupBCR * ABCR = new ArupBCR();
      ABCR->randomize(arup_threshold_activation, false);
      int index = add_Seeder(ABCR);
      ana << "   ... Generated BCR " << n << " index= " << index << "("
          << getSequence(index)->size << "):" << getSequence(index)->print() << endl;
   }
   ana << "... done.\n";

   // Opening the file streams for output during simulation - think of calling close_files(); later
   // to finish the job at the end of simulation

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

// logdata[4].open("tbeta.out");                    logdata[4] << "! Beta fuer die Summe der
// praesentierten Epitope.\n";
// logdata[5].open("tomega.out");                   logdata[5] << "! Omega fuer die Summe der
// praesentierten Epitope.\n";
// logdata[6].open("tbpersum.out");                 logdata[6] << "! % Centroblasten an den
// praesentierten Epitopen.\n";
// logdata[7].open("topersum.out");                 logdata[7] << "! % Output-Zellen an den
// praesentierten Epitopen.\n";
// logdata[8].open("tfb.out");                      logdata[8] << "! % der selektierten
// Centroblasten.\n";
// logdata[9].open("tfo.out");                      logdata[9] << "! % der selektierten
// Output-Zellen.\n";

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
   ana << "arupProtein space READY! (YEAAAAH !!)\n";
}
arupSpace::~arupSpace() {
   for (int i = 0; i < (int) poolSequences.size(); ++i) {
      if (poolSequences[i]) { delete poolSequences[i]; }
   }
}
// common with ss.cpp
void arupSpace::close_files() {
   cout << "Close files in shape space ... ";
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) { logdata[i].close(); }
   }
   logmeanaff.close();
   loghighaff.close();
   logdiversity.close();
   log048aff.close();

   ofstream outSeqSpace;
   outSeqSpace.open("seqspace.out");
   outSeqSpace << printarupProteins();
   outSeqSpace.close();
   cout << "done.\n";
}
/// the arupProtein should already have a type (BCR, TCR or AG), but arupProtein only will raise an
// error !!!
long arupSpace::index_adding_arupProtein(arupProtein * toAdd, long int ID_of_mother_arupProtein) {
   if ((int) poolSequences.size()
       != n_Sequences) {
      cerr
         << "ERR: arupSpace::index_adding_arupProtein, n_Sequences management error !!" << endl;
   }
   if (toAdd
       == NULL) {
      cerr << "ERR: arupSpace::index_adding_arupProtein(..., arupProtein* = NULL !!!)"
           << endl;
   }
   if (toAdd->getType()
       == numberTypesArup) {
      cerr
         <<
         "ERR: arupSpace::index_adding_arupProtein(), forbidden to add a arupProtein without type (should be TCR, Antigen or BCR)"
         << endl;
      return -1;
   }

   poolSequences.push_back(toAdd);
   if (toAdd->getType()
       != typeArupAG) {
      toAdd->max_affinity_to_antigens = maxAffinityarupProteinToAntigens(toAdd);
   } else { toAdd->max_affinity_to_antigens = 0.0; }
   indexParentSeq.push_back(ID_of_mother_arupProtein);
   n_Sequences++;

   // NOTE : when an antigen is added, the affinities of all BCRs have to be recomputed...
   //       (not used) and the bins of this new AG has to be filled with all the producers depending
   // on their affinities
   if (toAdd->getType() == typeArupAG) {
      for (int pos = 0; pos < n_Sequences; ++pos) {
         if (poolSequences[pos]->getType() == typeArupBCR) {
            poolSequences[pos]->max_affinity_to_antigens = maxAffinityarupProteinToAntigens(
               poolSequences[pos]);
            // (not used) toAdd->add_producer(arupProtein::seq_affinity(poolSequences[pos], toAdd,
            // R_affinity, maxSizeClusters, typeAffinityFunction),
            // poolSequences[pos]->nb_AB_producers());
         }
      }
   }
   // if you create a new BCR, nothing change because there is no yet cells with this BCR
   // if you create a new TCR, well, nobody cares !
   return n_Sequences - 1;
}
long arupSpace::getMutation(long from_position) {
   if ((from_position < 0)
       || (from_position
           >= n_Sequences)) {
      cerr << "arupProteinSpqce::getMutation(id=" << from_position
           << "), wrong id (max " << n_Sequences << " arupProteins" << endl;
      return -1;
   }
   arupProtein * theSeq = poolSequences[from_position];
   if (theSeq
       == NULL) {
      cerr << "arupSpace::getMutation, poolSequences[" << from_position << "] is NULL ! \n";
      return -1;
   }
   if (theSeq->getType()
       != typeArupBCR) {
      cerr
         << "arupSpace::getMutation, you are trying to mutate something else than a BCR" << endl;
      return -1;
   }
   ArupBCR * newBCR = new ArupBCR((ArupBCR*) theSeq);
   double deltaE = arupProtein::lawMutations->getRandValue();
   newBCR->mutate(deltaE);
   return index_adding_arupProtein(newBCR, from_position);
}
double arupSpace::affinity(long int seqId1, long int seqId2) {
   if ((seqId1 < 0)
       || (seqId2 < 0)) { cerr << "ERR : arupSpace::Affinity, negative indices" << endl; return 0; }
   if ((seqId1 >= n_Sequences)
       || (seqId2
           >= n_Sequences)) {
      cerr << "ERR : arupSpace::Affinity(seqId1, seqId2), an index is out of scope\n";
      return 0;
   }
   if (getSequence(seqId1)->getType()
       != typeArupBCR) {
      cerr
         << "ERR: ArupSpace::affinity(BCR_id, ANTIGEN_id) is not symmetric. Please use BCR first \n";
      return 0;
   }
   if (getSequence(seqId2)->getType()
       != typeArupAG) {
      cerr
         << "ERR: ArupSpace::affinity(BCR_id, ANTIGEN_id) is not symmetric. Please use AG second \n";
      return 0;
   }
   return amplitude * arupProtein::arup_affinity(getSequence(seqId1), getSequence(seqId2));
}
double arupSpace::affinity_norm(long int seqId1, long int seqId2) {
   return affinity(seqId1, seqId2) / (amplitude + 1e-9);
}
double arupSpace::affinity(long int seqId1, long int seqId2, double &tr) {
   double v = affinity(seqId1, seqId2);
   if ((tr < 0)
       || (tr
           >= 1)) {
      cerr << "arupSpace::affinity(..., ..., Threshold=" << tr
           << "), the threshold should be in [0;1]\n";
      return 0.0;
   }
   return (v - tr) / (1 - tr);
}
void arupSpace::correct_average_affinity(cells celltyp, long &pos, double &average) {
   // Never call this routine before the first cell of type celltyp was created in SS!
   if (sum_cell[celltyp] == 1) { average = 0; }
   // cout<<"average vorher="<<average<<"; add aff="<<best_affinity(pos)<<"; ";
   average = ((sum_cell[celltyp] - 1) * average + best_affinity_norm(pos)) / (sum_cell[celltyp]);
   // cout<<"sum="<<sum_cell[celltyp]<<"; threshold="<<average<<"\n";
}
double arupSpace::maxAffinityarupProteinToAntigens(arupProtein * seq) {
   if (!seq) {
      cerr << "arupSpace::maxAffinityarupProteinToAntigens , arupProtein is NULLLL\n";
      return 0;
   }
   double aff = -1e9;
   for (int i = 0; i < n_Antigen; ++i) {
      aff = max(aff, amplitude * arupProtein::arup_affinity(poolSequences[Antigens[i]], seq));
   }
   return aff;
}
double arupSpace::best_affinity(long pos) {
   if ((pos < 0)
       || (pos
           >= n_Sequences)) {
      cerr << "ERR: arupSpace::bestAffinity, invalid index: " << pos
           << ". Total nb of arupProteins :" << n_Sequences << endl;
      return 0;
   }
   if (poolSequences[pos]->max_affinity_to_antigens <= -1) {
      poolSequences[pos]->max_affinity_to_antigens
         = maxAffinityarupProteinToAntigens(poolSequences[pos]);
   }
   return poolSequences[pos]->max_affinity_to_antigens;
}
double arupSpace::best_affinity_norm(long pos) {
   if ((pos < 0)
       || (pos
           >= n_Sequences)) {
      cerr << "ERR: arupSpace::bestAffinityNorm, arupProtein index unknown : "
           << pos << ". Total nb of arupProteins :" << n_Sequences << endl;
      return 0;
   }
   return best_affinity(pos) / (amplitude + 1e-9);
}
int arupSpace::get_nearest_Antigen(long n) {
   double aff = -1;
   int IDclosest = -1;
   for (int i = 0; i < n_Antigen; ++i) {
      double newAff = arupProtein::arup_affinity(poolSequences[Antigens[i]], getSequence(n));
      if (newAff > aff) { aff = newAff; IDclosest = i; }
   }
   if (IDclosest
       < 0) {
      cerr
         <<
      "arupSpace::get_nearest_antigen(), could not find a closest antigen. Note : n_Antigen = "
         << n_Antigen << endl;
   }
   return IDclosest;                     // Might be more consistent to change it to return
                                         // get_antigen(IDclosest) !!!
}
double arupSpace::Abstandquad(long int &seqId1, long int &seqId2) {
   if ((seqId1 < 0)
       || (seqId2 < 0)) { cerr << "ERR : arupSpace::Abstandquad, negative indices" << endl;
                          return 0; }
   if ((seqId1 >= n_Sequences)
       || (seqId2
           >= n_Sequences)) {
      cerr << "ERR : arupSpace::Abstandquad(seqId1, seqId2), an index is out of scope\n";
      return 0;
   }
   return std::pow(arupProtein::continuoushamming(getSequence(seqId1), getSequence(seqId2)), 2.0);
}
long int arupSpace::add_Antigen(string newarupProtein) {
   ArupAntigen * newAG = new ArupAntigen(newarupProtein);
   return add_Antigen(newAG);
}
long int arupSpace::add_Antigen(ArupAntigen * newAG) {
   Antigens.push_back(index_adding_arupProtein((arupProtein*) newAG, -1));
   n_Antigen++;
   return Antigens[n_Antigen - 1];
}
bool arupSpace::add_new_Antigen() { return true; }               // antigen chosen by AffinitySpace
bool arupSpace::add_new_Antigen(long ag_index) { return true; }   // antigen index provided by the
                                                                  // calling routine

long int arupSpace::add_Seeder(string newarupProtein) {
   ArupBCR * newBCR = new ArupBCR(newarupProtein);
   return add_Seeder(newBCR);
}
long int arupSpace::add_Seeder(ArupBCR * newBCR) {
   Seeders.push_back(index_adding_arupProtein((arupProtein*) newBCR, -1));
   n_Seeder++;
   return Seeders[n_Seeder - 1];
}
// ============================================================
// Dateien schreiben
// ============================================================

double arupSpace::mean_affinity(double * affinities) {
   long i;
   double aff,all;
   for (i = 0; i < 15; i++) {
      affinities[i] = 0.;
   }
   for (i = n_Antigen; i < n_Sequences; i++) {
      // Find nearest antigen-type
      aff = best_affinity_norm(i);
      if (aff > 0.3) {
         affinities[3] += double (getSequence(i)->n_cell[sCB]);
         affinities[4] += double (getSequence(i)->n_cell[sCC]);
         affinities[5] += double (getSequence(i)->n_cell[sout]);
      }
      if (aff < 0.4) {
         affinities[6] += double (getSequence(i)->n_cell[sCB]);
         affinities[7] += double (getSequence(i)->n_cell[sCC]);
         affinities[8] += double (getSequence(i)->n_cell[sout]);
      }
      if ((aff >= 0.4) && (aff < 0.8)) {
         affinities[9] += double (getSequence(i)->n_cell[sCB]);
         affinities[10] += double (getSequence(i)->n_cell[sCC]);
         affinities[11] += double (getSequence(i)->n_cell[sout]);
      }
      if (aff >= 0.8) {
         affinities[12] += double (getSequence(i)->n_cell[sCB]);
         affinities[13] += double (getSequence(i)->n_cell[sCC]);
         affinities[14] += double (getSequence(i)->n_cell[sout]);
      }
      // add corresponding affinity with weight #cb and #cc to ccaff and cbaff
      affinities[0] += aff * double (getSequence(i)->n_cell[sCB]);
      affinities[1] += aff * double (getSequence(i)->n_cell[sCC]);
      affinities[2] += aff * double (getSequence(i)->n_cell[sout]);
   }
   // Divide through total cell numbers:
   all = affinities[0] + affinities[1];
   if (sum_cell[sCB] > 0) {
      affinities[0] /= sum_cell[sCB];
      affinities[3] /= sum_cell[sCB];
      affinities[6] /= sum_cell[sCB];
      affinities[9] /= sum_cell[sCB];
      affinities[12] /= sum_cell[sCB];
   } else { affinities[0] = 0.; }
   if (sum_cell[sCC] > 0) {
      affinities[1] /= sum_cell[sCC];
      affinities[4] /= sum_cell[sCC];
      affinities[7] /= sum_cell[sCC];
      affinities[10] /= sum_cell[sCC];
      affinities[13] /= sum_cell[sCC];
   } else { affinities[1] = 0.; }
   if (sum_cell[sout] > 0) {
      // cout<<outs<<" "<<ha_out<<" "<<sum_cell[sout]<<"\n";
      affinities[2] /= sum_cell[sout];
      affinities[5] /= sum_cell[sout];
      affinities[8] /= sum_cell[sout];
      affinities[11] /= sum_cell[sout];
      affinities[14] /= sum_cell[sout];
      // cout<<outs<<" "<<ha_out<<"\n";
   } else { affinities[2] = 0.; }
   if (sum_cell[sCB] + sum_cell[sCC] > 0) {
      all /= (sum_cell[sCB] + sum_cell[sCC]);
   } else { all = 0.; }
   return all;
}
void arupSpace::get_diversity(double time) {
   logdiversity << time << "   ";
   int diversity[number];
   for (int j = 0; j < number; j++) {
      diversity[j] = 0;
   }
   for (long i = 0; i < n_Sequences; i++) {
      if (getSequence(i)->n_cell[sCB] > 0.5) { ++diversity[sCB]; }
      if (getSequence(i)->n_cell[sCC] > 0.5) { ++diversity[sCC]; }
      if (getSequence(i)->n_cell[sout] > 0.5) { ++diversity[sout]; }
      if (getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC] + getSequence(i)->n_cell[sout]
          > 0.5) { ++diversity[total]; }
      if (getSequence(i)->n_cell[soutext] + getSequence(i)->n_cell[soutextproduce]
          > 0.5) { ++diversity[soutext]; }
   }
   logdiversity << diversity[total] << "   "
                << diversity[sCB] << "   "
                << diversity[sCC] << "   "
                << diversity[sout] << "   "
                << diversity[soutext] << "\n";
}
void arupSpace::mean_affinities_ag(double * mean_affinities, int ag_index) {
   double aff, nBCs, nOUT;
   double totalBC[8];
   // intialise affinities
   for (int i = 0; i < 8; i++) {
      mean_affinities[i] = 0.;
      totalBC[i] = 0.;
   }
   // go through all points of the AffinitySpace
   for (long i = 0; i < n_Sequences; i++) {
      // load affinities[] with corresponding values
      /* [0]: mean affinity of CB+CC to Antigen[ag_index]
       *  [1]: as [0] but only including cells with aff>0.1
       *  [2]: as [0] but only including cells with aff>0.2
       *  [3-5]: as [0-2] for accummulated output
       *  [6]: fraction of CB+CC cells with aff>0.3
       *  [7]: fraction of OUT cells with aff>0.3
       */
      // determine the affinity of this point to the antigen denoted by ag_index
      aff = affinity_norm(i, get_Antigen(ag_index));
      // get the total number of BCs at this point
      nBCs = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // get the total number of accumulated output at this point
      nOUT = getSequence(i)->n_cell[sout];
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
   }
   // normalise to get the mean affinities
   for (int i = 0; i < 6; i++) {
      if (totalBC[i] > 0) {
         mean_affinities[i] /= totalBC[i];
      } else {
         mean_affinities[i] = 0.;
      }
   }
   // get the fraction of aff>0.3 for BC
   if (totalBC[0] > 0) {
      mean_affinities[6] = totalBC[6] / totalBC[0];
   } else {
      mean_affinities[6] = 0.;
   }
   // ... and for accumulated output
   if (totalBC[3] > 0) {
      mean_affinities[7] = totalBC[7] / totalBC[3];
   } else {
      mean_affinities[7] = 0.;
   }
}
void arupSpace::to_multiag_files(double time) {
   ofstream
      saffin_ags, saffin_t10, saffin_t20,
      saffin_out_ags, saffin_out_t10, saffin_out_t20,
      shaffin_ag, shaffin_out_ag,
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
   shaffin_out_ag.open("shaffin_out_ags.out", ofstream::app);
   shaffin_out_ag << time << "  ";
   cross_gcr.open("cross_reactivity_gcr.out", ofstream::app);
   cross_gcr << time << "  ";
   cross_out.open("cross_reactivity_out.out", ofstream::app);
   cross_out << time << "  ";
   double fracsum = 0., fracOUTsum = 0.;
   // variables for cross reactivity
   double cross_reactivity[6];
   for (int i = 0; i < 6; i++) {
      cross_reactivity[i] = 0.;
   }
   for (int a = 0; a < n_Antigen; a++) {
      double mean_affinities[8];
      mean_affinities_ag(mean_affinities, a);
      saffin_ags << mean_affinities[0] << "  ";
      saffin_t10 << mean_affinities[1] << "  ";
      saffin_t20 << mean_affinities[2] << "  ";
      saffin_out_ags << mean_affinities[3] << "  ";
      saffin_out_t10 << mean_affinities[4] << "  ";
      saffin_out_t20 << mean_affinities[5] << "  ";
      shaffin_ag << mean_affinities[6] << "  ";
      shaffin_out_ag << mean_affinities[7] << "  ";
      fracsum += mean_affinities[6];
      fracOUTsum += mean_affinities[7];
      // make the sum of all antigens to get cross-reactivity (with thresholds) for BC and out
      for (int i = 0; i < 6; i++) {
         cross_reactivity[i] += mean_affinities[i];
      }
   }
   // add the sum of the fractions, which is a measure of cross-reactivity, to the fraction-file
   shaffin_ag << fracsum << "\n";
   shaffin_ag.close();
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
//MS
void arupSpace::meanSD_out_affinity(double * aff_mean, double * aff_sd) {cerr<<"notcoded";exit(1);}
double arupSpace::best_distance(long pos) {cerr<<"notcoded";exit(1);}

void arupSpace::to_ssfiles(double time) {
   // sum_check();
   for (short int i = 0; i < SSlogs; i++) {
      if ((i != sFDC) && (i != sTcell)) {
         if (i != soutdec) { logdata[i] << time << "     " << sum_cell[i]; }
         if (i == sout) {
            logdata[sout] << "    " << sum_cell[sout] - oldsum_cell[sout];
         } else if (i == soutdec) {
            logdata[sout] << "   " << sum_cell[soutdec] << "   ";
            if (sum_cell[sout] > 0) {
               logdata[sout] << sum_cell[soutdec] / sum_cell[sout] << "\n";
            } else { logdata[sout] << "0\n"; }
         } else if (i == sallapoptosis) {
            logdata[sallapoptosis] << "    " << sum_cell[sallapoptosis]
               - oldsum_cell[sallapoptosis] << "\n";
         } else { logdata[i] << "\n"; }
      }
      oldsum_cell[i] = sum_cell[i];
   }

   double affinities[15];
   double a;
   a = mean_affinity(affinities);
   logmeanaff << time << "   " << a << "   " << affinities[0] << "   "
              << affinities[1] << "   " << affinities[2] << "\n";
   loghighaff << time << "   " << affinities[3] << "   "
              << affinities[4] << "   " << affinities[5] << "\n";
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
/// Philippe 17.05.16 adapting the output functions for affinity histograms

int arupSpace::intN_mutation(arupProtein * x, arupProtein * y) {
   return arupProtein::continuoushamming(x,y);
}
int arupSpace::best_hamming_distance(long n) {
   int distance = 1.0e+09;  // stupid large number
   for (int nAg = 0; nAg < get_n_Antigen(); nAg++) {
      distance = min(distance, intN_mutation(getSequence(n),getSequence(get_Antigen(nAg))));
   }
   return distance;
}
double arupSpace::mean_hamming_distance(long n) {
   int sumdistance = 0;
   for (int nAg = 0; nAg < get_n_Antigen(); nAg++) {
      sumdistance += intN_mutation(getSequence(n),getSequence(get_Antigen(nAg)));
   }
   if (get_n_Antigen() > 0) {
      return double (sumdistance) / double (get_n_Antigen());
   }
   return 0;
}
double arupSpace::N_mutation(arupProtein * k, arupProtein * l) {
   double hamm = arupProtein::continuoushamming(k, l);
   return hamm * hamm;
}
double arupSpace::mean_affinity_norm(long pos) {
   double tmp = 0.;
   for (int j = 0; j < get_n_Antigen(); j++) {
      tmp += affinity_norm(pos, get_Antigen(j));
   }
   tmp /= double (get_n_Antigen());
   return tmp;
}
void arupSpace::write_gcbc_hamming(double time) {
   int max_hamming = arupProtein::arup_length_sequences;
   vector<int> nbc(max_hamming + 1, 0);          // histogram counters for each Hamming distance,
   vector<int> nbc_mean_ag(max_hamming + 1, 0);  // .. and for the mean Hamming distance for all
                                                 // antigens/epitopes

   // remove this line in the general case:
   /// Philippe : don't understand this warning
   cout << "WARNING! Modified wrong code in SS::best_hamming_distance(long&)\n";

   int nbchere = 0;
   int hamming = 0;
   double hamming_mean = 0.;
   for (long i = 0; i < n_Sequences; i++) {
      // get the total number of GC-BCs at this point
      nbchere = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // determine the Hamming distance of this point to the nearest antigen
      hamming = best_hamming_distance(i);
      hamming_mean = mean_hamming_distance(i);
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      if (hamming >= max_hamming) {
         nbc[max_hamming] += nbchere;
      } else {
         nbc[hamming] += nbchere;
      }
      // now do rounding of the mean hamming distance of sequence i:
      hamming = int (hamming_mean + 0.5);
      // ... and attribute to the
      if (hamming >= max_hamming) {
         nbc_mean_ag[max_hamming] += nbchere;
      } else {
         nbc_mean_ag[hamming] += nbchere;
      }
   }
   // open the output file (append)
   ofstream gcbc_hamming;
   gcbc_hamming.open("gcbc_hamming.out", ofstream::app);
   for (int n = 0; n <= max_hamming; n++) {
      gcbc_hamming << time / 24. << "   " << time << "   " << n
                   << "   " << nbc[n] << "   " << nbc_mean_ag[n] << "\n";
   }
   gcbc_hamming.close();
}
void arupSpace::write_gcbc_affinity(double time) {
   // +++++++++++++++ OPTION ++++++++++++++++
   double logfactor = 3.0;
   int n_bins = 10;
   double affarray[n_bins];
   affarray[0] = 0;
   affarray[1] = 0.0002 / logfactor;
   for (int b = 2; b < n_bins; b++) {
      affarray[b] = logfactor * affarray[b - 1];
   }
   /* for (int b = 0; b < n_bins; b++) {
    *  cout << affarray[b] << ",";
    *  } cout << "\n";  */
   // +++++++++++ end OPTION ++++++++++++++++

   // double check that the whole range of binding probabilities is covered
   if (affarray[n_bins - 1] * logfactor < 1.0) {
      cerr << "In sequenceSpace::write_gcbc_affinity(double): largest affinity is smaller than 1.\n"
           << "Abort simulation.\n";
      exit(1);
   }
   // define a set of counters for each Hamming distance
   vector<int> nbc(n_bins, 0);
   // .. and for the mean Hamming distance for all antigens/epitopes
   vector<int> nbc_mean_ag(n_bins, 0);
   int nbchere = 0;
   double affbest = 0;
   double affmean = 0.;
   // Go through all points of the shape space and add to the counters:
   for (long i = 0; i < n_Sequences; i++) {
      // get the total number of GC-BCs at this point
      nbchere = getSequence(i)->n_cell[sCB] + getSequence(i)->n_cell[sCC];
      // determine the Hamming distance of this point to the nearest antigen
      affbest = best_affinity_norm(i);
      affmean = mean_affinity_norm(i);
      // Add the number of GC-BC found to the counter of cells with Hamming-distance <hamming>
      int binbest = 0, binmean = 0;
      while (binbest<n_bins&&affbest> affarray[binbest]) {
         ++binbest;
      }   // now binbest is one too far, as affbest>0, binbest>0 as well
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
      gcbc_affinity << time / 24. << "   " << time << "   " << affarray[n] << "   "
                    << nbc[n] << "   " << nbc_mean_ag[n] << "\n";
   }
   gcbc_affinity.close();
}
string arupSpace::printarupProteins(bool showarupProteinsWithNoAliveCells, bool showTree) {
   stringstream ss;
   ss << n_Antigen << " arupProteins of antigens\n";
   for (int i = 0; i < n_Antigen; ++i) {
      ss << "Ag nr " << i << " ,ID=\t" << get_Antigen(i) << "\t("
         << getSequence(get_Antigen(i))->size << "\t" << getSequence(get_Antigen(i))->print()
         << "\n";
   }
   ss << n_Seeder << " arupProteins of seeders\n";
   for (int i = 0; i < n_Seeder; ++i) {
      ss << "Seeder BCR nr " << i << " ,ID=\t" << get_Seeder(i) << "\t(" << getSequence(get_Seeder(
                                                                                           i))->size
         << "\tAff=\t" << getSequence(get_Seeder(i))->max_affinity_to_antigens << "\t"
         << getSequence(get_Seeder(i))->print() << "\n";
   }
   ss << n_TCRs << " arupProteins of TCRs\n";
   for (int i = 0; i < n_TCRs; ++i) {
      ss << "TCR nr " << i << " ,ID=\t" << get_TCR(i) << "\t(" << getSequence(get_TCR(i))->size
         << "\t" << getSequence(get_TCR(i))->print() << "\n";
   }

   vector<vector<long> > tree;
   tree.resize(n_Sequences);
   for (int i = 0; i < n_Sequences; ++i) {
      if (indexParentSeq[i] >= 0) { tree[indexParentSeq[i]].push_back(i); }
   }

   ss << n_Sequences - n_Antigen - n_TCRs << " arupProteins of BCRs (seeders or not)\n";
   if (!showarupProteinsWithNoAliveCells) {
      ss
         << "List of arupProteins with alive cells carrying them:\n";
   }
   for (int i = 0; i < n_Sequences; ++i) {
      if (getSequence(i)->getType() == typeArupBCR) {
         if (showarupProteinsWithNoAliveCells
             || ((getSequence(i)->n_cell[sCB]) + (getSequence(i)->n_cell[sCC])
                 + (getSequence(i)->n_cell[sout]) > 0)) {
            arupProtein * seq = getSequence(i);    // doesn't need the cast actually
            ss << "ID=\t" << i << "\tNbCellsTot(BCO)=\t" << getSequence(i)->n_cell[sCB]
               + getSequence(i)->n_cell[sCC] + getSequence(i)->n_cell[sout] << "\t("
               << getSequence(i)->n_cell[sCB] << "-" << getSequence(i)->n_cell[sCC] << "-"
               << getSequence(i)->n_cell[sout] << ")\t(" << seq->size << "),affinity=\t"
               << best_affinity(i) << "\t" << seq->print() << "\n";
            if (showTree) {
               if (tree[i].size() > 0) {
                  ss << "\tDaughters are :\n";
                  for (int j = 0; j < (int) tree[i].size(); ++j) {
                     ss << "\t-->ID=\t" << tree[i][j] << "\t,nCC=\t"
                        << getSequence(tree[i][j])->n_cell[sCC] << "\taffinity=\t" << best_affinity(
                        tree[i][j]) << "\tImproving=\t" << best_affinity(tree[i][j])
                        - best_affinity(i) << "\n";
                  }
               }
            }
         }
      }
   }

   ss << external_cells.size()
      << " arupProteins IDs where there are external cells (producing or not) \n";
   for (int i = 0; i < (int) external_cells.size(); ++i) {
      ss << "ExtCell\t" << i << " ,ID=\t" << external_cells[i] << "\t";
      if (getSequence(external_cells[i])->getType()
          != typeArupBCR) {
         ss << "ERR: this is of type "
            << getSequence(external_cells[i])->typeInString() << endl;
      } else {
         ArupBCR * bseq = (ArupBCR*) getSequence(external_cells[i]);
         ss << "n_soutext=\t" << bseq->n_cell[soutext] << "\tn_soutextproduce:\t"
            << bseq->n_cell[soutextproduce] << endl;
      }
   }

   ss << ab_producers.size() << " arupProteins IDs where there are producing cells \n";
   for (int i = 0; i < (int) ab_producers.size(); ++i) {
      ss << "ProducerCel\t" << i << " ,ID=\t" << ab_producers[i] << "\t";
      if (getSequence(ab_producers[i])->getType()
          != typeArupBCR) {
         ss << "ERR: this is of type "
            << getSequence(ab_producers[i])->typeInString() << endl;
      } else {
         ArupBCR * bseq = (ArupBCR*) getSequence(ab_producers[i]);
         ss << "n_soutext=\t" << bseq->n_cell[soutext] << "\tn_soutextproduce:\t"
            << bseq->n_cell[soutextproduce] << endl;
      }
   }

   return ss.str();
}
bool comparupProteins(pair<double, arupProtein*> a, pair<double, arupProtein*> b) {
   return a.first > b.first;
}
string arupProtein::testeAffinityFunctions() {
/*
 *   stringstream out;
 *   out << "Testing the properties ot the affinity function for the following parameters : \n";
 *   out << "   ->     L= " << L << "\t(Size of arupProteins)" << endl;
 *   out << "   ->     R= " << R << "\t(specificity parameter)" << endl;
 *   out << "   -> maxCl= " << maxClusters << "\t(cluster size scale)" << endl;
 *   switch(typeAffinityFunction){
 *       case seqAff:{out << "   -> Using standard affinity (Saham's)\n"; break;}
 *       case seqAffNorm:{out << "   -> Using standard affinity normalized by maxCl^r\n"; break;}
 *       case seqAffWindow:{out << "   -> Using the maximum affinity of a sliding window\n"; break;}
 *   }
 *
 *   out << "==== Part 1 : enumerates all (if possible), or a lot of arupProteins and sort them by
 * affinity to get the best ones : ====" << endl;
 *
 #define resolutiondistrib 100
 #define maxarupProteinsToEnumeate 1100000
 *
 *   vector< pair<double, arupProtein*> > store;
 *   arupProtein* ref = new arupProtein(L);
 *   vector<double> distribution;
 *   distribution.resize(resolutiondistrib + 1);
 *
 *   int total = 0;
 *   int maxim = pow(2, L);
 *   if(L > 26) maxim = maxarupProteinsToEnumeate + 1;  // to avoid it to become negative ...
 *   bool enumerateAll = (maxim < maxarupProteinsToEnumeate);
 *   if(enumerateAll){
 *       for(int i = 0; i < maxim; ++i){
 *           arupProtein* a = new arupProtein(L, (long) i);
 *           double affi = arupProtein::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
 *           distribution[(int) (((double) resolutiondistrib)*affi)] += 1.0;     // put into the
 * histogram classes for affinity
 *           store.push_back( pair<double, arupProtein*>(affi,a));
 *           total++;
 *       }
 *   } else {
 *       for(int i = 0; i < maxarupProteinsToEnumeate; ++i){
 *           arupProtein* a = new arupProtein(L);
 *           a->randomize();
 *           double affi = arupProtein::seq_affinity(a, ref, R, maxClusters, typeAffinityFunction);
 *           distribution[(int) (((double) resolutiondistrib)*affi)] += 1.0;
 *           store.push_back( pair<double, arupProtein*>(affi,a));
 *           total++;
 *       }
 *   }
 *
 *   out << "Distribution of affinities\n";
 *   for(int i = 0; i < (int) distribution.size(); ++i){
 *       distribution[i] /= (double) total;
 *       out << i << "\t" << distribution[i] << "\t" << double(i) * (1.0 / (double)
 * resolutiondistrib) << "\t" << double(i + 1) * (1.0 / (double) resolutiondistrib) << endl;
 *   }
 *
 *   out << "\narupProteins and affinity, " << ((enumerateAll) ? " in the order of ID\n" : "
 * randomly generated\n");
 *   for(int i = 0; i < 200; ++i){
 *       out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
 *   }
 *
 *   out << "\narupProteins, sorted from the best, out of the " << min(maxim, (int)
 * maxarupProteinsToEnumeate) << " evaluated arupProteins\n";
 *   std::sort(store.begin(), store.end(), comparupProteins);
 *   for(int i = 0; i < 200; ++i){
 *       out << i << "\t" << store[i].second->print() << "\t" << store[i].first << "\n";
 *   }
 *   if(enumerateAll){
 *   out << "\nAffinity of arupProteins taken randomly\n";
 *       for(int i = 0; i < 100; ++i){
 *           arupProtein* seqtmp = new arupProtein(L);
 *           seqtmp->randomize();
 *           out << i << "\t" << arupProtein::seq_affinity(seqtmp, ref, R, maxClusters,
 * typeAffinityFunction) << "\t" << seqtmp->print() << "\n";
 *       }
 *   }
 *   for(int i = 0; i < (int) store.size(); ++i)
 *       delete store[i].second;
 *
 *
 *
 *   out << "==== Part 2 : Evaluating cross-reactivity in the system : ====" << endl;
 *
 *   int nbAntigens = 10;
 *   out << "Generating randomly " << nbAntigens << " antigens " << endl;
 *   vector<arupProtein*> ags;
 *   for(int i = 0; i < nbAntigens; ++i){
 *       arupProtein* seq = new arupProtein(L);
 *       seq->randomize();
 *       ags.push_back(seq);
 *       out << "\tAg nr " << i << "\t" << seq->print() << endl;
 *   }
 *   out << "\nNumber of antigens recognized by randomly generated arupProteins, based on
 * threshold\n";
 *
 *   out << "  -> (for the first 100 arupProteins : ) In the case of random arupProteins" << endl;
 *   total = 0;
 #define thresholdRecoAg 0.1
 *   int nbDiscardedSeq = 0; // arupProteins that don't recognize anything
 *   int countprint = 0;
 *   for(int k = 0; k < min(maxim, (int) maxarupProteinsToEnumeate); ++k){
 *       if(k == 100) out << "  -> (for the remaining arupProteins) for arupProteins recognizing at
 * least an antigen with affinity 0.1" << endl;
 *       total++;
 *
 *       // for each arupProtein,
 *       bool recoAtLeastOne = false;
 *       vector<double> nbRecoDepThresh(10, 0.0);
 *       vector<double> affinityEach(nbAntigens, 0.0);
 *       arupProtein* seqtmp = new arupProtein(L);
 *       seqtmp->randomize();
 *       for(int j = 0; j < nbAntigens; ++j){
 *           double thisAff = arupProtein::seq_affinity(seqtmp, ags[j], R, maxClusters,
 * typeAffinityFunction);
 *           if((thisAff > thresholdRecoAg) || (k < 100)) recoAtLeastOne = true; else nbDiscardedSeq
 * ++;
 *           affinityEach[j] = thisAff;
 *           for(int i = 0; i <= (int) (9.99 * thisAff); ++i){
 *               if(i < 10) nbRecoDepThresh[i] ++;
 *           }
 *       }
 *       if(recoAtLeastOne && (countprint < 5000)){
 *           countprint++;
 *           out << "RandSeq " << k << ", " << seqtmp->print() << " ";
 *           out << "nbAgPerThreshold:";
 *           for(int i = 0; i < 10; ++i){
 *               out << "\t" << nbRecoDepThresh[i];
 *           }
 *           out << "\taffPerAg:";
 *           for(int i = 0; i < nbAntigens; ++i){
 *               out << "\t" << affinityEach[i];
 *           }
 *           out << endl;
 *       }
 *       delete seqtmp;
 *   }
 *   out << "   ... Nb of arupProteins analyzed: " << total << endl;
 *   out << "   ... Nb of arupProteins discarded: " << nbDiscardedSeq << "(except the 100 first
 * ones, i.e. among the :" << total - 100 << " remaining)"<< endl;
 *
 *
 *
 *   out << "==== Part 3 : Evaluating the effect of mutations : ====" << endl;
 *
 *   arupProtein* start = new arupProtein(L);   // starting by '0000' : the best arupProtein
 *   out << "NbMut\tarupProtein\taffinity\n";
 *   for(int i = 0; i < 2*L; ++i){
 *       out << i << "\t" << start->print() << "\t" << arupProtein::seq_affinity(start, ref, R,
 * maxClusters, typeAffinityFunction) << endl;
 *       start->mutateOnePosition();
 *   }
 *
 *   out << "\tReaching a good affinity\t" << endl;
 *   start->randomize();
 *   double prevaff = arupProtein::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
 *
 *   bool stop = false;
 *   for(int i = 0; (i < L) && (!stop); ++i){
 *       out << "arupProtein : " << start->print() << "\tAff:\t" << prevaff << "\t";
 *       out << "PossibleMut:";
 *       vector<int> posGoodMutations;
 *       for(int i = 0; i < L; ++i){
 *           arupProtein stmp = arupProtein(start);
 *           stmp.content[i] = !stmp.content[i];
 *           double newaff = arupProtein::seq_affinity(&stmp, ref, R, maxClusters,
 * typeAffinityFunction);
 *           out << "\t" << newaff;
 *           if(newaff > prevaff) posGoodMutations.push_back(i);
 *       }
 *       out << endl;
 *       if(posGoodMutations.size() > 0){
 *           int nextmut = irandom(posGoodMutations.size()-1);
 *           start->content[posGoodMutations[nextmut]] = !
 * start->content[posGoodMutations[nextmut]];
 *           prevaff = arupProtein::seq_affinity(start, ref, R, maxClusters, typeAffinityFunction);
 *       } else {
 *           stop = true;
 *       }
 *   }
 *
 *   return out.str();*/
   return string("");
}
void arupSpace::testarupSpace() { }
/*
 *   //use cerr to catch errors and avoid delays between cout/cerr
 #define OutputTest cerr
 *   OutputTest << "============ Testing arupProteins ... ============\n" << endl;
 *   arupProtein* a = new arupProtein(10);
 *   OutputTest << "Empty arupProtein  a :" << a->print() << endl;
 *   a->mutateOnePosition();
 *   OutputTest << "One Mutation    a :" << a->print() << endl;
 *   a->mutateOnePosition();
 *   OutputTest << "Another one     a :" << a->print() << endl;
 *   a->randomize();
 *   OutputTest << "Randomized :    a :" << a->print() << endl;
 *   arupProtein* b = new arupProtein(string("01111111111110"));
 *   OutputTest << "New arupProtein    b :" << b->print() << endl;
 *   b->mutate(0.5);
 *   OutputTest << "Mutate 50%/base b:" << b->print() << endl;
 *   arupProtein* c = new arupProtein(b);
 *   OutputTest << "New arupProtein  c=b :" << c->print() << endl;
 *   OutputTest << ((arupProtein::compSeq(b,c)) ? "c equals b" : "problem : c != b") << endl;
 *   OutputTest << "affinity b-c (r=2): " << arupProtein::seq_affinity(b, c, 2.0, -1, seqAff) <<
 * endl;
 *   arupProtein* d = new arupProtein(string("1111100000"));
 *   arupProtein* e = new arupProtein(string("1111011111"));
 *   OutputTest << "affinity d-e (r=3): " << arupProtein::seq_affinity(d, e, 3.0, -1, seqAff) << "
 * between " << d->print() << "\t" << e->print() << endl;
 *   OutputTest << "hamming(d,e) = " << arupProtein::hamming(d,e) << endl;
 *
 *   //Antigen::number_of_bins = 10;
 *   OutputTest << "Getting type of a arupProtein : in enum, " << e->getType() << " and as string: "
 * << e->typeInString() << "\t" << e->print() << endl;
 *   BCR* s1 = new BCR(10);
 *   OutputTest << "Getting type of a BCR      : in enum, " << s1->getType() << " and as string: "
 * << s1->typeInString() << "\t" << s1->print() << endl;
 *   TCR* s2 = new TCR(10);
 *   OutputTest << "Getting type of a TCR      : in enum, " << s2->getType() << " and as string: "
 * << s2->typeInString() << "\t" << s2->print() << endl;
 *   Antigen* s3 = new Antigen(10);
 *   OutputTest << "Getting type of an Antigen : in enum, " << s3->getType() << " and as string: "
 * << s3->typeInString() << "\t" << s3->print() << endl;
 *
 *   s1->add_producer(0.15, 4);
 *   s1->add_producer(0.35, 5);
 *   s1->rem_producer(0.31, 2);
 *   s1->rem_producer(0.0, 1);
 *   s1->add_producer(0.0, 2);
 *   s1->printNbCells();
 *   s1->print();*/

/*
 *
 *   OutputTest << "============ Testing arupSpace ... ============\n" << endl;
 *
 *   ofstream ana("testFile.out");
 *   Parameter par;
 *   par.Value.size_arupProteins = 10;
 *   par.Value.init_antigen_Sequences = 8;
 *   par.Value.initAntigenSeqs.resize(100, string(""));
 *   par.Value.initAntigenSeqs[0] = string("0000000001");
 *   par.Value.initAntigenSeqs[0] = string("-1");
 *   par.Value.max_hamming_antigens = 2;
 *   par.Value.min_hamming_antigens = 1;
 *   par.Value.totalBss = 5; // nb seeder cells
 *   par.Value.initBCRSeqs.resize(100, string(""));
 *   par.Value.initBCRSeqs[0] = string("1111111110");
 *   par.Value.initBCRSeqs[0] = string("-1");
 *   par.Value.max_hamming_BCRs = 2;
 *   par.Value.min_initial_affinity_BCRs = 0.1;
 *   par.Value.max_initial_affinity_BCRs = 1.0;
 *   par.Value.initTCRSeqs.resize(100, string(""));
 *   par.Value.initTCRSeqs[0] = string("1111111110");
 *   par.Value.initTCRSeqs[0] = string("-1");
 *   par.Value.max_hamming_TCRs = 2;
 *   par.Value.min_initial_affinity_TCRs = 0.1;
 *   par.Value.max_initial_affinity_TCRs = 1.0;
 *   par.Value.R_affinity = 2.0;
 *   par.Value.pm_differentiation_time = 0.5;
 *   ana << "hi" << endl;
 *   arupSpace sp(par, ana);
 *
 *
 *   OutputTest << "arupProtein space at initialisation : "<< endl;
 *   OutputTest << sp.printarupProteins() << endl;
 *
 *   arupProtein* newseq = new arupProtein("0001111110");
 *   BCR* aBCR = new BCR(newseq);
 *   BCR* bBCR = new BCR("1100111000");
 *
 *   OutputTest << "when trying to add a arupProtein (without type) to the arupSpace, should raise
 * an error. " << endl;
 *   sp.index_adding_arupProtein(newseq);  // should raise an error
 *
 *   OutputTest << endl;
 *   long newIda = sp.index_adding_arupProtein(aBCR);  // should raise an error
 *   long newIdb = sp.index_adding_arupProtein(bBCR);  // should raise an error
 *   OutputTest << "inserting two BCRs to the space, and got the IDs, and best affinity to antigens
 * (automatically updated :" << endl;
 *   OutputTest << sp.getSequence(newIda)->print() << " with ID " << newIda << "\t" <<
 * sp.best_affinity(newIda) << endl;
 *   OutputTest << sp.getSequence(newIdb)->print() << " with ID " << newIdb << "\t" <<
 * sp.best_affinity(newIdb) << endl;
 *
 *   sp.add_cell(soutext, newIda);
 *   sp.add_cell(soutext, newIda);
 *   sp.add_cell(soutext, newIda);
 *   OutputTest << "After adding three cell to the first BCR arupProtein, now list of cells (n_cell)
 * for this arupProtein:" << endl;
 *   OutputTest << sp.getSequence(newIda)->printNbCells() << endl;
 *   cout << "differentiating the three cells for dt = 0.5 with pm_differentiation_time = 0.5\n";
 *   sp.PM_differentiate(soutext, soutextproduce, 0.5);
 *   OutputTest << "new state of the arupProtein space : " << endl;
 *   OutputTest << sp.printarupProteins(true) << endl;
 *
 *   Antigen* newA = new Antigen(string("01000001010"));
 *   Antigen* newB = new Antigen(string("01111001010"));
 *   long idNewA = sp.add_Antigen(newA);
 *   long idNewB = sp.add_Antigen(newB);
 *   OutputTest << "Now, adding two antigens : " << newA->print() << " ID=" << idNewA << "\t" <<
 * newB->print() << " ID=" << idNewB << endl;
 *   OutputTest << sp.printarupProteins(true) << endl;
 *
 *  }*/
