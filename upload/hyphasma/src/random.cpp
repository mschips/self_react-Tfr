#include "random.h"
#include <math.h>
// #include <iostream.h>

int irandom(int bis) {
  //cerr << "RAND_MAX="<<RAND_MAX << "; ";
  // return int(double(rand())*bis/RAND_MAXp1);
  int result = int (double (rand()) * bis / (double (RAND_MAX) + double (1)));
  //cerr << result << "; ";
  return result;
}
/* old version:
 * int irandom(int bis) {
 * // cout << "RAND_MAX="<<RAND_MAX << "; ";
 * int tmp=rand();
 * // For Linux:
 * double dtmp=double(tmp)*double(bis)/RAND_MAX;
 * // For Windows:
 * // double dtmp=double(bis)*double(tmp)/(RAND_MAX+1);
 * //  cout << "in irandom=" << dtmp << "," << int(dtmp) << "\n";
 *
 * // Sicherheit:
 * //if (int(dtmp)==bis) {
 * //cout << tmp<<"="<<bis<<"=result="<<int(dtmp)<<"\n";
 * //exit(1);
 * //}
 * return int(dtmp);
 * }
 */

long int lrandom(long int bis) {
   // cout << "Not created! See Numerical Recipes in C p.279.\n";
   return bis;
}
double drandom(double bis) {
   // return double(rand())*bis/RAND_MAXp1;
   return double (rand()) * bis / (double (RAND_MAX) + double (1));
}
double drandom() {
   // return double(rand())/RAND_MAXp1;
   return double (rand()) / (double (RAND_MAX) + double (1));
}
void random_sequence(dynarray<long> &m, long int max) {
   // Define a sorted sequence
   // dynarray<long> k(max,0,1);
   long k[max];
   long i,j;
   for (i = 0; i < max; i++) {
      k[i] = i;
   }
   for (i = 0; i < max; i++) {
      // Get random number element of [0,max-i-1]
      // ### Achtung unter WINDOWS!!! Hier gehts nur bis max=32768: lrandom machen
      // Unter LINUX o.k.
      j = irandom(int (max - i));
      // Save k[j] at the next place in array m
      m[i] = k[j];
      // Delete value from k[]
      // k.erase(j);
      // Shift last value to place j
      k[j] = k[max - i - 1];
   }
}
void random2_sequence(long * m, long max) {
   // Define a sorted sequence
   int k[max];
   int i;
   for (i = 0; i < max; i++) {
      k[i] = i;
   }
   while (max > 0) {
      // Get random number element of [0,max-1]
      i = irandom(int (max));
      --max;
      // Save k[i] at the next place in array m from the back
      m[max] = k[i];
      // Shift last value to place i in k
      k[i] = k[max];
   }
}
// =======================================================================
/* Generation of samples from normal distributions (global as used in different context)
 */

double inverse_erf(double x) {
   const double pi = 3.141592654;
   const double a = 8. * (pi - 3) / (3. * pi * (4. - pi));
   double z = 0;
   double kl = (2. / (pi * a) + 0.5 * log(1. - x * x));
   //  cout<<x<<"-->"<<kl<<">"<<sqrt(kl-log(1.-x*x)/a)<<"\n";
   z = sqrt(sqrt(kl * kl - log(1. - x * x) / a) - kl);
   if (x < 0) { z *= -1.; }
   return z;
}
double get_sample_from_normal(const double &mean, const double &width) {
   double xnew,r;
   r = -1;
   while (r < 0) {
      r = drandom();
   }
   double arg = 2. * r - 1.;
   xnew = mean + sqrt(2.) * width * inverse_erf(arg);
   return xnew;
}
double get_positive_sample_from_normal(const double &mean, const double &width) {
   double xnew,r;
   xnew = -1;
   while (xnew <= 0 || xnew >= 2 * mean) {
      // the second condition is needed for keeping the mean
      r = -1;
      while (r < 0) {
         r = drandom();
      }
      double arg = 2. * r - 1.;
      xnew = mean + sqrt(2.) * width * inverse_erf(arg);
   }
   return xnew;
}
// =============================================================
// =============================================================
// ============= Class pre_randomize ===========================
// =============================================================

pre_randomize::pre_randomize() { }
pre_randomize::pre_randomize(short dim2) {
   // Hier wird ein Feld erzeugt, in dem alle Permutation der Zahlen
   // 0 bis 2dim-1 stehen. Man kann dabei an alle moeglichen Sequenzen
   // von Nachbarn denken.

   // Berechne die Zahl der Permutationen
   permut = 1;
   for (long i = 2; i <= dim2; i++) {
      permut *= i;
   }

   // Merke Zahl der Felder pro set
   sl = dim2;

   // Definiere das zugehoerige Feld[0..permut-1][0..sl-1]:
   field = new short*[permut];
   for (long i = 0; i < permut; i++) {
      field[i] = new short[sl];
   }
   // for (long i=0; i<permut; i++) for (short j=0; j<sl; j++) field[i][j]=0;

   // Oeffne ein output-file
   out.open("permutat.out");
   out << "Alle Permutationen von " << sl << " Zahlen:\n";
   out << "Generiert mit pre_randomize in random.h.C\n";

   // Berechne die Permutationen:
   long wdh = permut / sl;
   short newnumber = 0;
   short nextnumber = 1;
   short count = 0;
   // cout<<"sl="<<sl<<"   wdh="<<wdh<<"   permut="<<permut<<"\n";
   for (short m = 0; m < sl; m++) {
      // cout<<"m="<<m<<": \n";
      newnumber = 0;
      // Durchlaufe alle Saetze:
      for (long a = 0; a < permut; a++) {
         // cout<<"  a="<<a<<": ";
         // Checke, ob newnumber schon vorkam -> erhoehen
         nextnumber = 1;
         while (nextnumber == 1) {
            nextnumber = 0;
            for (short j = 0; j < m; j++) {
               if (field[a][j] == newnumber) { nextnumber = 1; }
            }
            if (nextnumber == 1) { add_number(newnumber); }
         }
         // cout<<"write "<<newnumber<<"\n";
         // Ins Feld schreiben
         field[a][m] = newnumber;
         // Zahl der geschriebenen merken
         count++;
         // Nach wdh mal schreiben naechste Zahl
         if (count == wdh) {
            count = 0;
            add_number(newnumber);
         }
      }
      if (wdh > 1) { wdh /= (sl - m - 1); }
   }

   // Schreibe das Ergebnis in das file
   dynarray<short> tmp(sl,0,1);
   for (long i = 0; i < permut; i++) {
      get_set(i,tmp);
      for (short j = 0; j < sl; j++) {
         out << tmp[j] << "  ";
      }
      out << " = Satz " << i << "\n";
   }
   /*out<<"\n\n\nSchreibe nun einige zufaellig raus:\n";
    * for (long i=0; i<5000; i++) {
    * get_rand_set(tmp);
    * for (short j=0; j<sl; j++) out<<tmp[j]<<"  ";
    * out<<" = Versuch "<<i<<"\n";
    * }
    */
   out.close();
}
pre_randomize::~pre_randomize() { }

void pre_randomize::get_set(long n, dynarray<short> &xx) {
   for (short i = 0; i < sl; i++) {
      xx[i] = field[n][i];
   }
}
void pre_randomize::get_rand_set(dynarray<short> &xx) {
   get_set(irandom(permut),xx);
}
short* pre_randomize::get_set(long n) {
   return field[n];
}
short* pre_randomize::get_rand_set() {
   return get_set(irandom(permut));
}
void pre_randomize::add_number(short &n) {
   n++;
   if (n == sl) { n = 0; }
}
// ============================================
// Class Gauss-randomize
// ============================================

gauss_randomize::gauss_randomize() {
   arraydim = 100;
   // Define the array
   field = new double[arraydim];
}
gauss_randomize::gauss_randomize(short dataset) {
  if (dataset == 1) { gauss_initialize(1.05,1.05,181,0,3.141592654); } else if (dataset == 2) {
      cyster_initialize();
   } else { field = new double[10]; }
}
double gauss_randomize::gauss(double &x, double &x0, double &width) {
   return exp(-1. * (x - x0) * (x - x0) / (width * width));
}
void gauss_randomize::gauss_initialize(double x0, double width, int Nx, double xmin, double xmax) {
   /* Values between xmin and xmax at precision (xmax-xmin)/(Nx-1) are saved in an array.
    * [N-1 in the denominator ensures that there are really Nx possible values including
    * both interval limits [xmin,xmax]]
    * The frequency of occurrence corresponds to a Gaussian distribution
    * centered at x0 with width "width".
    * The dimension of array is a result of the minimum ymin of the Gaussian in the
    * considered interval [xmin,xmax]. The vertical resolution is then given by
    * ymin/Ny, where Ny=10 in standard situation.
    */
   // Calculate the x-interval
   double dx = (xmax - xmin) / (double (Nx - 1));
   // Find the minimum of the Gaussian (ymin)
   double ymin = 2.;
   for (int i = 0; i < Nx; i++) {
      double x = xmin + double (i) * dx;
      double y = gauss(x,x0,width);
      if (y < ymin) { ymin = y; }
   }
   // cout<<"ymin="<<ymin<<"\n";
   if (ymin < 0.001) { ymin = 0.001; }
   // Define the vertical resolution (Ny is a local constant double)
   double dy = ymin / Ny;
   // cout<<"dy="<<dy<<"\n";
   // Calculate the total number of necessary entries in the array
   arraydim = 0;
   for (int i = 0; i < Nx; i++) {
      double x = xmin + double (i) * dx;
      arraydim += int (gauss(x,x0,width) / dy + 0.5);
   }
   // cout<<"arraydim="<<arraydim<<"\n";
   // Define the array
   field = new double[arraydim];
   // Fill in the entries
   int ind = 0;
   for (int i = 0; i < Nx; i++) {
      // cout<<"make i="<<i<<" ... ";
      double x = xmin + double (i) * dx;
      // double x=int((xmin+double(i)*dx)*100.0/dx+0.5)*dx/100.0;
      // cout<<" x="<<x<<" ... ";
      int n = int (gauss(x,x0,width) / dy + 0.5);
      // cout<<" n="<<n<<" ...  ind="<<ind<<" ... ";
      for (int j = 0; j < n; j++) {
         field[ind] = x;
         ++ind;
      }
      // cout<<" done. ind="<<ind<<"\n";
   }
   cout << "Gauss-turning-angle-distribution initialised.\n";
}
double gauss_randomize::cyster07angle_wt(int angle) {
   switch (angle) {
      case 5:
         return 23;
         break;

      case 15:
         return 61;
         break;

      case 25:
         return 85;
         break;

      case 35:
         return 97;
         break;

      case 45:
         return 100;
         break;

      case 55:
         return 82;
         break;

      case 65:
         return 79;
         break;

      case 75:
         return 62;
         break;

      case 85:
         return 53;
         break;

      case 95:
         return 44;
         break;

      case 105:
         return 40;
         break;

      case 115:
         return 35;
         break;

      case 125:
         return 30;
         break;

      case 135:
         return 25;
         break;

      case 145:
         return 20;
         break;

      case 155:
         return 17;
         break;

      case 165:
         return 9;
         break;

      case 175:
         return 4;
         break;
   }
   return -1;
}
void gauss_randomize::cyster_initialize() {
   /* Does the same as the Gauss-constructor but using a fixed dataset.
    * Different datasets may be stored here and may be called by the variable dataset.
    */
   // Use Allen 2007 Figure S5 C (green) for the turning angle distribution
   // define the number of x-values
   int Nx = 18;
   // define the minimum x-value
   int xmin = 5;
   // Calculate the x-interval
   int dx = 10;
   // Find the minimum occuring value
   double ymin = 4.;
   // Define the vertical resolution (Ny is a local constant double)
   double dy = ymin / Ny;
   // Calculate the total number of necessary entries in the array
   arraydim = 0;
   for (int i = 0; i < Nx; i++) {
      int x = xmin + i * dx;
      arraydim += int (cyster07angle_wt(x) / dy + 0.5);
   }
   // Define the array
   field = new double[arraydim];
   // Fill in the entries
   int ind = 0;
   for (int i = 0; i < Nx; i++) {
      int x = xmin + i * dx;
      int n = int (cyster07angle_wt(x) / dy + 0.5);
      double xrad = double (x) * 3.141592654 / 180.;
      for (int j = 0; j < n; j++) {
         field[ind] = xrad;
         ++ind;
      }
   }
   cout << "Allen et al. 2007 turning-angle-distribution for WT initialised.\n";
}
gauss_randomize::~gauss_randomize() {
   delete[] field;
}
double gauss_randomize::get_distribution_value() {
   // cout<<"get irandom(arraydim+1) --> ";
   // int ind=irandom(arraydim+1);
   // cout<<ind<<"\n";
   // return field[ind];
   return field[irandom(arraydim + 1)];
}
