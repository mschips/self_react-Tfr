#include "kinetics.h"

GCkinetics::GCkinetics(const double &par_dt) {
   max_GCvolume = 0;
   n_vol_values = 0;
   dt = par_dt;
   for (int a = 0; a < n_vol_max; a++) {
      GCvolume_t[a] = 0;
      t_GCvolume[a] = 0.;
   }
   n_chi2_values = 0;
   chi2_sum = 0.;
}
double GCkinetics::liu91(const int &day) {
   switch (day) {
      case 1:
         return 0.05;
         break;

      case 2:
         return 0.16;
         break;

      case 3:
         return 0.9;
         break;

      case 4:
         return 1.0;
         break;

      case 6:
         return 0.54;
         break;

      case 7:
         return 0.41;
         break;

      case 14:
         return 0.22;
         break;

      case 21:
         return 0.06;
         break;
   }
   return -1.;
}
double GCkinetics::hol92(const int &day) {
   switch (day) {
      case 3:
         return 0.31;
         break;

      case 5:
         return 1.0;
         break;

      case 10:
         return 0.42;
         break;

      case 17:
         return 0.04;
         break;
   }
   return -1.;
}
double GCkinetics::get_chi2(const int &a, const short &liu, const int &day) {
   cout << "get_chi_volume(liu=" << liu << ", day=" << day << ")\n";
   double chi;
   if (liu == 1) {
      chi = (liu91(day) - double (GCvolume_t[a]) / double (max_GCvolume)) / liu91(day);
      // cout<<"CHI! "<<time<<" liu="<<liu91(day)*maxGCvolume<<" vol="<<get_GCvolume()<<"
      // chi="<<chi;
   } else {
      chi = (hol92(day) - double (GCvolume_t[a]) / double (max_GCvolume)) / hol92(day);
      // cout<<"CHI! "<<time<<" hol="<<hol92(day)*maxGCvolume<<" vol="<<get_GCvolume()<<"
      // chi="<<chi;
   }
   chi *= chi;
   // cout<<" chi2="<<chi<<"\n";
   ++n_chi2_values;
   chi2_sum += chi;
   return chi;
}
double GCkinetics::get_chi2(const int &a) {
   // Go through the experimental points:
   double chi2 = 0.;
   if ((t_GCvolume[a] >= 24.) && (t_GCvolume[a] < 24. + dt)) {
      chi2 = get_chi2(a,1,1);
   } else if ((t_GCvolume[a] >= 48.) && (t_GCvolume[a] < 48. + dt)) {
      chi2 = get_chi2(a,1,2);
   } else if ((t_GCvolume[a] >= 72.) && (t_GCvolume[a] < 72. + dt)) {
      chi2 = get_chi2(a,1,3);
      chi2 += get_chi2(a,0,3);
   } else if ((t_GCvolume[a] >= 96.) && (t_GCvolume[a] < 96. + dt)) {
      chi2 = get_chi2(a,1,4);
   } else if ((t_GCvolume[a] >= 120.) && (t_GCvolume[a] < 120. + dt)) {
      chi2 = get_chi2(a,0,5);
   } else if ((t_GCvolume[a] >= 144.) && (t_GCvolume[a] < 144. + dt)) {
      chi2 = get_chi2(a,1,6);
   } else if ((t_GCvolume[a] >= 168.) && (t_GCvolume[a] < 168. + dt)) {
      chi2 = get_chi2(a,1,7);
   } else if ((t_GCvolume[a] >= 240.) && (t_GCvolume[a] < 240. + dt)) {
      chi2 = get_chi2(a,0,10);
   } else if ((t_GCvolume[a] >= 336.) && (t_GCvolume[a] < 336. + dt)) {
      chi2 = get_chi2(a,1,14);
   } else if ((t_GCvolume[a] >= 408.) && (t_GCvolume[a] < 408. + dt)) {
      chi2 = get_chi2(a,0,17);
   } else if ((t_GCvolume[a] > 504. - dt) && (t_GCvolume[a] <= 504.)) {
      chi2 = get_chi2(a,1,21);
   }
   return chi2;
}
void GCkinetics::save_GCvolume(const long &vol, const double &t) {
   cout << "In save_GCvolume(...)...\n";
   cout << "vol=" << vol << "; t=" << t << "\n";
   GCvolume_t[n_vol_values] = vol;
   t_GCvolume[n_vol_values] = t;  // save time in hours
   ++n_vol_values;
}
void GCkinetics::check_GCvolume(const long &vol, const double &t) {
   // Check for a larger maximum GC volume:
   if (vol > max_GCvolume) { max_GCvolume = vol; }
   // Save all those for which an experimental value exists:
   if (((t >= 24.) && (t < 24. + dt))
       || ((t >= 48.) && (t < 48. + dt))
       || ((t >= 72.) && (t < 72. + dt))
       || ((t >= 96.) && (t < 96. + dt))
       || ((t >= 120.) && (t < 120. + dt))
       || ((t >= 144.) && (t < 144. + dt))
       || ((t >= 168.) && (t < 168. + dt))
       || ((t >= 240.) && (t < 240. + dt))
       || ((t >= 336.) && (t < 336. + dt))
       || ((t >= 408.) && (t < 408. + dt))
       || ((t > 504. - dt) && (t <= 504.))
       ) {
      save_GCvolume(vol,t);
   }
}
// Returns the total sigma-squared for the kinetics
double GCkinetics::get_sigma2_GCvolume() {
   ofstream chi_volume;
   chi_volume.open("vol_chi2.out");
   chi_volume << "! time  :  # exp-values  :   chi^2  :  sum_of_chi^2  :  ";
   if (use_sigma_n_1 == 0.) { chi_volume << "sigma_(n)\n"; } else { chi_volume << "sigma_(n-1)\n"; }

   for (int a = 0; a < n_vol_values; a++) {
      double chi2 = get_chi2(a);
      cout << "chi2(" << GCvolume_t[a] << "," << t_GCvolume[a] << ")=" << chi2 << "\n";
      // Write to file
      chi_volume << t_GCvolume[a] << "   " << n_chi2_values << "   " << chi2 << "   "
                 << chi2_sum << "   ";
      if (n_chi2_values - use_sigma_n_1 > 0.) {
         chi_volume << chi2_sum / (double (n_chi2_values) - use_sigma_n_1);
      } else {
         chi_volume << "*";
      }
      chi_volume << "\n";
   }

   chi_volume.close();
   if (n_chi2_values - use_sigma_n_1 > 0.) {
      return chi2_sum / (double (n_chi2_values) - use_sigma_n_1);
   }
   return -1.;
}
