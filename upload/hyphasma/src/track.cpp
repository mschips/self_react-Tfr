#include "track.h"
#include <string.h>
#include <math.h>

TRACK::TRACK() {
   wrote_data = false;
}
void TRACK::init(const double &xr_ext, const double &dt_ext,
                 const short &d, const int * prodim, const int &objno,
                 const double &from, const double &until, const double &t_interval,
                 const int &v_resol, const double &delta_v,
                 const int &alpha_resol, const double &delta_alpha,
                 const int &s_resol, const double &delta_s) {
   if (N_cells > 10) {
      cout << "WARNING! The number of cell-types is larger than 10.\n"
           << "         This might cause problems in the automatically\n"
           << "         generated gle-files because of two needed digits\n"
           << "         in the file names.\n"
           << "WARNING ended.\n";
   }
   x_resolution = xr_ext;
   dt = dt_ext;
   dim = d;
   for (short i = 0; i < 3; i++) {
      if (i < dim) {
         PPDIM[i] = prodim[i];
      } else {
         PPDIM[i] = 0;
      }
   }
   N_OBJECTS = objno;
   max_N_OBJECTS = factor_N_OBJECTS * N_OBJECTS;
   if (objno > 93) {
      cout << "WARNING! The number of tracked objects is larger than 93.\n"
           << "         This might cause problems in the automatically\n"
           << "         generated gle-files.\n"
           << "WARNING ended.\n";
   }
   movement_list = new dynarray<track_data>[max_N_OBJECTS];
   time = -1.0;  // will be provided in every call
   TRACK::TRACKFROM = from;
   TRACK::TRACKUNTIL = until;
   TRACK::DELTA_T = t_interval;
   TRACK::V_RESOLUTION = v_resol;
   TRACK::DELTA_V = delta_v;
   TRACK::ALPHA_RESOLUTION = alpha_resol;
   TRACK::DELTA_ALPHA = delta_alpha;
   TRACK::S_RESOLUTION = s_resol;
   TRACK::DELTA_S = delta_s;

   set_nhess(dim);
}
void TRACK::set_nhess(short &d) {
   // ++++++++++++ OPTION +++++++++++++++++++++++
   // direction of plane for zone transition analysis
   TRACK::nhess[0] = 0.;
   TRACK::nhess[1] = 0.;
   TRACK::nhess[2] = 0.;
   if (dim == 2) { TRACK::nhess[1] = 1.; } else if (dim == 3) {
      TRACK::nhess[2] = 1.;
   }
   // (0,1,0) in 2D and (0,0,1) in 3D.
   // ++++++++ end OPTION +++++++++++++++++++++++
}
void TRACK::set_N_OBJECTS(int n) {
   if (n < max_N_OBJECTS) { N_OBJECTS = n; } else {
      cout << "ERROR: It is not possible to reset the max-number of tracked objects\n"
           << "       to a number larger than the "
           << factor_N_OBJECTS << " times the number of initialised objects.\n"
           << "       If you want to have more entries change factor_N_OBJECTS in track.h.\n";
      exit(1);
   }
}
void TRACK::get_trackraw() {
   ifstream fff("trackraw.out");
   // get the saved number of tracked objects:
   fff >> N_OBJECTS
   >> TRACKFROM
   >> TRACKUNTIL
   >> DELTA_T
   >> x_resolution
   >> dt
   >> PPDIM[0]
   >> PPDIM[1]
   >> PPDIM[2]
   >> dim
   >> V_RESOLUTION
   >> DELTA_V
   >> ALPHA_RESOLUTION
   >> DELTA_ALPHA
   >> S_RESOLUTION
   >> DELTA_S;

   set_nhess(dim);

   // Show the read track-initialising data:
   cout << "Data as read from trackraw.out:\n";
   cout << "N_OBJECTS=" << N_OBJECTS << "\n"
        << "TRACKFROM=" << TRACKFROM << "\n"
        << "TRACKUNTIL=" << TRACKUNTIL << "\n"
        << "DELTA_T=" << DELTA_T << "\n"
        << "V_RESOLUTION=" << V_RESOLUTION << "\n"
        << "DELTA_V=" << DELTA_V << "\n"
        << "ALPHA_RESOLUTION=" << ALPHA_RESOLUTION << "\n"
        << "DELTA_ALPHA=" << DELTA_ALPHA << "\n"
        << "S_RESOLUTION=" << S_RESOLUTION << "\n"
        << "DELTA_S=" << DELTA_S << "\n"
	<< "INCLUDE_INCONTACT=" << INCLUDE_INCONTACT << "\n"
	<< "x_resolution=" << x_resolution << "\n"
        << "dt=" << dt << "\n"
        << "PPDIM[0]=" << PPDIM[0] << "\n"
        << "PPDIM[1]=" << PPDIM[1] << "\n"
        << "PPDIM[2]=" << PPDIM[2] << "\n"
        << "dim=" << dim << "\n";

   movement_list = new dynarray<track_data>[N_OBJECTS];
   // no reserved space as data analysis from file
   // initialise all N_OBJECTS movement-lists:
   for (int xx = 0; xx < N_OBJECTS; xx++) {
      movement_list[xx].delete_all();
   }
   // Get the data until the end of the file:
   int old_a = -1;
   while (fff.eof() == 0) {
      wrote_data = true;
      track_data tmp;
      int a; 
      int s; 
      int act;
      fff >> tmp.t
      >> a                                         // object index
      >> act                                       // type of action (see enum in track.h)
      >> tmp.r[0] >> tmp.r[1] >> tmp.r[2]
      >> tmp.pol[0] >> tmp.pol[1] >> tmp.pol[2]
      >> tmp.elongation >> tmp.l2s_axis
      >> tmp.fdc_contact_time
      >> s;                                        // object status
      if (fff.fail() == false) {
         tmp.status = states(s);
         tmp.action = action_types(act);
         if (a != old_a) { cout << a << "."; old_a = a; }
         movement_list[a].add(tmp);
      }
   }
   fff.close();
   time = -1.0;
   if (N_cells > 10) {
      cout << "WARNING! The number of cell-types is larger than 10.\n"
           << "         This might cause problems in the automatically\n"
           << "         generated gle-files because of two needed digits\n"
           << "         in the file names.\n"
           << "WARNING ended.\n";
   }
   if (N_OBJECTS > 93) {
      cout << "WARNING! The number of tracked objects is larger than 93.\n"
           << "         This might cause problems in the automatically\n"
           << "         generated gle-files.\n"
           << "WARNING ended.\n";
   }
}
void TRACK::get_trackraw(int Npertype) {
   if (Npertype < 0) { get_trackraw(); } else {
      ifstream fff("trackraw.out");
      // get the saved number of tracked objects:
      fff >> N_OBJECTS
      >> TRACKFROM
      >> TRACKUNTIL
      >> DELTA_T
      >> x_resolution
      >> dt
      >> PPDIM[0]
      >> PPDIM[1]
      >> PPDIM[2]
      >> dim
      >> V_RESOLUTION
      >> DELTA_V
      >> ALPHA_RESOLUTION
      >> DELTA_ALPHA
      >> S_RESOLUTION
      >> DELTA_S;

      // Show the read track-initialising data:
      cout << "Data as read from trackraw.out:\n";
      cout << "N_OBJECTS=" << N_OBJECTS << "\n"
           << "TRACKFROM=" << TRACKFROM << "\n"
           << "TRACKUNTIL=" << TRACKUNTIL << "\n"
           << "DELTA_T=" << DELTA_T << "\n"
           << "V_RESOLUTION=" << V_RESOLUTION << "\n"
           << "DELTA_V=" << DELTA_V << "\n"
           << "ALPHA_RESOLUTION=" << ALPHA_RESOLUTION << "\n"
           << "DELTA_ALPHA=" << DELTA_ALPHA << "\n"
           << "S_RESOLUTION=" << S_RESOLUTION << "\n"
           << "DELTA_S=" << DELTA_S << "\n"
	   << "INCLUDE_INCONTACT=" << INCLUDE_INCONTACT << "\n"
           << "x_resolution=" << x_resolution << "\n"
           << "dt=" << dt << "\n"
           << "PPDIM[0]=" << PPDIM[0] << "\n"
           << "PPDIM[1]=" << PPDIM[1] << "\n"
           << "PPDIM[2]=" << PPDIM[2] << "\n"
           << "dim=" << dim << "\n";

      set_nhess(dim);

      movement_list = new dynarray<track_data>[N_OBJECTS];
      // initialise all N_OBJECTS movement-lists:
      for (int xx = 0; xx < N_OBJECTS; xx++) {
         movement_list[xx].delete_all();
      }

      // introduce an array for counting the objects of each type
      int nfilled[N_cells];
      for (int xx = 0; xx < N_cells; xx++) {
         nfilled[xx] = 0;
      }
      int nread[N_cells];
      for (int xx = 0; xx < N_cells; xx++) {
         nread[xx] = 0;
      }

      cout << "Read " << Npertype << " objects per type ...\n";
      N_OBJECTS = 0;
      // Get the data until the end of the file:
      int old_a = -1;   // , old_s=-1;
      bool acceptedobject = true;
      // +++++++++++++++++++ OPTION ++++++++++++++++++++++++++
      // Put a value larger than zero in order to skip the first
      // offset cells of each type.
      int offset = 0;
      // standard value is 0
      if (offset > 0) { cout << "\nWARNING!!! offset>0 in TRACK::get_track_raw(int)!\n"; }
      // ++++++++++++++++END OPTION ++++++++++++++++++++++++++
      while (fff.eof() == 0) {
         wrote_data = true;
         track_data tmp;
         int a;
         int s;
         int act;
         fff >> tmp.t
         >> a
         >> act
         >> tmp.r[0] >> tmp.r[1] >> tmp.r[2]
         >> tmp.pol[0] >> tmp.pol[1] >> tmp.pol[2]
         >> tmp.elongation >> tmp.l2s_axis
         >> tmp.fdc_contact_time
         >> s;
         if (fff.fail() == false) {
            if ((acceptedobject && (a == old_a))
                || ((a != old_a) && (nread[s] >= offset) && (nfilled[s] < Npertype))) {
               if (a != old_a) {
                  cerr << a << ".";
                  old_a = a;
                  //	  old_s=s;
                  acceptedobject = true;
                  ++nfilled[s];
                  ++N_OBJECTS;
               }
               tmp.status = states(s);
               tmp.action = action_types(act);
               movement_list[N_OBJECTS - 1].add(tmp);
            } else {
               if (a != old_a) { ++nread[s]; }
               old_a = a;
               acceptedobject = false;
            }
         }
      }
      fff.close();
      cout << "\nTotal number of tracked objects = " << N_OBJECTS << "\n";
      cout << "Number of tracked cells for each type: ";
      for (int i = 0; i < N_cells; i++) {
         cout << nfilled[i] << ",";
         if ((nfilled[i] != 0) && (nfilled[i] < Npertype)) {
            cout << "\n\nWARNING! The requested number " << Npertype
                 << " of objects per type was not found!\n"
                 << "         Type " << i << " has only " << nfilled[i] << " objects.\n";
         }
      }
      cout << "\n";
      time = -1.0;
      if (N_cells > 10) {
         cout << "WARNING! The number of cell-types is larger than 10.\n"
              << "         This might cause problems in the automatically\n"
              << "         generated gle-files because of two needed digits\n"
              << "         in the file names.\n"
              << "WARNING ended.\n";
      }
      if (N_OBJECTS > 93) {
         cout << "WARNING! The number of tracked objects is larger than 93.\n"
              << "         This might cause problems in the automatically\n"
              << "         generated gle-files.\n"
              << "WARNING ended.\n";
      }
   }
}
TRACK::~TRACK() {
   delete[] movement_list;
}
// statics:
double TRACK::TRACKFROM = 0.0;
double TRACK::TRACKUNTIL = 1.0; // hours
double TRACK::DELTA_T = 0.01666666666666666666666666667; // 1 minute
int TRACK::V_RESOLUTION = 100;
double TRACK::DELTA_V = 2.0; // microns/min
int TRACK::ALPHA_RESOLUTION = 36;
double TRACK::DELTA_ALPHA = 5.0; // microns/min
int TRACK::S_RESOLUTION = 100;
double TRACK::DELTA_S = 0.1; // axis-ratio
double TRACK::nhess[3] = { 0,0,0 };
bool TRACK::INCLUDE_INCONTACT = true;

// =================================================================================

void TRACK::Write_movement(long i, states s, double t, double * r, double * pol,
                           double elong, double l2s, double fct,
                           action_types action_type) {
   track_data tmp;
   tmp.status = s;
   tmp.t = t;
   for (short a = 0; a < dim; a++) {
      tmp.r[a] = r[a];
      tmp.pol[a] = pol[a];
   }
   if (dim == 2) {
      tmp.r[2] = 0;
      tmp.pol[2] = 0;
   }
   tmp.elongation = elong;
   tmp.l2s_axis = l2s;
   tmp.fdc_contact_time = fct;
   tmp.action = action_type;
   movement_list[i].add(tmp);
   wrote_data = true;
}
// =================================================================================

long TRACK::add_new_track(states s, double t, double * r, double * pol,
                          double elong, double l2s) {
   // A new track is added to the list of tracks and initialised with one
   // movement of the object with global index N_OBJECTS with position etc.
   long i = N_OBJECTS;
   set_N_OBJECTS(N_OBJECTS + 1);
   Write_movement(i,s,t,r,pol,elong,l2s,-1,trackini);
   cout << "Dynamically added track number " << N_OBJECTS << ".\n";
   return i;
}
// =================================================================================

// A movement of the object with global index i is inserted into the motility data:
// The new position at time t is (x,y,z).
void TRACK::Write_movement(long i, states s, double t, double * r, double * pol,
                           double elong, double l2s, action_types action_type) {
   Write_movement(i,s,t,r,pol,elong,l2s,-1,action_type);
}
// void TRACK::Write_movement(long i, states s, double t, double*r) {
//  Write_movement(i,s,t,r,1,1);
// }

// =================================================================================

void TRACK::Write_movement(long i, states s, double t, double *r, double * pol,
                           action_types action_type) {
   Write_movement(i,s,t,r,pol,1,1,-1,action_type);
}
// =================================================================================

void TRACK::Stop_movement(long i, states s, double t, double * r, double * pol) {
   track_data tmp;
   tmp.status = s;
   tmp.t = t;
   for (short a = 0; a < dim; a++) {
      tmp.r[a] = r[a];
      tmp.pol[a] = pol[a];
   }
   if (dim == 2) { tmp.r[2] = 0; tmp.pol[2] = 0; }
   tmp.elongation = 0;
   tmp.l2s_axis = 0;
   tmp.fdc_contact_time = -1;
   tmp.action = deathORend;
   movement_list[i].add(tmp);
}
void TRACK::increment_fileno(fileno &tmp) {
   // suffix declared in setparam.h
   int i = 3;
   char weiter = 1;
   while (i >= 0 && weiter == 1) {
      if (tmp[i] != '9') {
         ++tmp[i];
         weiter = 0;
      } else {
         tmp[i] = '0';
         i--;
      }
   }
   if (weiter == 1) {
      cerr << "Too large number in increment_fileno(..) !\n";
   }
}
double TRACK::get_deltar(double * a, double * b) {
   // a and b contain the arrays of indices on the lattice.
   // They correspond to the positions of objects in dim-dimensions.
   // The euklidean distance between both points a and b is returned,
   // where the distance is rescaled with the lattice resolution in microns.
   // Thus, the returned value is also in microns.
   double dr = 0;
   for (short d = 0; d < dim; d++) {
      // get the components for the distance
      dr += ((a[d] - b[d]) * (a[d] - b[d]) * x_resolution * x_resolution);
   }
   return sqrt(dr);
}
double TRACK::get_deltar(const int &i, const int &ab, const int &bis) {
   // Like get_deltar(double*, double*)
   // Calculates the path including all sub-paths.
   // ab contains the index in movement_list[i] of the start-position
   // bis contains the last index in the same list.
   // The summed euclidean path between the position at ab and bis is returned.
   // The distance is rescaled with the lattice resolution in microns.
   // Thus, the returned value is also in microns.
   double dr = 0.;
   double drtmp;
   for (int x = ab; x < bis; x++) {
      drtmp = 0.;
      for (short d = 0; d < dim; d++) {
         // get the components for the distance
         drtmp += ((movement_list[i][x + 1].r[d] - movement_list[i][x].r[d])
                   * (movement_list[i][x + 1].r[d] - movement_list[i][x].r[d])
                   * x_resolution * x_resolution);
      }
      dr += (sqrt(drtmp));
   }
   return dr;
}
void TRACK::Show_tracks() {
   cout << "\nSave raw track-data (y corrected for orientation) "
        << "in trackdata.out and track_raw.out:\n";

   /*
    * ofstream fff("trackdata.out");
    * for (int a=0; a<N_OBJECTS; a++) {
    * fff<<"Object "<<a<<" of type "<<movement_list[a][0].status<<":\n ";
    * for (int b=0; b<movement_list[a].benutzt(); b++)
    *  fff<<"t="<<movement_list[a][b].t
    * <<"   action="<<movement_list[a][b].action
    * <<",  x="<<movement_list[a][b].r[0]
    * <<", y="<<PPDIM[1]-movement_list[a][b].r[1]-1
    * <<", z="<<movement_list[a][b].r[2]
    * <<",  x="<<movement_list[a][b].pol[0]
    * <<", y="<<movement_list[a][b].pol[1]
    * <<", z="<<movement_list[a][b].pol[2]
    * <<",  elonga="<<movement_list[a][b].elongation
    * <<", l2s_axis="<<movement_list[a][b].l2s_axis
    * <<",  fdc_contact_time="<<movement_list[a][b].fdc_contact_time
    * <<";\n";
    * fff<<"\n";
    * }
    * fff.close();
    */

   ofstream ff2("trackraw.out");
   ff2 << N_OBJECTS << "  "
       << TRACKFROM << "  "
       << TRACKUNTIL << "  "
       << DELTA_T << "  "
       << x_resolution << "  "  // external
       << dt << "  "  // external
       << PPDIM[0] << "  "
       << PPDIM[1] << "  "
       << PPDIM[2] << "  "
       << dim << "\n"
       << V_RESOLUTION << "  "
       << DELTA_V << "  "
       << ALPHA_RESOLUTION << "  "
       << DELTA_ALPHA << "  "
       << S_RESOLUTION << "  "
       << DELTA_S << "  "
       << "\n";
   for (int a = 0; a < N_OBJECTS; a++) {
      for (int b = 0; b < movement_list[a].benutzt(); b++) {
         ff2 << movement_list[a][b].t << "  "
             << a << "   "
             << movement_list[a][b].action << "   "
             << movement_list[a][b].r[0] << "  "
             << movement_list[a][b].r[1] << "  "
             << movement_list[a][b].r[2] << "   "
             << movement_list[a][b].pol[0] << "  "
             << movement_list[a][b].pol[1] << "  "
             << movement_list[a][b].pol[2] << "   "
             << movement_list[a][b].elongation << "  "
             << movement_list[a][b].l2s_axis << "   "
             << movement_list[a][b].fdc_contact_time << "  "
             << movement_list[a][b].status << "\n";
      }
   }
   ff2.close();
}
void TRACK::get_color(int i, char color[10]) {
   int coli = i - int (i / 10) * 10;
   switch (coli) {
      case  0:
         strcat(color,"black");
         break;

      case  1:
         strcat(color,"blue");
         break;

      case  2:
         strcat(color,"red");
         break;

      case  3:
         strcat(color,"green");
         break;

      case  4:
         strcat(color,"yellow");
         break;

      case  5:
         strcat(color,"cyan");
         break;

      case  6:
         strcat(color,"magenta");
         break;

      case  7:
         strcat(color,"orange");
         break;

      case  8:
         strcat(color,"brown");
         break;

      case  9:
         strcat(color,"grey");
         break;
   }
}
void TRACK::get_marker(short c, char marker[10]) {
   if (c == 1) {
      strcat(marker,"triangle");       // CB
   } else if (c == 2) {
      strcat(marker,"circle");            // CC
   } else if (c == 4) {
      strcat(marker,"diamond");            // OUT
   } else if (c == 5) {
      strcat(marker,"square");            // TC
   } else { strcat(marker,"star"); }
}
void TRACK::get_cellname(short c, char cc[10]) {
   switch (c) {
      case  0:
         strcat(cc,"nocell");
         break;

      case  1:
         strcat(cc,"CB");
         break;

      case  2:
         strcat(cc,"CC");
         break;

      case  3:
         strcat(cc,"FDC");
         break;

      case  4:
         strcat(cc,"OUT");
         break;

      case  5:
         strcat(cc,"TC");
         break;

      case  6:
         strcat(cc,"blast1");
         break;

      case  7:
         strcat(cc,"blast2");
         break;
   }
}
void TRACK::make_trace_types_a() {
   // cout<<"enter... ";
   fileno typeno = "0000";
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      // Find number of entries for that type
      int nc = 0;
      // cout<<"count "<<c<<"... ";
      for (int i = 0; i < N_OBJECTS; i++) {
         if (movement_list[i][get_first_index(i)].status == c) { ++nc; }
      }
      if (nc > 0) {
         // cout<<"found objects... ";
         char datname[30] = "trace_a";
         strcat(datname,typeno);
         strcat(datname,".gle");
         ofstream x(datname);
         char cellname[10] = "";
         get_cellname(c,cellname);
         x << "size 15 15\n" << "amove 2 2\n" << "set font TEXCMR\n" << "begin graph\n"
           << "  size 12 12\n" << "  fullsize\n"
           << "  xtitle \" ";
         if (dim == 2) { x << "x"; } else if (dim == 3) {
            x << "y";
         }
         x << " [microns]\" hei .7 font TEXCMR\n"
           << "  x2title \" " << cellname << "-tracking in GC \" hei .7 font TEXCMR\n"
           << "  ytitle \" ";
         if (dim == 2) { x << "y"; } else if (dim == 3) {
            x << "z";
         }
         x << " [microns]\" hei .7 font TEXCMR\n"
           << "  xaxis font TEXCMR\n" << "  yaxis font TEXCMR\n"
           << "  xaxis min 0 max ";
         if (dim == 2) { x << int (PPDIM[0] * x_resolution + 0.5); } else if (dim == 3) {
            x << int (PPDIM[1] * x_resolution + 0.5);
         }
         x << "\n" << "  yaxis min 0 max ";
         if (dim == 2) { x << int (PPDIM[1] * x_resolution + 0.5); } else if (dim == 3) {
            x << int (PPDIM[2] * x_resolution + 0.5);
         }
         x << "\n";
         // cout<<"header done. ";

         // put the beginning and ends:
         int d = N_OBJECTS + 1;
         char marker[10] = "";
         get_marker(c,marker);
         short cc = 0;
         if (dim == 2) { cc = 0; } else if (dim == 3) {
            cc = 1;
         }
         // x<<"  data  trace_end000"<<c<<".out d"<<d<<"=c1,c2\n";
         x << "  data  trace_end000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2 + cc
           << "\n";
         x << "  d" << d << " marker f" << marker << " msize 0.3 color black\n";
         ++d;
         // x<<"  data  trace_dead000"<<c<<".out d"<<d<<"=c1,c2\n";
         x << "  data  trace_dead000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2
            + cc << "\n";
         x << "  d" << d << " marker " << marker << " msize 0.3 color red\n";
         ++d;
         x << "  data  trace_begin000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2
            + cc << "\n";
         x << "  d" << d << " marker f" << marker << " msize 0.3 color green\n";
         // cout<<"beg&ends done.";

         // Go through the groups:
         int ngroups = int (N_OBJECTS / MAX_COLUMNS) + 1;
         char group = 'a';
         int data = 1;
         for (int g = 0; g < ngroups; g++) {
            // cout<<"in group "<<g<<"... ";
            int from = g * MAX_COLUMNS;
            int upto = from + MAX_COLUMNS;
            if (upto > N_OBJECTS) { upto = N_OBJECTS; }
            int col = 2;
            for (int i = from; i < upto; i++) {
               // cout<<"do object "<<i<<"... ";
               // if there is an object of type c between from and upto write the read-command
               if (movement_list[i][get_first_index(i)].status == c) {
                  x << "  data  trace_a_" << group << "000" << c << ".out d"
                    << data << "=c" << col + cc << ",c" << col + cc + 1 << "\n";
                  char color[10] = "";
                  get_color(i,color);
                  x << "  d" << data << " lstyle 1 lwidth 0.06 color " << color << "\n";
                  ++data;
                  col += 3;
               }
               // cout<<" done.";
            }
            ++group;
         }
         // cout<<" end of groups. ";
         x << "end graph\n" << "begin key\n" << "  hei .4\n" << "  position tr\n";
         x << "  text \" " << cellname << "-begin \" marker f" << marker
           << " msize 0.5 color green\n";
         x << "  text \" " << cellname << "-end \" marker f" << marker
           << " msize 0.5 color black\n";
         x << "  text \" dead " << cellname << " \" marker " << marker
           << " msize 0.5 color red\n";
         x << "end key\n";
         x.close();
         // cout<<" wrote "<<c<<"\n";
      }
      // cout<<" end of type "<<c<<"\n";
   }
}
void TRACK::make_traces_r() {
   ofstream x("traces_r.gle");
   x << "size 15 15\n" << "amove 2 2\n" << "set font TEXCMR\n" << "begin graph\n"
     << "  size 12 12\n" << "  fullsize\n"
     << "  xtitle \" ";
   if (dim == 2) { x << "x"; } else if (dim == 3) {
      x << "y";
   }
   x << " [microns]\" hei .7 font TEXCMR\n"
     << "  x2title \" cell-traces starting at (0,0) \" hei .7 font TEXCMR\n"
     << "  ytitle \" ";
   if (dim == 2) { x << "y"; } else if (dim == 3) {
      x << "z";
   }
   x << " [microns]\" hei .7 font TEXCMR\n"
     << "  xaxis font TEXCMR\n" << "  yaxis font TEXCMR\n"
     << "  xaxis min -150 max 150\n" << "  yaxis min -150 max 150\n";
   fileno typeno = "0000";
   short cc = 0;
   if (dim == 2) { cc = 0; } else if (dim == 3) {
      cc = 1;
   }
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      if (used_cell_type[c] == true) {
         int d = N_OBJECTS + c;
         char marker[10] = "";
         get_marker(c,marker);
         x << "  data  trace_end000" << c << ".out d" << d << "=c" << 4 + cc << ",c" << 5 + cc
           << "\n";
         x << "  d" << d << " marker f" << marker << " msize 0.4 color black\n";
         d += (N_cells - 1);
         x << "  data  trace_dead000" << c << ".out d" << d << "=c" << 4 + cc << ",c" << 5
            + cc << "\n";
         x << "  d" << d << " marker " << marker << " msize 0.4 color red\n";
      }
   }
   for (int i = 0; i < N_OBJECTS; i++) {
      int whichgroup = int (i / MAX_COLUMNS);
      int col = 3 * (i - whichgroup * MAX_COLUMNS) + 2;
      char color[10] = "";
      char group = 'a';
      for (int bla = 0; bla < whichgroup; bla++) {
         ++group;
      }
      get_color(i,color);
      x << "  data  traces_r_" << group << ".out  d" << i + 1 << "=c" << col + cc << ",c" << col
         + cc + 1 << "\n"
        << "  d" << i + 1 << " lstyle 1 lwidth 0.06 color " << color << "\n";
   }
   x << "end graph\n" << "begin key\n" << "  hei .4\n" << "  position tr\n";
   for (short c = CB; c < N_cells; c++) {
      if (used_cell_type[c] == true) {
         char cellname[10] = "";
         get_cellname(c,cellname);
         char marker[10] = "";
         get_marker(c,marker);
         x << "  text \" " << cellname << "-end \" marker f" << marker
           << " msize 0.5 color black\n";
         x << "  text \" dead " << cellname << " \" marker " << marker
           << " msize 0.5 color red\n";
      }
   }
   x << "end key\n";
   x.close();
}
void TRACK::make_traces_a() {
   ofstream x("traces_a.gle");
   x << "size 15 15\n" << "amove 2 2\n" << "set font TEXCMR\n" << "begin graph\n"
     << "  size 12 12\n" << "  fullsize\n"
     << "  xtitle \" ";
   if (dim == 2) { x << "x"; } else if (dim == 3) {
      x << "y";
   }
   x << " [microns]\" hei .7 font TEXCMR\n"
     << "  x2title \" cell-tracking in the GC \" hei .7 font TEXCMR\n"
     << "  ytitle \" ";
   if (dim == 2) { x << "y"; } else if (dim == 3) {
      x << "z";
   }
   x << " [microns]\" hei .7 font TEXCMR\n"
     << "  xaxis font TEXCMR\n" << "  yaxis font TEXCMR\n"
     << "  xaxis min 0 max 450\n" << "  yaxis min 0 max 450\n";
   fileno typeno = "0000";
   short cc = 0;
   if (dim == 2) { cc = 0; } else if (dim == 3) {
      cc = 1;
   }
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      if (used_cell_type[c] == true) {
         int d = N_OBJECTS + c;
         char marker[10] = "";
         get_marker(c,marker);
         x << "  data  trace_end000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2 + cc
           << "\n";
         x << "  d" << d << " marker f" << marker << " msize 0.3 color black\n";
         d += (N_cells - 1);
         x << "  data  trace_dead000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2
            + cc << "\n";
         x << "  d" << d << " marker " << marker << " msize 0.3 color red\n";
         d += (N_cells - 1);
         x << "  data  trace_begin000" << c << ".out d" << d << "=c" << 1 + cc << ",c" << 2
            + cc << "\n";
         x << "  d" << d << " marker f" << marker << " msize 0.3 color green\n";
      }
   }
   for (int i = 0; i < N_OBJECTS; i++) {
      int whichgroup = int (i / MAX_COLUMNS);
      int col = 3 * (i - whichgroup * MAX_COLUMNS) + 2;
      char color[10] = "";
      char group = 'a';
      for (int bla = 0; bla < whichgroup; bla++) {
         ++group;
      }
      get_color(i,color);
      x << "  data  traces_a_" << group << ".out  d" << i + 1 << "=c" << col + cc << ",c" << col
         + cc + 1 << "\n"
        << "  d" << i + 1 << " lstyle 1 lwidth 0.06 color " << color << "\n";
   }
   x << "end graph\n" << "begin key\n" << "  hei .4\n" << "  position tr\n";
   for (short c = CB; c < N_cells; c++) {
      if (used_cell_type[c] == true) {
         char cellname[10] = "";
         get_cellname(c,cellname);
         char marker[10] = "";
         get_marker(c,marker);
         x << "  text \" " << cellname << "-begin \" marker f" << marker
           << " msize 0.5 color green\n";
         x << "  text \" " << cellname << "-end \" marker f" << marker
           << " msize 0.5 color black\n";
         x << "  text \" dead " << cellname << " \" marker " << marker
           << " msize 0.5 color red\n";
      }
   }
   x << "end key\n";
   x.close();
}
void TRACK::make_reach_dist_x() {
   ofstream x("reach_dist_x.gle");
   x << "size 15 15\n" << "amove 2 2\n" << "set font TEXCMR\n" << "begin graph\n"
     << "  size 12 12\n" << "  fullsize\n"
     << "  xtitle \" sqrt(time) [sqrt(min)]\" hei .7 font TEXCMR\n"
     << "  x2title \" Reached distances of individual cells \" hei .7 font TEXCMR\n"
     << "  ytitle \" Reached distance [microns]\" hei .7 font TEXCMR\n"
     << "  xaxis font TEXCMR\n" << "  yaxis font TEXCMR\n"
     << "  xaxis min 0\n" << "  yaxis min 0\n";
   for (short i = 0; i < N_OBJECTS; i++) {
      if (movement_list[i][get_first_index(i)].status != nocell) {
         char marker[10] = "";
         get_marker(i,marker);
         x << "  data  reach_dist_x.out d" << i + 1 << "=c2,c" << i + 3 << "\n";
         short lstyle = int (i / 10) + 1;
         char color[10] = "";
         get_color(i,color);
         x << "  d" << i + 1 << " lstyle " << lstyle << " lwidth 0.05 color " << color << "\n";
      }
   }
   x << "end graph\n";
   x.close();
}
void TRACK::make_gle() {
   make_traces_r();
   make_traces_a();
   make_trace_types_a();
   make_reach_dist_x();
}
void TRACK::Write_files() {
   Write_files(true);
}
void TRACK::Write_files(bool generate_raw) {
   if (wrote_data) {
      if (generate_raw) { Show_tracks(); } else {
         cout << "\nTracking analysis ...\n";

         // Find the tracked cell types
         for (short c = CB; c < N_cells; c++) {
            used_cell_type[c] = false;
            for (int i = 0; i < N_OBJECTS; i++) {
               int j = get_first_index(i);
               if (j < 0) { cerr << "j negative\n"; }
               if ((j >= 0) && (movement_list[i][j].status == c)) { used_cell_type[c] = true; }
            }
         }

         cout << "  write trace_beginXXXX.out ...\n";
         char name0[20] = "trace_begin";
         Write_trace_begin(name0);

         // Decompose the set of objects into several file in order to limit the number of columns
         int ngroups = int (N_OBJECTS / MAX_COLUMNS) + 1;
         char group[2] = "";

         if (ngroups > 24) {
            cout << "Too many tracked objects for traces! --> Ignore traces.\n";
         } else {
            cout << "  write traces_r.out ...\n";
            group[0] = 'a';
            for (short a = 0; a < ngroups; a++) {
               char name1[20] = "traces_r_";
               strcat(name1,group);
               Write_tracks(name1,a * MAX_COLUMNS,true);
               ++group[0];
            }

            cout << "  write trace_rXXXX.out ...\n";
            group[0] = 'a';
            for (short a = 0; a < ngroups; a++) {
               char name2[20] = "trace_r_";
               strcat(name2,group);
               Write_track_types(name2,a * MAX_COLUMNS,true);
               ++group[0];
            }

            cout << "  write traces_a.out ...\n";
            group[0] = 'a';
            for (short a = 0; a < ngroups; a++) {
               char name3[20] = "traces_a_";
               strcat(name3,group);
               Write_tracks(name3,a * MAX_COLUMNS,false);
               ++group[0];
            }

            cout << "  write trace_aXXXX.out ...\n";
            group[0] = 'a';
            for (short a = 0; a < ngroups; a++) {
               char name4[20] = "trace_a_";
               strcat(name4,group);
               Write_track_types(name4,a * MAX_COLUMNS,false);
               ++group[0];
            }
         }

         cout << "  write trace_endXXXX.out ...\n";
         char name5[20] = "trace_end";
         Write_trace_end(name5);

         cout << "  write trace_deadXXXX.out ...\n";
         char name6[20] = "trace_dead";
         Write_trace_dead(name6);

         cout << "  write transzone.out ...\n";
         for (short a = 0; a < 10; a++) {
            char name7[20] = "transzone";
            switch (a) {
               case 0:
                  break;

               case 1:
                  strcat(name7,"1");
                  break;

               case 2:
                  strcat(name7,"2");
                  break;

               case 3:
                  strcat(name7,"3");
                  break;

               case 4:
                  strcat(name7,"4");
                  break;

               case 5:
                  strcat(name7,"5");
                  break;

               case 6:
                  strcat(name7,"6");
                  break;

               case 7:
                  strcat(name7,"7");
                  break;

               case 8:
                  strcat(name7,"8");
                  break;

               case 9:
                  strcat(name7,"9");
                  break;
            }
            Write_trans_zone(name7,nhess,double (10 * a));
         }

         cout << "  write ncellzone.out ...\n";
         char name8[20] = "ncellzone";
         Write_cellnumber_in_zones(name8,nhess,160.);

         make_gle();

         cout << "  write reach_dist.out ...\n";
         Write_reached_distance();

         cout << "  write turning_angle.out ...\n";
         Write_turning_angle();

         cout << "  write fdccontact.out ...\n";
         Write_fdc_contact_times();

         // ##### Make loops as above here
         cout << "  write v_t.out vmean_t.out vhisto.out ...\n";
         Write_speed();

         cout << "  write axis....out ... \n";
         Write_shape();
         // ##### until here.
         cout << "end of tracking analysis.\n\n";
      }   // end if raw data not written
   } else { cout << "No tracking data saved!\n"; }
}
int TRACK::get_last_index(int obj) {
   // return movement_list[i].benutzt()-1; // former version of it for constant TRACKUNTIL

   int ind = 0;
   while (ind < movement_list[obj].benutzt()
          && movement_list[obj][ind].t <= TRACKUNTIL) {
      ++ind;
   }
   --ind;
   return ind;
}
int TRACK::get_first_index(int obj) {
   int ind = 0;
   while (ind < movement_list[obj].benutzt()
          && movement_list[obj][ind].t < TRACKFROM + DELTA_T * 1.e-08) {
      ++ind;
   }
   --ind;
   if (ind >= movement_list[obj].benutzt()) { --ind; }
   if (ind < 0) { ind = 0; }
   return ind;
}
void TRACK::Write_trace_begin(char prename[20]) {
   fileno typeno = "0000";
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      // Find number of entries for that type
      int nc = 0;
      for (int i = 0; i < N_OBJECTS; i++) {
         if (movement_list[i][get_first_index(i)].status == c) { ++nc; }
      }
      // Start writing only if there are corresponding cells:
      if (nc > 0) {
         char datname[30] = "";
         strcat(datname,prename);
         strcat(datname,typeno);
         strcat(datname,".out");
         ofstream fff(datname);
         fff << "! Cell-trace start-points for " << nc << " cells of type " << c << ":\n"
             << "! Start at t=" << TRACKFROM << " hours\n";
         fff << "! x-   y-   z-position in microns\n";

         for (int i = 0; i < N_OBJECTS; i++) {
            int firsti = get_first_index(i);
            if (movement_list[i][firsti].status == c) {
               for (short d = 0; d < dim; d++) {
                  // if (d!=1)
                  if ((d == 0) || ((dim == 3) && (d == 1))) {
                     fff << movement_list[i][firsti].r[d] * x_resolution << "  ";
                  } else {
                     fff << (PPDIM[d] - movement_list[i][firsti].r[d] - 1) * x_resolution
                         << "  ";
                  }
               }
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
               fff << "\n";
            }
         }
         fff.close();
      }   // if nc>0
   }  // for (cell-types)
}
void TRACK::Write_trace_end(char prename[20]) {
   fileno typeno = "0000";
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      // Find number of entries for that type
      int nc = 0;
      for (int i = 0; i < N_OBJECTS; i++) {
         if (movement_list[i][get_first_index(i)].status == states(c)) { ++nc; }
      }
      // cerr<<"run type "<<c<<" ...\n";
      // Start writing only if there are corresponding cells:
      if (nc > 0) {
         // cerr<<"found "<<nc<<"entries.\n";
         char datname[30] = "";
         strcat(datname,prename);
         strcat(datname,typeno);
         strcat(datname,".out");
         // cerr<<"open file "<<datname<<"\n";
         ofstream fff(datname);
         fff << "! Cell-trace end-points for still living cells of type " << c << ":\n"
             << "! Position at t=" << TRACKUNTIL << " hours\n";
         fff << "! x-   y-   z-position in microns: first 3 absolute, last 3 relative\n";

         for (int i = 0; i < N_OBJECTS; i++) {
            int lasti = get_last_index(i);
            int firsti = get_first_index(i);
            if ((movement_list[i][firsti].status == states(c))
                && (movement_list[i][lasti].status != nocell)) {
               // cerr<<"write object "<<i<<"\n";
               // Write absolute values
               for (short d = 0; d < dim; d++) {
                  // if (d!=1)
                  if ((d == 0) || ((dim == 3) && (d == 1))) {
                     fff << movement_list[i][lasti].r[d] * x_resolution << "  ";
                  } else {
                     fff << (PPDIM[d] - movement_list[i][lasti].r[d] - 1) * x_resolution
                         << "  ";
                  }
               }
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
               fff << "  ";
               // Write relative values
               for (short d = 0; d < dim; d++) {
                  //	    if (d!=1)
                  if ((d == 0) || ((dim == 3) && (d == 1))) {
                     fff << (movement_list[i][lasti].r[d] - movement_list[i][firsti].r[d])
                        * x_resolution << "  ";
                  } else {
                     fff << (movement_list[i][firsti].r[d] - movement_list[i][lasti].r[d])
                        * x_resolution << "  ";
                  }
               }
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
               fff << "\n";
            }
         }
         fff.close();
      }   // if nc>0
   }  // for (cell-types)
}
void TRACK::Write_trace_dead(char prename[20]) {
   fileno typeno = "0000";
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      // Find number of entries for that type
      int nc = 0;
      for (int i = 0; i < N_OBJECTS; i++) {
         if (movement_list[i][get_first_index(i)].status == states(c)) { ++nc; }
      }

      // Start writing only if there are corresponding cells:
      if (nc > 0) {
         char datname[30] = "";
         strcat(datname,prename);
         strcat(datname,typeno);
         strcat(datname,".out");
         ofstream fff(datname);
         fff << "! Cell-trace end-points for dead cells of type " << c << ":\n"
             << "! Position at t=" << TRACKUNTIL << " hours\n";
         fff << "! x-   y-   z-position in microns: first 3 absolute, last 3 relative\n";

         for (int i = 0; i < N_OBJECTS; i++) {
            int lasti = get_last_index(i);
            int firsti = get_first_index(i);
            if ((movement_list[i][firsti].status == states(c))
                && (movement_list[i][lasti].status == nocell)) {
               // Write absolute values
               for (short d = 0; d < dim; d++) {
                  //	    if (d!=1)
                  if ((d == 0) || ((dim == 3) && (d == 1))) {
                     fff << movement_list[i][lasti].r[d] * x_resolution << "  ";
                  } else {
                     fff << (PPDIM[d] - movement_list[i][lasti].r[d] - 1) * x_resolution
                         << "  ";
                  }
               }
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
               fff << "  ";
               // Write relative values
               for (short d = 0; d < dim; d++) {
                  //	    if (d!=1)
                  if ((d == 0) || ((dim == 3) && (d == 1))) {
                     fff << (movement_list[i][lasti].r[d] - movement_list[i][firsti].r[d])
                        * x_resolution << "  ";
                  } else {
                     fff << (movement_list[i][firsti].r[d] - movement_list[i][lasti].r[d])
                        * x_resolution << "  ";
                  }
               }
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
               fff << "\n";
            }
         }
         fff.close();
      }   // if nc>0
   }  // for (cell-types)
}
int TRACK::get_min_a(int obj) {
   // returns the index of the minimum position in direction a
   // where a is y in 2D and z in 3D
   int firsti = get_first_index(obj);
   int lasti = get_last_index(obj);
   double amin = movement_list[obj][firsti].r[dim - 1];
   int ai = firsti;
   for (int a = firsti; a <= lasti; a++) {
      if (movement_list[obj][a].r[dim - 1] < amin) { ai = a; }
   }
   return ai;
}
int TRACK::get_max_a(int obj) {
   // returns the index of the maximum position in direction a
   // where a is y in 2D and z in 3D
   int firsti = get_first_index(obj);
   int lasti = get_last_index(obj);
   double amax = movement_list[obj][firsti].r[dim - 1];
   int ai = firsti;
   for (int a = firsti; a <= lasti; a++) {
      if (movement_list[obj][a].r[dim - 1] > amax) { ai = a; }
   }
   return ai;
}
// written 22.6.2007 mmh
void TRACK::Write_trans_zone(char prename[20], double * hess, double thick) {
   /* Writes a file containing the number of cells that crossed the zone limits.
    * The time is the beginning of the interval.
    * The duration is the tracking duration.
    * (0,0,0) is the lower left front corner of the reaction volume.
    * The zone is defined a normed vector hess defining a line or a plane (2D vs 3D).
    * The position of the border line/plane is determined by
    * a z-value z0 [(0,z0) in 2D and (0,0,z0) in 3D]
    * which is varied to cover the whole reaction volume.
    *
    * The equation of the line/plane is Hesssche-Normalform
    *           \vec{hess}\cdot(\vec{x}-\vec{z0})=0
    * A point \vec{r} is below that line/plane when
    *           (\vec{r}-\vec{z0})\cdot\vec{hess}<0
    * and above if >0.
    *
    * thick impose a certain domain that has to be passed in order to count as transmigration.
    * It is given in micron. Thus in lattics constants the shift of the plane
    * is thick/x_resolution. The y-value defining the plane has to be shifted
    * by thick/n_y in order to get a parallel plane at distance thick.
    * Depending on the initial point being below or above the plane the y-value
    * has to be shifted up or down, respectively, by thick/n_y.
    */

   // At first normalise the hess vector which is provided as direction only:
   double hessnorm = 0.;
   for (short i = 0; i < dim; i++) {
      hessnorm += hess[i] * hess[i];
   }
   if (hessnorm < 1.e-08) { cout << "Error in TRACK::Write_transzone(...)!\n"; exit(1); }
   hessnorm = sqrt(hessnorm);
   for (short i = 0; i < dim; i++) {
      hess[i] /= hessnorm;
   }

   /* Find the range of z0 values that guarantees that the lines/planes with orientation
    * hess cover the whole reaction volume:
    * Calculate the crossing points of the line/plane with the y (2D) or z (3D) axis
    * for the line/plane passing through all corners of the reaction volume and
    * take the minimum value of these z0-values.
    * In 2D the two corner are (0,0) and (xmax,0).
    * In 3D the four corner are (0,0,0), (xmax,0,0), (0,0,zmax) and (xmax,0,zmax).
    * Here xmax and ymax are the dimension of the reaction volume.
    * Solve for the corner \vec{c}=(a,0,b) the equation
    *             \vec{hess}\cdot(\vec{c}-\vec{z_0}) = 0
    * (see readme.tex for details).
    */

   // dimension of the reaction volume:
   int max[dim];
   for (short i = 0; i < dim; i++) {
      max[i] = PPDIM[i] - 1;
   }
   double a0min = 0.;
   double tmp;
   // note that a0min is z0min in 3D and y0min in 2D
   if (dim == 2) {
      tmp = hess[0] * max[0] / hess[1];
      if (tmp < a0min) { a0min = tmp; }
   } else if (dim == 3) {
      tmp = hess[0] * max[0] / hess[2];
      if (tmp < a0min) { a0min = tmp; }
      tmp = hess[1] * max[1] / hess[2];
      if (tmp < a0min) { a0min = tmp; }
      tmp = (hess[0] * max[0] + hess[1] * max[1]) / hess[2];
      if (tmp < a0min) { a0min = tmp; }
      /*
       * tmp=hess[0]*max[0]/hess[1];
       * if (tmp<y0min) y0min=tmp;
       * tmp=hess[2]*max[2]/hess[1];
       * if (tmp<y0min) y0min=tmp;
       * tmp=(hess[0]*max[0]+hess[2]*max[2])/hess[1];
       * if (tmp<y0min) y0min=tmp;
       */
   }
   // Now a0min is the starting value in vertical direction!

   // The number of planes to be considered is fixed here:
   // ++++++++++++ OPTION +++++++++++++++++++++++
   int numplanes = 15;
   // ++++++++ end OPTION +++++++++++++++++++++++
   double plane_steps = double (PPDIM[dim - 1] - 1) / double (numplanes - 1);
   // note that this uses parallel planes in 1-direction for 2D and 2-direction in 3D
   // these can also be not parallel to the 0- and 1- direction in 3D!

   // Open the file
   ofstream fff(strcat(prename,".out"));
   fff << "! Number of zone crossing for each cell type and in each direction.\n"
       << "! y0(micron) active-length time duration(hours)"
       << "  sum-up sum-down sum-all   type1up type1down type1all  "
       << "  type2up type2down type2all ...\n";

   // Go through all planes
   for (int j = 0; j < numplanes; j++) {
      double a0;
      a0 = a0min + double (j) * plane_steps;
      a0 = PPDIM[dim - 1] - a0 - 1;
      double da[3];   // defines thickness for trans-migration
      double va0[3];   // defines position vector of plane
      if (dim == 2) {
         va0[0] = 0.;
         va0[1] = a0;
         va0[2] = 0.;
         // Shift in 1(y)-direction to get a plane at distance thick
         da[0] = 0.;
         if (nhess[1] == 0) { da[1] = 0.; } else {
            da[1] = thick / (2. * x_resolution * nhess[1]);
         }
         da[2] = 0.;
      } else if (dim == 3) {
         va0[0] = 0.;
         va0[1] = 0.;
         va0[2] = a0;
         // Shift in 2(z)-direction to get a plane at distance thick
         da[0] = 0.;
         da[1] = 0.;
         if (nhess[2] == 0) { da[2] = 0.; } else {
            da[2] = thick / (2. * x_resolution * nhess[2]);
         }
      }
      fff << (PPDIM[dim - 1] - a0 - 1) * x_resolution << "  0  " << TRACKFROM << "  "
          << TRACKUNTIL - TRACKFROM << "    ";
      /* #### Note that the length of the path (size of the plane) in the reaction volume
       * is not yet determined! This length also depends on the shape of the reaction
       * volume (square or sphere). One might consider the size of the plane to be crossed
       * to be a suitable weight for the expected number of observed crossings.
       */

      int up[N_cells];
      int down[N_cells];
      int all[N_cells];
      for (int i = 0; i < N_cells; i++) {
         up[i] = 0;
         down[i] = 0;
         all[i] = 0;
      }
      // Go through all objects
      for (int i = 0; i < N_OBJECTS; i++) {
         int firsti = get_first_index(i);
         if (movement_list[i][firsti].status != nocell) {
            // Add this cell to all
            ++all[int (movement_list[i][firsti].status)];

            // get the min and max values for a-direction
            // #### does assume hess in (dim-1)-direction!!!
            int amini = get_min_a(i);
            int amaxi = get_max_a(i);
            double hess_product_min = 0.;
            for (short k = 0; k < dim; k++) {
               hess_product_min += (movement_list[i][amini].r[k] - va0[k] + da[k]) * hess[k];
            }
            double hess_product_max = 0.;
            for (short k = 0; k < dim; k++) {
               hess_product_max += (movement_list[i][amaxi].r[k] - va0[k] - da[k]) * hess[k];
            }
            if ((hess_product_min < 0.) && (hess_product_max >= 0.)) {
               if (amaxi > amini) {
                  ++down[int (movement_list[i][firsti].status)];
               } else { ++up[int (movement_list[i][firsti].status)]; }
            }

            // get the position of this object and determine on which side of the plane it is
            /*
             * double hess_product=0.;
             * for (short k=0; k<dim; k++)
             * hess_product+=(movement_list[i][firsti].r[k]-vy0[k]+dy[k])*hess[k];
             * short below=2;
             * // below is meant in coordinates, thus,
             * // below=true corresponds to above in the reaction volume with inverted y-axis
             * if (hess_product<0.) below=1;
             * else {
             * hess_product=0.;
             * for (short k=0; k<dim; k++)
             *  hess_product+=(movement_list[i][firsti].r[k]-vy0[k]-dy[k])*hess[k];
             * if (hess_product>=0.) below=0;
             * else {
             *  below=2;
             *  --all[int(movement_list[i][firsti].status)];
             * }
             * }
             * // find the last entry
             * int lasti=get_last_index(i);
             * hess_product=0.;
             * if (below==1) {
             * for (short k=0; k<dim; k++)
             *  hess_product+=(movement_list[i][lasti].r[k]-vy0[k]-dy[k])*hess[k];
             * if (hess_product>=0.) ++down[int(movement_list[i][firsti].status)];
             * }
             * else if (below==0) {
             * for (short k=0; k<dim; k++)
             *  hess_product+=(movement_list[i][lasti].r[k]-vy0[k]+dy[k])*hess[k];
             * if (hess_product<0.) ++up[int(movement_list[i][firsti].status)];
             * }
             */
         }
      }   // end for i
          // sum up all cell types
      for (int i = 1; i < N_cells; i++) {
         // position 0 corresponds to state "nocell" and is thus free to be used:
         up[0] += up[i];
         down[0] += down[i];
         all[0] += all[i];
      }
      // Write the result to the file (index 0 contains the sum!)
      for (int i = 0; i < N_cells; i++) {
         fff << up[i] << "  " << down[i] << "  " << all[i] << "    ";
      }
      fff << "\n";
   }

   fff.close();
}
// written 28.2.2010 mmh
void TRACK::Write_cellnumber_in_zones(char prename[20], double * hess, double zpos0) {
   /* Writes a file containing the number of cells on both sides of the
    * plane defined by hess.
    * The zone is defined a normed vector hess defining a line or a plane (2D vs 3D).
    * The position of the border line/plane is determined by
    * a z-value zpos0 [(0,zpos0) in 2D and (0,0,zpos0) in 3D].
    * The equation of the line/plane is Hesssche-Normalform
    *           \vec{hess}\cdot(\vec{x}-\vec{zpos0})=0
    * A point \vec{r} is below that line/plane when
    *           (\vec{r}-\vec{zpos0})\cdot\vec{hess}<0
    * and above if >0.
    *
    * The time is the beginning of the interval.
    * The duration is the tracking duration.
    * (0,0,0) is the lower left front corner of the reaction volume.
    */

   // At first normalise the hess vector which is provided as direction only:
   double hessnorm = 0.;
   for (short i = 0; i < dim; i++) {
      hessnorm += hess[i] * hess[i];
   }
   if (hessnorm < 1.e-08) { cout << "Error in TRACK::Write_transzone(...)!\n"; exit(1); }
   hessnorm = sqrt(hessnorm);
   for (short i = 0; i < dim; i++) {
      hess[i] /= hessnorm;
   }

   /* Find the range of z0 values that guarantees that the lines/planes with orientation
    * hess cover the whole reaction volume:
    * Calculate the crossing points of the line/plane with the y (2D) or z (3D) axis
    * for the line/plane passing through all corners of the reaction volume and
    * take the minimum value of these z0-values.
    * In 2D the two corner are (0,0) and (xmax,0).
    * In 3D the four corner are (0,0,0), (xmax,0,0), (0,0,zmax) and (xmax,0,zmax).
    * Here xmax and ymax are the dimension of the reaction volume.
    * Solve for the corner \vec{c}=(a,0,b) the equation
    *             \vec{hess}\cdot(\vec{c}-\vec{z_0}) = 0
    * (see readme.tex for details).
    */

   // Open the file
   ofstream fff(strcat(prename,".out"));
   fff << "! Number of cells above and below the plane through " << zpos0
       << " microns (from top) defined by hess vector (" << hess[0] << "," << hess[1];
   if (dim == 3) { fff << "," << hess[2]; }
   fff << ").\n"
       << "! time(hours) all-cells cells-above cells-below fraction-above fraction-below\n";

   // rescale position of border line/plane with lattice constant
   zpos0 = long (zpos0 / x_resolution + 0.5) - 1;

   // double a0=PPDIM[dim-1]-zpos0-1; // position from below
   double va0[3];  // defines position vector of plane
   if (dim == 2) {
      va0[0] = 0.;
      va0[1] = zpos0;
      va0[2] = 0.;
   } else if (dim == 3) {
      va0[0] = 0.;
      va0[1] = 0.;
      va0[2] = zpos0;
   }

   int up[N_cells];
   int down[N_cells];
   int all[N_cells];

   // Run this for all times in the time interval:
   double t = TRACKFROM;

   // cerr<<"Go through time steps ...\n";
   // Now successively increase the evaluated time
   while (t <= TRACKUNTIL + DELTA_T) {
      // cerr<<"t="<<t<<" ...\n";
      fff << 60. * (t - TRACKFROM) << "    ";
      for (int i = 0; i < N_cells; i++) {
         up[i] = 0;
         down[i] = 0;
         all[i] = 0;
      }

      int m[N_OBJECTS];
      for (int a = 0; a < N_OBJECTS; a++) {
         m[a] = 0;
      }
      for (int i = 0; i < N_OBJECTS; i++) {
         // cerr<<"Work on i="<<i<<"...\n";
         int firsti = get_first_index(i);
         if (firsti
             < 0) { cerr << "negative get_first_index(" << i << ") in Write_tracks(...)!\n"; }
         if (movement_list[i][firsti].status != nocell) {
            // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
            // cerr<<"i="<<i<<"m[i]="<<m[i]<<"; t="<<t<<"; mov-list[i][m[i]].t="
            // <<movement_list[i][m[i]].t<<"; benutzt="<<movement_list[i].benutzt()<<"\n";
            while (m[i] < movement_list[i].benutzt()
                   && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
               // cerr<<" in while: mov-list["<<i<<"]["<<m[i]<<"].t="
               // <<movement_list[i][m[i]].t<<" compares to t+="<<t+DELTA_T*1.e-08<<"\n";
               ++m[i];
            }
            // Nun steht der Index eins zu weit!
            --m[i];
            /* Same reflection as in Write_tracks(); */
            if (m[i] < 0) { m[i] = 0; }
            // cerr<<"result: i="<<i<<", m[i]="<<m[i]<<"; \n";
            else
            // If the cell still exists do the counting:
            if (movement_list[i][m[i]].status != nocell) {
               ++all[int (movement_list[i][firsti].status)];
               double hess_product = 0.;
               for (short k = 0; k < dim; k++) {
                  hess_product += (movement_list[i][m[i]].r[k] - va0[k]) * hess[k];
               }
               if (hess_product > 0.) {
                  ++down[int (movement_list[i][firsti].status)];
               } else { ++up[int (movement_list[i][firsti].status)]; }
            }
         }    // if movement_list[i]>0
      }   // for i
          // sum up all cell types
      for (int i = 1; i < N_cells; i++) {
         // position 0 corresponds to state "nocell" and is thus free to be used:
         up[0] += up[i];
         down[0] += down[i];
         all[0] += all[i];
      }
      // Write the result to the file (index 0 contains the sum!)
      //  for (int i=0; i<N_cells; i++) fff<<up[i]<<"  "<<down[i]<<"  "<<all[i]<<"    ";
      fff << all[0] << "  " << up[0] << "  " << down[0] << "  ";
      if (all[0]
          > 0) {
         fff << double (up[0]) / double (all[0]) << "  " << double (down[0]) / double (all[0]);
      } else { fff << "0  0"; }
      fff << "\n";
      t += DELTA_T;
   }  // while
   fff.close();
}
void TRACK::Write_tracks(char prename[20], int from, bool relative) {
   /* The tracks with the initial position of every object in the centre:
    * The time between two plotted points is DELTA_T.
    */
   // find out the limits
   int upto = from + MAX_COLUMNS;
   if (upto > N_OBJECTS) { upto = N_OBJECTS; }
   // cerr<<"1:"<<upto<<","<<MAX_COLUMNS<<"\n";

   ofstream fff(strcat(prename,".out"));

   if (relative) { fff << "! Relative"; } else { fff << "! Absolute"; }
   fff << " cell-traces\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   fff << "! time(min)   relative position of the cell, each 3 columns one cell\n";

   // save reference time:
   double t = TRACKFROM;
   // Introduce an index for every tracked cell
   // cerr<<"N_OBJECTS="<<N_OBJECTS<<"\n";
   int m[N_OBJECTS];
   for (int a = 0; a < N_OBJECTS; a++) {
      m[a] = 0;
   }

   // cerr<<"Write initial line ...\n";
   // Write initial line
   fff << "0    ";
   // for (int i=0; i<N_OBJECTS; i++) {
   for (int i = from; i < upto; i++) {
      // cerr<<"Work on object "<<i<<"...\n";
      int firsti = get_first_index(i);
      if (firsti < 0) { cerr << "negative get_first_index(" << i << ") in Write_tracks(...)!\n"; }
      if (movement_list[i][firsti].status != nocell) {
         if (relative) { fff << "0  0  0  "; } else {
            fff << movement_list[i][firsti].r[0] * x_resolution << "  ";
            if (dim == 2) {
               fff << (PPDIM[1] - movement_list[i][firsti].r[1] - 1) * x_resolution << "  ";
            } else if (dim == 3) {
               fff << movement_list[i][firsti].r[1] * x_resolution << "  ";
            } else { fff << "0  "; }
            if (dim == 3) {
               fff << (PPDIM[2] - movement_list[i][firsti].r[2] - 1) * x_resolution << "  ";
            } else { fff << "0  "; }
         }
      } else { fff << "*  *  *  "; }
      fff << "  ";
   }
   fff << "\n";

   // cerr<<"Go through time steps ...\n";
   // Now successively increase the evaluated time
   while (t <= TRACKUNTIL) {
      t += DELTA_T;
      // cerr<<"t="<<t<<" ...\n";
      fff << 60. * (t - TRACKFROM) << "    ";
      // for (int i=0; i<N_OBJECTS; i++) {
      for (int i = from; i < upto; i++) {
         // cerr<<"Work on i="<<i<<"...\n";
         int firsti = get_first_index(i);
         if (firsti
             < 0) { cerr << "negative get_first_index(" << i << ") in Write_tracks(...)!\n"; }
         if (movement_list[i][firsti].status != nocell) {
            // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
            // cerr<<"i="<<i<<"m[i]="<<m[i]<<"; t="<<t<<";
            // mov-list[i][m[i]].t="<<movement_list[i][m[i]].t<<";
            // benutzt="<<movement_list[i].benutzt()<<"\n";
            while (m[i] < movement_list[i].benutzt()
                   && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
               // cerr<<" in while: mov-list["<<i<<"]["<<m[i]<<"].t="<<movement_list[i][m[i]].t<<"
               // compares to t+="<<t+DELTA_T*1.e-08<<"\n";
               ++m[i];
            }
            // Nun steht der Index eins zu weit!
            --m[i];
            /* ### I have reinserted the following line (was outcommented before) because in the
             * case
             * of an object that started tracking later only, the while loop is not entered and thus
             * m[i]=-1 results from the previous line. */
            // cerr<<"result: i="<<i<<", m[i]="<<m[i]<<"; \n";
            if (m[i] < 0) { m[i] = 0; } else {
               // Die Position schreiben
               for (short d = 0; d < dim; d++) {
                  if (relative) {
                     //	    if (d!=1)
                     if ((d == 0) || ((dim == 3) && (d == 1))) {
                        fff
                           << (movement_list[i][m[i]].r[d] - movement_list[i][firsti].r[d])
                           * x_resolution << "  ";
                     } else {
                        fff
                           << (movement_list[i][firsti].r[d] - movement_list[i][m[i]].r[d])
                           * x_resolution << "  ";
                     }
                  } else {
                     //	    if (d!=1)
                     if ((d == 0) || ((dim == 3) && (d == 1))) {
                        fff << movement_list[i][m[i]].r[d] * x_resolution << "  ";
                     } else {
                        fff << (PPDIM[d] - movement_list[i][m[i]].r[d] - 1)
                           * x_resolution << "  ";
                     }
                  }
               }
               // cerr<<"Wrote i="<<i<<"\n";
               for (short d = dim; d < 3; d++) {
                  fff << "0  ";
               }
            }
         }    // if movement_list[i]>0
         else { fff << "*  *  *  "; }
         fff << "  ";
      }   // for i
      fff << "\n";
   }

   fff.close();
}
void TRACK::Write_track_types(char prename[20], int from, bool relative) {
   /* The tracks with the initial position of every object in the centre:
    * The time between two plotted points is DELTA_T.
    *
    * Note that the objects are distributed in different files. The number
    * of objects of each type per file is not constant. It is the number that
    * appears in the set of objects between object from and from+MAX_COLUMNS.
    * This can be zero or MAX_COLUMNS or anything in between.
    ### If later the corresponding gle-file is generated automatically
    ##################this has to enter the generation routine (the analysis has to be repeated)
    ##################in order to read the right file to find the corresponding object. ###
    */

   // find out the limits
   int upto = from + MAX_COLUMNS;
   if (upto > N_OBJECTS) { upto = N_OBJECTS; }

   // Introduce an index for every tracked cell
   int m[N_OBJECTS];
   for (int a = 0; a < N_OBJECTS; a++) {
      m[a] = 0;
   }

   fileno typeno = "0000";
   for (short c = CB; c < N_cells; c++) {
      increment_fileno(typeno);
      // Find number of entries for that type
      int nc = 0;
      for (int i = from; i < upto; i++) {
         if (movement_list[i][get_first_index(i)].status == c) { ++nc; }
      }
      // Start writing only if there are corresponding cells:
      if (nc > 0) {
         char datname[30] = "";
         strcat(datname,prename);
         strcat(datname,typeno);
         strcat(datname,".out");
         ofstream fff(datname);
         if (relative) { fff << "! Relative"; } else { fff << "! Absolute"; }
         fff << " cell-traces for " << nc << " cells of type " << c << ":\n"
             << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
             << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
         fff << "! time(min)   relative position of the cell, each 3 columns one cell\n";

         // starting time:
         double t = TRACKFROM;
         // Write initial line
         fff << 60. * (t - TRACKFROM) << "    ";
         for (int i = from; i < upto; i++) {
            int firsti = get_first_index(i);
            if (movement_list[i][firsti].status == c) {
               if (relative) { fff << "0  0  0    "; } else {
                  for (short d = 0; d < dim; d++) {
                     //	      if (d!=1)
                     if ((d == 0) || ((dim == 3) && (d == 1))) {
                        fff << movement_list[i][firsti].r[d] * x_resolution << "  ";
                     } else {
                        fff << (PPDIM[d] - movement_list[i][firsti].r[d] - 1)
                           * x_resolution << "  ";
                     }
                  }
                  for (short d = dim; d < 3; d++) {
                     fff << "0  ";
                  }
                  fff << "  ";
               }
            }
         }
         fff << "\n";

         // Now successively increase the evaluated time
         while (t <= TRACKUNTIL) {
            t += DELTA_T;
            fff << t - TRACKFROM << "    ";
            for (int i = from; i < upto; i++) {
               int firsti = get_first_index(i);
               if (movement_list[i][firsti].status == c) {
                  // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
                  while (m[i] < movement_list[i].benutzt()
                         && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
                     ++m[i];
                  }
                  // Nun steht der Index eins zu weit!
                  --m[i];
                  /* Same reflection as in Write_tracks(); */
                  if (m[i] < 0) { m[i] = 0; } else {
                     // Die Position schreiben
                     for (short d = 0; d < dim; d++) {
                        if (relative) {
                           //		if (d!=1)
                           if ((d == 0) || ((dim == 3) && (d == 1))) {
                              fff
                                 << (movement_list[i][m[i]].r[d]
                                     - movement_list[i][firsti].r[d]) * x_resolution << "  ";
                           } else {
                              fff
                                 << (movement_list[i][firsti].r[d]
                                     - movement_list[i][m[i]].r[d]) * x_resolution << "  ";
                           }
                        } else {
                           //		if (d!=1)
                           if ((d == 0) || ((dim == 3) && (d == 1))) {
                              fff << movement_list[i][m[i]].r[d] * x_resolution << "  ";
                           } else {
                              fff << (PPDIM[d] - movement_list[i][m[i]].r[d] - 1)
                                 * x_resolution << "  ";
                           }
                        }
                     }
                     for (short d = dim; d < 3; d++) {
                        fff << "0  ";
                     }
                     fff << "  ";
                  }
               }
            }
            fff << "\n";
         }    // end while()
         fff.close();
      }   // end if (nc>0)
   }  // end for (states ..)
}
void TRACK::Write_reached_distance() {
   /* The reached distance of the objects is plotted agains t and sqrt(t):
    * The time between two plotted distances is DELTA_T.
    */
   ofstream fff("reach_dist.out");
   ofstream ffx("reach_dist_x.out");

   fff << "! Reached distance of tracked cells:\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   fff << "! time(min) : sqrt(time) : #cell : mean_delta_r "
       << ": standard_deviation : sqrt(mean_delta_r^2) : standard_deviation\n";
   ffx << "! Reached distance of tracked cells:\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   ffx << "! time(min) sqrt(time) single_cell_delta_r ...\n";

   // startint time:
   double t = TRACKFROM;
   // Introduce an index for every tracked cell
   int m[N_OBJECTS];
   double dr[N_OBJECTS];
   for (int a = 0; a < N_OBJECTS; a++) {
      m[a] = get_first_index(a);
      dr[a] = 0.;
   }

   // Write initial line
   fff << "0.0  0.0   0.0  0.0  0.0   0.0  0.0\n";
   ffx << "0.0  0.0   ";
   for (int i = 0; i < N_OBJECTS; i++) {
      if (movement_list[i][get_first_index(i)].status != nocell) {
         ffx << "0  ";
      } else { ffx << "*  "; }
   }
   ffx << "\n";

   // Now successively increase the evaluated time
   while (t <= TRACKUNTIL) {
      t += DELTA_T;
      // cout<<"t="<<t<<" ...\n";
      fff << 60. * (t - TRACKFROM) << "  " << sqrt(60. * (t - TRACKFROM)) << "    ";
      ffx << 60. * (t - TRACKFROM) << "  " << sqrt(60. * (t - TRACKFROM)) << "    ";
      double obj_pool = 0;
      // Calculate reached distance for all objects at this time
      for (int i = 0; i < N_OBJECTS; i++) {
         int firsti = get_first_index(i);
         if (movement_list[i][firsti].status != nocell) {
            // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
            while (m[i] < movement_list[i].benutzt()
                   && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
               ++m[i];
            }
            // Nun steht der Index eins zu weit!
            --m[i];
            // cout<<"result: i="<<i<<", m[i]="<<m[i]<<"; \n";
            if (m[i] < 0) { m[i] = 0; } else {
               // If the cell is still alive
               if (movement_list[i][m[i]].status != nocell) {
                  // calculate the distance
                  /*dr[i]=0;
                   * for (short d=0; d<dim; d++)
                   * dr[i]+=((movement_list[i][m[i]].r[d]-movement_list[i][0].r[d])
                   *(movement_list[i][m[i]].r[d]-movement_list[i][0].r[d])
                   * x_resolution*x_resolution);
                   * dr[i]=sqrt(dr[i]); */
                  dr[i] = get_deltar(movement_list[i][m[i]].r,movement_list[i][firsti].r);
                  ++obj_pool;
               }
            }
         }    // if movement_list[i]>0
      }   // for i

      // Calculate average and standard deviation
      double meandr = 0., meandr2 = 0.;
      if (obj_pool > 0) {
         for (int i = 0; i < N_OBJECTS; i++) {
            if (movement_list[i][m[i]].status != nocell) {
               meandr += dr[i];
               meandr2 += dr[i] * dr[i];
            }
         }
         meandr /= double (obj_pool);
         meandr2 /= double (obj_pool);
         meandr2 = sqrt(meandr2);
      } else { meandr = -1; meandr2 = -1; }
      double sigma = 0., sigma2 = 0.;
      if (obj_pool > 1) {
         for (int i = 0; i < N_OBJECTS; i++) {
            if (movement_list[i][m[i]].status != nocell) {
               sigma += (dr[i] - meandr) * (dr[i] - meandr);
               sigma2 += (dr[i] - meandr2) * (dr[i] - meandr2);
            }
         }
         sigma /= (obj_pool - 1);
         sigma = sqrt(sigma);
         sigma2 /= (obj_pool - 1);
         sigma2 = sqrt(sigma2);
      } else { sigma = 0; sigma2 = 0; }

      // Write it to the file (one line)
      fff << obj_pool << "    "
          << meandr << "  " << sigma << "    "
          << meandr2 << "  " << sigma2 << "\n";
      for (int i = 0; i < N_OBJECTS; i++) {
         if (movement_list[i][m[i]].status != nocell) {
            ffx << dr[i] << "  ";
         } else { ffx << "*  "; }
      }
      ffx << "\n";
   }

   fff.close();
   ffx.close();
}
double TRACK::get_scalarproduct(const double * a, const double * b) {
   double sum = 0.0;
   for (short d = 0; d < dim; d++) {
      sum += (a[d] * b[d]);
   }
   // This is used for the scalar product of normed vectors only!
   // Values larger than 1 lead to an error in acos()
   if (sum > 1.0) { sum = 1.0; }
   if (sum < -1.0) { sum = -1.0; }
   return sum;
}
double TRACK::get_2norm(const double * k, const double * l) {
   double distance = 0.;
   for (short i = 0; i < dim; i++) {
      double x = k[i] - l[i];
      distance += (x * x);
   }
   return sqrt(distance);
}
void TRACK::Write_turning_angle() {
   long alpha[ALPHA_RESOLUTION + 1];
   for (int a = 0; a < ALPHA_RESOLUTION + 1; a++) {
      alpha[a] = 0;
   }
   long total_turns = 0;
   double mean_angle = 0., mean_histo_angle = 0.;
   double old_pol[3];
   //  double old_pos[3];

   for (int i = 0; i < N_OBJECTS; i++) {
     //      cerr<<"i="<<i<<": ";
      int firsti = get_first_index(i);
      int lasti = get_last_index(i);
      //      cerr << firsti << ","<<lasti<<"; ";
      // if (movement_list[i].benutzt()>0)
      if (lasti > firsti) {
         for (short a = 0; a < dim; a++) {
            old_pol[a] = movement_list[i][firsti].pol[a];
            // old_pos[a]=
	    //            movement_list[i][firsti].r[a];
         }
      }
      // cout<<"Object no "<<i<<" with "<<movement_list[i].benutzt()<<" entries:\n";
      // for (long j=1; j<movement_list[i].benutzt(); j++) {
      for (long j = firsti + 1; j <= lasti; j++) {
	//cerr<<j<<".t="<<movement_list[i][j].t<<": ";
         // Get the angle between subsequent entries
         /*
          * cout<<"old_pol=("<<old_pol[0]<<","<<old_pol[1]<<old_pol[2]<<") ";
          * cout<<"mov_pol=("
          * <<movement_list[i][j].pol[0]<<","
          * <<movement_list[i][j].pol[1]<<","
          * <<movement_list[i][j].pol[2]<<") ";
          */
         // Only include in the statistics when the turning angle is larger than zero.
         // Note that this ignores the rare event that the turning angle was re-chosen
         // but with the same value!
         // if (alph>1.e-03) {
         if (movement_list[i][j].action == polarisation) {
            // cout<<"include!";
            double alph = acos(get_scalarproduct(old_pol,movement_list[i][j].pol));
            // cout<<alph<<" -> ";
            alph *= (180. / 3.141592654);     // transform RAD into DEG
            // cout<<alph<<": ";
            // save alph in the histogram:
            int aindex = int (alph / DELTA_ALPHA);
            if (aindex > ALPHA_RESOLUTION) { aindex = ALPHA_RESOLUTION; }
            ++alpha[aindex];
            // if (aindex==0) cout<<"aindex=0: alpha[0]="<<alpha[aindex]<<"; ";
            ++total_turns;
            mean_angle += alph;
            mean_histo_angle += (double (aindex) + 0.5) * DELTA_ALPHA;
            // save new vars as old
            for (short a = 0; a < dim; a++) {
               old_pol[a] = movement_list[i][j].pol[a];
               // old_pos[a]=
               //movement_list[i][j].r[a];
            }
         }    // end if
              // cout<<"\n";
      }   // end for j
   }  // end for i

   // Average turning angle:
   if (total_turns > 0) {
      mean_angle /= total_turns;
      mean_histo_angle /= total_turns;
   } else { mean_angle = 0; mean_histo_angle = 0; }

   ofstream fff;
   fff.open("turning_angle.out");
   fff << "! Turning angle histogram:\n"
       << "! alpha(grad) : #occurrence : N_of_turns : #perNturns : "
       << "N_OBJECTS : #per_OBJ\n";
   for (int i = 0; i <= ALPHA_RESOLUTION; i++) {
      // fff<<i*DELTA_ALPHA<<"   "
      fff << (double (i) + 0.5) * DELTA_ALPHA << "   "
          << alpha[i] << "    "
          << total_turns << "  "
          << double (alpha[i]) / double (total_turns) << "    "
          << N_OBJECTS << "  "
          << double (alpha[i]) / double (N_OBJECTS)
          << "\n";
   }
   fff.close();

   fff.open("turning_mean.out");
   fff << "! average turning angle of cells in grad\n"
       << "! alphamean : alpha-histo-mean : # turns : 0.0\n";
   fff << mean_angle << "  " << mean_histo_angle << "  " << total_turns << "   0.0\n";
   fff.close();
}
void TRACK::Write_speed() {
   /* The time course of the object speed:
    * Write in N_OBJECT columns (t speed).
    * The speed is calculated on the basis of time intervals of DELTA_T
    * and of the initial and final position in these intervals.
    *
    * Frequency of occurence of objects speeds:
    * Use V_RESOLUTION different speeds with DELTA_V interval.
    * The smallest speed is 0.
    * The speed is calculated on the basis of time intervals of DELTA_T
    * and of the initial and final position in these intervals.
    * In an extra column the result based on real pathlengths is given.
    */

   ofstream fff("v_t.out");
   fff << "! Time course of speed of tracked cells:\n"
       << "! Based on dx observed in the time interval DELTA_T (1).\n"
       << "! Based on the path run in the interval DELTA_T (2).\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   fff << "! time(min) speed of cells (1)  of cells (2)\n";

   ofstream ffx("vmean_t.out");
   ffx << "! Time course of mean speed of all tracked cells:\n"
       << "! Based on dx observed in the time interval DELTA_T (1).\n"
       << "! Based on the path run in the interval DELTA_T (2).\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   ffx << "! time(min) : mean speed of cells (1) : standard deviation :"
       << "! mean speed of cells (2) : standard deviation\n";

   // cerr<<"initialise histo array ...\n";

   long speed_path[V_RESOLUTION + 1];  // Based on pathlength
   long speed_dx[V_RESOLUTION + 1];  // Based on reached distance in interval DELTA_T
   //  long speed_euklid[V_RESOLUTION+1]; // ???
   for (int a = 0; a < V_RESOLUTION + 1; a++) {
      speed_path[a] = 0;
      speed_dx[a] = 0;
      // speed_euklid[a]=0;
   }

   // cerr<<"initialise indices ...\n";

   // save reference time:
   double t = TRACKFROM;
   // Introduce an index for every tracked cell
   int m[N_OBJECTS], mlast[N_OBJECTS];
   bool immobile[N_OBJECTS];
   double dr[N_OBJECTS], drpath[N_OBJECTS];
   for (int a = 0; a < N_OBJECTS; a++) {
      m[a] = get_first_index(a);
      mlast[a] = m[a];
      dr[a] = 0.;
      drpath[a] = 0.;
      immobile[a] = false;
   }
   double lastr[N_OBJECTS][dim];
   for (int i = 0; i < N_OBJECTS; i++) {
      int firsti = m[i];
      if (movement_list[i][firsti].status != nocell) {
         for (short d = 0; d < 3; d++) {
            lastr[i][d] = movement_list[i][firsti].r[d];
         }
      } else {
         for (short d = 0; d < 3; d++) {
            lastr[i][d] = -99999.;
         }
      }
   }

   // cerr<<"initialise single cell histo arrays ...\n";

   //  long speed1_dx[N_OBJECTS][V_RESOLUTION+1];
   long * *speed1_dx;
   speed1_dx = new long*[N_OBJECTS];
   for (int i = 0; i < N_OBJECTS; i++) {
      speed1_dx[i] = new long[V_RESOLUTION + 1];
   }
   //  long speed1_path[N_OBJECTS][V_RESOLUTION+1];
   long * *speed1_path;
   speed1_path = new long*[N_OBJECTS];
   for (int i = 0; i < N_OBJECTS; i++) {
      speed1_path[i] = new long[V_RESOLUTION + 1];
   }
   for (int i = 0; i < N_OBJECTS; i++) {
      for (int j = 0; j <= V_RESOLUTION; j++) {
         speed1_dx[i][j] = 0;
         speed1_path[i][j] = 0;
      }
   }

   double mean_obj_pool = 0.;
   int n_timesteps = 0;
   short switch_immobile_to = 2;
   //  cerr<<"start while ...\n";

   // Now successively increase the evaluated time
   while (t <= TRACKUNTIL - DELTA_T) {
      t += DELTA_T;
      ++n_timesteps;
      // cerr<<"t="<<t<<" ...\n";
      fff << 60. * (t - TRACKFROM) << "  ";
      ffx << 60. * (t - TRACKFROM) << "  ";
      // variables for mean speed
      double obj_pool = 0.;
      double sumdr = 0., sumdrpath = 0.;
      // Calculate speeds for all objects at this time
      for (int i = 0; i < N_OBJECTS; i++) {
         int firsti = get_first_index(i);
         // cerr<<"  work on object "<<i<<" with first-index="<<firsti<<"...\n";
         if (movement_list[i][firsti].status != nocell) {
            // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
            while (m[i] < movement_list[i].benutzt()
                   && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
	      if (INCLUDE_INCONTACT == false) {
		if (movement_list[i][m[i]].action == incontact) { switch_immobile_to = 1; }
		if (movement_list[i][m[i]].action == offcontact) { switch_immobile_to = 0; }
		// MMH2Marta: added this to also set immobile for TFR
                ///Marta2MMH yes, but this is non functioning until we set tracking in cellman right?
                /// (TODO!)
		if (movement_list[i][m[i]].action == inTFRcontact) { switch_immobile_to = 1; }
		if (movement_list[i][m[i]].action == offTFRcontact) { switch_immobile_to = 0; }
	      }
               ++m[i];
            }
            // Nun steht der Index eins zu weit!
            --m[i];
            // cerr<<"result: i="<<i<<", m[i]="<<m[i]<<"; \n";
            if (m[i] < 0) { m[i] = 0; } else {
	      /*
	      if (i == 9 && t >= 146.0 && t < 148.0) {
		cerr<<"t="<<t<<", m["<<i<<"]="<<m[i]<<", action="<<movement_list[i][m[i]].action
		    <<", immobile["<<i<<"]="<<immobile[i]<<":\n";
	      }
	      */
              // If the cell is still alive and has a past
	      if (movement_list[i][m[i]].status != nocell && not(immobile[i])) {
		 /* note that we must not check that the status is movement because
		    the cell has to be counted in order to get the right mean speed.
		 */
                  // cerr<<"  cell is still alive: do analysis ...\n";
                  // calculate the run distance during the last time step DELTA_T.
                  dr[i] = get_deltar(movement_list[i][m[i]].r,lastr[i]);
                  drpath[i] = get_deltar(i,mlast[i],m[i]);

                  // save the new position
                  for (short d = 0; d < dim; d++) {
                     lastr[i][d] = movement_list[i][m[i]].r[d];
                  }
                  mlast[i] = m[i];
                  // add to the sum of all objects;
                  sumdr += dr[i];
                  sumdrpath += drpath[i];
                  // increment the number of treated objects
                  obj_pool += 1.0;
                  // Calculate the speed
                  double v = dr[i] / (60. * DELTA_T);
                  double vpath = drpath[i] / (60. * DELTA_T);
                  // write to file
                  fff << v << "  " << vpath << "  ";       // speed in microns/min
                  // Save this speed in the histogram
                  int vi = int (v / DELTA_V + 0.5);
                  //int vi = int (v / DELTA_V);
                  if (vi > V_RESOLUTION) { vi = V_RESOLUTION; }
                  ++speed_dx[vi];
                  ++speed1_dx[i][vi];
                  int vipath = int (vpath / DELTA_V + 0.5);
                  // int vipath = int (vpath / DELTA_V);
                  if (vipath > V_RESOLUTION) { vipath = V_RESOLUTION; }
                  ++speed_path[vipath];
                  ++speed1_path[i][vipath];
		  /*
		  if (i == 9 && t >= 146.0 && t < 148.0) {
		    cout<<"dr[i]="<<dr[i]<<", v="<<v<<", vi="<<vi<<"\n";
		  }
		  */

               } else { fff << "*  *  "; }
            }
	    if (switch_immobile_to < 2) {
	      immobile[i] = bool(switch_immobile_to);
	      switch_immobile_to = 2;
	    }
         }    // if movement_list[i]>0
         else { fff << "*  *  "; }
      }   // for i
      fff << "\n";

      // cerr<<"get mean speed ...\n";

      // Calculate the mean speed from the summed displacements
      double meanv = -1.0, meanvpath = -1.0;
      mean_obj_pool += obj_pool;
      if (obj_pool > 0) {
         // Note that dead cells are not included in obj_pool.
         // Thus, the average and the sd are based on the number of living cells at each t.

         meanv = sumdr / (60. * DELTA_T * obj_pool);
         meanvpath = sumdrpath / (60. * DELTA_T * obj_pool);

         // Calculate the corresponding standard deviation
         double sigma = 0., sigmapath = 0.;
         if (obj_pool > 1) {
            for (int i = 0; i < N_OBJECTS; i++) {
               if (movement_list[i][m[i]].status != nocell) {
                  sigma += (dr[i] / (60. * DELTA_T) - meanv)
                           * (dr[i] / (60. * DELTA_T) - meanv);
                  sigmapath += (drpath[i] / (60. * DELTA_T) - meanvpath)
                               * (drpath[i] / (60. * DELTA_T) - meanvpath);
               }
            }
            sigma /= (obj_pool - 1);
            sigma = sqrt(sigma);
            sigmapath /= (obj_pool - 1);
            sigmapath = sqrt(sigmapath);
         } else { sigma = 0; sigmapath = 0; }

         // Write to file
         ffx << meanv << "  " << sigma << "    " << meanvpath << "  " << sigmapath << "\n";
      } else { ffx << "*   0   *   0"; }
   }
   fff.close();
   ffx.close();

   // cerr<<"make histogram ...\n";

   // Calculate the average of tracked cell numbers per timestep
   if (n_timesteps > 0) { mean_obj_pool /= double (n_timesteps); }

   long total_dx_measurements = 0, total_path_measurements = 0;
   for (int i = 0; i <= V_RESOLUTION; i++) {
      total_dx_measurements += speed_dx[i];
      total_path_measurements += speed_path[i];
   }
   fff.open("vhisto.out");
   fff << "! Speed histogram:\n"
       << "! (1) based on positions at DELTA_T\n"
       << "! (2) based on real pathlength\n"
       << "! v(micron/min) : #(1) : #(2) : mean_cell# : #percell(1) : #percell(2) : "
       << "N_OBJECTS : #perN(1) : #perN(2) : #perN_sd(1) : #perN_sd(2) : "
       << "total(1)measurements : % (1) occurrence : "
       << "total(2)measurements : % (2) occurrence\n";
   for (int i = 0; i <= V_RESOLUTION; i++) {
      // cerr<<"  write v-step "<<i<<"...\n";
      // fff<<i*DELTA_V<<"   "
      // fff << (double (i) + 0.5) * DELTA_V << "   "
      fff << (double (i)) * DELTA_V << "   "
          << speed_dx[i] << "  "
          << speed_path[i] << "    "
          << mean_obj_pool << "    "
          << double (speed_dx[i]) / mean_obj_pool << "  "
          << double (speed_path[i]) / mean_obj_pool << "    "
          << N_OBJECTS << "    ";
      // Note that the standard deviation is ill-defined if based on mean_obj_pool.
      // Thus, the averages are given with respect to N_OBJECTS as well, irrespective
      // of the number of cells died during the time window of tracking:
      // Get averages of numbers of occurrences of speeds:
      double n_v_percell_dx = double (speed_dx[i]) / double (N_OBJECTS);
      double n_v_percell_path = double (speed_path[i]) / double (N_OBJECTS);
      // Calculate standard deviations:
      double sigma = 0., sigma_path = 0.;
      if (N_OBJECTS > 1) {
         for (int j = 0; j < N_OBJECTS; j++) {
            sigma += ((double (speed1_dx[j][i]) - n_v_percell_dx)
                      * (double (speed1_dx[j][i]) - n_v_percell_dx));
            sigma_path += ((double (speed1_path[j][i]) - n_v_percell_path)
                           * (double (speed1_path[j][i]) - n_v_percell_path));
         }
         sigma /= (N_OBJECTS - 1);
         sigma = sqrt(sigma);
         sigma_path /= (N_OBJECTS - 1);
         sigma_path = sqrt(sigma_path);
      }
      // Write to file:
      fff << n_v_percell_dx << "  "
          << n_v_percell_path << "    "
          << sigma << "  "
          << sigma_path << "   "
          << total_dx_measurements << "  "
          << double (speed_dx[i]) / double (total_dx_measurements) << "   "
          << total_path_measurements << "  "
          << double (speed_path[i]) / double (total_path_measurements)
          << "\n";
   }
   fff.close();

   /*
    * cerr<<"deliberate memory\n";
    * for (int i=0; i<=V_RESOLUTION; i++) {
    * delete [] speed1_dx[i];
    * delete [] speed1_path[i];
    * }
    * delete [] speed1_dx;
    * delete [] speed1_path;
    */

   // cerr<<"write vmean_all.out ...\n";
   fff.open("vmean_all.out");
   fff << "! average speed of cells in micron/minute\n"
       << "! vmean (interval) : # moves : vmean (path in interval) : # moves : 0.0\n";
   int n = 0, npath = 0;
   double vm = 0., vmpath = 0.;
   for (int i = 0; i < V_RESOLUTION; i++) {
      n += speed_dx[i];
      npath += speed_path[i];
      // vm+=speed_dx[i]*double(i)*DELTA_V;
      // vmpath+=speed_path[i]*double(i)*DELTA_V;
      vm += double (speed_dx[i]) * (double (i) + 0.5) * DELTA_V;
      vmpath += double (speed_path[i]) * (double (i) + 0.5) * DELTA_V;
   }
   if (n > 0) { vm /= double (n); }
   if (npath > 0) { vmpath /= double (npath); }
   fff << vm << "  " << n << "    " << vmpath << "  " << npath << "     0.0\n";
   fff.close();
}
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================

void TRACK::Write_fdc_contact_times() {
   double DC = 1.;  // min
   int C_RESOLUTION = int (60. * double (TRACKUNTIL - TRACKFROM) / DC + 0.5);
   long contact[C_RESOLUTION + 1];
   for (int a = 0; a < C_RESOLUTION + 1; a++) {
      contact[a] = 0;
   }
   for (int i = 0; i < N_OBJECTS; i++) {
      // cout<<"object no="<<i<<": ";
      int firsti = get_first_index(i);
      int lasti = get_last_index(i);
      // cout<<"first="<<firsti<<", last="<<lasti<<"; ";
      for (int j = firsti; j <= lasti; j++) {
         if (movement_list[i][j].action == fdcdetachment) {
            // cout<<"j="<<j<<", contact="<<movement_list[i][j].fdc_contact_time<<";; ";
            int ci = int (60. * movement_list[i][j].fdc_contact_time / DC);
            if (ci > C_RESOLUTION) { ci = C_RESOLUTION; }
            ++contact[ci];
         }
      }
      // cout<<"\n";
   }

   long total_measurements = 0;
   for (int i = 0; i <= C_RESOLUTION; i++) {
      total_measurements += contact[i];
   }

   ofstream fff;
   fff.open("fdccontact.out");
   fff << "! CC-FDC-contact-time histogram:\n"
       << "! contact time min : #occurrence : total#measurement : %occurrence\n";
   for (int i = 0; i <= C_RESOLUTION; i++) {
      // fff<<i*DELTA_V<<"   "
      fff << (double (i) + 0.5) * DC << "   "
          << contact[i] << "  "
          << total_measurements << "   "
          << double (contact[i]) / double (total_measurements)
          << "\n";
   }
   fff.close();

   fff.open("fdccmean.out");
   fff << "! average CC-FDC-contact time in minutes\n"
       << "! mean-contact : #measurements :  0.0\n";
   int n = 0;
   double cm = 0.;
   for (int i = 0; i < C_RESOLUTION; i++) {
      n += contact[i];
      cm += double (contact[i]) * (double (i) + 0.5) * DC;
   }
   if (n > 0) { cm /= double (n); }
   fff << cm << "  " << n << "    0.0\n";
   fff.close();
}
// ====================================================================
// ====================================================================
// ====================================================================
// ====================================================================

void TRACK::Write_shape() {
   /* The time course of the object shape index:
    * Long axis to radius is saved in .elongation
    * Long to short axis is saved in .l2s_axis
    * Measurements are done at intervals DELTA_T.
    *
    * Frequency of occurence of object shape:
    * Minimum shape index is 1.
    * Use S_RESOLUTION different speeds with DELTA_S interval.
    */
   // full analogy to Write_speed()
   // replace "v" by "si"

   ofstream fff("axis_t.out");
   fff << "! Time course of elongation of tracked cells:\n"
       << "! Based on long2short-axis (1).\n"
       << "! Based on long2radius (2).\n"
       << "! The long axis is calculated in direction of polarity.\n"
       << "! The short axis is calculated perpendicular to the polarity.\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   fff << "! time(min) : elongation (1) :  elongation (2)\n";

   ofstream ffx("axismean_t.out");
   ffx << "! Time course of mean elongation of all tracked cells:\n"
       << "! Based on long2short-axis (1).\n"
       << "! Based on long2radius (2).\n"
       << "! The long axis is calculated in direction of polarity.\n"
       << "! The short axis is calculated perpendicular to the polarity.\n"
       << "! Start at " << TRACKFROM << " until " << TRACKUNTIL << " hours.\n"
       << "! Tracking time interval " << 60. * DELTA_T << " minutes.\n";
   ffx << "! time(min) : mean elongation (1) : standard deviation :"
       << "! mean elongation (2) : standard deviation\n";

   long axis_l2s[S_RESOLUTION + 1];  // Based on long2short axis
   long axis_l2r[S_RESOLUTION + 1];  // Based on long axis to radius
   for (int a = 0; a < S_RESOLUTION + 1; a++) {
      axis_l2s[a] = 0;
      axis_l2r[a] = 0;
   }

   // save reference time:
   double t = TRACKFROM;
   // Introduce an index for every tracked cell
   int m[N_OBJECTS];
   // Save current elongation values somewhere:
   double l2s[N_OBJECTS], l2r[N_OBJECTS];
   for (int a = 0; a < N_OBJECTS; a++) {
      m[a] = 0;
      l2s[a] = 0.;
      l2r[a] = 0.;
   }

   long axis1_l2s[N_OBJECTS][S_RESOLUTION + 1];
   long axis1_l2r[N_OBJECTS][S_RESOLUTION + 1];
   for (int i = 0; i < N_OBJECTS; i++) {
      for (int j = 0; j < S_RESOLUTION; j++) {
         axis1_l2s[i][j] = 0;
         axis1_l2r[i][j] = 0;
      }
   }

   double mean_obj_pool = 0.;
   int n_timesteps = 0;

   // Now successively increase the evaluated time
   while (t <= TRACKUNTIL - DELTA_T) {
      t += DELTA_T;
      ++n_timesteps;
      // cerr<<"t="<<t<<" ...\n";
      fff << 60. * (t - TRACKFROM) << "  ";
      ffx << 60. * (t - TRACKFROM) << "  ";
      // variables for mean speed
      double obj_pool = 0.;
      double suml2s = 0., suml2r = 0.;
      // Calculate speeds for all objects at this time
      for (int i = 0; i < N_OBJECTS; i++) {
         int firsti = get_first_index(i);
         if (movement_list[i][firsti].status != nocell) {
            // Gehe in Zelle i solange vor bis die aktuelle Zeit t erreicht ist
            while (m[i] < movement_list[i].benutzt()
                   && movement_list[i][m[i]].t < t + DELTA_T * 1.e-08) {
               ++m[i];
            }
            // Nun steht der Index eins zu weit!
            --m[i];
            // cout<<"result: i="<<i<<", m[i]="<<m[i]<<"; \n";
            if (m[i] < 0) { m[i] = 0; } else {
               // If the cell is still alive and has a past
               if (movement_list[i][m[i]].status != nocell) {
                  l2s[i] = movement_list[i][m[i]].l2s_axis;
                  l2r[i] = movement_list[i][m[i]].elongation;
                  // add to the sum of all objects;
                  suml2s += l2s[i];
                  suml2r += l2r[i];
                  // increment the number of treated objects
                  obj_pool += 1.0;
                  // write to file
                  fff << l2s[i] << "  " << l2r[i] << "  ";       // ratios

                  // Save these elongations in the histogram
                  int ax_i = int (l2s[i] / DELTA_S + 0.5);
                  if (ax_i > S_RESOLUTION) { ax_i = S_RESOLUTION; }
                  ++axis_l2s[ax_i];
                  ++axis1_l2s[i][ax_i];
                  ax_i = int (l2r[i] / DELTA_S + 0.5);
                  if (ax_i > S_RESOLUTION) { ax_i = S_RESOLUTION; }
                  ++axis_l2r[ax_i];
                  ++axis1_l2r[i][ax_i];
               } else { fff << "*  *  "; }
            }
         }    // if movement_list[i]>0
         else { fff << "*  *  "; }
      }   // for i
      fff << "\n";

      // Calculate the mean speed from the summed displacements
      double meanl2s = -1.0, meanl2r = -1.0;
      mean_obj_pool += obj_pool;
      if (obj_pool > 0) {
         meanl2s = suml2s / obj_pool;
         meanl2r = suml2r / obj_pool;

         // Calculate the corresponding standard deviation
         double sigmal2s = 0., sigmal2r = 0.;
         if (obj_pool > 1) {
            for (int i = 0; i < N_OBJECTS; i++) {
               if (movement_list[i][m[i]].status != nocell) {
                  sigmal2s += (l2s[i] - meanl2s) * (l2s[i] - meanl2s);
                  sigmal2r += (l2r[i] - meanl2r) * (l2r[i] - meanl2r);
               }
            }
            sigmal2s /= (obj_pool - 1);
            sigmal2s = sqrt(sigmal2s);
            sigmal2r /= (obj_pool - 1);
            sigmal2r = sqrt(sigmal2r);
         }

         // Write to file
         ffx << meanl2s << "  " << sigmal2s << "    " << meanl2r << "  " << sigmal2r << "\n";
      } else { ffx << "*   0   *   0"; }
   }
   fff.close();
   ffx.close();

   if (n_timesteps > 0) { mean_obj_pool /= double (n_timesteps); }

   fff.open("axishisto.out");
   fff << "! Elongation histogram:\n"
       << "! (1) based on long to short axis\n"
       << "! (2) based on long axis to radius\n"
       << "! elongation : #(1) : #(2) : mean_cell# : #percell(1) : #percell(2) : "
       << "N_OBJECTS : #perN(1) : #perN(2) : #perN_sd(1) : #perN_sd(2)\n";
   for (int i = 0; i < S_RESOLUTION; i++) {
      fff << i * DELTA_S << "   "
          << axis_l2s[i] << "  "
          << axis_l2r[i] << "    "
          << mean_obj_pool << "    "
          << double (axis_l2s[i]) / mean_obj_pool << "  "
          << double (axis_l2r[i]) / mean_obj_pool << "    "
          << N_OBJECTS << "    ";
      // Note that the standard deviation is ill-defined if based on mean_obj_pool.
      // Thus, the averages are given with respect to N_OBJECTS as well, irrespective
      // of the number of cells died during the time window of tracking:
      // Get averages of numbers of occurrences of speeds:
      double n_l2s_percell = double (axis_l2s[i]) / double (N_OBJECTS);
      double n_l2r_percell = double (axis_l2r[i]) / double (N_OBJECTS);
      // Calculate standard deviations:
      double sigma_l2s = 0., sigma_l2r = 0.;
      if (N_OBJECTS > 1) {
         for (int j = 0; j < N_OBJECTS; j++) {
            sigma_l2s += ((double (axis1_l2s[j][i]) - n_l2s_percell)
                          * (double (axis1_l2s[j][i]) - n_l2s_percell));
            sigma_l2r += ((double (axis1_l2r[j][i]) - n_l2r_percell)
                          * (double (axis1_l2r[j][i]) - n_l2r_percell));
         }
         sigma_l2s /= (N_OBJECTS - 1);
         sigma_l2s = sqrt(sigma_l2s);
         sigma_l2r /= (N_OBJECTS - 1);
         sigma_l2r = sqrt(sigma_l2r);
      }
      // Write to file:
      fff << double (axis_l2s[i]) / double (N_OBJECTS) << "  "
          << double (axis_l2r[i]) / double (N_OBJECTS) << "    "
          << sigma_l2s << "  "
          << sigma_l2r
          << "\n";
   }
   fff.close();

   fff.open("axismean_all.out");
   fff << "! average elongation of cells\n"
       << "! long2short axis : # measurements : longaxis2radius : # measurements : 0.0\n";
   int nl2s = 0, nl2r = 0;
   double axisml2s = 0., axisml2r = 0.;
   for (int i = 0; i < S_RESOLUTION; i++) {
      nl2s += axis_l2s[i];
      nl2r += axis_l2r[i];
      axisml2s += axis_l2s[i] * i * DELTA_S;
      axisml2r += axis_l2r[i] * i * DELTA_S;
   }
   if (nl2s > 0) { axisml2s /= nl2s; }
   if (nl2r > 0) { axisml2r /= nl2r; }
   fff << axisml2s << "  " << nl2s << "    " << axisml2r << "  " << nl2r << "     0.0\n";
   fff.close();
}
