#include "space.h"
// #include "signals.h"
#include <math.h>
#include <string.h>

/*
 * states operator++(const states& a) {
 * states b;
 * if (a==nocell) b=CB;
 * else if (a==CB) b=CC;
 * else if (a==CC) b=FDC;
 * else if (a==FDC) b=TC;
 * else if (a==TC) b=out;
 * else if (a==out) b=blast1;
 * else if (a==blast1) b=blast2;
 * else b=N_cells;
 * return b;
 * }
 */

spacepoint::spacepoint() {
   cell = nocell;
   listi = -1;
   FDClisti = -1;
}
spacepoint::spacepoint(long int &i, states &c, long &listindex, long &FDClistindex) {
   cell = c;
   listi = listindex;
   FDClisti = FDClistindex;
}
spacepoint::~spacepoint() { }
char operator ==(const spacepoint &a, const spacepoint &b) {
   if ((a.cell != b.cell) || (a.listi != b.listi) || (a.FDClisti != b.FDClisti)) { return 0; }
   return 1;
}
char operator !=(const spacepoint &a, const spacepoint &b) {
   return !((a == b) == 1);
}
// ### flag Neumann oder Moore einfuehren

space::space()
   : grid() {
   cellknot = new spacepoint[pointnumber];
}
space::space(Parameter &par, ofstream &ana)
   : grid(cellspace,
          par.Value.system,
          par.Value.DimSpace,
          par.Value.dx,
          par.Value.GC_radius,
          par.Value.gridsize,
          par.Value.vol_shape,
          par.Value.obstacles,
          par.Value.wall_level,
          par.Value.collagen_density,
          par.Value.collagen_cluster,
          par.Value.wall_width,
          par.Value.slit_number,
          par.Value.slit_width,
          par.Value.use_specific_turning_angles,
          ana) {
   // grid(cellspace,par,ana);
   cellknot = new spacepoint[pointnumber];
   for (long n = 0; n < pointnumber; n++) {
      if (knot[n].status == external) { cellknot[n].cell = N_cells; }
   }
   zone_separator = get_zone_separation_index();
   ana << "zone_separation_lattice_index = " << zone_separator << "\n";
}
space::~space() {
   // cout<<"in ~space()\n";
   delete[] cellknot;
}
void space::set_knot(const long &i, const states &s, const long &li) {
   // cout<<"i="<<i<<"; s="<<s<<"; li="<<li<<"\n";
   knot[i].status = object;
   cellknot[i].cell = s;
   cellknot[i].listi = li;
   // cout<<"end set_knot.\n";
}
void space::clear_knot(const long &i) {
   knot[i].status = nothing;
   cellknot[i].cell = nocell;
   cellknot[i].listi = -1;
}
long space::get_random_empty_position() {
  /* Returns the index of an empty position on the grid within a reaction volume.
     If the grid is full, -1 is returned.
  */
  long getpos = -1;
  long k[dim];
  bool foundaplace = false;
  long trials = 0;
  // OPTION ++++++++++++++++++++++++++++++++++++++++++++++++
  // long factor = 1; 
  long max_trials = pointnumber; // or factor*pointnumber;
  // end OPTION ++++++++++++++++++++++++++++++++++++++++++++
  /* max_trials is a smooth condition: when the grid runs full, new divisions are
   * already prevented when it becomes unlikely to find an empty space within a
   * reasonable number of random trials. Thus, the system will rarely reach the
   * state of really being full.
   */
  while (not (foundaplace) && trials < max_trials) {
    for (int i = 0; i < dim; i++) {
      k[i] = irandom(prodimvec[i]);
    }
    getpos = Index(k);
    // As the <states> value <nocell> is only possible for <grid_states> unequal <external>,
    // it is sufficient to check <states>==<nocell> (even though the random values above
    // might generate points outside the reaction volume):
    if ((getpos != -1) && (cellknot[getpos].cell == nocell)) {
      foundaplace = true;
    }
    ++trials;
  }
  return getpos;
}
long space::get_zone_separation_index() {
  /* The following definition divides the GC volume into DZ and LZ of equal size.
   * The first half of the indices on the lattice defines the LZ
   * The routine returns the first index on the grid attributed to the DZ.
   * Note that the central plane of the reaction volume is counted as LZ! 
   */
  long max = long (prodimvec[dim - 1] / 2) + 1;
  for (short d = 0; d < dim - 1; d++) {
    max *= prodimvec[d];
  }
  return max;
}
short space::object_border(const long &i) {
   // return 0 if all neighbors of i belong to the same object, 1 otherwise
   // uses Neumann neighborhood only
   for (short d = 0; d < dim2; d++) {
      if ((knot[i].near_n[d] == -1)
          || (cellknot[i].cell != cellknot[knot[i].near_n[d]].cell)
          || (cellknot[i].listi != cellknot[knot[i].near_n[d]].listi)) {
         return 1;
      }
   }
   return 0;
}
long space::check_neighbors_for_border(const long i, const long &ref, double &min_distance,
                                       dynarray<long> &done) {
   // cout<<"start check_neighbors("<<i<<") ... ";
   // go through all neighbors
   long nn[dim2];
   dynarray<long> donehere(dim2,1,0);
   short n = 0;
   short checked = 0;
   double distance;
   for (short d = 0; d < dim2; d++) {
      // Bedingung: nur probieren falls nicht schon geschehen
      if (done.find(knot[i].near_n[d]) == -1) {
         ++checked;
         // check if one of them is an object_border-point
         if (object_border(knot[i].near_n[d]) == 1) {
            nn[n] = knot[i].near_n[d];
            distance = Abstand(nn[n],ref);
            // der Fall "gleicher Abstand" kann zwar durch Rundung in der
            // folgenden Bedingung falsch beurteilt werden, da aber ohnehin
            // eine zufaellige Nachauswahl getroffen wird ist die zufaellige
            // Auswahl durch Rundungsfehler ebensogut.
            if ((distance < min_distance + 1.e-08) && (distance > min_distance - 1.e-08)) {
               ++n;
            } else if (distance < min_distance) {
               min_distance = distance;
               nn[0] = knot[i].near_n[d];
               n = 1;
            }
         } else {
            donehere.add(knot[i].near_n[d]);
            // in donehere werden alle gespeichert, die hier getestet wurden
            // und keine Randpunkte sind. Diese werden unten weiteruntersucht.
            // ### kann es sein, dass ein Randpunkt, der naeher zu ref liegt
            // ### nur ueber einen anderen Randpunkt gefunden werden kann?
            // ### Ich glaube nicht, bin aber nicht sicher.
         }
         done.add(knot[i].near_n[d]);
      }
   }
   // done.show();
   // Falls checked=0 wurde kein neuer Nachbar mehr gefunden: Abbruch mit -1
   if (checked == 0) { return -1; }
   // Falls n=0 wurde kein border-point gefunden
   // Falls n=1 wurde der border-point gefunden, der am naehesten zu ref liegt
   // Falls n>1 wurden mehrere zur ref gleichweit entfernte border-points gefunden
   if ((n == 0) && (Abstand(long (i),long (ref)) <= min_distance + 1)) {
      // cout<<"  go through neighbors("<<i<<") ... ";
      for (short d = 0; d < dim2; d++) {
         if ((done.find(knot[i].near_n[d]) == -1) || (donehere.find(knot[i].near_n[d]) != -1)) {
            // mache nur falls der Punkt nicht schon frueher untersucht wurde,
            // mache trotzdem, falls er nur eben oben als Nachbar untersucht wurde.
            nn[n] = check_neighbors_for_border(knot[i].near_n[d],ref,min_distance,done);
            if (nn[n] != -1) {
               // es wurde was gefunden
               // distance=Abstand(nn[n],long(ref));
               // cout<<"  distance="<<distance<<" min_distance="<<min_distance<<"\n";
               // if (distance==min_distance) ++n;
               // else if (distance<min_distance) {
               // min_distance=distance;
               ++n;
               // }
            }
         }
      }
      // Behalte nur die Punkte, die einen kleineren Abstand bedeuten:
      // n enthaelt die Zahl der gefundenen Punkte:
      short d = 0;
      // starte ab 1 gefundenem Punkt (d.h. n=1)
      while (d < n) {
         distance = Abstand(nn[d],long (ref));
         // cout<<"(n,d)=("<<n<<","<<d<<"): dist(nn[d])="<<distance<<"\n";
         if (distance > min_distance + 1.e-08) {
            --n;
            // cout<<"kopiere nn["<<n<<"] auf nn["<<d<<"], d.h. "<<nn[n]<<" auf "<<nn[d]<<"\n";
            nn[d] = nn[n];
         } else { ++d; }
      }
   }
   if (n == 0) {
      // cout<<i<<" (n=0) returns -1\n";
      return -1;
   }
   int n1 = 0;
   if (n > 1) { n1 = irandom(n); }
   // liefert eine Zahl aus 0,1,...,n-1 (und n mit quasi-Null Wahrscheinlichkeit)
   if (n1 == n) { n1 = n - 1; }
   // cout<<i<<" returns "<<nn[n1]<<"\n";
   // char c; cin>>c;
   return nn[n1];
}
long space::next_border(const long &i) {
   /* returns the next point belonging to the same object as i and being
    * at the border of this object. This routine makes sense only if called
    * from a point i that is an internal point of the object without direct
    * contact to other objects.
    */
   if (object_border(i) == 1) { return i; } else {
      dynarray<long> done(dim2 * dim2,1,0);
      done.add(i);
      double min_distance = 2.0 * prodim2;
      return check_neighbors_for_border(i,i,min_distance,done);
   }
}
short space::self(const long &i, const long &li, const states &s) {
   if ((i != -1) && (cellknot[i].cell == s) && (cellknot[i].listi == li)) {
      return 1;
   } else { return 0; }
}
short space::get_n_self(const long &i, const long &li,
                        const states &s, long excluded) {
   /* Calculates the number of neighbors to lattice point i that belong to
    * the cell of type "status" denoted by the cell-list-index li. If a specific point should
    * not be counted it can be excluded by handing over its index to excluded
    */
   short nself = 0;
   long k;
   if (excluded == -1) { excluded = -2; }
   for (short j = 0; j < dim2; j++) {
      k = knot[i].near_n[j];
      if ((k != excluded) && (self(k,li,s) == 1)) { ++nself; }
   }
   return nself;
}
long space::end_of_tail(const long &start, long back,
                        const long &object, const states &s) {
   /*
    * Follows up the neighbours of "start" up to a neighbour which has more
    * than two or only one neighbour belonging to the same "object" in the
    * same "s". "back" is excluded from consideration and defines the
    * direction of search.
    * Returns:
    * -1*index is returned if a neighbour at "index" with three or more neighbours is found.
    * The index of the neighbour is returned, if this neighbour has only one neighbour.
    * If "0" is returned this means that it is an end-of-tail, as at "index=0"
    *    three neighbours don't exist.
    *
    * This routine tests, if "start" is part of a tail of the "object" which
    * extends in the direction opposite of "back". If so the end of the tail
    * is returned, if not -1 is returned.
    */
   short nself;
   long k,newstart,nowstart;
   nowstart = start;
   newstart = -1;
   short begin = 1;
   while (nowstart != start || begin == 1) {
      begin = 0;
      nself = 0;
      for (short j = 0; j < dim2; j++) {
         k = knot[nowstart].near_n[j];
         if ((k != back) && (self(k,object,s) == 1)) {
            ++nself;
            newstart = k;
         }
      }
      if (nself == 0) {
         // i.e. no further neighbour --> end of tail!
         return nowstart;
      }
      if (nself > 1) {
         // i.e. more than two neighbours (including "back") --> no tail!
         return -1 * nowstart;
      }
      // else: only one neighbour found in addition to "back" -->
      // take the actual point as "back" point
      back = nowstart;
      // and the re-start point as its neighbour
      nowstart = newstart;
      // and look for the next neighbours ...
   }
   // Running in a circle. --> Tail closing in itself --> return the starting point!
   return start;
}
long space::decompose(const long &i, const long &li, const states &s) {
   /*
    * Checks if for all neighbors x of i on the lattice that belong to the
    * cell li with state s a neighbor y of i exists that has a common
    * neighbor with x belonging to the same cell li in the same state.
    *
    * Returns -1 if the movement of i is forbidden, i.e. the cell would be
    *    decomposed.
    * Returns the index of a fragment if this fragment can be moved without risk.
    *    This is normally "i", but a different fragment may also be moved instead
    *    (see problem of tail shortening).
    */
   // #### ist es noetig, dass man alle x durchlaeuft, bevor abgebrochen wird??? CPU-time
   long xx,y,xn,yn;
   short k,l,m,found;
   short totalok = 0;
   // Find the number of next neighbors belonging to the same cell:
   short n_self = get_n_self(i,li,s,-1);
   // cout<<"\n ... in decompose: n_self="<<n_self;
   // If only one neighbor exists, the movement of the fragment at index i is allowed:
   if (n_self < 2) { return i; } else {
      // Otherwise, go through all neighbors of i:
      for (short n = 0; n < dim2; n++) {
         found = 0;
         // Get the index of the neighbor n:
         xx = knot[i].near_n[n];
         // Test only if xx belongs to cell li and if xx is inside the reaction volume:
         if (self(xx,li,s) == 1) {
            m = 0;
            // Go through all other neighbors of point i:
            while (m < dim2 && found == 0) {
               if (n == m) {
                  ++m;       // the neighbor under consideration is not reconsidered here
               } else {
                  // Get index of one of the other neighbors of point i:
                  y = knot[i].near_n[m];
                  /* Consider this neighbor y if not outside the reaction volume and
                   * if y belongs to the same cell */
                  if (self(y,li,s) == 1) {
                     // Look if both, the neighbors n and m are connected to each other:
                     k = 0;
                     // To this end, go through all neighbors of xx and ...
                     while (k < dim2 && found == 0) {
                        // ... get its index and ...
                        xn = knot[xx].near_n[k];
                        l = 0;
                        // ... check if xn is a neighbor of y:
                        if ((xn != i) && (self(xn,li,s) == 1)) {
                           while (l < dim2 && found == 0) {
                              yn = knot[y].near_n[l];
                              /* If, indeed, xn is a neighbor of y, yn==xn, and a connection
                               * of xx and y (vs. n and m) exists. Save it in found=1: */
                              if (yn == xn) { found = 1; } else { ++l; }
                           }
                           ++k;
                        } else { ++k; }
                     }
                     ++m;
                  } else { ++m; }
               }
            }
         } else { ++totalok; }
         totalok += found;
         /*
          * Each neighbor of i is either not belonging to cell li or a connection
          * to another neighbor of i has been found. If not in at least one case,
          * found will be zero in this case and consequently totalok will be smaller
          * than the total number of neighboring places, i.e. dim2. This is the
          * test criterion.
          */
      }
   }
   /*
    * Problem: For ring structures, the cell will seem to be disconnected by
    * the movement of any ring elements without being it. This is related to the
    * fact, that the present criterion is a local one. Far connection are not
    * detected. This error may lead to immobilized cells as in this example of
    * a 10 fragment cell:    XXXX
    *                        X  O
    *            YYYY
    * Ways out:
    * 1) Prevent such ring structures from being built. How can a ring be
    *    locally detected? A barycenter which is either empty or occupied
    * by a different object is not conclusive to identify a ring. But
    * reading point 3), every ring in question here is necessarily symmetric.
    * Otherwise fragments would exist that can be moved. Therefore, indeed,
    * all stable and immobile ring structures will have the property that
    * the barycenter is not occupied by the cell itself. One promising
    * solution is, to add additional tests for the movement if this criterion
    * appears to apply. This procedure would return then: o.k. you can move
    * i. In the calling routine the barycenter criterion is tested and
    * according to the result the target place for i is more carefully
    * tested because of risks to build a stable ring structure.
    * 2) Remember that the ring when it has been built is not really connected
    *    at the last fragment, i.e. O in the above example has moved in last,
    * and came from the bottom row of Y's. Then the fragment O has to remember
    * that it is connected to Y but not to X. Then the right X can be moved
    * as it has only one connection to the cell:   XXX
    *                                             XX O
    *                                 YYYY
    *    This may lead to other problems: Another X coming down in between X and O
    * will not be connected to O (according to the same criterion) neither to X.
    * We have to look for a criterion that remains local concerning connectivity
    * but avoids stable and immobile ring structures.
    * 3) Note, that such a stable ring structure can only appear for the total number
    *    of fragments being such that a ring with width of one fragment can be built.
    * Thicker rings always allow single fragments to move into the center (if not
    * another cell is enclosed by the ring). Therefore, it may be possible to
    * make an extra test for some specific fragment numbers that allow stable
    * ring states.
    * 4) It is always possible to replace this test routine here, by a general
    *    connectivity test as provided in check_connectivity(...). It is simply
    * a question of CPU time! Then the ring structure could be resolved as the
    * movement of every fragment in the ring will not destroy the connectivity
    * of the cell which is then correctly tested, using all cell fragments.
    * This is the always possible solution but it is non-local and therefore
    * time consuming.
    * Took solution 1) now!!!
    */
   if (totalok == dim2) {
      // cout<<"This move is allowed!!!! totalok="<<totalok<<" !!!!!!!!!!!!!!\n";
      // cout<<" -> allow "<<i<<" to move.\n";
      return i;
   } else {
      /*
       * Reconsider the case of two neighbours belonging to the same cell:
       * This may be a tail of one fragment width. Such tail develop from
       * the movement of the cell and are rather stable because all fragments
       * but the last are forbidden to be moved. The shortening of the tail
       * becomes rather unprobable. Therefore, here a criterion is checked
       * in order to decide if this is a fragment within a cell-tail and if
       * yes, the movement is allowed for the last fragment of the tail
       * instead.
       */
      if (n_self == 2) {
         /* One can be sure that both neighbours of "i" are not directly connected by one neighbour,
          * as this has been excluded in the former test routine. */
         // Test the neighbours again:
         long threeneighbours = 0;
         for (short n = 0; n < dim2; n++) {
            // Get the index of the neighbor n:
            xx = knot[i].near_n[n];
            // If xx is one of the two neighbours continue:
            if (self(xx,li,s) == 1) {
               xn = end_of_tail(xx,i,li,s);
               if (xn >= 0) {
                  // cout<<"move "<<xn<<" instead of "<<i<<"!\n";
                  // cout<<get_n_self(xn,li,s,-2)<<"; ";
                  return xn;
                  /* Note, that if the cell is composed of a line of fragments only,
                   * i.e. that two ends of tail exist, one of these ends would be
                   * chosen here by chance. This means that in this case always the
                   * neighbour (opening for a tail) which is situated at the beginning
                   * of the list would be taken for movement, implying an anisotropy.
                   * However, this more pathological constellation will be rather seldomly,
                   * such that it may be neglected (###).
                   */
               } else {
                  if (threeneighbours == 0) {
                     threeneighbours = xn;
                  } else if (threeneighbours == xn) {
                     return i;
                  }
                  /* This last case allows the originally intended move of i:
                   * If we are here, "i" has two self-neighbours.
                   * The tails starting from "i" have both been tested.
                   * The first tail ended in a point with 3 neighbours that has
                   * been saved in "3neighbours" (return of end_of_tail()).
                   * The second tail ends on the same point. It follows that
                   * the underlying fragment arrangement is a ring with widht
                   * of one fragment. The move is obligatory in order to dissolve
                   * the ring and respects the compactness of the cell.
                   */
               }
            }
         }
      }
      // cout<<"This move is forbidden!!!! totalok="<<totalok<<" !!!!!!!!!!!!!!\n";
      // cout<<" -> forbid\n\n";
      return -1;
   }
}
