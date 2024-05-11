
#include "fix_wall_reflect_cyli.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "error.h"
#include "input.h"
#include "lattice.h"
#include "modify.h"
#include "update.h"
#include "variable.h"
#include <cstring>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

FixWallReflectCyli::FixWallReflectCyli(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 8) error->all(FLERR,"Illegal fix wall/reflect/cyli command");

  dynamic_group_allow = 1;

  // parse args
  int scaleflag = 1;

  rad = utils::numeric(FLERR,arg[3],false,lmp);
  bottom = utils::numeric(FLERR,arg[4],false,lmp);
  lid = utils::numeric(FLERR,arg[5],false,lmp);

  if (arg[4]>arg[5])
      error->all(FLERR,"Cannot be upside down");

  if (strcmp(arg[6],"units") == 0) {
    if (strcmp(arg[7],"box") == 0) scaleflag = 0;
    else if (strcmp(arg[7],"lattice") == 0) scaleflag = 1;
    }
  else error->all(FLERR,"Illegal fix wall/reflect/cyli command");

  // error check
  if (domain->xperiodic)
    error->all(FLERR,"Cannot use fix wall/reflect/cyli in periodic dimension");
  if (domain->yperiodic)
    error->all(FLERR,"Cannot use fix wall/reflect/cyli in periodic dimension");
  if (domain->dimension == 2)
    error->all(FLERR,"Cannot use fix wall/reflect/cyli for a 2d simulation");

  if (scaleflag) {
      xscale = domain->lattice->xlattice;
      yscale = domain->lattice->ylattice;
      zscale = domain->lattice->zlattice;
  }
  else xscale = yscale = zscale = 1.0;
}

/* ---------------------------------------------------------------------- */

FixWallReflectCyli::~FixWallReflectCyli()
{
  if (copymode) return;
}

/* ---------------------------------------------------------------------- */

int FixWallReflectCyli::setmask()
{
  int mask = 0;
  mask |= POST_INTEGRATE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixWallReflectCyli::init()
{
  int nrigid = 0;
  for (int i = 0; i < modify->nfix; i++)
    if (modify->fix[i]->rigid_flag) nrigid++;

  if (nrigid && comm->me == 0)
    error->warning(FLERR,"Should not allow rigid bodies to bounce off "
                   "relecting walls");
}

/* ---------------------------------------------------------------------- */

void FixWallReflectCyli::post_integrate()
{
  wall_particle(rad,bottom,lid);
}

/* ----------------------------------------------------------------------
   this method may be overwritten by a child class
------------------------------------------------------------------------- */

void FixWallReflectCyli::wall_particle(double rad, double bottom, double lid)
{
  int i;
  double delta_x, delta_y, delta_z;

  double **x = atom->x;
  double **v = atom->v;
  double **f = atom->f;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;
  double dt_step = update->dt;

  for (i = 0; i < nlocal; i++) {
    if (mask[i] & groupbit) {

      double dist_cyli = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);

      if (dist_cyli > rad) {
        double init_x = x[i][0] - v[i][0]*dt_step;
        double init_y = x[i][1] - v[i][1]*dt_step;
        double init_z = x[i][2] - v[i][2]*dt_step;

        double delx = x[i][0] - init_x;
        double dely = x[i][1] - init_y;
        double delz = x[i][2] - init_z;

        int p = update->ntimestep;
        int q = atom->tag[i];
 
        double t1 = (-2*(init_x*delx+init_y*dely)+sqrt(4*(init_x*delx+init_y*dely)*(init_x*delx+init_y*dely)-4*(delx*delx+dely*dely)*((init_x*init_x)+(init_y*init_y)-rad*rad)))/2*(delx*delx+dely*dely);
        double t2 = (-2*(init_x*delx+init_y*dely)-sqrt(4*(init_x*delx+init_y*dely)*(init_x*delx+init_y*dely)-4*(delx*delx+dely*dely)*((init_x*init_x)+(init_y*init_y)-rad*rad)))/2*(delx*delx+dely*dely);

        if ((p == 1260) && (q == 31543)) {
            char str1[512];
            sprintf(str1,"Reflecting particle info 1" BIGINT_FORMAT " "
                    TAGINT_FORMAT " "" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g",
                    update->ntimestep,atom->tag[i],x[i][0],x[i][1],x[i][2],init_x,init_y,init_z,delx,dely,delz,t1,t2);
            error->warning(FLERR,str1,0);
        }

        // Intersection point
        double inter1_x = init_x + delx*t1;
        double inter1_y = init_y + dely*t1;
        double inter1_z = init_z + delz*t1;

        double inter2_x = init_x + delx*t2;
        double inter2_y = init_y + dely*t2;
        double inter2_z = init_z + delz*t2;

        double dist1 = sqrt((x[i][0] - inter1_x)*(x[i][0] - inter1_x) + (x[i][1] - inter1_y)*(x[i][1] - inter1_y) + (x[i][2] - inter1_z)*(x[i][2] - inter1_z));
        double dist2 = sqrt((x[i][0] - inter2_x)*(x[i][0] - inter2_x) + (x[i][1] - inter2_y)*(x[i][1] - inter2_y) + (x[i][2] - inter2_z)*(x[i][2] - inter2_z));

        if (dist1<dist2) {
          delta_x = x[i][0] - inter1_x;
          delta_y = x[i][1] - inter1_y;
          delta_z = x[i][2] - inter1_z;

          x[i][0] = inter1_x - delta_x;
          x[i][1] = inter1_y - delta_y;
          x[i][2] = inter1_z - delta_z;      
        }

        else {
          delta_x = x[i][0] - inter2_x;
          delta_y = x[i][1] - inter2_y;
          delta_z = x[i][2] - inter2_z;

          x[i][0] = inter2_x - delta_x;
          x[i][1] = inter2_y - delta_y;
          x[i][2] = inter2_z - delta_z;
        }  
        
        if ((p == 1260) && (q == 31543)) {
            char str1[512];
            sprintf(str1,"Reflecting particle info 2" BIGINT_FORMAT " "
                    TAGINT_FORMAT " "" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g",
                    update->ntimestep,atom->tag[i],x[i][0],x[i][1],x[i][2],delta_x,delta_y,delta_z,dist1,dist2,inter1_x,inter1_y,inter1_z,inter2_x,inter2_y,inter2_z);
            error->warning(FLERR,str1,0);
        }

        //NAN Check      
        for (int k = 0; k < 3; k++) { 
            if (!std::isfinite(x[i][k])) {
            char str3[128];
            sprintf(str3," NAN Check " BIGINT_FORMAT " "
                    TAGINT_FORMAT " "" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g"" %g",
                    update->ntimestep,atom->tag[i],init_x,init_y,init_z,x[i][0],x[i][1],x[i][2],dist_cyli,t1,t2);
            error->warning(FLERR,str3,0);
            error->one(FLERR,"Non-numeric atom coords - simulation unstable");
            }
        }
        
        // flip velocities and forces
        v[i][0] = -v[i][0];
        v[i][1] = -v[i][1];
        v[i][2] = -v[i][2];

        f[i][0] = -f[i][0];
        f[i][1] = -f[i][1];
        f[i][2] = -f[i][2];       
      }
    }
  }
}
