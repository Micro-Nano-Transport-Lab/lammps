#ifdef FIX_CLASS

FixStyle(wall/reflect/cyli,FixWallReflectCyli)

#else

#ifndef LMP_FIX_WALL_REFLECT_CYLI_H
#define LMP_FIX_WALL_REFLECT_CYLI_H

#include "fix.h"

namespace LAMMPS_NS {

class FixWallReflectCyli : public Fix {
 public:
  enum{NONE=0,CONSTANT};

  FixWallReflectCyli(class LAMMPS *, int, char **);
  virtual ~FixWallReflectCyli();
  int setmask();
  void init();
  void post_integrate();

 protected:
  double rad, bottom, lid;
  double xscale,yscale,zscale;

  virtual void wall_particle(double rad, double bottom, double lid);
};

}

#endif
#endif
