
#include <iomanip>
#include "iPic3D.h"


// init static HDF5 groups
//#if defined(FTI_CKPT) && defined(FTI_HDF5)
//FTIT_H5Group topoGroup;
//FTIT_H5Group fieldGroup;
//FTIT_H5Group BxGroup;
//FTIT_H5Group ByGroup;
//FTIT_H5Group BzGroup;
//FTIT_H5Group ExGroup;
//FTIT_H5Group EyGroup;
//FTIT_H5Group EzGroup;
//#endif
//int dummy = 1;

using namespace iPic3D;

int main(int argc, char **argv) {
  
  iPic3D::c_Solver KCode;
  bool b_err = false;

  /* ------------------------------ */
  /* 0- Initialize the solver class */
  /* ------------------------------ */

  KCode.Init(argc, argv);
  KCode.InjectBoundaryParticles();
  KCode.GatherMoments();

  /* ------------ */
  /* 1- Main loop */
  /* ------------ */

  for (int i = KCode.FirstCycle(); i <= KCode.LastCycle(); i++) {

    if (KCode.get_myrank() == 0) cout << " ======= Cycle " << i << " ======= " << endl;

    /* ----------------------------------------------------- */
    /* 2- Calculate fields and move particles                */
    /*    Exit if there is a memory issue with the particles */
    /* ----------------------------------------------------- */

    KCode.UpdateCycleInfo(i);
    KCode.CalculateField();

    b_err = KCode.ParticlesMover();

    if (!b_err) KCode.CalculateBField();
    if (!b_err) KCode.GatherMoments();
    if ( b_err) i = KCode.LastCycle() + 1;

    /* --------------- */
    /* 3- Output files */
    /* --------------- */

    KCode.WriteOutput(i);
    KCode.WriteConserved(i);
    KCode.WriteRestart(i);

  }

  KCode.Finalize();

  return 0;
}
