#ifndef NEK5000_BRIDGE_H
#define NEK5000_BRIDGE_H

#include <mpi.h>

#ifdef __cplusplus
extern "C" {
#endif

#ifdef UPCASE
#  define FORTRAN_NAME(low,up) up
#else
#ifdef UNDERSCORE
#  define FORTRAN_NAME(low,up) low##_
#else
#  define FORTRAN_NAME(low,up) low
#endif
#endif

#define sensei_bridge_initialize FORTRAN_NAME(sensei_bridge_initialize, SENSEI_BRIDGE_INITIALIZE)
#define sensei_bridge_update FORTRAN_NAME(sensei_bridge_update, SENSEI_BRIDGE_UPDATE)
#define sensei_bridge_finalize FORTRAN_NAME(sensei_bridge_finalize, SENSEI_BRIDGE_FINALIZE)

void sensei_bridge_initialize(MPI_Comm* comm, char* casename, int* elems,
      int* vx_dim, int* vy_dim, int* vz_dim, double* vmesh_x, double* vmesh_y, double* vmesh_z,
      int* px_dim, int* py_dim, int* pz_dim, double* pmesh_x, double* pmesh_y, double* pmesh_z,
      double* vel_x, double* vel_y, double* vel_z, double* pressure, double* temp, double* jacobian, int* t_dim);
void sensei_bridge_update(int* tstep, double* time);
void sensei_bridge_finalize();

#ifdef __cplusplus
} // extern "C"
#endif

#endif
