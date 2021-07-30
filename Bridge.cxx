#include "Bridge.h"
#include "DataAdaptor.h"

#include <ConfigurableAnalysis.h>
#include <vtkNew.h>
#include <vtkDataObject.h>
#include <vtkDoubleArray.h>
#include <vtkSOADataArrayTemplate.h>
#include <Profiler.h>

namespace BridgeInternals
{
  static vtkSmartPointer<nek_sensei::DataAdaptor> GlobalDataAdaptor;
  static vtkSmartPointer<sensei::ConfigurableAnalysis> GlobalAnalysisAdaptor;
}

//-----------------------------------------------------------------------------
void sensei_bridge_initialize(MPI_Comm* comm, char* casename, int* elems,
      int* vx_dim, int* vy_dim, int* vz_dim, double* vmesh_x, double* vmesh_y, double* vmesh_z,
      int* px_dim, int* py_dim, int* pz_dim, double* pmesh_x, double* pmesh_y, double* pmesh_z,
      double* vel_x, double* vel_y, double* vel_z, double* pressure, double* temp, double* jacobian, int* t_dim,
      double* x_min, double* x_max, double* y_min, double* y_max, double* z_min, double* z_max)
{
  sensei::Profiler::Initialize();
  sensei::Profiler::StartEvent("nek::sensei_bridge::initialize");
  int arrayLen = (*vx_dim) * (*vy_dim) * (*vz_dim) * (*elems);
  BridgeInternals::GlobalDataAdaptor = vtkSmartPointer<nek_sensei::DataAdaptor>::New();
  BridgeInternals::GlobalDataAdaptor->SetCommunicator(*comm);
  BridgeInternals::GlobalDataAdaptor->Initialize(0,
        *vx_dim, *vy_dim, *vz_dim, *elems, vmesh_x, vmesh_y, vmesh_z, x_min, x_max, y_min, y_max, z_min, z_max,
        vel_x, vel_y, vel_z, pressure, arrayLen);
  if(*vx_dim != *px_dim){
    BridgeInternals::GlobalDataAdaptor->Initialize(1,
        *px_dim, *py_dim, *pz_dim, *elems, pmesh_x, pmesh_y, pmesh_z, x_min, x_max, y_min, y_max, z_min, z_max,
        vel_x, vel_y, vel_z, pressure, arrayLen);
  }

/*
  vtkSOADataArrayTemplate<double> *velocity = vtkSOADataArrayTemplate<double>::New();
  velocity->SetNumberOfComponents(3);
  velocity->SetArray(0, vel_x, arrayLen, true, true);
  velocity->SetArray(1, vel_y, arrayLen, true, true);
  velocity->SetArray(2, vel_z, arrayLen, true, true);
  velocity->SetName("velocity");
*/
  //BridgeInternals::GlobalDataAdaptor->AddArray(0, "velocity", velocity);

  /*
  BridgeInternals::GlobalDataAdaptor->AddArray(1, "pressure", pressure);
  for(int i = 0; i < *t_dim; i++)
    BridgeInternals::GlobalDataAdaptor->AddArray(0, "t"+i, temp + arrayLen*i);

  BridgeInternals::GlobalDataAdaptor->AddArray(1, "jacobian", jacobian);
*/
  BridgeInternals::GlobalAnalysisAdaptor = vtkSmartPointer<sensei::ConfigurableAnalysis>::New();
  BridgeInternals::GlobalAnalysisAdaptor->SetCommunicator(*comm);
  std::string caseStr = casename;
  BridgeInternals::GlobalAnalysisAdaptor->Initialize(caseStr + ".xml");
  sensei::Profiler::EndEvent("nek::sensei_bridge::initialize");
}

//-----------------------------------------------------------------------------
void sensei_bridge_update(int* tstep, double* time)
{
  sensei::Profiler::StartEvent("nek::sensei_bridge::update");
  BridgeInternals::GlobalDataAdaptor->SetDataTime(*time);
  BridgeInternals::GlobalDataAdaptor->SetDataTimeStep(*tstep);
  BridgeInternals::GlobalAnalysisAdaptor->Execute(BridgeInternals::GlobalDataAdaptor);
  BridgeInternals::GlobalDataAdaptor->ReleaseData();
  sensei::Profiler::EndEvent("nek::sensei_bridge::update");
}

//-----------------------------------------------------------------------------
void sensei_bridge_finalize()
{
  sensei::Profiler::StartEvent("nek::sensei_bridge::finalize");
  BridgeInternals::GlobalAnalysisAdaptor->Finalize();
  BridgeInternals::GlobalDataAdaptor = nullptr;
  BridgeInternals::GlobalAnalysisAdaptor = nullptr;
  sensei::Profiler::EndEvent("nek::sensei_bridge::finalize");
  sensei::Profiler::Finalize();
}
