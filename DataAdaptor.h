#ifndef NEK_DATAADAPTOR_H
#define NEK_DATAADAPTOR_H

#include <DataAdaptor.h>
#include <senseiConfig.h>

#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>

#include <map>
#include <string>

#include <sdiy/master.hpp>

namespace nek5000{

struct NekMesh{
  double *mesh_x, *mesh_y, *mesh_z;
  int length;
  int meshLen, x_dim, y_dim, z_dim, elems;
  vtkDataSet *mesh;
  int steps = 0;
  sdiy::DiscreteBounds DomainExtent;
};

class DataAdaptor : public sensei::DataAdaptor {
public:
  static DataAdaptor* New();
  senseiTypeMacro(DataAdaptor, sensei::DataAdaptor);

  /// Initialize a mesh
  void InitMesh(int index, int x_dim, int y_dim, int z_dim, int elems,
      double* mesh_x, double* mesh_y, double* mesh_z);

  /// Set the pointers to simulation memory.
  void AddArray(int index, const std::string& name, vtkAbstractArray* data);
  void AddArray(int index, const std::string& name, double* data);

  /// Finalize the data adaptor.
  void Finalize();

  // SENSEI API
  int GetNumberOfMeshes(unsigned int &numMeshes) override;

  int GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &md) override;

  int GetMesh(const std::string &meshName, bool structureOnly,
    vtkDataObject *&mesh) override;

  int AddArray(vtkDataObject* mesh, const std::string &meshName,
    int association, const std::string &arrayName) override;

  int ReleaseData() override;

protected:
  DataAdaptor(): meshLen(0){}
  ~DataAdaptor(){}

  int ValidateMeshName(const std::string &meshName);
  int ValidateAssociation(int association);
  int GetMeshIndex(const std::string &meshName);

  using ArraysType = std::map<std::string, vtkAbstractArray*>;
  ArraysType arrays[2];

  NekMesh meshes[2];
  int meshLen;

private:
  DataAdaptor(const DataAdaptor&) = delete;
  void operator=(const DataAdaptor&) = delete;
};

}

#endif
