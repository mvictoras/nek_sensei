#ifndef NEK_DATAADAPTOR_H
#define NEK_DATAADAPTOR_H

#include <DataAdaptor.h>

namespace nek_sensei
{

class DataAdaptor : public sensei::DataAdaptor 
{
public:
  static DataAdaptor* New();
  senseiTypeMacro(DataAdaptor, sensei::DataAdaptor);

  // @brief Initialize the data adaptor.
  ///
  /// This initializes the data adaptor. This must be called once per simulation run.
  /// @param nblocks is the total number of blocks in the simulation run.
  void Initialize(int index, int x_dim, int y_dim, int z_dim, int elems,
      double* mesh_x, double* mesh_y, double* mesh_z,
      double* x_min, double* x_max, double* y_min, double* y_max, double* z_min, double* z_max,
      double* vel_x, double* vel_y, double* vel_z, double* pressure, int vel_size);

  /// Set the extents for local blocks.
  void SetBlockExtent(int xmin, int xmax, int ymin,
      int ymax, int zmin, int zmax);
	
	void SetDomainExtent(int xmin, int xmax, int ymin,
   int ymax, int zmin, int zmax);

  /// Set the bounds for local blocks.
  void SetBlockBounds(double xmin, double xmax, double ymin,
      double ymax, double zmin, double zmax);
	
	void SetDomainBounds(double xmin, double xmax, double ymin,
   double ymax, double zmin, double zmax);


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
  DataAdaptor();
  ~DataAdaptor();

private:
  DataAdaptor(const DataAdaptor&) = delete;
  void operator=(const DataAdaptor&) = delete;

  struct InternalsType;
  InternalsType *Internals;
};

}

#endif
