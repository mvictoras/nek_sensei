#include "DataAdaptor.h"
#include "MeshMetadata.h"
#include "Error.h"

#include <vtkObjectFactory.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkSOADataArrayTemplate.h>

#include <cassert>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <fstream>

namespace nek5000{

  static
  vtkUnstructuredGrid *newUnstructuredBlock(int x_dim, int y_dim, int z_dim, int elems,
      double* mesh_x, double* mesh_y, double* mesh_z, bool structureOnly) {
    
    vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
    if(!structureOnly) {
      int arrayLen = x_dim * y_dim * z_dim * elems;
      vtkSOADataArrayTemplate<double> *pointsData = vtkSOADataArrayTemplate<double>::New();
      pointsData->SetNumberOfComponents(3);
      pointsData->SetArray(0, mesh_x, arrayLen,  true, true);
      pointsData->SetArray(1, mesh_y, arrayLen, false, true);
      pointsData->SetArray(2, mesh_z, arrayLen, false, true);
      
      //vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
      vtkPoints *points = vtkPoints::New();
      points->SetDataTypeToDouble();
      points->SetNumberOfPoints(arrayLen * 3);
      points->SetData(pointsData);

      ug->SetPoints(points);
      pointsData->Delete();

      // calculate cell ids based off of Nek5000's format
      vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
      int x_len = x_dim-1;
      int y_len = y_dim-1;
      int cells_type;
      if(z_dim == 1){
        int numCells = x_len*y_len*elems;
        cells_type = VTK_QUAD;
#ifdef VTK_CELL_ARRAY_V2
        vtkIdTypeArray *nlist = vtkIdTypeArray::New();
        nlist->SetNumberOfValues(numCells * 5);
        vtkIdType *nl = nlist->GetPointer(0);
        for(int elem=0; elem<elems; ++elem){
          for(int y=0; y < y_len; ++y){
            for(int x=0; x < x_len; ++x){
              nl[0] = 4;
              nl[1] = x + x_dim*(y + y_dim*elem);
              nl[2] = nl[1] + x_dim;
              nl[3] = nl[2] + 1;
              nl[4] = nl[1] + 1;
              nl += 5;
            }
          }
        }
        cells->SetCells(numCells, nlist);
        nlist->Delete();
#else
        // 2D mesh using quads
        vtkIdType* cellData = cells->WritePointer(numCells, numCells*5);
        for(int elem=0; elem<elems; ++elem){
          for(int y=0; y < y_len; ++y){
            for(int x=0; x < x_len; ++x){
              cellData[0] = 4;
              cellData[1] = x + x_dim*(y + y_dim*elem);
              cellData[2] = cellData[1] + x_dim;
              cellData[3] = cellData[2] + 1;
              cellData[4] = cellData[1] + 1;
              cellData += 5;
            }
          }
        }
#endif
      }else{
        // 3D mesh using hexahedrons
        int z_len = z_dim-1;
        cells_type = VTK_HEXAHEDRON;
        int numCells = x_len*y_len*z_len*elems;
#ifdef VTK_CELL_ARRAY_V2
        vtkIdTypeArray *nlist = vtkIdTypeArray::New();
        nlist->SetNumberOfValues(numCells * 9);
        vtkIdType *nl = nlist->GetPointer(0);
        for(int elem=0; elem<elems; ++elem){
          for(int z=0; z < z_len; ++z){
            for(int y=0; y < y_len; ++y){
              for(int x=0; x < x_len; ++x){
                nl[0] = 8;
                nl[1] = x + x_dim*(y + y_dim*(z + z_dim*elem));
                nl[2] = nl[1] + x_dim*y_dim;
                nl[3] = nl[2] + x_dim;
                nl[4] = nl[1] + x_dim;
                nl[5] = nl[1] + 1;
                nl[6] = nl[2] + 1;
                nl[7] = nl[3] + 1;
                nl[8] = nl[4] + 1;
                nl += 9;
              }
            }
          }
        }
        cells->SetCells(numCells, nlist);
        nlist->Delete();
#else
        vtkIdType* cellData = cells->WritePointer(numCells, numCells*9);
        for(int elem=0; elem<elems; ++elem){
          for(int z=0; z < z_len; ++z){
            for(int y=0; y < y_len; ++y){
              for(int x=0; x < x_len; ++x){
                cellData[0] = 8;
                cellData[1] = x + x_dim*(y + y_dim*(z + z_dim*elem));
                cellData[2] = cellData[1] + x_dim*y_dim;
                cellData[3] = cellData[2] + x_dim;
                cellData[4] = cellData[1] + x_dim;
                cellData[5] = cellData[1] + 1;
                cellData[6] = cellData[2] + 1;
                cellData[7] = cellData[3] + 1;
                cellData[8] = cellData[4] + 1;
                cellData += 9;
              }
            }
          }
        }
#endif
      }

      // set unstructured grid
      ug->SetCells(cells_type, cells);
    }
    return ug;
  }

  //-----------------------------------------------------------------------------
  senseiNewMacro(DataAdaptor);

  //-----------------------------------------------------------------------------
  void DataAdaptor::InitMesh(int index, int x_dim, int y_dim, int z_dim, int elems,
                double* mesh_x, double* mesh_y, double* mesh_z){
    // update last mesh
    if(meshLen <= index)
      meshLen = index + 1;

    // set points pointers
    meshes[index].mesh_x = mesh_x;
    meshes[index].mesh_y = mesh_y;
    meshes[index].mesh_z = mesh_z;
    meshes[index].x_dim = x_dim;
    meshes[index].y_dim = y_dim;
    meshes[index].z_dim = z_dim;
    meshes[index].elems = elems;

    // set length
    meshes[index].length = x_dim * y_dim * z_dim * elems;
  }

  //-----------------------------------------------------------------------------
  void DataAdaptor::Finalize(){
    int nRanks = 1;
    int rank = 0;

    MPI_Comm_rank(this->GetCommunicator(), &rank);
    MPI_Comm_size(this->GetCommunicator(), &nRanks);

    if(rank == 0) {
      std::string pvdFileName = "mesh.pvd";
      std::ofstream pvdFile(pvdFileName);

      if (!pvdFile)
      {
        SENSEI_ERROR("Failed to open " << pvdFileName << " for writing")
      }

      pvdFile << "<?xml version=\"1.0\"?>" << endl
        << "<VTKFile type=\"Collection\" version=\"0.1\""
           " byte_order=\"LittleEndian\" compressor=\"\">" << endl
        << "<Collection>" << endl;

      for (long i = 0; i < 2; ++i)
        {
        for(long k = 0; k < nRanks; ++k) {
          for (long j = 0; j < meshes[i].steps; ++j)
            {
              std::ostringstream oss;
              oss << "mesh_"
                << std::setw(6) << std::setfill('0') << k << "_"
                << std::setw(6) << std::setfill('0') << j << ".vtu";
              std::string fileName = oss.str();

            pvdFile << "<DataSet timestep=\"" << j
              << "\" group=\"\" part=\"\" file=\"" << fileName
              << "\"/>" << endl;
            }
          }
        }

      pvdFile << "</Collection>" << endl
        << "</VTKFile>" << endl;



    }

    meshes[0].mesh_x = NULL;
    meshes[0].mesh_y = NULL;
    meshes[0].mesh_z = NULL;
    meshes[0].mesh_x = NULL;
    meshes[0].mesh_y = NULL;
    meshes[1].mesh_z = NULL;
    arrays[0].clear();
    arrays[1].clear();
  }

  //-----------------------------------------------------------------------------
  void DataAdaptor::AddArray(int index, const std::string& name, vtkAbstractArray* data){
    // ensure that the first mesh is used when meshes are deduplicated
    arrays[index%meshLen][name] = data;
  }

  //-----------------------------------------------------------------------------
  void DataAdaptor::AddArray(int index, const std::string& name, double* data){
    vtkDoubleArray* vtkArray = vtkDoubleArray::New();
    vtkArray->SetName(name.c_str());
    vtkArray->SetArray(data, meshes[index%meshLen].length, true);
    // ensure that the first mesh is used when meshes are deduplicated
    arrays[index%meshLen][name] = vtkArray;
  }


  //-----------------------------------------------------------------------------
  int DataAdaptor::ValidateMeshName(const std::string &meshName){
    if(GetMeshIndex(meshName) < 0) {
      SENSEI_ERROR("No mesh \"" << meshName << "\"")
      return 0;
    }
    return 1;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::ValidateAssociation(int association){
    /* Not worrying about this now
    if(association != vtkDataObject::FIELD_ASSOCIATION_POINTS){
      SENSEI_ERROR("Only point data on mesh")
      return 0;
    }
    */
    return 1;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::GetMeshIndex(const std::string &meshName){
    if(meshLen == 1){
      if(meshName == "mesh")
        return 0;
    }else if(meshLen == 2){
      if(meshName == "vmesh")
        return 0;
      else if(meshName == "pmesh")
        return 1;
    }
    return -1;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes){
    numMeshes = meshLen;
    return 0;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) {
    if(id >= meshLen){
      if(meshLen == 1)
        SENSEI_ERROR("Mesh id is out of range. 1 mesh available.")
      else
        SENSEI_ERROR("Mesh id is out of range. 2 meshes available.")
      return -1;
    }

    int rank = 0;
    int nRanks = 1;

    MPI_Comm_rank(this->GetCommunicator(), &rank);
    MPI_Comm_size(this->GetCommunicator(), &nRanks);

    if(meshLen == 1){
      metadata->MeshName = "mesh";
    }else if(meshLen == 2){
      metadata->MeshName = (id==0? "vmesh": "pmesh");
    }
    
    metadata->MeshType = VTK_UNSTRUCTURED_GRID;
    metadata->CoordinateType = VTK_DOUBLE;
    metadata->NumArrays = 2;
    metadata->ArrayName = {"pressure", "velocity"};
    metadata->ArrayCentering = {vtkDataObject::CELL, vtkDataObject::CELL};
    metadata->ArrayComponents = {1, 3};
    metadata->ArrayType = {VTK_DOUBLE, VTK_DOUBLE};
    metadata->StaticMesh = 0;

    return 0;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::GetMesh(const std::string &meshName, bool structureOnly,
        vtkDataObject *&mesh){
    if(!ValidateMeshName(meshName)) return -1;

    int meshIndex = GetMeshIndex(meshName);

    mesh = newUnstructuredBlock(meshes[meshIndex].x_dim, meshes[meshIndex].y_dim, meshes[meshIndex].z_dim, meshes[meshIndex].elems, 
          meshes[meshIndex].mesh_x, meshes[meshIndex].mesh_y, meshes[meshIndex].mesh_z, structureOnly);

    meshes[meshIndex].steps++;

    return 0;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName,
        int association, const std::string& arrayName){
    if(!ValidateMeshName(meshName)) return -1;
    if(!ValidateAssociation(association)) return -1;

    int meshIndex = GetMeshIndex(meshName);
    ArraysType::iterator iter = arrays[meshIndex].find(arrayName);
    if(iter == arrays[meshIndex].end()){
      SENSEI_ERROR("No array named \"" << arrayName << "\"")
      return -1;
    }

    if(iter->second->GetNumberOfComponents() > 0) {
      // add vtkDoubleArray to mesh point data
      vtkUnstructuredGrid* grid = vtkUnstructuredGrid::SafeDownCast(mesh);
      grid->GetPointData()->AddArray(iter->second);
    }

    return 0;
  }

  //-----------------------------------------------------------------------------
  int DataAdaptor::ReleaseData(){
    return 0;
  }

}
