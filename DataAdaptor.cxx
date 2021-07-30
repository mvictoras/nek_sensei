#include "DataAdaptor.h"
#include "MeshMetadata.h"
#include "MPIUtils.h"
#include "VTKUtils.h"
#include "Error.h"

#include <vtkObjectFactory.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkSOADataArrayTemplate.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkUnstructuredGrid.h>

#include <sdiy/master.hpp>

#include <cassert>
#include <algorithm>
#include <ostream>
#include <sstream>
#include <fstream>

#include <Profiler.h>

typedef struct BlockType {
	int size;
	std::array<int, 3> dim;
	std::array<double*, 3> mesh;
	double *velocity_x;
	double *velocity_y;
	double *velocity_z;
	double *pressure;

  double vel_size;
} Block;


static
unsigned long getBlockNumCells(Block* block)
{

  if (block->dim[2] == 1) 
    return (block->dim[0] - 1) * (block->dim[1] - 1) * block->size;
  else 
    return (block->dim[0] - 1) * (block->dim[1] - 1) * (block->dim[2] - 1) * block->size;
}

static
long getBlockNumPoints(Block* block)
{
  return block->dim[0] * block->dim[1] * block->dim[2] * block->size * 3;
}

static
std::array<double,2> getArrayRange(unsigned long nSize, double* data)
{
	double bmin = std::numeric_limits<double>::max();
	double bmax = std::numeric_limits<double>::lowest();
	for (unsigned long i = 0; i < nSize; ++i)
	{
		bmin = std::min(bmin, data[i]);
		bmax = std::max(bmax, data[i]);
	}
	return {bmin,bmax};
}


/// @param gmin global min of this array
/// @param gmax global max of this array
static
std::array<double,2> getArrayRange(unsigned long nSize, double* data, std::array<double,2> &gRange)
{
	std::array<double,2> range = getArrayRange(nSize, data);	
	gRange[0] = std::min(gRange[0], range[0]);
	gRange[1] = std::max(gRange[1], range[1]);
	return range;
}

static
vtkUnstructuredGrid *newUnstructuredBlock(Block *block, bool structureOnly) 
{
  
  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
  if(!structureOnly) 
  {
    int arrayLen = block->dim[0] * block->dim[1] * block->dim[2] * block->size;
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    points->SetNumberOfPoints(arrayLen);
    double* mesh_x = block->mesh[0];
    double* mesh_y = block->mesh[1];
    double* mesh_z = block->mesh[2];
    //for(int i=0; i<arrayLen; ++i)
    //  points->SetPoint(i, mesh_x[i], mesh_y[i], mesh_z[i]);

    sensei::Profiler::StartEvent("DataAdaptor::newUnstructuredBlock::vtkSOAArrayTemplate");

    vtkSOADataArrayTemplate<double> *pointsData = vtkSOADataArrayTemplate<double>::New();
    pointsData->SetNumberOfComponents(3);
    pointsData->SetArray(0, block->mesh[0], arrayLen,  true, true);
    pointsData->SetArray(1, block->mesh[1], arrayLen, false, true);
    pointsData->SetArray(2, block->mesh[2], arrayLen, false, true);
    
    vtkPoints *points = vtkPoints::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(arrayLen * 3);
    points->SetData(pointsData);
    sensei::Profiler::EndEvent("DataAdaptor::newUnstructuredBlock::vtkSOAArrayTemplate");

    ug->SetPoints(points);
    //pointsData->Delete();

    // calculate cell ids based off of Nek5000's format
    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    int x_len = block->dim[1]-1;
    int y_len = block->dim[2]-1;
    int cells_type;
    if(block->dim[2] == 1)
    {
      int numCells = x_len*y_len*block->size;
      cells_type = VTK_QUAD;
#ifdef VTK_CELL_ARRAY_V2
      vtkIdTypeArray *nlist = vtkIdTypeArray::New();
      nlist->SetNumberOfValues(numCells * 5);
      vtkIdType *nl = nlist->GetPointer(0);
      for(int elem=0; elem<size; ++elem)
      {
        for(int y=0; y < y_len; ++y)
        {
          for(int x=0; x < x_len; ++x)
          {
            nl[0] = 4;
            nl[1] = x + block->dim[0]*(y + block->dim[1]*elem);
            nl[2] = nl[1] + block->dim[0];
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
      for(int elem=0; elem<block->size; ++elem)
      {
        for(int y=0; y < y_len; ++y)
        {
          for(int x=0; x < x_len; ++x)
          {
            cellData[0] = 4;
            cellData[1] = x + block->dim[0]*(y + block->dim[1]*elem);
            cellData[2] = cellData[1] + block->dim[0];
            cellData[3] = cellData[2] + 1;
            cellData[4] = cellData[1] + 1;
            cellData += 5;
          }
        }
      }
#endif
    }
    else
    {
      // 3D mesh using hexahedrons
      int z_len = block->dim[2]-1;
      cells_type = VTK_HEXAHEDRON;
      int numCells = x_len*y_len*z_len*block->size;
#ifdef VTK_CELL_ARRAY_V2
      vtkIdTypeArray *nlist = vtkIdTypeArray::New();
      nlist->SetNumberOfValues(numCells * 9);
      vtkIdType *nl = nlist->GetPointer(0);
      for(int elem=0; elem<size; ++elem)
      {
        for(int z=0; z < z_len; ++z)
        {
          for(int y=0; y < y_len; ++y)
          {
            for(int x=0; x < x_len; ++x)
            {
              nl[0] = 8;
              nl[1] = x + block->dim[0]*(y + block->dim[1]*(z + block->dim[2]*elem));
              nl[2] = nl[1] + block->dim[0]*block->dim[1];
              nl[3] = nl[2] + block->dim[0];
              nl[4] = nl[1] + block->dim[0];
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
      for(int elem=0; elem<block->size; ++elem)
      {
        for(int z=0; z < z_len; ++z)
        {
          for(int y=0; y < y_len; ++y)
          {
            for(int x=0; x < x_len; ++x)
            {
              cellData[0] = 8;
              cellData[1] = x + block->dim[0]*(y + block->dim[1]*(z + block->dim[2]*elem));
              cellData[2] = cellData[1] + block->dim[0]*block->dim[1];
              cellData[3] = cellData[2] + block->dim[0];
              cellData[4] = cellData[1] + block->dim[0];
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
static 
void getBounds(const sdiy::DiscreteBounds &db, double *ext) 
{
  ext[0] = db.min[0];
  ext[1] = db.max[0];
  ext[2] = db.min[1];
  ext[3] = db.max[1];
  ext[4] = db.min[2];
  ext[5] = db.max[2];
}

static
void getBlockExtent(const sdiy::DiscreteBounds &db, int *ext) 
{
  // converts from DIY layout to VTK
  ext[0] = db.min[0];
  ext[1] = db.max[0];
  ext[2] = db.min[1];
  ext[3] = db.max[1];
  ext[4] = db.min[2];
  ext[5] = db.max[2];
}

namespace nek_sensei
{


struct DataAdaptor::InternalsType
{
  InternalsType() : NumBlocks(0), Origin{0,0,0}, 
    Spacing{1,1,1}, NumGhostCells(0) {}

	long NumBlocks;                                   // total number of blocks on all ranks
  sdiy::DiscreteBounds DomainExtent;                // global index space
	sdiy::DiscreteBounds BlockExtents;                // local block extents, indexed by global block id
  sdiy::DiscreteBounds DomainBounds;                // global index space
	sdiy::DiscreteBounds BlockBounds;                 // local block extents, indexed by global block id
  Block* BlockData;                                 // local data array, indexed by block id

  double Origin[3];                                 // lower left corner of simulation domain
  double Spacing[3];                                // mesh spacing

  int NumGhostCells;                                // number of ghost cells
};

//-----------------------------------------------------------------------------
senseiNewMacro(DataAdaptor);

//-----------------------------------------------------------------------------
DataAdaptor::DataAdaptor() :
  Internals(new DataAdaptor::InternalsType())
{
}

//-----------------------------------------------------------------------------
DataAdaptor::~DataAdaptor()
{
  delete this->Internals;
}

//-----------------------------------------------------------------------------
void DataAdaptor::Initialize(int index, int x_dim, int y_dim, int z_dim, int elems,
              double* mesh_x, double* mesh_y, double* mesh_z,
              double *x_min, double* x_max, double* y_min, double* y_max, double* z_min, double* z_max,
              double* vel_x, double* vel_y, double* vel_z, double* pressure, int vel_size)
{

	int nRanks = 1;

  MPI_Comm_size(this->GetCommunicator(), &nRanks);

  this->Internals->NumBlocks = nRanks;

	std::array<double,2> x_range = getArrayRange(elems, mesh_x);
	std::array<double,2> y_range = getArrayRange(elems, mesh_y);
	std::array<double,2> z_range = getArrayRange(elems, mesh_z);

	this->SetDomainExtent(
      0, x_dim - 1, 
      0, y_dim - 1, 
      0, z_dim - 1);

  this->SetBlockExtent(
      0, x_dim - 1, 
      0, y_dim - 1,
      0, z_dim - 1);

  this->SetBlockBounds(
      x_range[0], x_range[1],
      y_range[0], y_range[1],
      z_range[0], z_range[1]);

  this->SetDomainBounds(
      *x_min, *x_max,
      *y_min, *y_max,
      *z_min, *z_max);

	this->Internals->BlockData = new Block();
	this->Internals->BlockData->size = elems;
	this->Internals->BlockData->dim[0] = x_dim;
	this->Internals->BlockData->dim[1] = y_dim;
	this->Internals->BlockData->dim[2] = z_dim;
	this->Internals->BlockData->mesh[0] = mesh_x;
	this->Internals->BlockData->mesh[1] = mesh_y;
	this->Internals->BlockData->mesh[2] = mesh_z;
  this->Internals->BlockData->velocity_x = vel_x;
  this->Internals->BlockData->velocity_y = vel_y;
  this->Internals->BlockData->velocity_z = vel_z;
  this->Internals->BlockData->vel_size = vel_size;
  this->Internals->BlockData->pressure = pressure;
}

//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------
void DataAdaptor::SetBlockExtent(int xmin, int xmax, int ymin,
   int ymax, int zmin, int zmax)
{
  this->Internals->BlockExtents.min[0] = xmin;
  this->Internals->BlockExtents.min[1] = ymin;
  this->Internals->BlockExtents.min[2] = zmin;

  this->Internals->BlockExtents.max[0] = xmax;
  this->Internals->BlockExtents.max[1] = ymax;
  this->Internals->BlockExtents.max[2] = zmax;
}

//-----------------------------------------------------------------------------
void DataAdaptor::SetDomainExtent(int xmin, int xmax, int ymin,
   int ymax, int zmin, int zmax)
{
  this->Internals->DomainExtent.min[0] = xmin;
  this->Internals->DomainExtent.min[1] = ymin;
  this->Internals->DomainExtent.min[2] = zmin;

  this->Internals->DomainExtent.max[0] = xmax;
  this->Internals->DomainExtent.max[1] = ymax;
  this->Internals->DomainExtent.max[2] = zmax;
}

//-----------------------------------------------------------------------------
void DataAdaptor::SetDomainBounds(double xmin, double xmax, double ymin,
   double ymax, double zmin, double zmax)
{
  this->Internals->DomainBounds.min[0] = xmin;
  this->Internals->DomainBounds.min[1] = ymin;
  this->Internals->DomainBounds.min[2] = zmin;

  this->Internals->DomainBounds.max[0] = xmax;
  this->Internals->DomainBounds.max[1] = ymax;
  this->Internals->DomainBounds.max[2] = zmax;
}

//-----------------------------------------------------------------------------
void DataAdaptor::SetBlockBounds(double xmin, double xmax, double ymin,
   double ymax, double zmin, double zmax)
{
  this->Internals->BlockBounds.min[0] = xmin;
  this->Internals->BlockBounds.min[1] = ymin;
  this->Internals->BlockBounds.min[2] = zmin;

  this->Internals->BlockBounds.max[0] = xmax;
  this->Internals->BlockBounds.max[1] = ymax;
  this->Internals->BlockBounds.max[2] = zmax;
}


//-----------------------------------------------------------------------------
int DataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes){
  numMeshes = 1;
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) 
{

  sensei::Profiler::StartEvent("DataAdaptor::GetMeshMetadata");
	if (id > 2)
	{
		SENSEI_ERROR("invalid mesh id " << id)
    return -1;
	}

  int rank = 0;
  int nRanks = 1;

  MPI_Comm_rank(this->GetCommunicator(), &rank);
  MPI_Comm_size(this->GetCommunicator(), &nRanks);

	// Three types of meshes: mesh, vmesh, pmesh
	// XXX - How to distinguish between them?
	metadata->MeshName = (id == 0 ? "mesh" : "vmesh");	

  /// There are 1 local blocks per rank
	int nBlocks = 1;
  metadata->MeshType = VTK_MULTIBLOCK_DATA_SET;
  metadata->BlockType = VTK_UNSTRUCTURED_GRID;
  metadata->CoordinateType = VTK_DOUBLE;
  metadata->NumBlocks = this->Internals->NumBlocks;
  metadata->NumBlocksLocal = {nBlocks};
  metadata->NumGhostCells = this->Internals->NumGhostCells;
  metadata->NumArrays = 2;
  metadata->ArrayName = {"pressure", "velocity"};
  metadata->ArrayCentering = {vtkDataObject::POINT, vtkDataObject::POINT};
  metadata->ArrayComponents = {1, 3};
  metadata->ArrayType = {VTK_DOUBLE, VTK_DOUBLE};
  metadata->StaticMesh = 0;
 
  if (metadata->Flags.BlockExtentsSet()) 
	{
    std::array<int,6> ext;
    getBlockExtent(this->Internals->DomainExtent, ext.data());
    metadata->Extent = std::move(ext);

    metadata->BlockExtents.reserve(nBlocks);

    getBlockExtent(this->Internals->BlockExtents, ext.data());
    metadata->BlockExtents.emplace_back(std::move(ext));
  }

  if (metadata->Flags.BlockBoundsSet()) 
	{
    std::array<double,6> bounds;
		getBounds(this->Internals->DomainBounds, bounds.data());
    metadata->Bounds = std::move(bounds);

    metadata->BlockBounds.reserve(nBlocks);
    
    getBounds(this->Internals->BlockBounds, bounds.data());
    metadata->BlockBounds.emplace_back(std::move(bounds));
  }

  if (metadata->Flags.BlockSizeSet()) 
	{
    long nCells = getBlockNumCells(this->Internals->BlockData);
    long nPts = getBlockNumPoints(this->Internals->BlockData);
    metadata->BlockNumCells.push_back(nCells);
    metadata->BlockNumPoints.push_back(nPts);
    metadata->BlockCellArraySize.push_back(9*nCells);
	}
	  
	if (metadata->Flags.BlockDecompSet()) {
    metadata->BlockOwner.push_back(rank);
    metadata->BlockIds.push_back(rank);
  }

  if (metadata->Flags.BlockArrayRangeSet())
	{
    std::array<double,2> gvRange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
		std::array<double,2> gpRange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
    
    unsigned long nCells = this->Internals->BlockData->vel_size;
    std::array<double,2> vel_x_range = getArrayRange(nCells, this->Internals->BlockData->velocity_x, gvRange);
    std::array<double,2> vel_y_range = getArrayRange(nCells, this->Internals->BlockData->velocity_y, gvRange);
    std::array<double,2> vel_z_range = getArrayRange(nCells, this->Internals->BlockData->velocity_z, gvRange);
    std::array<double,2> pr_range = getArrayRange(nCells, this->Internals->BlockData->pressure, gpRange);
    std::array<double, 2> vel_range = {std::min({vel_x_range[0], vel_y_range[0], vel_z_range[0]}), std::max({vel_x_range[1], vel_y_range[1], vel_z_range[1]})};
    metadata->BlockArrayRange.push_back({pr_range, vel_range});
		metadata->ArrayRange.push_back(gpRange); 
		metadata->ArrayRange.push_back(gvRange); 
	}
	sensei::Profiler::EndEvent("DataAdaptor::GetMeshMetadata");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::GetMesh(const std::string &meshName, bool structureOnly,
	vtkDataObject *&mesh)
{
  sensei::Profiler::StartEvent("DataAdaptor::GetMesh");

	mesh = nullptr;

  if ((meshName != "mesh") && (meshName != "pmesh") && (meshName != "vmesh"))
	{
    SENSEI_ERROR("nek5000 provides meshes named \"mesh\", \"pmesh\","
      " and \"vmesh\". you requested \"" << meshName << "\"")
    return -1;
	}
  int rank = 0;

  MPI_Comm_rank(this->GetCommunicator(), &rank);

	vtkMultiBlockDataSet *mb = vtkMultiBlockDataSet::New();
  mb->SetNumberOfBlocks(this->Internals->NumBlocks);

  vtkUnstructuredGrid *ug = newUnstructuredBlock(this->Internals->BlockData, structureOnly);
  mb->SetBlock(rank, ug);
  ug->Delete();

	mesh = mb;
  sensei::Profiler::EndEvent("DataAdaptor::GetMesh");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName,
      int association, const std::string& arrayName)
{
  sensei::Profiler::StartEvent("DataAdaptor::AddArray");

	vtkMultiBlockDataSet *mb = dynamic_cast<vtkMultiBlockDataSet*>(mesh);
  if (!mb)
	{
		SENSEI_ERROR("unexpected mesh type "
      << (mesh ? mesh->GetClassName() : "nullptr"))
    return -1;
	}

	int rank = 0;

  MPI_Comm_rank(this->GetCommunicator(), &rank);

	vtkDataObject *blk = mb->GetBlock(rank);
	if (!blk)
	{
		SENSEI_ERROR("encountered empty block at index " << rank)
		return -1;
	}
       
	vtkDoubleArray *da = nullptr;
	vtkDataSetAttributes *dsa = nullptr;

	if( arrayName == "pressure" ) 
	{
		dsa = blk->GetAttributes(vtkDataObject::POINT);
		da = vtkDoubleArray::New();
		da->SetName("pressure");
		da->SetArray(this->Internals->BlockData->pressure, this->Internals->BlockData->vel_size, 1);	
	}
	
	if (dsa && da) {
		dsa->AddArray(da);
		da->Delete();
	}
	sensei::Profiler::EndEvent("DataAdaptor::AddArray");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::ReleaseData(){
  return 0;
}
}
