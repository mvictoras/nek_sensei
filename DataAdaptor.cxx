#include "DataAdaptor.h"
#include "MeshMetadata.h"
#include "MPIUtils.h"
#include "Error.h"

#include <svtkObjectFactory.h>
#include <svtkDoubleArray.h>
#include <svtkPointData.h>
#include <svtkCellArray.h>
#include <svtkSOADataArrayTemplate.h>
#include <svtkMultiBlockDataSet.h>
#include <svtkUnsignedCharArray.h>
#include <svtkUnstructuredGrid.h>

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
  double *vorticity_x;
  double *vorticity_y;
  double *vorticity_z;
	double *pressure;
  double *temperature;
  double *jacobian;

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
  return block->dim[0] * block->dim[1] * block->dim[2] * block->size;
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
svtkUnstructuredGrid *newUnstructuredBlock(Block *block, bool structureOnly) 
{
  svtkUnstructuredGrid *ug = svtkUnstructuredGrid::New();
  if(!structureOnly) 
  {
    int arrayLen = block->dim[0] * block->dim[1] * block->dim[2] * block->size;
    sensei::Profiler::StartEvent("nek::DataAdaptor::newUnstructuredBlock::svtkSOAArrayTemplate");

    svtkSOADataArrayTemplate<double> *pointsData = svtkSOADataArrayTemplate<double>::New();
    pointsData->SetNumberOfComponents(3);
    pointsData->SetArray(0, block->mesh[0], arrayLen,  true, true);
    pointsData->SetArray(1, block->mesh[1], arrayLen, false, true);
    pointsData->SetArray(2, block->mesh[2], arrayLen, false, true);
    
    svtkPoints *points = svtkPoints::New();
    points->SetDataTypeToDouble();
    points->SetNumberOfPoints(arrayLen * 3);
    points->SetData(pointsData);
    sensei::Profiler::EndEvent("nek::DataAdaptor::newUnstructuredBlock::svtkSOAArrayTemplate");

    ug->SetPoints(points);
    points->Delete();
    pointsData->Delete();

    // calculate cell ids based off of Nek5000's format
    int x_len = block->dim[1]-1;
    int y_len = block->dim[2]-1;
    if(block->dim[2] == 1)
    {
      svtkIdType ncells = x_len*y_len*block->size;
      
      svtkIdTypeArray *nlist = svtkIdTypeArray::New();
      nlist->SetNumberOfValues(ncells * 4);

      svtkUnsignedCharArray *cellTypes = svtkUnsignedCharArray::New();
      cellTypes->SetNumberOfValues(ncells);

      svtkIdTypeArray *cellLocations = svtkIdTypeArray::New();
      cellLocations->SetNumberOfValues(ncells + 1);

      svtkIdType *nl = nlist->GetPointer(0);
      unsigned char *ct = cellTypes->GetPointer(0);
      svtkIdType *cl = cellLocations->GetPointer(0);
      int offset = 0;
      for(int elem=0; elem<block->size; ++elem)
        for(int y=0; y < y_len; ++y)
            for(int x=0; x < x_len; ++x)
          {
           
            *ct++ = SVTK_QUAD;

            *cl++ = offset;
            offset += 4;

            nl[0] = x + block->dim[0]*(y + block->dim[1]*elem);
            nl[1] = nl[0] + block->dim[0];
            nl[2] = nl[1] + 1;
            nl[3] = nl[0] + 1;

            nl += 4;
          }

      *cl = offset;

      svtkCellArray *cells = svtkCellArray::New();
      cells->SetData(cellLocations, nlist);

      ug->SetCells(cellTypes, cells);

      nlist->Delete();
      cellTypes->Delete();
      cellLocations->Delete();
      cells->Delete();
    }
    else
    {
      // 3D mesh using hexahedrons
      int z_len = block->dim[2]-1;

      svtkIdType ncells = x_len*y_len*z_len*block->size;

      svtkIdTypeArray *nlist = svtkIdTypeArray::New();
      nlist->SetNumberOfValues(ncells * 8);

      svtkUnsignedCharArray *cellTypes = svtkUnsignedCharArray::New();
      cellTypes->SetNumberOfValues(ncells);

      svtkIdTypeArray *cellLocations = svtkIdTypeArray::New();
      cellLocations->SetNumberOfValues(ncells + 1);

      svtkIdType *nl = nlist->GetPointer(0);
      unsigned char *ct = cellTypes->GetPointer(0);
      svtkIdType *cl = cellLocations->GetPointer(0);
      int offset = 0;
      for(int elem=0; elem<block->size; ++elem)
        for(int z=0; z < z_len; ++z)
          for(int y=0; y < y_len; ++y)
            for(int x=0; x < x_len; ++x)
            {

              *ct++ = SVTK_HEXAHEDRON;

              *cl++ = offset;
              offset += 8;

              nl[0] = x + block->dim[0]*(y + block->dim[1]*(z + block->dim[2]*elem));
              nl[1] = nl[0] + block->dim[0]*block->dim[1];
              nl[2] = nl[1] + block->dim[0];
              nl[3] = nl[0] + block->dim[0];
              nl[4] = nl[0] + 1;
              nl[5] = nl[1] + 1;
              nl[6] = nl[2] + 1;
              nl[7] = nl[3] + 1;

              nl += 8;
            }
      *cl = offset;
      
      svtkCellArray *cells = svtkCellArray::New();
      cells->SetData(cellLocations, nlist);

      ug->SetCells(cellTypes, cells);

      nlist->Delete();
      cellTypes->Delete();
      cellLocations->Delete();
      cells->Delete();
    }
    
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
  sdiy::DiscreteBounds VelocityArrayRange;          // Range of velocity values
  sdiy::DiscreteBounds VorticityArrayRange;         // Range of vorticity values
  sdiy::DiscreteBounds PressureArrayRange;          // Range of pressure values
  sdiy::DiscreteBounds TemperatureArrayRange;       // Range of temperature values
  sdiy::DiscreteBounds JacobianArrayRange;          // Range of jacobian values
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
  delete this->Internals->BlockData;
  delete this->Internals;
}

//-----------------------------------------------------------------------------
void DataAdaptor::Initialize(int index, int x_dim, int y_dim, int z_dim, int elems,
              double* mesh_x, double* mesh_y, double* mesh_z,
              double* vel_x, double* vel_y, double* vel_z, 
              double* vort_x, double* vort_y, double* vort_z,
              double* pressure, double *temp, double *jacobian,
              double *x_min, double* x_max, double* y_min, double* y_max, double* z_min, double* z_max,
              double* vel_x_min, double* vel_x_max, double* vel_y_min, double* vel_y_max, double* vel_z_min, double* vel_z_max,
              double* vort_x_min, double* vort_x_max, double* vort_y_min, double* vort_y_max, double* vort_z_min, double* vort_z_max,
              double* pr_min, double* pr_max, double* temp_min, double* temp_max, double* jac_min, double* jac_max, int vel_size)

{

	int nRanks = 1;
  int i =0 ;
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

  this->SetArrayRange(
      this->Internals->VelocityArrayRange,
      *vel_x_min, *vel_x_max,
      *vel_y_min, *vel_y_max,
      *vel_z_min, *vel_z_max
      );
   
  this->SetArrayRange(
      this->Internals->VorticityArrayRange,
      *vel_x_min, *vel_x_max,
      *vel_y_min, *vel_y_max,
      *vel_z_min, *vel_z_max
      );

  //this->SetArrayRange(this->Internals->PressureArrayRange, *pr_min, *pr_max );
  //this->SetArrayRange(this->Internals->TemperatureArrayRange, *temp_min, *temp_max );
  //this->SetArrayRange(this->Internals->JacobianArrayRange, *jac_min, *jac_max );

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
  this->Internals->BlockData->vorticity_x = vort_x;
  this->Internals->BlockData->vorticity_y = vort_y;
  this->Internals->BlockData->vorticity_z = vort_z;
  this->Internals->BlockData->vel_size = vel_size;
  this->Internals->BlockData->pressure = pressure;
  this->Internals->BlockData->temperature = temp;
  this->Internals->BlockData->jacobian = jacobian;
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
void DataAdaptor::SetArrayRange(sdiy::DiscreteBounds &bounds, double xmin, double xmax, double ymin,
   double ymax, double zmin, double zmax)
{
  bounds.min[0] = xmin;
  bounds.min[1] = ymin;
  bounds.min[2] = zmin;

  bounds.max[0] = xmax;
  bounds.max[1] = ymax;
  bounds.max[2] = zmax;
}

//-----------------------------------------------------------------------------
void DataAdaptor::SetArrayRange(sdiy::DiscreteBounds &bounds, double min, double max) {
  bounds.min[0] = min;
  bounds.max[0] = max;
}

//-----------------------------------------------------------------------------
int DataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes){
  numMeshes = 1;
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::GetMeshMetadata(unsigned int id, sensei::MeshMetadataPtr &metadata) 
{
  sensei::Profiler::StartEvent("nek::DataAdaptor::GetMeshMetadata");
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
  metadata->MeshType = SVTK_MULTIBLOCK_DATA_SET;
  metadata->BlockType = SVTK_UNSTRUCTURED_GRID;
  metadata->CoordinateType = SVTK_DOUBLE;
  metadata->NumBlocks = this->Internals->NumBlocks;
  metadata->NumBlocksLocal = {nBlocks};
  metadata->NumGhostCells = this->Internals->NumGhostCells;
  metadata->NumArrays = 9;
  metadata->ArrayName = {"velocity_x", "velocity_y", "velocity_z", "vorticity_x", "vorticity_y", "vorticity_z", "pressure", "temperature", "jacobian"};
  metadata->ArrayCentering = {svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT, svtkDataObject::POINT};
  metadata->ArrayComponents = {1, 1, 1, 1, 1, 1, 1, 1, 1};
  metadata->ArrayType = {SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE, SVTK_DOUBLE};
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

  //if (metadata->Flags.BlockSizeSet()) 
	//{
    long nCells = getBlockNumCells(this->Internals->BlockData);
    long nPts = getBlockNumPoints(this->Internals->BlockData);
    metadata->BlockNumCells.push_back(nCells);
    metadata->BlockNumPoints.push_back(nPts);
    metadata->BlockCellArraySize.push_back(8*nCells);
	//}
	  
	if (metadata->Flags.BlockDecompSet()) {
    metadata->BlockOwner.push_back(rank);
    metadata->BlockIds.push_back(rank);
  }

  if (metadata->Flags.BlockArrayRangeSet())
	{
    std::array<double,2> gvRange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
		std::array<double,2> gpRange = {std::numeric_limits<double>::max(), std::numeric_limits<double>::lowest()};
    
    unsigned long nCells = this->Internals->BlockData->vel_size;
    std::array<double,2> vel_x_BlockRange = getArrayRange(nCells, this->Internals->BlockData->velocity_x, gvRange);
    std::array<double,2> vel_y_BlockRange = getArrayRange(nCells, this->Internals->BlockData->velocity_y, gvRange);
    std::array<double,2> vel_z_BlockRange = getArrayRange(nCells, this->Internals->BlockData->velocity_z, gvRange);
    std::array<double,2> vort_x_BlockRange = getArrayRange(nCells, this->Internals->BlockData->vorticity_x, gvRange);
    std::array<double,2> vort_y_BlockRange = getArrayRange(nCells, this->Internals->BlockData->vorticity_y, gvRange);
    std::array<double,2> vort_z_BlockRange = getArrayRange(nCells, this->Internals->BlockData->vorticity_z, gvRange);
    std::array<double,2> pr_BlockRange = getArrayRange(nCells, this->Internals->BlockData->pressure, gpRange);
    std::array<double,2> temp_BlockRange = getArrayRange(nCells, this->Internals->BlockData->temperature, gpRange);
    std::array<double,2> jac_BlockRange = getArrayRange(nCells, this->Internals->BlockData->jacobian, gpRange);
    metadata->BlockArrayRange.push_back({vel_x_BlockRange, vel_y_BlockRange, vel_z_BlockRange, 
        vort_x_BlockRange, vort_y_BlockRange, vort_z_BlockRange, 
        pr_BlockRange, temp_BlockRange, jac_BlockRange});

		std::array<double,2> vel_x_range = { this->Internals->VelocityArrayRange.min[0], this->Internals->VelocityArrayRange.max[0] };
    std::array<double,2> vel_y_range = { this->Internals->VelocityArrayRange.min[1], this->Internals->VelocityArrayRange.max[1] };
    std::array<double,2> vel_z_range = { this->Internals->VelocityArrayRange.min[2], this->Internals->VelocityArrayRange.max[2] };
    std::array<double,2> vort_x_range = { this->Internals->VorticityArrayRange.min[0], this->Internals->VorticityArrayRange.max[0] };
    std::array<double,2> vort_y_range = { this->Internals->VorticityArrayRange.min[1], this->Internals->VorticityArrayRange.max[1] };
    std::array<double,2> vort_z_range = { this->Internals->VorticityArrayRange.min[2], this->Internals->VorticityArrayRange.max[2] };
    std::array<double,2> pr_range = { this->Internals->PressureArrayRange.min[0], this->Internals->PressureArrayRange.max[0] };
    std::array<double,2> temp_range = { this->Internals->TemperatureArrayRange.min[0], this->Internals->TemperatureArrayRange.max[0] };
    std::array<double,2> jac_range = { this->Internals->JacobianArrayRange.min[0], this->Internals->JacobianArrayRange.max[0] };
    metadata->ArrayRange.push_back(vel_x_range);
    metadata->ArrayRange.push_back(vel_y_range);
    metadata->ArrayRange.push_back(vel_z_range);
    metadata->ArrayRange.push_back(vort_x_range);
    metadata->ArrayRange.push_back(vort_y_range);
    metadata->ArrayRange.push_back(vort_z_range);
    metadata->ArrayRange.push_back(pr_range);
    metadata->ArrayRange.push_back(temp_range);
    metadata->ArrayRange.push_back(jac_range); 
	  
  }
	sensei::Profiler::EndEvent("nek::DataAdaptor::GetMeshMetadata");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::GetMesh(const std::string &meshName, bool structureOnly,
	svtkDataObject *&mesh)
{

  sensei::Profiler::StartEvent("nek::DataAdaptor::GetMesh");
	mesh = nullptr;

  if ((meshName != "mesh") && (meshName != "pmesh") && (meshName != "vmesh"))
	{
    SENSEI_ERROR("nek5000 provides meshes named \"mesh\", \"pmesh\","
      " and \"vmesh\". you requested \"" << meshName << "\"")
    return -1;
	}
  int rank = 0;

  MPI_Comm_rank(this->GetCommunicator(), &rank);

	svtkMultiBlockDataSet *mb = svtkMultiBlockDataSet::New();
  mb->SetNumberOfBlocks(this->Internals->NumBlocks);

  svtkUnstructuredGrid *ug = newUnstructuredBlock(this->Internals->BlockData, structureOnly);
  mb->SetBlock(rank, ug);
  ug->Delete();

	mesh = mb;
  sensei::Profiler::EndEvent("nek::DataAdaptor::GetMesh");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::AddArray(svtkDataObject* mesh, const std::string &meshName,
      int association, const std::string& arrayName)
{
  sensei::Profiler::StartEvent("nek::DataAdaptor::AddArray");
	svtkMultiBlockDataSet *mb = dynamic_cast<svtkMultiBlockDataSet*>(mesh);
  if (!mb)
	{
		SENSEI_ERROR("unexpected mesh type "
      << (mesh ? mesh->GetClassName() : "nullptr"))
    return -1;
	}

	int rank = 0;

  MPI_Comm_rank(this->GetCommunicator(), &rank);

	svtkDataObject *blk = mb->GetBlock(rank);
	if (!blk)
	{
		SENSEI_ERROR("encountered empty block at index " << rank)
		return -1;
	}
       
	svtkDoubleArray *da = nullptr;
	svtkDataSetAttributes *dsa = nullptr;

	if( arrayName == "velocity_x" && this->Internals->BlockData->velocity_x != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("velocity_x");
		da->SetArray(this->Internals->BlockData->velocity_x, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "velocity_y" && this->Internals->BlockData->velocity_y != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("velocity_y");
		da->SetArray(this->Internals->BlockData->velocity_y, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "velocity_z" && this->Internals->BlockData->velocity_z != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("velocity_z");
		da->SetArray(this->Internals->BlockData->velocity_z, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "vorticity_x" && this->Internals->BlockData->vorticity_x != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("vorticity_x");
		da->SetArray(this->Internals->BlockData->vorticity_x, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "vorticity_y" && this->Internals->BlockData->vorticity_y != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("vorticity_y");
		da->SetArray(this->Internals->BlockData->vorticity_y, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "vorticity_z" && this->Internals->BlockData->vorticity_z != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("vorticity_z");
		da->SetArray(this->Internals->BlockData->vorticity_z, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "pressure" && this->Internals->BlockData->pressure != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("pressure");
		da->SetArray(this->Internals->BlockData->pressure, this->Internals->BlockData->vel_size, 1);	
	}
  else if( arrayName == "temperature" && this->Internals->BlockData->temperature != nullptr ) 
	{
		dsa = blk->GetAttributes(svtkDataObject::POINT);
		da = svtkDoubleArray::New();
		da->SetName("temperature");
		da->SetArray(this->Internals->BlockData->temperature, this->Internals->BlockData->vel_size, 1);	
	} /*
  else if( arrayName == "jacobian" && this->Internals->BlockData->velocity_y != nullptr ) 
	{
		//dsa = blk->GetAttributes(svtkDataObject::POINT);
		//da = svtkDoubleArray::New();
		//da->SetName("jacobian");
		//da->SetArray(this->Internals->BlockData->pressure, this->Internals->BlockData->vel_size, 1);	
	}
  */
	if (dsa && da) {
		dsa->AddArray(da);
		da->Delete();
	}
	sensei::Profiler::EndEvent("nek::DataAdaptor::AddArray");
  return 0;
}

//-----------------------------------------------------------------------------
int DataAdaptor::ReleaseData(){
  return 0;
}

}
