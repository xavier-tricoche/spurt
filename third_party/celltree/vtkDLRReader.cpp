#include <vtkDataArraySelection.h>
#include <vtkErrorCode.h>
#include <vtkUnstructuredGrid.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkObjectFactory.h>
#include <vtkFieldData.h>
#include <vtkPointData.h>
#include <vtkByteSwap.h>
#include <vtkFloatArray.h>
#include <vtkIntArray.h>
#include <vtkCellType.h>
#include <vtkCellArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkStdString.h>

#include <netcdf.h>
#include <libgen.h>
#include <map>

#include "vtkDLRReader.hpp"

vtkCxxRevisionMacro(vtkDLRReader, "$Revision: 1.9 $");
vtkStandardNewMacro(vtkDLRReader);

// -------------------------------------------------------------------------

class nc_exception 
{
};

// -------------------------------------------------------------------------

static void nccheck( int result )
{
    if( result != NC_NOERR )
    {
        // vtkErrorMacro( "Could not read from NetCDF file\n" );
        throw nc_exception();
    }
}

// -------------------------------------------------------------------------

vtkDLRReader::vtkDLRReader()
{
    this->SetNumberOfInputPorts(0);
    this->FileName          = NULL;
    this->NumberOfVariables = 0;
    this->VariableSelection = vtkDataArraySelection::New();
}

//----------------------------------------------------------------------------

vtkDLRReader::~vtkDLRReader()
{
    if (this->FileName)
    {
        delete [] this->FileName;
    }
    
    this->VariableSelection->Delete();
    
    delete[] this->ComponentNumber;
    delete[] this->VariableName;
}

//----------------------------------------------------------------------------

void vtkDLRReader::PrintSelf( ostream& os, vtkIndent indent )
{
    this->Superclass::PrintSelf(os,indent);

    os << indent << "File Name: " 
       << (this->FileName ? this->FileName : "(none)") << endl;

    os << indent << "Number Of Variables: " << this->NumberOfVariables << endl;
    
    for (int i=0; i < this->NumberOfVariables; i++)
    {
        os << "\tVariableName[" << i << "] = " 
           << this->VariableName[i] << endl;
        os << "\tComponentNumber[" << i << "] = " 
           << this->ComponentNumber[i] << endl;
        os << "\tVariableSelection->GetArraySetting(" << i << ") = " 
           << (this->VariableSelection->GetArraySetting(i) ? "ENABLED" : "DISABLED") << endl;
        os << endl;
    }
}

//----------------------------------------------------------------------------

int vtkDLRReader::RequestInformation( vtkInformation *vtkNotUsed(request),
                                      vtkInformationVector **vtkNotUsed(inputVector),
                                      vtkInformationVector *vtkNotUsed(outputVector))
{
    if( !this->FileName )
    {
        vtkErrorMacro("No filename specified");
        return 0;
    }

    std::ifstream in( this->FileName );
    in >> GridFileName >> PvalFileName;

    if( !in )
    {
        this->SetErrorCode(vtkErrorCode::UnrecognizedFileTypeError);
        vtkErrorMacro( "Unable to read DLR file" );
        return 0;
    }

    in.close();

    vtkStdString path;
    
    {
        char* tmp = strdup( this->FileName );
        path = dirname( tmp );
        free( tmp );
    }

    GridFileName = path + '/' + GridFileName;
    PvalFileName = path + '/' + PvalFileName;

    int ncid, result;
    
    if( NC_NOERR != (result = nc_open( PvalFileName, NC_NOWRITE, &ncid ) ) )
    {
		std::cerr << "Grid file name = " << GridFileName << '\n';
		std::cerr << "Pval file name = " << PvalFileName << '\n';
		std::cerr << "result = " << result << '\n';
        this->SetErrorCode( vtkErrorCode::CannotOpenFileError );
        vtkErrorMacro( "Unable to open NetCDF file" );
        return 0;
    }
                      
    int nvars = 0;                                  
    nc_inq_nvars( ncid, &nvars );
    
    std::map< vtkStdString, int > vars;
    
    for( int varid=0; varid<nvars; ++varid )
    {
        int realdim = 1;
        char name[NC_MAX_NAME+1], *realname = name;
        nc_inq_varname( ncid, varid, name );
        
        if( !strncmp( "y_", name, 2 ) || !strncmp( "z_", name, 2 ) )
            continue;
                
        if( !strncmp( "x_", name, 2 ) )
        {
            realdim   = 3;
            realname += 2;
        }
        
        vars.insert( std::make_pair( realname, realdim ) );
    }

    nc_close( ncid );

    // Fields associated with each particle point: velocity, mass, tag
    NumberOfVariables = vars.size();

    VariableName    = new vtkStdString[NumberOfVariables];
    ComponentNumber = new vtkIdType[NumberOfVariables];

    std::map< vtkStdString, int >::iterator vi = vars.begin();
    
    for( unsigned int i=0; i<vars.size(); ++i, ++vi )
    {
        VariableName[i] = vi->first;
        ComponentNumber[i] = vi->second;

        VariableSelection->AddArray( VariableName[i].c_str() );
    }
                                                                                
    return 1;
}

//----------------------------------------------------------------------------
int vtkDLRReader::RequestData( vtkInformation *vtkNotUsed(request),
                               vtkInformationVector **vtkNotUsed(inputVector),
                               vtkInformationVector *outputVector )
{
    // get the info object
    vtkInformation *outInfo = outputVector->GetInformationObject(0);
                                                                            
    // get the output
    vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));
                                                                            
    vtkDebugMacro( << "Reading DLR file");
                                                                            
    // If RequestInformation() failed the FileStream will be NULL
    // if( this->FileStream == NULL )
    // {
    //     return 0;
    // }
                                                                            
    // Read the file into the output unstructured grid
    this->ReadFile(output);
    return 1;
}

//----------------------------------------------------------------------------
void vtkDLRReader::ReadFile( vtkUnstructuredGrid *output )
{
    this->SetErrorCode( vtkErrorCode::NoError );

    int ncid, result;
    
    if( NC_NOERR != (result = nc_open( GridFileName, NC_NOWRITE, &ncid ) ) )
    {
        this->SetErrorCode( vtkErrorCode::CannotOpenFileError );
        vtkErrorMacro( "Specified filename not found" );
        return;
    }

    ReadPoints( ncid, output );
    ReadCells( ncid, output );
    
    nc_close( ncid );

    if( NC_NOERR != (result = nc_open( PvalFileName, NC_NOWRITE, &ncid ) ) )
    {
        this->SetErrorCode( vtkErrorCode::CannotOpenFileError );
        vtkErrorMacro( "Specified filename not found" );
        return;
    }

    ReadVelocity( ncid, output );
    
#if 0    
    
    vtkFloatArray *velocity = vtkFloatArray::New();
  vtkFloatArray *mass     = vtkFloatArray::New();
  vtkIntArray *tag        = vtkIntArray::New();

  // Allocate space in the unstructured grid for all nodes
  output->Allocate(this->NumberOfNodes, this->NumberOfNodes);
  output->SetPoints(points);

  // Allocate velocity array if requested, add to point and cell data
  if (this->PointDataArraySelection->GetArraySetting(USE_VELOCITY))
    {
    velocity->SetName("velocity");
    velocity->SetNumberOfComponents(DIMENSION);
    velocity->SetNumberOfTuples(this->NumberOfNodes);
    output->GetPointData()->AddArray(velocity);
    if (!output->GetPointData()->GetVectors())
      {
      output->GetPointData()->SetVectors(velocity);
      }
    }
  
  // Allocate mass array if requested, add to point and cell data
  if (this->PointDataArraySelection->GetArraySetting(USE_MASS))
    {
    mass->SetName("mass");
    mass->SetNumberOfComponents(1);
    mass->SetNumberOfTuples(this->NumberOfNodes);
    output->GetPointData()->AddArray(mass);
    if (!output->GetPointData()->GetScalars())
      {
      output->GetPointData()->SetScalars(mass);
      }
    }

  // Allocate tag array if requested, add to point and cell data
  if (this->PointDataArraySelection->GetArraySetting(USE_TAG))
    {
    tag->SetName("tag");
    tag->SetNumberOfComponents(1);
    tag->SetNumberOfTuples(this->NumberOfNodes);
    output->GetPointData()->AddArray(tag);
    if (!output->GetPointData()->GetScalars())
      {
      output->GetPointData()->SetScalars(tag);
      }
    }

  const int numFloats = 7;
  const int numInts = 1;
  float block[numFloats]; // x,xvel,y,yvel,z,zvel,mass
  int iBlock[numInts]; // id
  int j = 0;
  double min[DIMENSION], max[DIMENSION];
  bool firstTime = true;

  for (int i = 0; i < DIMENSION; i++)
    {
    min[i] = 0;
    max[i] = -1;
    }

  // Loop to read all particle data
  for (vtkIdType i = this->PositionRange[0]; 
       i <= this->PositionRange[1]; 
       i += this->Stride)
    {
    j++;
    double progress = double(j) / double(this->NumberOfNodes);
    if (int(100 * progress) % 5 == 0)
      {
      this->UpdateProgress(progress);
      }

    // If stride > 1 we use seek to position to read the data record
    if (this->Stride > 1)
      {
      vtkIdType position = NUMBER_OF_DATA * i * BYTES_PER_DATA;
      this->FileStream->seekg(position, ios::beg);
      }

    // Read the floating point part of the data
    this->FileStream->read((char*) block, numFloats * sizeof(float));

    int returnValue = this->FileStream->gcount();
    if (returnValue != numFloats * (int)sizeof(float))
      {
      cerr << "vtkDLRReader Warning: read " << returnValue 
           << " floats" << endl;
      this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
      continue;
      }

    // Read the integer part of the data
    this->FileStream->read((char *)iBlock, numInts * sizeof(int));
    returnValue = this->FileStream->gcount();
    if (returnValue != numInts * (int)sizeof(int))
      {
      cerr << "vtkDLRReader Warning: read " << returnValue << " ints" << endl;
      this->SetErrorCode(vtkErrorCode::PrematureEndOfFileError);
      continue;
      }

    // These files are always little-endian
    vtkByteSwap::Swap4LERange(block, numFloats);

    // Negative value is an error so wraparound if it occurs
    if (block[X] < 0.0) 
      {
      block[X] = this->BoxSize + block[X];
      }
    if (block[Y] < 0.0) 
      {
      block[Y] = this->BoxSize + block[Y];
      }
    if (block[Z] < 0.0) 
      {
      block[Z] = this->BoxSize - block[Z];
      }

    // Insert the location into the point array
    vtkIdType vtkPointID = 
      points->InsertNextPoint(block[X], block[Y], block[Z]);
    if (this->MakeCells)
      {
      output->InsertNextCell(1, 1, &vtkPointID);
      }

    // Collect extents of positions
    if (firstTime == true)
      {
      min[0] = max[0] = block[X];
      min[1] = max[1] = block[Y];
      min[2] = max[2] = block[Z];
      firstTime = false;
      }
    else
      {
      if (min[0] > block[X]) 
        {
        min[0] = block[X];
        }
      if (max[0] < block[X]) 
        {
        max[0] = block[X];
        }
      if (min[1] > block[Y]) 
        {
        min[1] = block[Y];
        }
      if (max[1] < block[Y]) 
        {
        max[1] = block[Y];
        }
      if (min[2] > block[Z]) 
        {
        min[2] = block[Z];
        }
      if (max[2] < block[Z]) 
        {
        max[2] = block[Z];
        }
      }

    // Store velocity data if requested
    if (this->PointDataArraySelection->GetArraySetting(USE_VELOCITY))
      {
      velocity->SetComponent(vtkPointID, 0, block[X_VELOCITY]);
      velocity->SetComponent(vtkPointID, 1, block[Y_VELOCITY]);
      velocity->SetComponent(vtkPointID, 2, block[Z_VELOCITY]);
      }

    // Store mass data if requested
    if (this->PointDataArraySelection->GetArraySetting(USE_MASS))
      {
      mass->SetComponent(vtkPointID, 0, block[MASS]);
      }

    // Store tag data if requested
    if (this->PointDataArraySelection->GetArraySetting(USE_TAG))
      {
      tag->SetComponent(vtkPointID, 0, iBlock[0]);
      }
    } // end loop over PositionRange

  // Set the point extents on the output data
  GetOutput(0)->SetWholeExtent((int)floor(min[0]), (int)ceil(max[0]),
                               (int)floor(min[1]), (int)ceil(max[1]),
                               (int)floor(min[2]), (int)ceil(max[2]));
  GetOutput(0)->SetWholeBoundingBox(0.0, this->BoxSize,
                                    0.0, this->BoxSize,
                                    0.0, this->BoxSize);

  // Clean up internal storage
  velocity->Delete();
  mass->Delete();
  tag->Delete();
  points->Delete();
  output->Squeeze();
 
  // Close the file stream just read
  delete this->FileStream;
  this->FileStream = NULL;
#endif

  output->Squeeze();
}


//----------------------------------------------------------------------------
int vtkDLRReader::GetNumberOfVariables()
{
    return this->VariableSelection->GetNumberOfArrays();
}

//----------------------------------------------------------------------------
void vtkDLRReader::EnableAllVariables()
{
    this->VariableSelection->EnableAllArrays();
}

//----------------------------------------------------------------------------
void vtkDLRReader::DisableAllVariables()
{
    this->VariableSelection->DisableAllArrays();
}

//----------------------------------------------------------------------------
const char* vtkDLRReader::GetVariableName(int index)
{
    return this->VariableName[index].c_str();
}

//----------------------------------------------------------------------------
int vtkDLRReader::GetVariableStatus(const char* name)
{
    return this->VariableSelection->ArrayIsEnabled(name);
}

//----------------------------------------------------------------------------
void vtkDLRReader::SetVariableStatus(const char* name, int status)
{
    if( status )
        this->VariableSelection->EnableArray(name);
    else
        this->VariableSelection->DisableArray(name);
}

//----------------------------------------------------------------------------
vtkUnstructuredGrid* vtkDLRReader::GetOutput()
{
    return this->GetOutput(0);
}

//----------------------------------------------------------------------------
vtkUnstructuredGrid* vtkDLRReader::GetOutput(int idx)
{
    if(idx)
        return NULL;
    else
        return vtkUnstructuredGrid::SafeDownCast( this->GetOutputDataObject(idx) );
}
//----------------------------------------------------------------------------
bool vtkDLRReader::ReadPoints( int ncid, vtkUnstructuredGrid* grid )
{
    const char* varname[] = { "points_xc", "points_yc", "points_zc" };

    vtkFloatArray* coords  = vtkFloatArray::New();
    float* tmp = 0;
    
    try
    {
        int dimid;
        nccheck( nc_inq_dimid( ncid, "no_of_points", &dimid ) );

        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, dimid, &dimlen ) );

        int varid[3];
        for( unsigned int i=0; i<3; ++i )
        {
            int ndims, vdimids[NC_MAX_VAR_DIMS];
        
            nccheck( nc_inq_varid( ncid, varname[i], varid+i ) );
            nccheck( nc_inq_varndims( ncid, varid[i], &ndims ) );
            nccheck( nc_inq_vardimid( ncid, varid[i], vdimids ) );
        
            if( ndims > 1 || vdimids[0] != dimid )
                nccheck( NC_EDIMSIZE );
        }
    
        coords->SetNumberOfComponents( 3 );
        coords->SetNumberOfTuples( dimlen );
    
        tmp = new float[dimlen];
    
        for( unsigned int i=0; i<3; ++i )
        {
            nccheck( nc_get_var_float( ncid, varid[i], tmp ) );

            float* wptr = coords->WritePointer( 0, 3*dimlen ) + i;

            for( unsigned int i=0; i<dimlen; ++i, wptr += 3 )
                *wptr = tmp[i];      
        }
        
        delete[] tmp;
        
        vtkPoints* points = vtkPoints::New();
        points->SetData( coords );
        
        grid->SetPoints( points );
    }
    catch( nc_exception )
    {
        delete[] tmp;
        coords->Delete();
        coords = 0;
    }
    
    return coords;
}
//----------------------------------------------------------------------------
bool vtkDLRReader::ReadCells( int ncid, vtkUnstructuredGrid* grid )
{
    vtkCellArray*         cells     = vtkCellArray::New();
    vtkIdTypeArray*       locations = vtkIdTypeArray::New();
    vtkUnsignedCharArray* types     = vtkUnsignedCharArray::New();

    const struct {
        unsigned char type;
        unsigned int  size;
        const char*   varname;
    } ctypes[4] = {
        { VTK_TETRA,       4, "points_of_tetraeders" },
        { VTK_HEXAHEDRON,  8, "points_of_hexaeders" },
        { VTK_PYRAMID,     5, "points_of_pyramids" },
        { VTK_WEDGE,       6, "points_of_prisms" },
    };

    int* tmp;

    try
    {
        int varid[4];
        size_t size[4], numcells = 0, numindices = 0, maxindices = 0;

        for( int i=0; i<4; ++i )
        { 
            int result = nc_inq_varid( ncid, ctypes[i].varname, &varid[i] );
            
            if( result == NC_ENOTVAR )
            {
                varid[i] = -1;
                continue;
            }
            else
                nccheck( result );

            int ndims;
            nccheck( nc_inq_varndims( ncid, varid[i], &ndims ) );
        
            if( ndims!=2 )
                nccheck( NC_EVARSIZE );

            int dimids[2];
            nccheck( nc_inq_vardimid( ncid, varid[i], dimids ) );
            
            size_t dimlen[2];
            nccheck( nc_inq_dimlen( ncid, dimids[0], &dimlen[0] ) );
            nccheck( nc_inq_dimlen( ncid, dimids[1], &dimlen[1] ) );
        
            if( dimlen[1] != ctypes[i].size )
                nccheck( NC_EVARSIZE );
            
            size[i]     = dimlen[0];
            numcells   += dimlen[0];
            numindices += dimlen[0]*dimlen[1];

            if( dimlen[0]*dimlen[1] > maxindices )
                maxindices = dimlen[0]*dimlen[1];
        }
                
        int* tmp = new int[maxindices];
        
        vtkIdType*     iptr = cells->WritePointer( numcells, numcells + numindices );
        vtkIdType*     lptr = locations->WritePointer( 0, numcells );
        unsigned char* tptr = types->WritePointer( 0, numcells );
        
        for( unsigned int i=0, ci=0, ii=0; i<4; ++i )
        {
            if( varid[i] < 0 )
                continue;
                
            nccheck( nc_get_var_int( ncid, varid[i], tmp ) );

            const int* tmpptr = tmp;
            for( unsigned int j=0; j<size[i]; ++j, ++ci )
            {
                tptr[ci] = ctypes[i].type;
                lptr[ci] = ii;

                iptr[ii++] = ctypes[i].size;
                
                for( unsigned int k=0; k<ctypes[i].size; ++k )
                    iptr[ii++] = *(tmpptr++);
            }
        }
        
        delete[] tmp;
        
        grid->SetCells( types, locations, cells );
    }
    catch( nc_exception )
    {
        cells->Delete();
        types->Delete();
        locations->Delete();
        return false;
    }
    
    return true;
}

// -------------------------------------------------------------------------

bool vtkDLRReader::ReadVelocity( int ncid, vtkUnstructuredGrid* grid )
{
    const char* varname[] = { "x_velocity", "y_velocity", "z_velocity" };

    vtkFloatArray* vel  = vtkFloatArray::New();
    float* tmp = 0;
    
    try
    {
        int dimid;
        nccheck( nc_inq_dimid( ncid, "no_of_points", &dimid ) );

        size_t dimlen;
        nccheck( nc_inq_dimlen( ncid, dimid, &dimlen ) );

        int varid[3];
        for( unsigned int i=0; i<3; ++i )
        {
            int ndims, vdimids[NC_MAX_VAR_DIMS];
        
            nccheck( nc_inq_varid( ncid, varname[i], varid+i ) );
            nccheck( nc_inq_varndims( ncid, varid[i], &ndims ) );
            nccheck( nc_inq_vardimid( ncid, varid[i], vdimids ) );
        
            if( ndims > 1 || vdimids[0] != dimid )
                nccheck( NC_EDIMSIZE );
        }
    
        vel->SetNumberOfComponents( 3 );
        vel->SetNumberOfTuples( dimlen );
        vel->SetName( "velocity" );
    
        tmp = new float[dimlen];
    
        for( unsigned int i=0; i<3; ++i )
        {
            nccheck( nc_get_var_float( ncid, varid[i], tmp ) );

            float* wptr = vel->WritePointer( 0, 3*dimlen ) + i;

            for( unsigned int i=0; i<dimlen; ++i, wptr += 3 )
                *wptr = tmp[i];      
        }
        
        delete[] tmp;
        
        grid->GetPointData()->SetVectors( vel );
    }
    catch( nc_exception )
    {
        delete[] tmp;
        vel->Delete();
        vel = 0;
    }
    
    return vel;
}
