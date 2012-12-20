/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#include "itkMINCImageIO.h"
#include <stdio.h>
#include "vnl/vnl_vector.h"

namespace itk
{


bool MINCImageIO::CanReadFile(const char *file)
{
  if ( *file == 0 )
    {
    itkDebugMacro(<< "No filename specified.");
    return false;
    }

  std::string::size_type mncPos = filename.rfind(".mnc");
  if ( (mncPos != std::string::npos)
	&& (mncPos == filename.length() - 4) )
    {
    return true;
    }

  mncPos = filename.rfind(".MNC");
  if ( (mncPos != std::string::npos)
	&& (mncPos == filename.length() - 4) )
    {
    return true;
    }

  mncPos = filename.rfind(".mnc2");
  if ( (mncPos != std::string::npos)
	&& (mncPos == filename.length() - 5) )
    {
    return true;
    }

  mncPos = filename.rfind(".MNC2");
  if ( (mncPos != std::string::npos)
	&& (mncPos == filename.length() - 5) )
    {
    return true;
    }

//   mihandle_t volume;
//
//   if ( miopen_volume(file, MI2_OPEN_READ, &volume) < 0 )
//     {
//     itkDebugMacro(<< " Can not open File:" << file << "\n");
//     return false;
//     }
//   if ( miclose_volume(volume) < 0 )
//     {
//     itkDebugMacro(<< " Can not close File:" << file << "\n");
//     return false;
//     }

  return true;
}

void MINCImageIO::Read(void *buffer)
{
  mihandle_t volume;

  // call to minc2.0 function to open the file
  if ( miopen_volume(m_FileName.c_str(), MI2_OPEN_READ, &volume) < 0 )
    {
    // Error opening the volume
    itkDebugMacro("Could not open file \"" << m_FileName.c_str() << "\".");
    return;
    }
  mitype_t volume_data_type;
  if ( miget_data_type(volume, &volume_data_type) < 0 )
    {
    itkDebugMacro(" Can not get volume data type!!\n");
    }
  int slice_scaling_flag;
  // find out whether the data has slice scaling
  if ( miget_slice_scaling_flag(volume, &slice_scaling_flag) < 0 )
    {
    itkDebugMacro(" Can not get slice scaling flag!!\n");
    }
  double valid_max, valid_min;
  // find the data valid range
  if ( miget_volume_valid_range(volume, &valid_max, &valid_min) < 0 )
    {
    itkDebugMacro(" Can not get valid max or min!!\n");
    }

  miclass_t volume_data_class;
  if ( miget_data_class(volume, &volume_data_class) < 0 )
    {
    itkDebugMacro(" Can not get volume data class!!\n");
    }

  misize_t start[MINC_MAXDIM + 1];
  misize_t count[MINC_MAXDIM + 1];

  // figure out how many dimensions out of the total NDims
  // are used by this class
  size_t usefulDimensions = 0;
  size_t i;
  for ( i = 0; i < MINC_MAXUSE; i++ )
    {
    if ( this->m_DimensionIndices[i] != -1 )
      {
      usefulDimensions++;
      }
    }

  // fill out the array of dimension handles,"regularly sampled"
  // the dimensions will be retrieved in file order
  midimhandle_t *hdims = new midimhandle_t[usefulDimensions];
  if ( miget_volume_dimensions(volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, this->m_NDims,
                               hdims) < 0 )
    {
    itkDebugMacro(" Can not get dimension handles!!\n");
    return;
    }

  midimhandle_t *apparent_order = new  midimhandle_t[usefulDimensions];
  // order of dim_indices x,y,z,t,vector-dimension
  // apparent order vector-dimension,t,z,y,x
  size_t j = 0;
  for ( i = 0; i < 5; i++ )
    {
    if ( this->m_DimensionIndices[i] != -1 )
      {
      apparent_order[j] = hdims[this->m_DimensionIndices[i]];
      j++;
      }
    }

  //check to see if app order same as file order
  for ( i = 0; i < usefulDimensions; i++ )
    {
    if ( this->m_DimensionIndices[i] != static_cast< int >( i ) )
      {
      // set apparent order of dimensions so data can be accesed in that order
      if ( miset_apparent_dimension_order(volume, usefulDimensions, apparent_order) < 0 )
        {
        itkDebugMacro(" Can not get apparent dimension order!!\n");
        }
      break;
      }
    }
  // clean dynamic space allocated for dimension handles
  delete[] hdims;
  delete[] apparent_order;
  //set the unused dimension to start 0 and offset 1 if ANY
  i = 0;
  for ( i = 0; i < ( this->m_NDims - usefulDimensions ); i++ )
    {
    start[i] = 0;
    count[i] = 1;
    }

  start[i] = 0;
  count[i] = this->GetDimensions(2);
  i++;
  start[i] = 0;
  count[i] = this->GetDimensions(1);
  i++;
  start[i] = 0;
  count[i] = this->GetDimensions(0);
  i++;

  //now take care of vector or time dimension
  // z,y,x, t, vector_dimension
  if ( usefulDimensions == 5 )
    {
    start[i] = 0;
    misize_t icount = count[i];
    if ( miget_dimension_size(apparent_order[3], &icount) < 0 )
      {
      itkDebugMacro(" Can not get dimension size \n");
      }
    count[i] = icount;
    i++;
    start[i] = 0;
    icount = count[i];
    if ( miget_dimension_size(apparent_order[4], &icount) < 0 )
      {
      itkDebugMacro(" Can not get dimension size \n");
      }
    count[i] = icount;
    }
  if ( usefulDimensions == 4 )
    {
    start[i] = 0;
    misize_t icount = count[i];
    if ( miget_dimension_size(apparent_order[3], &icount) < 0 )
      {
      itkDebugMacro(" Can not get dimension size \n");
      }
    count[i] = icount;
    }

  // set the data class for the file
  switch ( volume_data_class )
    {
    case MI_CLASS_REAL:
    case MI_CLASS_COMPLEX:
    case MI_CLASS_INT:
      if ( slice_scaling_flag )
        {
        if ( miget_real_value_hyperslab(volume, volume_data_type, start, count, buffer) < 0 )
          {
          itkDebugMacro(" Can not get real value hyperslab!!\n");
          }
        }
      else
        {
        if ( miget_voxel_value_hyperslab(volume, volume_data_type, start, count, buffer) < 0 )
          {
          itkDebugMacro(" Can not get voxel value hyperslabs!!\n");
          }
        }
      break;
    case MI_CLASS_LABEL:
      if ( miget_voxel_value_hyperslab(volume, volume_data_type, start, count, buffer) < 0 )
        {
        itkDebugMacro(" Can not get voxel value hyperslabs!!\n");
        }
      break;
    case MI_CLASS_UNIFORM_RECORD:
      {
      itkDebugMacro(" Leave this until minc2.0 support it complete!!\n");
      }
    default:
      return;
    }

  if ( miclose_volume(volume) < 0 )
    {
    itkDebugMacro(" Can not close volume!\n");
    }
}


void MINCImageIO::CleanupDimensions(void)
{
  for ( int i = 0; i < this->m_NDims; i++ )
    {
     if(this->m_DimensionName[i])
       free(this->m_DimensionName[i]);
     this->m_DimensionName[i]=NULL;
    }

  if(this->m_DimensionName) delete[] this->m_DimensionName;
  if(this->m_DimensionSize) delete[] this->m_DimensionSize;
  if(this->m_DimensionStart) delete[] this->m_DimensionStart;
  if(this->m_DimensionStep) delete[] this->m_DimensionStep;
  if(this->m_DimensionIndices) delete[] this->m_DimensionIndices;

  this->m_DimensionName  = NULL;
  this->m_DimensionSize  = NULL;
  this->m_DimensionStart = NULL;
  this->m_DimensionStep  = NULL;
  this->m_DimensionIndices = NULL;

}


void MINCImageIO::AllocateDimensions(int nDims)
{
  this->CleanupDimensions();

  m_NDims=nDims;

  this->m_DimensionName  = new const char*[m_NDims];
  this->m_DimensionSize  = new misize_t[m_NDims];
  this->m_DimensionStart = new double[m_NDims];
  this->m_DimensionStep  = new double[m_NDims];

  for ( int i = 0; i < this->m_NDims; i++ )
    {
    this->m_DimensionName[i]  = NULL;
    this->m_DimensionSize[i]  = 0;
    this->m_DimensionStart[i] = 0.0;
    this->m_DimensionStep[i]  = 0.0;
    }

  for ( int i = 0; i < sizeof(m_DimensionIndices) ; i++ )
  {
    this->m_DimensionIndices[i] = -1;
  }

}

// close existing volume, cleanup internal structures
void MINCImageIO::CloseVolume(void)
{
  if( m_volume )
    miclose_volume( m_volume );
  m_volume = NULL;

  this->CleanupDimensions();
}

MINCImageIO::MINCImageIO()
{
  this->m_NDims = 0;
  this->m_DimensionName  = NULL;
  this->m_DimensionSize  = NULL;
  this->m_DimensionStart = NULL;
  this->m_DimensionStep  = NULL;

  m_Complex = 0;
  m_volume = NULL;

  for ( int i = 0; i < sizeof(m_DimensionIndices) ; i++ )
  {
    this->m_DimensionIndices[i] = -1;
  }
}

MINCImageIO::~MINCImageIO()
{
  CloseVolume();
}

void MINCImageIO::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "NDims: " << m_NDims << std::endl;
}

void MINCImageIO::ReadImageInformation()
{
  this->CloseVolume();
  // call to minc2.0 function to open the file
  if ( miopen_volume(m_FileName.c_str(), MI2_OPEN_READ, &m_volume) < 0 )
    {
    // Error opening the volume
    itkDebugMacro("Could not open file \"" << m_FileName.c_str() << "\".");
    return;
    }

  // find out how many dimensions are there regularly sampled
  // dimensions only
  int ndims;
  if ( miget_volume_dimension_count(m_volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, &ndims) < 0 )
    {
    itkDebugMacro("Could not get the number of dimensions in the volume!");
    return;
    }
  this->AllocateDimensions(ndims);


  // get dimension handles in FILE ORDER (i.e, the order as they are
  // submitted to file)
  midimhandle_t *hdims = new midimhandle_t[m_NDims];
  if ( miget_volume_dimensions(volume, MI_DIMCLASS_ANY, MI_DIMATTR_ALL, MI_DIMORDER_FILE, m_NDims,
                               hdims) < 0 )
    {
    itkDebugMacro("Could not get dimension handles!");
    return;
    }

  for (int i = 0; i < m_NDims; i++ )
    {
    char *name;
    if ( miget_dimension_name(hdims[i], &name) < 0 )
      {
      // Error getting dimension name
      itkDebugMacro("Could not get dimension name!");
      return;
      }

    this->m_DimensionName[i] = name;
    if(!strcmp(name,MIxspace) || !strcmp(name,MIxfrequency)) //this is X space
      m_DimensionIndices[1]=i;
    else if(!strcmp(name,MIyspace) || !strcmp(name,MIyfrequency)) //this is Y space
      m_DimensionIndices[2]=i;
    else if(!strcmp(name,MIzspace) || !strcmp(name,MIzfrequency)) //this is Z space
      m_DimensionIndices[3]=i;
    else if(!strcmp(name,MIvector_dimension) ) //this is vector space
      m_DimensionIndices[0]=i;
    else if(!strcmp(name,MItime) || !strcmp(name,MItfrequency)) //this is time space
      m_DimensionIndices[4]=i;
    else
      {
      itkDebugMacro(<<"Unsupported MINC dimension:"<<name);
      return;
      }
    }

  // fill the DimensionSize by calling the following MINC2.0 function
  if ( miget_dimension_sizes(hdims, m_NDims, m_DimensionSize) < 0 )
    {
    // Error getting dimension sizes
    itkDebugMacro("Could not get dimension sizes!");
    return;
    }

  if ( miget_dimension_separations(hdims, MI_ORDER_FILE, m_NDims, m_DimensionStep) < 0 )
    {
    itkDebugMacro(<<" Could not dimension sizes");
    return;
    }

  if ( miget_dimension_starts(hdims, MI_ORDER_FILE, this->m_NDims, m_DimensionStart) < 0 )
    {
    itkDebugMacro(<<" Could not dimension sizes");
    return;
    }

  // analyze dimension order and decide how to restructure the data

  mitype_t volume_data_type;
  if ( miget_data_type(volume, &volume_data_type) < 0 )
    {
    itkDebugMacro(" Can not get volume data type!!\n");
    }

  // find out whether the data has slice scaling
  miboolean_t slice_scaling_flag=0;
  miboolean_t global_scaling_flag=0;

  if ( miget_slice_scaling_flag(volume, &slice_scaling_flag) < 0 )
    {
    itkDebugMacro(" Can not get slice scaling flag!!\n");
    }

  //voxel valid range
  double valid_min,valid_max;
  //get the voxel valid range
  if(miget_volume_valid_range(m_volume,&valid_max,&valid_min) < 0 )
  {
    itkDebugMacro(" Can not get volume valid range!!\n");
  }

  //real volume range, only awailable when slice scaling is off
  double volume_min=0.0,volume_max=1.0;
  if( !slice_scaling_flag )
  {
    if( miget_volume_range(m_volume,&valid_max,&valid_min) < 0 )
    {
      itkDebugMacro(" Can not get volume range!!\n");
    }
    global_scaling_flag=(volume_min==valid_min && volume_max==valid_max);
  }

  miclass_t volume_data_class;

  if ( miget_data_class(volume, &volume_data_class) < 0 )
    {
    itkDebugMacro(" Could not get data class");
    return;
    }

  // set the file data type
  if(slice_scaling_flag || global_scaling_flag)
  {
    switch ( volume_data_type )
      {
      case MI_TYPE_FLOAT:
        this->SetComponentType(FLOAT);
        break;
      case MI_TYPE_DOUBLE:
        this->SetComponentType(DOUBLE);
        break;
      case MI_TYPE_FCOMPLEX:
        this->SetComponentType(FLOAT);
        break;
      case MI_TYPE_DCOMPLEX:
        this->SetComponentType(DOUBLE);
        break;
      default:
        this->SetComponentType(FLOAT);
        break;
      } //end of switch
    //file will have do
  } else {
    switch ( volume_data_type )
      {
      case MI_TYPE_BYTE:
        this->SetComponentType(CHAR);
        break;
      case MI_TYPE_UBYTE:
        this->SetComponentType(UCHAR);
        break;
      case MI_TYPE_SHORT:
        this->SetComponentType(SHORT);
        break;
      case MI_TYPE_USHORT:
        this->SetComponentType(USHORT);
        break;
      case MI_TYPE_INT:
        this->SetComponentType(INT);
        break;
      case MI_TYPE_UINT:
        this->SetComponentType(UINT);
        break;
      case MI_TYPE_FLOAT:
        this->SetComponentType(FLOAT);
        break;
      case MI_TYPE_DOUBLE:
        this->SetComponentType(DOUBLE);
        break;
      case MI_TYPE_SCOMPLEX:
        this->SetComponentType(SHORT);
        break;
      case MI_TYPE_ICOMPLEX:
        this->SetComponentType(INT);
        break;
      case MI_TYPE_FCOMPLEX:
        this->SetComponentType(FLOAT);
        break;
      case MI_TYPE_DCOMPLEX:
        this->SetComponentType(DOUBLE);
        break;
      default:
        itkDebugMacro("Bad data type ");
        return;
      } //end of switch
  }

  this->ComputeStrides();

  switch ( volume_data_class )
    {
    case MI_CLASS_REAL:
      this->SetPixelType(SCALAR);
      break;
    case MI_CLASS_INT:
      this->SetPixelType(SCALAR);
      break;
    case MI_CLASS_LABEL:
      this->SetPixelType(SCALAR);
      // create an array of label names and values
      // not sure how to pass this to itk yet!
      break;
    case MI_CLASS_COMPLEX:
      m_Complex = 1;
      this->SetPixelType(COMPLEX);
      numberOfComponents *= 2;
      break;
    default:
      itkDebugMacro("Bad data class ");
      return;
    } //end of switch


  this->SetNumberOfDimensions(m_NDims);
  this->SetNumberOfComponents(numberOfComponents);

  if ( miclose_volume(volume) < 0 )
    {
    // Error closing the volume
    itkDebugMacro("Could not close file \"" << m_FileName.c_str() << "\".");
    return;
    }

  //clean dynamic allocation
  delete[] sizes;
  delete[] hdims;
}

void MINCImageIO::SetDimensionName(size_t i, char *name)
{
  if ( name )
    {
    this->Modified();
    this->m_DimensionName[i] = name;
    }
  else
    {
    return;
    }
}

bool MINCImageIO::CanWriteFile(const char *name)
{
  std::string filename = name;

  // transform filename to lower case to make checks case-insensitive
  std::transform(filename.begin(), filename.end(), filename.begin(), ( int ( * )(int) ) std::tolower);

  if (  filename == "" )
    {
    itkDebugMacro(<< "No filename specified.");
    return false;
    }

  std::string::size_type mncPos = filename.rfind(".mnc");
  if ( ( mncPos != std::string::npos )
       && ( mncPos > 0 )
       && ( mncPos == filename.length() - 4 ) )
    {
    return true;
    }

  mncPos = filename.rfind(".mnc2");
  if ( ( mncPos != std::string::npos )
       && ( mncPos > 0 )
       && ( mncPos == filename.length() - 5 ) )
    {
    return true;
    }

  return false;
}

/*
 * fill out the appropriate header information
*/
void MINCImageIO::WriteImageInformation(void)
{
  std::cout << "WriteImageInformation" << std::endl;
  // FIXME: implement this!
}

template< class TBuffer >
void MINCComputeScalarRange(int itkNotUsed(
                              Strides)[3], int Sizes[3], int nComponents, double & maxval, double & minval,
                            TBuffer *buffer)
{
// This differs with ITK Journal version in that
// no longer skipping ahead by strides, but just plodding through
// buffer one point at a time.

// FIXME: this difference is probably something that should be thought
// through more clearly.
// Rupert Brooks, August 23, 2007
  double tmpminval = 1000.0;
  double tmpmaxval = 0.0;
  size_t idX, idY, idZ;

  for ( idZ = 0; idZ < Sizes[2]; idZ++ )
    {
    for ( idY = 0; idY < Sizes[1]; idY++ )
      {
      for ( idX = 0; idX < Sizes[0] * nComponents; idX++ )
        {
        TBuffer val = *buffer++; //Strides[0];
        if ( val > tmpmaxval )
          {
          tmpmaxval = val;
          }
        if ( val < tmpminval )
          {
          tmpminval = val;
          }
        }
      //buffer += Strides[1];
      }
    //buffer += Strides[2];
    }
  maxval = tmpmaxval;
  minval = tmpminval;
}

//void MINCImageIO::
template< class TBuffer >
void MINCWriteHyperSlab(mihandle_t volume,
                        size_t ndims,
                        mitype_t minctype,
                        const misize_t start[],
                        const misize_t count[],
                        TBuffer *buffer)
{
  size_t i, j;

  misize_t tmpstart[MINC_MAXDIM + 1];
  misize_t tmpcount[MINC_MAXDIM + 1];
  // calculate the number of voxels per slice
  size_t size = 1;

  for ( i = 1; i < ndims; i++ )
    {
    tmpstart[i] = start[i];

    tmpcount[i] = count[i];

    size = size * count[i];
    }

  // allocate memory for a slice
  TBuffer *tmpbuffer = new TBuffer[size];
  for ( i = 0; i < count[0]; i++ )
    {
    // set start to the current slice
    tmpstart[0] = i;
    tmpcount[0] = 1;

    // copy the slice
    for ( j = 0; j < size; j++ )
      {
      tmpbuffer[j] = buffer[j];
      }

    // write the slice
    miset_voxel_value_hyperslab(volume, minctype,
                                tmpstart, tmpcount, tmpbuffer);

    // move on to next slice
    buffer += size;
    }

  delete[] tmpbuffer;
}

void MINCImageIO::Write(const void *buffer)
{
  size_t i, j;
  size_t ncomp;

// in general we cannot assume that m_original start or m_DirectionCosines exist
// have to recompute them

  // FIXME: i use a vnl_matrix here, because i know how to invert it, and the
  // itk::Matrix
  // used for the m_DirectionCosines I dont know how.  However, maybe this
  // matrix
  // should be a different type anyway, so its variable size.
  // This inelegant solution is Rupert being too lazy to RTFM, not for any good
  // reason
  vnl_matrix< double > dircosmatrix(3, 3);
  for ( i = 0; i < 3; i++ )
    {
    for ( j = 0; j < 3; j++ )
      {
      m_DirectionCosines[i][j] = this->m_Direction[j][i];
      dircosmatrix[i][j] = m_DirectionCosines[i][j];
      }
    }
  vnl_matrix< double > inverseDirectionCosines = vnl_matrix_inverse< double >(dircosmatrix);
  /*
  for( i = 0; i < 3; i++ )
    {
    m_OriginalStart[i]=0;
    for( j = 0; j < 3; j++ )
      {
      m_OriginalStart[i] += inverseDirectionCosines[i][j] * this->GetOrigin(j);
      }
    }
  */
  for ( i = 0; i < 3; i++ )
    {
    m_OriginalStart[i] = this->GetOrigin(i);
    }
  ncomp = this->GetNumberOfComponents();

  // ensure that the type is valid
  if ( m_Complex )
    {
    if ( this->GetComponentType() == CHAR || this->GetComponentType() == UCHAR )
      {
      itkDebugMacro("MINC does not support 8-bit complex types");
      return;
      }
    else if ( this->GetComponentType() == USHORT || this->GetComponentType() == UINT )
      {
      itkDebugMacro("MINC does not support unsigned complex types");
      return;
      }
    }

  // get the dimension order set by the user (copy it because
  // CheckDimensionOrder will modify it)
  char userdimorder[MINC_MAXDIM + 1];
  if ( this->GetDimensionOrder() != 0 )
    {
    strncpy(userdimorder, this->GetDimensionOrder(), MINC_MAXDIM);
    userdimorder[MINC_MAXDIM] = '\0';
    }
  else
    {
    strncpy(userdimorder, "zyx", MINC_MAXDIM);
    }

  // check the dimension order, add any dimensions that
  //  the users have set sizes for but which aren't in DimensionOrder
  if ( this->CheckDimensionOrder(userdimorder) == 0 )
    {
    return;
    }

  // three dimensions (x,y,z) are always present
  size_t ndims = 3;
  size_t vsize = 0; // vector_dimension size
  size_t usize = 0; // produce of all user defined dimensions
  size_t tsize = 0; // time dimension size
  size_t xsize = this->GetDimensions(0);
  size_t ysize = this->GetDimensions(1);
  size_t zsize = this->GetDimensions(2);

  const char *vname = "vector_dimension";
  const char *tname = "tspace";
  const char *xname = "xspace";
  const char *yname = "yspace";
  const char *zname = "zspace";

  // update the names of dimensions if the user has specified
  // e.g. xfrequency instead of xspace, and also get the sizes
  // for all non-spatial dimensions
  for ( i = 0; i <= MINC_MAXDIM; i++ )
    {
    if ( this->m_DimensionName[i] )
      {
      // the first char in the name
      int dimchar = this->m_DimensionName[i][0];

      if ( dimchar == 'v' )
        { // add v dimension
        vsize = this->m_DimensionSize[i];
        vname = this->m_DimensionName[i];
        if ( vsize > 0 )
          {
          ndims++;
          }
        }
      else if ( dimchar == 't' )
        { // add t dimension
        tsize = this->m_DimensionSize[i];
        tname = this->m_DimensionName[i];
        if ( tsize > 0 )
          {
          ndims++;
          }
        }
      else if ( dimchar == 'x' )
        {
        // xsize is calculated from extent
        xname = this->m_DimensionName[i];
        }
      else if ( dimchar == 'y' )
        {
        // ysize is calculated from extent
        yname = this->m_DimensionName[i];
        }
      else if ( dimchar == 'z' )
        {
        // zsize is calculated from extent
        zname = this->m_DimensionName[i];
        }
      else
        { // add other dimensions as vector dimensions
        if ( this->m_DimensionSize[i] > 0 )
          {
          if ( usize == 0 )
            {
            usize = this->m_DimensionSize[i];
            }
          else
            {
            usize *= this->m_DimensionSize[i];
            }
          ndims++;
          }
        }
      }
    }

  // csize depends on whether the data is complex
  size_t csize = ( m_Complex ? 2 : 1 );
  // check the number of components: the product of the sizes of
  // the complex, vector, and time dimensions should equal the
  // number of scalar components in the vtkImageData
  if ( ( tsize == 0 ? 1 : tsize ) * ( vsize == 0 ? 1 : vsize )
       * ( usize == 0 ? 1 : usize ) * csize != ncomp )
    {
    // the number of components does not match extra dimension sizes,
    // so try to make a vector dimension to account for this.

    // calculate the product of the dimension sizes except for vsize
    size_t dimprod = ( usize == 0 ? 1 : usize ) * ( tsize == 0 ? 1 : tsize ) * csize;

    if ( ncomp % dimprod != 0 )
      {
      itkDebugMacro(
        "Write: Product of non-spatial dimension sizes does not match number of scalar components: " << tsize << "*"
                                                                                                     << usize
        * csize << "!=" << ncomp);
      return;
      }

    // add the vector dimension if it was missing
    if ( vsize == 0 )
      {
      ndims++;
      // add to userdimorder unless it was already there
      char *cp;
      for ( cp = userdimorder; *cp != '\0'; cp++ )
        {
        if ( *cp == 'v' )
          {
          break;
          }
        }
      if ( *cp == '\0' ) // if didn't find 'v'
        {
        *cp++ = 'v';
        *cp++ = '\0';
        }
      }

    // calcuate the proper size of the vector dimension
    vsize = ncomp / dimprod;
    }

  // remove any unused dimensions from userdimorder
  i = 0;
  char *cp;
  for ( cp = userdimorder; *cp != '\0'; cp++ )
    {
    int dimchar = *cp;
    if ( dimchar == xname[0] || dimchar == yname[0] || dimchar == zname[0] )
      {
      userdimorder[i++] = dimchar;
      }
    else if ( ( dimchar == tname[0] && tsize != 0 )
              || ( dimchar == vname[0] && vsize != 0 ) )
      {
      userdimorder[i++] = dimchar;
      }
    else
      {
      for ( j = 0; j <= MINC_MAXDIM; j++ )
        {
        if ( this->m_DimensionName[j] && this->m_DimensionName[j][0] == dimchar )
          {
          userdimorder[i++] = dimchar;
          break;
          }
        }
      }
    }

  if ( i != ndims )
    {
    itkDebugMacro("Failed sanity check: i != ndims ("
                  << i << "!=" << ndims << ")");
    return;
    }

  userdimorder[i++] = '\0'; // terminate the string

  int           result = 0;
  mihandle_t    volume;
  midimhandle_t dim[MINC_MAXDIM + 1];
  misize_t offsets[MINC_MAXDIM + 1];
  misize_t counts[MINC_MAXDIM + 1];
  // create the MINC dimensions in the file order given by userdimorder,
  // remove any characters from userdimorder that don't match a used
  // dimension.
  for ( i = 0; i < ndims; i++ )
    {
    unsigned long dimsize = 0;
    midimclass_t  dimclass = MI_DIMCLASS_ANY;
    const char *  dimname = 0;
    int           dimchar = userdimorder[i];
    double        dimseparation = 1.0;
    double        dimstart = 0.0;

    if ( dimchar == xname[0] || dimchar == yname[0] || dimchar == zname[0] )
      { // spatial or spatial frequency
      if ( dimchar == xname[0] )
        {
        dimname = xname;
        dimsize = xsize;
        dimseparation = this->GetSpacing(0);
        dimstart = m_OriginalStart[0];
        }
      else if ( dimchar == yname[0] )
        {
        dimname = yname;
        dimsize = ysize;
        dimseparation = this->GetSpacing(1);
        dimstart = m_OriginalStart[1];
        }
      else /* if (dimchar == zname[0]) */
        {
        dimname = zname;
        dimsize = zsize;
        dimseparation = this->GetSpacing(2);
        dimstart = m_OriginalStart[2];
        }
      dimclass = MI_DIMCLASS_SPATIAL;
      if ( strcmp(&dimname[1], "frequency") == 0 )
        {
        dimclass = MI_DIMCLASS_SFREQUENCY;
        }
      }
    else if ( dimchar == tname[0] && tsize != 0 )
      { // time or tfrequency
      dimname = tname;
      dimsize = tsize;
      dimclass = MI_DIMCLASS_TIME;
      if ( strcmp(&dimname[1], "frequency") == 0 )
        {
        dimclass = MI_DIMCLASS_TFREQUENCY;
        }
      }
    else if ( dimchar == vname[0] && vsize != 0 )
      {                              // vector dimension
      dimname = vname;               // default is "vector_dimension"
      dimsize = vsize;               // default is all leftovers
      dimclass = MI_DIMCLASS_RECORD; // unknown
      }
    else
      { // other dimensions
        // search through user-defined dimensions
      for ( j = 0; j <= MINC_MAXDIM; j++ )
        {
        if ( this->m_DimensionName[j] && this->m_DimensionName[j][0] == dimchar )
          {
          dimname = this->m_DimensionName[j];
          dimsize = this->m_DimensionSize[j];
          dimclass = MI_DIMCLASS_USER; // unknown
          }
        }
      }

    //create MINC.0 file
    micreate_dimension(dimname, dimclass, MI_DIMATTR_REGULARLY_SAMPLED, dimsize, &dim[i]);

    // modify some parameters
    if ( dimclass == MI_DIMCLASS_SPATIAL || dimclass == MI_DIMCLASS_SFREQUENCY )
      {
      miset_dimension_units(dim[i], "mm");
      miset_dimension_start(dim[i], /* MI_ORDER_APPARENT, */ dimstart);
      miset_dimension_separation(dim[i], /* MI_ORDER_APPARENT,*/ dimseparation);

      double dircos[3];
      if ( dimname[0] == 'x' )
        {
        dircos[0] = m_DirectionCosines[0][0];
        dircos[1] = m_DirectionCosines[1][0];
        dircos[2] = m_DirectionCosines[2][0];
        }
      else if ( dimname[0] == 'y' )
        {
        dircos[0] = m_DirectionCosines[0][1];
        dircos[1] = m_DirectionCosines[1][1];
        dircos[2] = m_DirectionCosines[2][1];
        }
      else if ( dimname[0] == 'z' )
        {
        dircos[0] = m_DirectionCosines[0][2];
        dircos[1] = m_DirectionCosines[1][2];
        dircos[2] = m_DirectionCosines[2][2];
        }

      miset_dimension_cosines(dim[i], dircos);
      }
    else if ( dimclass == MI_DIMCLASS_TIME || dimclass == MI_DIMCLASS_TFREQUENCY )
      {
      miset_dimension_units(dim[i], "s");
      }
    }

  // set the file data type
  mitype_t minctype;
  switch ( this->GetComponentType() )
    {
    case CHAR:
      minctype = MI_TYPE_BYTE;
      break;
    case UCHAR:
      minctype = MI_TYPE_UBYTE;
      break;
    case SHORT:
      minctype = MI_TYPE_SHORT;
      break;
    case USHORT:
      minctype = MI_TYPE_USHORT;
      break;
    case INT:
      minctype = MI_TYPE_INT;
      break;
    case UINT:
      minctype = MI_TYPE_UINT;
      break;
    case FLOAT:
      minctype = MI_TYPE_FLOAT;
      break;
    case DOUBLE:
      minctype = MI_TYPE_DOUBLE;
      break;
    default:
      itkDebugMacro("Bad data type ");
      return;
    } //end of switch

  // find the class
  miclass_t mincclass = MI_CLASS_REAL;
  if ( m_Complex )
    {
    mincclass = MI_CLASS_COMPLEX;
    }

  result = micreate_volume(m_FileName.c_str(), ndims, dim, minctype, mincclass, NULL, &volume);
  if ( result >= 0 )
    {
    result = micreate_volume_image(volume);
    }
  // check for failed file open
  if ( result < 0 )
    {
    for ( i = 0; i < ndims; i++ )
      {
      mifree_dimension_handle(dim[i]);
      }
    itkDebugMacro( "Couldn't open minc file " << m_FileName.c_str() );
    return;
    }
  // set the apparent order before writing hyperslabs
  const char *dimnames[MINC_MAXDIM + 1];
  i = 0;
  dimnames[i] = zname; counts[i] = zsize; offsets[i++] = 0;
  dimnames[i] = yname; counts[i] = ysize; offsets[i++] = 0;
  dimnames[i] = xname; counts[i] = xsize; offsets[i++] = 0;
  if ( tsize != 0 )
    { // add time if it exists
    dimnames[i] = tname; counts[i] = tsize; offsets[i++] = 0;
    }
  // add other dimensions
  for ( j = 0; j <= MINC_MAXDIM; j++ )
    {
    if ( this->m_DimensionName[j]
         && strcmp(this->m_DimensionName[j], zname) != 0
         && strcmp(this->m_DimensionName[j], yname) != 0
         && strcmp(this->m_DimensionName[j], xname) != 0
         && strcmp(this->m_DimensionName[j], tname) != 0
         && strcmp(this->m_DimensionName[j], vname) != 0 )
      {
      dimnames[i] = this->m_DimensionName[j];
      counts[i] = this->m_DimensionSize[j];
      offsets[i++] = 0;
      }
    }
  // finally add a vector_dimension if we had to create one
  if ( vsize != 0 )
    {
    dimnames[i] = vname;
    counts[i] = vsize;
    offsets[i++] = 0;
    }
  // convert to string to compare to userdimorder
  char dimorder[MINC_MAXDIM + 1];
  for ( i = 0; i < ndims; i++ )
    {
    dimorder[i] = dimnames[i][0];
    }
  dimorder[ndims] = '\0';

  if ( strcmp(dimorder, userdimorder) != 0 )
    {
    miset_apparent_dimension_order_by_name(volume, ndims, (char **)dimnames);
    }

  // writing data in slice by slice
  switch ( this->GetComponentType() )
    {
    case CHAR:
      minctype = MI_TYPE_BYTE;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (char *)buffer);
      break;
    case UCHAR:
      minctype = MI_TYPE_UBYTE;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (unsigned char *)buffer);
      break;
    case SHORT:
      minctype = MI_TYPE_SHORT;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (short *)buffer);
      break;
    case USHORT:
      minctype = MI_TYPE_USHORT;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (unsigned short *)buffer);
      break;
    case INT:
      minctype = MI_TYPE_INT;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (int *)buffer);
      break;
    case UINT:
      minctype = MI_TYPE_UINT;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (unsigned int *)buffer);
      break;
    case FLOAT:
      minctype = MI_TYPE_FLOAT;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (float *)buffer);
      break;
    case DOUBLE:
      minctype = MI_TYPE_DOUBLE;
      MINCWriteHyperSlab(volume, ndims, minctype, offsets, counts, (double *)buffer);
      break;
    default:
      itkDebugMacro("Bad data type ");
      return;
    } //end of switch

  // set the min/max
  if ( mincclass == MI_CLASS_REAL
       && minctype != MI_TYPE_FLOAT && minctype != MI_TYPE_DOUBLE )
    {
    // need to calculate the min/max for the specified extent
    double minval, maxval;
    int    Strides[3];
    Strides[0] = m_Strides[0];
    Strides[1] = m_Strides[1];
    Strides[2] = m_Strides[2];
    int Sizes[3];
    /* returns empty !??
    Sizes[0] = this->m_DimensionSize[0];
    Sizes[1] = this->m_DimensionSize[1];
    Sizes[2] = this->m_DimensionSize[2];*/
    Sizes[0] = this->GetDimensions(0);
    Sizes[1] = this->GetDimensions(1);
    Sizes[2] = this->GetDimensions(2);
    switch ( minctype )
      {
      case MI_TYPE_BYTE:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (char *)buffer);
        break;
      case MI_TYPE_UBYTE:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (unsigned char *)buffer);
        break;
      case MI_TYPE_SHORT:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (short *)buffer);
        break;
      case MI_TYPE_USHORT:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (unsigned short *)buffer);
        break;
      case MI_TYPE_INT:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (int *)buffer);
        break;
      case MI_TYPE_UINT:
        MINCComputeScalarRange(Strides, Sizes, this->GetNumberOfComponents(), maxval, minval, (unsigned int *)buffer);
        break;
      default:
        itkDebugMacro("Bad data type ");
        return;
      } //end of switch
    miset_volume_valid_range(volume, maxval, minval);
    miset_volume_range(volume, maxval * m_Scale + m_Shift, minval * m_Scale + m_Shift);
    }
  miclose_volume(volume);
//   for ( i = 0; i < ndims; i++ )
//     {
//     mifree_dimension_handle(dim[i]);
//     }
}

// void MINCImageIO::XYZFromDirectionCosines(midimhandle_t *hdims, int *dim_indices, size_t *numberOfComponents)
// {
//   midimclass_t dim_class;
//   double       direction_cosines[3];
//   double       dircos[3][3] =
//               { { 1, 0, 0 },
//                 { 0, 1, 0 },
//                 { 0, 0, 1 } };
//
//   // figure out present dimension in the order of either
//   // xspace,yspace,zspace, time or xfrequency,yfrequency,zfrequency, tfrequency
//   // --> x,y,z,t and vector-dimension
//   size_t i = 0;
//   size_t counter = 0;
//   size_t counter2 = 5;
//
//   for ( i = 0; i < this->m_NDims; i++ )
//     {
//     if ( miget_dimension_class(hdims[i], &dim_class) < 0 )
//       {
//       // Error getting dimension class
//       itkDebugMacro("Could not get dim class\"" << m_FileName.c_str() << "\".");
//       return;
//       }
//
//     switch ( dim_class )
//       {
//       case MI_DIMCLASS_SPATIAL:
//       case MI_DIMCLASS_SFREQUENCY:
//         // if none of xspace,yspace or zspace
//         // use direction cosines to figure out which dimension is x,y,z
//         if ( miget_dimension_cosines(hdims[i], direction_cosines) < 0 )
//           {
//           // Error getting dimension direction cosine
//           itkDebugMacro("Could not getdirection cosines!");
//           return;
//           }
//         // fill the dircos array for now, order x,y,z
//         // figure out the order later
//         dircos[0][counter] = direction_cosines[0];
//         dircos[1][counter] = direction_cosines[1];
//         dircos[2][counter] = direction_cosines[2];
//         dim_indices[counter] = i;
//         counter++;
//         break;
//       case MI_DIMCLASS_TIME:
//       case MI_DIMCLASS_TFREQUENCY:
//         dim_indices[3] = i;
//         *numberOfComponents *= this->m_DimensionSize[i];
//         break;
//       // check for vector dimensions
//       case MI_DIMCLASS_RECORD:
//         dim_indices[4] = i;
//         *numberOfComponents *= this->m_DimensionSize[i];
//         break;
//       case MI_DIMCLASS_ANY:
//       case MI_DIMCLASS_USER:
//         dim_indices[counter2] = i;
//         counter2++;
//         break;
//         //default:
//         // any other dimension is ignored!!
//         //return;
//       } // end of switch
//     }   //end of for
//   // fill in the itk matrix for direction cosines
//
//   m_DirectionCosines.Fill(0.0);
//   m_DirectionCosines.SetIdentity();
//
//   //figure out from direction cosines which dimension is what
//   // largest z component (3 dimensions) --> zspace
//   // then largest y component (2 dimension) --> yspace
//   // last remaining dimension xspace.
//   //VF: sorry, this is not helpfull
//
//   int temp;
//   if ( counter == 3 ) // three spatial dimensions
//     {
//     if ( ( dircos[2][0] >= dircos[2][1] ) && ( dircos[2][0] >= dircos[2][2] ) )
//       {
//       // index 0 is the z dimension
//       m_DirectionCosines[0][2] = dircos[0][0];
//       m_DirectionCosines[1][2] = dircos[1][0];
//       m_DirectionCosines[2][2] = dircos[2][0];
//
//       if ( dircos[1][1] >= dircos[1][2] )
//         {
//         m_DirectionCosines[0][1] = dircos[0][1];
//         m_DirectionCosines[1][1] = dircos[1][1];
//         m_DirectionCosines[2][1] = dircos[2][1];
//         m_DirectionCosines[0][0] = dircos[0][2];
//         m_DirectionCosines[1][0] = dircos[1][2];
//         m_DirectionCosines[2][0] = dircos[2][2];
//         }
//       else
//         {
//         temp = dim_indices[1];
//         dim_indices[1] = dim_indices[2];
//         dim_indices[2] = temp;
//         m_DirectionCosines[0][1] =  dircos[0][2];
//         m_DirectionCosines[1][1] =  dircos[1][2];
//         m_DirectionCosines[2][1] =  dircos[2][2];
//         m_DirectionCosines[0][0] =  dircos[0][1];
//         m_DirectionCosines[1][0] =  dircos[1][1];
//         m_DirectionCosines[2][0] =  dircos[2][1];
//         }
//       }
//     else if ( ( dircos[2][1] >= dircos[2][0] ) && ( dircos[2][1] >= dircos[2][2] ) )
//       {
//       m_DirectionCosines[0][2] = dircos[0][1];
//       m_DirectionCosines[1][2] = dircos[1][1];
//       m_DirectionCosines[2][2] = dircos[2][1];
//       temp = dim_indices[0];
//       dim_indices[0] = dim_indices[1];
//       if ( dircos[1][0] >= dircos[1][2] )
//         {
//         dim_indices[1] = temp;
//         m_DirectionCosines[0][1] = dircos[0][0];
//         m_DirectionCosines[1][1] = dircos[1][0];
//         m_DirectionCosines[2][1] = dircos[2][0];
//         m_DirectionCosines[0][0] = dircos[0][2];
//         m_DirectionCosines[1][0] = dircos[1][2];
//         m_DirectionCosines[2][0] = dircos[2][2];
//         }
//       else
//         {
//         dim_indices[1] =  dim_indices[2];
//         dim_indices[2] = temp;
//         m_DirectionCosines[0][1] = dircos[0][2];
//         m_DirectionCosines[1][1] = dircos[1][2];
//         m_DirectionCosines[2][1] = dircos[2][2];
//         m_DirectionCosines[0][0] = dircos[0][0];
//         m_DirectionCosines[1][0] = dircos[1][0];
//         m_DirectionCosines[2][0] = dircos[2][0];
//         }
//       }
//     else
//       {
//       m_DirectionCosines[0][2] = dircos[0][2];
//       m_DirectionCosines[1][2] = dircos[1][2];
//       m_DirectionCosines[2][2] = dircos[2][2];
//       temp = dim_indices[0];
//       dim_indices[0] = dim_indices[2];
//       if ( dircos[1][0] >= dircos[1][1] )
//         {
//         dim_indices[2] = dim_indices[1];
//         dim_indices[1] = temp;
//         m_DirectionCosines[0][1] = dircos[0][0];
//         m_DirectionCosines[1][1] = dircos[1][0];
//         m_DirectionCosines[2][1] = dircos[2][0];
//         m_DirectionCosines[0][0] = dircos[0][1];
//         m_DirectionCosines[1][0] = dircos[1][1];
//         m_DirectionCosines[2][0] = dircos[2][1];
//         }
//       else
//         {
//         dim_indices[2] = temp;
//         m_DirectionCosines[0][0] = dircos[0][0];
//         m_DirectionCosines[1][0] = dircos[1][0];
//         m_DirectionCosines[2][0] = dircos[2][0];
//         m_DirectionCosines[0][1] = dircos[0][1];
//         m_DirectionCosines[1][1] = dircos[1][1];
//         m_DirectionCosines[2][1] = dircos[2][1];
//         }
//       }
//     }
// }
//
} // end namespace itk
