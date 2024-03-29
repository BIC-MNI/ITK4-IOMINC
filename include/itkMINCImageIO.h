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
/**
 *         The specification for this file format is taken from the
 *         web site http://www.bic.mni.mcgill.ca/ServicesSoftware/MINC
 * \author Vladimir S. FONOV
 *         Brain Imaging Center, Montreal Neurological Institute, McGill University, Montreal Canada 2012
 * \author Leila Baghdadi
 *         Mouse Imaging Centre, Toronto, Canada 2005.
 */

#ifndef __itkMINCImageIO_h
#define __itkMINCImageIO_h


#include "itkImageIOBase.h"

#include "itkMatrix.h"

#include <itk_minc2.h>

namespace itk
{
/** \class MINCImageIO
 *
 * \author Leila Baghdadi
 * \brief Class that defines how to read MINC file format. 
 * 
 * \ingroup ITKIOMINC
 *
 * Note, like ITK, MINC is N dimensional and dimensions
 * can be submitted in any arbitrary order. Here we make sure the
 * dimensions are ordered as xspace, yspace, zspace, time and
 * vector_dimension and so on or xfrequencey, yfrequency, zfrequency,
 * tfrequency and vector_dimension and so on NOTE** This class only
 * reads the regularly sampled dimensions as I am not sure how to deal
 * with "iregularly sampled" dimensions yet!
 *
 * This code was contributed in the Insight Journal paper:
 * "MINC2.0 IO Support for ITK"
 * by Baghdadi L.
 * http://hdl.handle.net/1926/191
 * http://www.insight-journal.org/browse/publication/88
 *
 * \ingroup IOFilters
 *
 */
class ITK_EXPORT MINCImageIO:public ImageIOBase
{
public:
  /** Standard class typedefs. */
  typedef MINCImageIO           Self;
  typedef ImageIOBase           Superclass;
  typedef SmartPointer< Self >  Pointer;
  typedef Matrix< float, 3, 3 > MatrixType;
  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(MINCImageIO, ImageIOBase);

  /*-------- This part of the interface deals with reading data. ------ */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanReadFile(const char *);

  /** Set the spacing and dimension information for the set filename. */
  virtual void ReadImageInformation();

  /** Reads the data from disk into the memory buffer provided. */
  virtual void Read(void *buffer);

  /*-------- This part of the interfaces deals with writing data. ----- */

  /** Determine the file type. Returns true if this ImageIO can read the
   * file specified. */
  virtual bool CanWriteFile(const char *);

  /** Writes the spacing and dimensions of the image.
   * Assumes SetFileName has been called with a valid file name. */
  virtual void WriteImageInformation();

  /** Writes the data to disk from the memory buffer provided. Make sure
   * that the IORegion has been set properly. */
  virtual void Write(const void *buffer);

protected:
  MINCImageIO();
  ~MINCImageIO();
  void PrintSelf(std::ostream & os, Indent indent) const;
  void WriteSlice(std::string & fileName, const void *buffer);

  int  m_NDims; /*Number of dimensions*/

  // dimension size and start and step, in FILE ORDER!
  
  const char  **m_DimensionName; 
  misize_t     *m_DimensionSize;
  double       *m_DimensionStart;
  double       *m_DimensionStep;
  int           m_DimensionIndices[5];
  midimhandle_t *m_MincFileDims;
  midimhandle_t *m_MincApparentDims;
  mitype_t      m_volume_type;
  miclass_t     m_volume_class;
  
  
  // MINC2 volume handle , currently opened
  mihandle_t   m_volume;

  MatrixType m_DirectionCosines;
  // complex type images, composed of complex numbers
  //int m_Complex;
  
  // will assign m_NDims and allocate all internal buffers to hold the information
  void AllocateDimensions(int nDims);
  
  // cleanup internal buffers
  void CleanupDimensions(void);
  
  // close existing volume, cleanup internal structures
  void CloseVolume(void);

private:
  MINCImageIO(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented

};
} // end namespace itk

#endif // __itkMINCImageIO_h
