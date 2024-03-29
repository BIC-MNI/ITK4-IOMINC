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
#ifndef __itkMINCImageIOFactory_h
#define __itkMINCImageIOFactory_h

#include "itkObjectFactoryBase.h"
#include "itkImageIOBase.h"

namespace itk
{
/** \class MINCImageIOFactory
 * \brief Create instances of MINCImageIO objects using an object factory.
 *
 * \ingroup ITKIOMINC
 *
 * This code was contributed in the Insight Journal paper:
 * "MINC.0 IO Support for ITK"
 * by Baghdadi L.
 * http://hdl.handle.net/1926/191
 * http://www.insight-journal.org/browse/publication/88
 *
 */
class ITK_EXPORT MINCImageIOFactory:public ObjectFactoryBase
{
public:
  /** Standard class typedefs. */
  typedef MINCImageIOFactory         Self;
  typedef ObjectFactoryBase          Superclass;
  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Class methods used to interface with the registered factories. */
  virtual const char * GetITKSourceVersion(void) const;

  virtual const char * GetDescription(void) const;

  /** Method for class instantiation. */
  itkFactorylessNewMacro(Self);
  static MINCImageIOFactory * FactoryNew() { return new MINCImageIOFactory; }
  /** Run-time type information (and related methods). */
  itkTypeMacro(MINCImageIOFactory, ObjectFactoryBase);

  /** Register one factory of this type  */
  static void RegisterOneFactory(void)
  {
    MINCImageIOFactory::Pointer MINCFactory = MINCImageIOFactory::New();

    ObjectFactoryBase::RegisterFactory(MINCFactory);
  }

protected:
  MINCImageIOFactory();
  ~MINCImageIOFactory();

private:
  MINCImageIOFactory(const Self &); //purposely not implemented
  void operator=(const Self &);      //purposely not implemented
};
} // end namespace itk

#endif
